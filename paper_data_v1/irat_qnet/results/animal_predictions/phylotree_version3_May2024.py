#!/usr/bin/env python
# coding: utf-8

# In[1]:


import copy
import pandas as pd
import numpy as np
import Levenshtein
from tqdm import tqdm


# In[2]:


mf=pd.read_csv('../../../../paper_data_v3/irat_enet/results/animal_predictions/animal_combined_with_irat.csv')#.assign(is_irat=0)
mf=mf.rename(columns={'name':'id','ha_sequence':'ha','na_sequence':'na','emergence':'emergence_risk'}).assign(impact_risk=0)
mf.to_csv('v3_animalrisk.csv')


# In[3]:


# parameters
threshold=6.25  # draw tree above this emergenec risk threshold
CONSTRUCT_PHYLO=False  # use constructed tree (tree construction takes hours)
OUTPUT_DIR='./'
PHYLO_TREE_DIR='./'
PHYLO_DIR='./'
COMBINED_RESULTS='v3_animalrisk.csv'
#COMBINED_RESULTS='/home/ishanu/Downloads/combined_results_irat.csv'
#CONSTRUCT_PHYLO=False
num_collapsed=19  # number of mutations within which leaves are collapsed
VERBOSE=False


# In[4]:


mf


# In[5]:


tmpf=pd.read_csv(COMBINED_RESULTS)
IRATSEQ=tmpf[tmpf.is_irat==1].id.values.astype(str)


# # class definition to hold multi sequence data

# In[6]:


class SeqInfo(object):
    """Holds information regarding the sequence.
    
    """
    def __init__(self, seq, 
                 protein,
                 accession,
                 subtype=None,
                 id=None,
                 name=None,
                 host=None, 
                 date=None, 
                 erisk=None,
                 irisk=None,
                 risk_flag=None,
                 country=None):
        self.name = name
        self.id = id
        self.protein=protein
        self.subtype=subtype        
        self.seq = seq
        self.accession = accession 
        self.host = host
        self.date = date
        self.erisk = erisk
        self.irisk = irisk
        self.risk_flag = risk_flag
        self.country = country
        
class MultipleSeqInfo(object):
    """Holds information regarding multiple sequences.
    
    Args:
        dataframe (pandas.DataFrame): list of records parsed from NCBI
        accessionname (str): column name for accession id
        proteinname (str): protein name 
        risk_threshold (float): emergence risk threshold to compute distance matrix
    """
    def __init__(self,
                 dataframe,
                 accessionname,
                 proteinname,
                 risk_threshold=6.2):
        
        self.seq_infos = {}
        self.risk_threshold = risk_threshold
        for i in np.arange(dataframe.index.size):
            record=dataframe.iloc[i,:]
            seqinfo = SeqInfo(
                name=record.id,
                seq=record[proteinname], 
                protein=proteinname,
                accession=record[accessionname],
                subtype=record.subtype,
                erisk=record.emergence_risk,
                irisk=record.impact_risk,
                risk_flag = record.emergence_risk > self.risk_threshold,
                host=None,
                date=None,
                country=None)
            #print(record.predicted_emergence_score > self.risk_threshold)
            self.seq_infos[seqinfo.accession] = seqinfo
            
    
    def compute_L_distance_matrix(self):
        highriskseq = pd.DataFrame.from_dict({key:val.seq 
                                              for (key,val) in self.seq_infos.items() 
                                              if val.risk_flag},orient='index',columns=['seq'])
        num=highriskseq.index.size
        d=np.zeros([num,num])
        for x in tqdm(np.arange(num*num)):
            j=x//num
            i=x-num*j
            if i > j:
                d[i,j] = Levenshtein.distance(highriskseq.seq.values[i],
                                                  highriskseq.seq.values[j])
        ds=pd.DataFrame(d)        
        ds=(ds+ds.transpose())
        ds.columns=highriskseq.index.values
        self.highriskdistancematrix=ds.copy()
        
        self.highriskdistancematrix.to_csv('dm'+str(self.risk_threshold)+'.csv',index=None)
        return 
    
    
    def accessions_to_subtype(self, accessions):
        """Create a dictionary mapping the accession to the host.
        """
        
        subtypes = []
        for accession in accessions:
            seqinfo = self.seq_infos[accession]
            subtypes.append(seqinfo.subtype)
            
        return subtypes

    def accessions_to_host(self, accessions):
        """Create a dictionary mapping the accession to the host.
        """
        
        hosts = []
        for accession in accessions:
            seqinfo = self.seq_infos[accession]
            hosts.append(seqinfo.host)
        return hosts
           

df=pd.read_csv(COMBINED_RESULTS,index_col=0).reset_index()
ALLinfoHA=MultipleSeqInfo(df.reset_index(),'ha_accession','ha',risk_threshold=threshold)


# In[8]:


df[df.emergence_risk>threshold].subtype.value_counts()


#ALLinfoHA.compute_L_distance_matrix()

import dill as pickle
# Pickling (serializing) the object
with open('ALLinfo.pkl', 'wb') as file:
    pickle.dump(ALLinfoHA, file)

# Unpickling (deserializing) the object
with open('ALLinfo.pkl', 'rb') as file:
    ALLinfoHA = pickle.load(file)


from Bio.Phylo import TreeConstruction
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import Entrez
from Bio import SeqIO


def load_dm(file_, upper_diag=True):
    """Load the distance matrix. 
    
    Also, do some preprocessing. 
    """
    
    df = pd.read_csv(file_)
    #df.set_index('Unnamed: 0', inplace=True)
    #assert np.all(df.columns == df.index)
    
    # drop duplicate columns after reading csv
    #df = df.loc[:, ~df.columns.str.replace("(\.\d+)$", "").duplicated()]
    
    if upper_diag:
        df = df + df.T
    return df

def save_tree(tree, file_name, save_type='xml'):
    """Saved the created phylogenetic tree."""
    
    if save_type == 'pickle':
        graph = Phylo.to_networkx(tree)
        save_pickled(graph, file_name)
    elif save_type == 'xml':
        Phylo.write(tree, file_name, 'phyloxml')
    else:
        raise ValueError('Not a correct save type.')
    
def pandas_dm_to_biopython_dm(dm):
    """Convert the pandas distance matrix to the biopython distance matrix.
    
    Returns:
        biopython distance matrix
    """
    
    accessions = dm.columns
    bio_dm = []
    for i, accession in enumerate(accessions):
        bio_dm.append(list(dm.iloc[i, :i+1].values))
        
    bio_dm = TreeConstruction._DistanceMatrix(
        list(dm.columns), 
        bio_dm)
    
    return bio_dm

def distance_matrix_to_phylo_tree(dm, outfile=None):
    """Create a phylogenetic tree from the distance matrix."""
    
    dm = pandas_dm_to_biopython_dm(dm)
    
    treeConstructor = TreeConstruction.DistanceTreeConstructor()
    tree = treeConstructor.nj(dm)
    
    if outfile is not None:
        save_tree(tree, outfile)


from ete3 import Tree, TreeStyle
from ete3 import Phyloxml
from ete3 import AttrFace, faces, Tree, NodeStyle, TreeStyle

def load_pickled(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f, encoding='latin')


def get_farthest_node(tree, sequence):
    return (tree&sequence).get_farthest_node()

def get_all_accessions_from_tree(tree):
    return [leaf_node.name for leaf_node in tree.get_leaves()]

def remove_certain_hosts_from_tree(tree, hosts):
    """Remove leaf nodes if the host of that leaf is in `hosts`"""
    
    tree = copy.deepcopy(tree)
    
    removed_accessions = []
    for leaf_node in tree.get_leaves():
        if leaf_node.host in hosts:
            leaf_node.detach()
            
    return tree

def set_midpoint_outgroup(tree):
    tree.set_outgroup(tree.get_midpoint_outgroup())


def load_tree(filename, type_='phyloxml'):
    """Load saved phylogenetic tree.
    """
    
    if type_ == 'phyloxml':
        project = Phyloxml()
        project.build_from_file(filename)

        for tree in project.get_phylogeny():
            break

        t=tree
        
    elif type_ == 'newick':
        t = Tree(filename, format=1)
    else:
        raise ValueError('Not a correct type.')
    
    return t



#CONSTRUCT_PHYLO=True
if CONSTRUCT_PHYLO:
    ALL_dm_ldistance = load_dm(
        OUTPUT_DIR + 'dm'+str(threshold)+'.csv', 
        upper_diag=False)
    
    distance_matrix_to_phylo_tree(
        ALL_dm_ldistance, PHYLO_TREE_DIR + 'ldistance'+str(threshold)+'.xml')


# # convert phyloxml tree to newick tree to manipulate trees

# In[13]:


Phylo.convert(
    PHYLO_DIR + 'ldistance'+str(threshold)+'.xml','phyloxml',
    PHYLO_DIR + 'ldistance'+str(threshold)+'.nhx','newick')

ltree = load_tree(
    PHYLO_DIR + 'ldistance'+str(threshold)+'.nhx',
    type_='newick')


# # label nodes in tree to add other attributes like subtype risk etc

# In[14]:


def bandify(val,min=6.05,max=6.71,flag=False):
    maptoten=int(np.ceil(((val-min)/(max-min))*10))
    L=' '+u'\u2580'*maptoten
    
    if flag:
        L=L+u'\u21DD'
    return L

def bandify(val,min=6,max=8,flag=False):
    maptoten=int(np.ceil(((val-min)/(max-min))*10))
    L=' '+u'\u2501'*maptoten
    if flag:
        L=L+u'\u2605'
    L=L+' '
    return L

def label_nodes(
        tree, 
        recordinfo):
    """Label the nodes of the tree.
    
    We label nodes on whether:
        it is covid19
    """
    
    tree = copy.deepcopy(tree)
    
    for node in tree:
        name = node.name      
        node.subtype = recordinfo.seq_infos[name].subtype
        node.erisk =recordinfo.seq_infos[name].erisk
        node.id = recordinfo.seq_infos[name].name + bandify(recordinfo.seq_infos[name].erisk,min=threshold)
        if VERBOSE:
            print(node.name,node.subtype,node.id,node.erisk)
    return tree


# # construct labelled tree

# In[15]:


labelled_tree=label_nodes(
    ltree, ALLinfoHA)


# # functions to collapse similar leaves

# In[16]:


def prune_nodes(t):
    # collapsed nodes are labeled, so you locate them and prune them
    for n in t.search_nodes(collapsed=True):
        for ch in n.get_children():
            ch.detach()
            
            
def mean(array):
    return sum(array)/float(len(array))

def cache_distances(tree):
    ''' precalculate distances of all nodes to the root''' 
    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist

def closest_node(node, node2tips, root_distance):
    """Find the closest node."""
    
    tips = []
    distances = []
    for tip in node2tips[node]:
        distances.append(root_distance[tip]-root_distance[node])
        tips.append(tip)
        #     index = np.argmin([root_distance[tip]-root_distance[node] for tip in node2tips[node]])
    index = np.argmin(distances)
    return tips[index]

def riskiest_node(node, node2tips):
    """Find the riskiest node."""
    
    tips = []
    risks = []
    for tip in node2tips[node]:
        risks.append(tip.erisk)
        tips.append(tip)
        #     index = np.argmin([root_distance[tip]-root_distance[node] for tip in node2tips[node]])
    index = np.argmax(risks)
    return tips[index]


def all_collapsed(node, node2tips,AllrecordInfo):
    """Find all nodes in collapsed set."""
    
    tips = []
    for tip in node2tips[node]:
        tips.append(AllrecordInfo.seq_infos[tip.name].name)
    return tips

def all_collapsed(node, node2tips,AllrecordInfo):
    """Find all nodes in collapsed set."""
    
    tips = []
    for tip in node2tips.get(node, []):  # use get() to avoid KeyError if node has no entries in node2tips
        tips.append(AllrecordInfo.seq_infos[tip.name].name)
        
    # if no tips were found and the node is a leaf, append the node itself
    if not tips and node.is_leaf():
        tips.append(AllrecordInfo.seq_infos[node.name].name)

    return tips


def collapse(tree, min_dist,AllrecordInfo):
    # cache the tip content of each node to reduce the number of times the tree is traversed
    
    tree = copy.deepcopy(tree)
    
    node2tips = tree.get_cached_content()
    root_distance = cache_distances(tree)

    for node in tree.get_descendants('preorder'):
        IRAT=False
        if not node.is_leaf():
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])
            if VERBOSE:
                print(avg_distance_to_tips)
            if avg_distance_to_tips <= min_dist:
                # do whatever, ete support node annotation, deletion, labeling, etc.
            
                #closest_name = closest_node(node, node2tips, root_distance).name
                
                # find if this name is in IRATSEQ ie any IRAT seq is in teh collapsed nosdees
                all_collapsed_nodes = all_collapsed(node, node2tips,AllrecordInfo)
                #print(all_collapsed_nodes)
                
                
                
                for i in IRATSEQ:
                    if i in all_collapsed_nodes:
                        IRAT=True
                        break
                                       
                closest_name = riskiest_node(node, node2tips).name
                node.subtype = AllrecordInfo.seq_infos[closest_name].subtype
                node.id = AllrecordInfo.seq_infos[closest_name].name + bandify(AllrecordInfo.seq_infos[closest_name].erisk,min=threshold,flag=IRAT) 
                node.name = '%s (%g)' %(closest_name,avg_distance_to_tips)
                
            
                node.add_features(collapsed=True)

                # set drawing attribute so they look collapsed when displayed with tree.show()
                node.img_style['draw_descendants'] = False
        else:
            all_collapsed_nodes = all_collapsed(node, node2tips,AllrecordInfo)
            for i in IRATSEQ:
                if i in all_collapsed_nodes:
                    IRAT=True
                    break
            node.id = AllrecordInfo.seq_infos[node.name].name + bandify(AllrecordInfo.seq_infos[node.name].erisk,min=threshold,flag=IRAT) 

    return tree


# In[17]:


# collapse leaved


# In[18]:


num_collapsed=20
ltree_collapsed = collapse(
    labelled_tree, 
    min_dist=num_collapsed, 
    AllrecordInfo=ALLinfoHA)

prune_nodes(ltree_collapsed)


# 
# # code for actual rendering

# In[19]:


# COLBAT='DarkRed'
# COLRAT='SteelBlue'
COLH3N2='#0077FF'
COLH1N1='#551177'
COLH1N2='#337733'
COLH5N1='#991133'
COLH7N9='#aa6622'
COLH9N2='#8877FF'
COLH3N3='#FF0000'
COLDEF='#aaaaaa'

FS=50
PW=10

def nodeAttribConstruct(color, node):
    N = AttrFace(
        "id", fsize=FS, 
        text_prefix=" ",penwidth=PW,ftype='Arial',
        fgcolor=color,fstyle='bold')
    faces.add_face_to_node(N, node, 1, position="branch-right")
    return N

def layout(node):
    if node.is_leaf():
        if  node.subtype == 'H1N1':
            N = nodeAttribConstruct(COLH1N1,node)
        elif node.subtype == 'H3N2':
            N = nodeAttribConstruct(COLH3N2,node)
        elif node.subtype == 'H3N3':
            N = nodeAttribConstruct(COLH3N3,node)
        elif node.subtype == 'H7N9':
            N = nodeAttribConstruct(COLH7N9,node)
        elif node.subtype == 'H9N2':
            N = nodeAttribConstruct(COLH9N2,node)
        elif node.subtype == 'H1N2':
            N = nodeAttribConstruct(COLH1N2,node)
        elif node.subtype == 'H5N1':
            N = nodeAttribConstruct(COLH5N1,node)
        else:
            N = nodeAttribConstruct(COLDEF,node)
            



            
def render_tree(tree, outfile):# all_seq_data, display_type='nearest_host'):
    """Render the tree inside the file to a circular 
    phylogenetic tree.
    
    NOTE: outfile should be in .pdf for best visuals
    Returns:
    """
    #tree = Tree(nwfile,format=1)

    ts = TreeStyle()
    ns = NodeStyle()
    ts.show_leaf_name = False
    #ts.rotation = 90
    ts.mode = "r"
    #ts.arc_start = -360 # 0 degrees = 3 o'clock
    #ts.arc_span = 360
    ts.scale=4
    ts.show_scale=False
    ts.branch_vertical_margin = .5 # 10 pixels between adjacent branches
    # ts.show_branch_length=True
    #ts.min_leaf_separation=10
    #ts.optimal_scale_level='full'
    #ts.branch_vertical_margin=0
    
    ns.hz_line_width=2
    ns.vt_line_width=1
    #ts.layout_fn = layout
    ns["vt_line_width"] = 16
    ns["hz_line_width"] = 16
    #     ns['fsize'] = 20
    for n in tree.traverse():
        n.set_style(ns)
        
    #all_accessions = all_seq_data['accessions'].values
    for n in tree:
        ts.layout_fn = layout

        
    tree.set_style(ns)
    tree.set_style(ts)
    
    ax=tree.render(
        outfile, 
        dpi=300, 
        w=500,
        tree_style=ts)

r=5
medrisknames=[]
for node in ltree_collapsed:
    if ALLinfoHA.seq_infos[node.name.split()[0]].erisk > r:
        medrisknames=np.append(medrisknames,
                                ALLinfoHA.seq_infos[node.name.split()[0]].name)
medrisknames=df[df.id.isin(medrisknames)][['id','subtype',
                                               'ha_accession',
                                               'na_accession',
                                               'impact_risk',
                                               'emergence_risk']].sort_values('emergence_risk',
                                                                                         ascending=False)
medrisknames.to_csv('allriskystrains_collapsed.csv',index=None)
medrisknames[medrisknames.id.str.contains('Ohio')]

ax=render_tree(
    ltree_collapsed,
    'riskyphylo'+str(threshold)+'_collapsed_'+str(num_collapsed)+'.pdf')


import dill as pickle
# Pickling (serializing) the object
with open('tree_collapsed.pkl', 'wb') as file:
    pickle.dump(ltree_collapsed, file)

# Unpickling (deserializing) the object
with open('tree_collapsed.pkl', 'rb') as file:
    ltree_collapsed = pickle.load(file)


r=6.057
r=6.5
highrisknames=[]
for node in ltree_collapsed:
    if ALLinfoHA.seq_infos[node.name.split()[0]].erisk > r:
        highrisknames=np.append(highrisknames,
                                ALLinfoHA.seq_infos[node.name.split()[0]].name)
highrisknamesdf=df[df.id.isin(highrisknames)][['id','subtype',
                                               'ha_accession',
                                               'na_accession',
                                               'impact_risk',
                                               'emergence_risk']].sort_values('emergence_risk',
                                                                                         ascending=False)
highrisknamesdf = highrisknamesdf.rename(columns={'id':'strain',
                                                  'ha_accession':'HA accession',
                                                  'na_accession':'NA accession',
                                                  'impact_risk':'predicted IRAT impact',
                                                  'emergence_risk':'predicted IRAT emergence'}).set_index('strain')
highrisknamesdf#.drop_duplicates()


# In[26]:


COLDICT={'H1N1':'colh1n1x!40',
         'H3N2':'colh3n2x!30',
         'H7N9':'colh7n9x!40','H9N2':'colh9n2x!30',
         'H1N2':'colh1n2x!30','H5N1':'colh5n1x!20',
         'H3N3':'colh3n3x!50'}
def rowcolor(row):
    return '\\rowcolor{' + COLDICT[row.subtype]+'}' + row['strain']

highrisknamesdf1 = highrisknamesdf.reset_index(drop=False)
highrisknamesdf1['strain']=highrisknamesdf1.apply(rowcolor,axis=1)
highrisknamesdf1=highrisknamesdf1.set_index('strain')

highrisknamesdf1=highrisknamesdf1.drop('predicted IRAT impact',axis=1)


from zedstat.textable import textable
textable(highrisknamesdf1.head(68),
         tabname='../../../../tex/overleaf3/Figures/tabdata/highrisk_v3.tex',
         FORMAT='%1.4f',INDEX=True,
         TABFORMAT='L{2.75in}|L{.45in}|L{.60in}|L{.6in}|C{1in}',LNTERM='\\\\\n')



#'H1N1', 'H1N2', 'H3N2', 'H3N3', 'H5N1', 'H7N9'
count=0
Subtype={'H1N1':0,'H3N2':0,'H7N9':0,'H9N2':0,'H5N1':0,'H1N2':0,'H3N3':0}
MaxriskStrain={'H1N1':None,'H3N2':None,'H7N9':None,'H9N2':None,'H5N1':None,'H1N2':None,'H3N3':None}
Subtype_strat={ 6.5:{'H1N1':0,'H3N2':0,'H7N9':0,'H9N2':0,'H5N1':0,'H1N2':0,'H3N3':0},
               7:{'H1N1':0,'H3N2':0,'H7N9':0,'H9N2':0,'H5N1':0,'H1N2':0,'H3N3':0},   
               7.5:{'H1N1':0,'H3N2':0,'H7N9':0,'H9N2':0,'H5N1':0,'H1N2':0,'H3N3':0},
               7.7:{'H1N1':0,'H3N2':0,'H7N9':0,'H9N2':0,'H5N1':0,'H1N2':0,'H3N3':0}}
for node in ltree_collapsed:
    Subtype[node.subtype]=Subtype[node.subtype]+1
    if MaxriskStrain[node.subtype] is None:
        MaxriskStrain[node.subtype]=(ALLinfoHA.seq_infos[node.name.split()[0]].name,
                                     ALLinfoHA.seq_infos[node.name.split()[0]].erisk)
    else:
        if MaxriskStrain[node.subtype][1]<ALLinfoHA.seq_infos[node.name.split()[0]].erisk:
            MaxriskStrain[node.subtype]=(ALLinfoHA.seq_infos[node.name.split()[0]].name,
                                         ALLinfoHA.seq_infos[node.name.split()[0]].erisk)
    for r in [6.5,7,7.5,7.7]:
        if r<ALLinfoHA.seq_infos[node.name.split()[0]].erisk:
            #print(node.subtype,ALLinfoHA.seq_infos[node.name.split()[0]].erisk)
            Subtype_strat[r][node.subtype]=Subtype_strat[r][node.subtype]+1
        
    count=count+1
maxriskdf = pd.DataFrame(MaxriskStrain)
Subtype_strat_df = pd.DataFrame(Subtype_strat)
print(count)

dcf=100*(Subtype_strat_df/Subtype_strat_df.sum(axis=0))
df_=Subtype_strat_df
dcfs=' ('+dcf.round(2).astype(str)+'\%)'
df__=df_.astype(str)+dcfs


# In[49]:


df__.columns=df__.columns.astype(str)
df__.index.name='subtype'
df__.columns=['score $>$'+x for x in df__.columns]

# In[55]:


from zedstat.textable import textable
#Subtype_strat_df.columns=[str(x) for x in Subtype_strat_df.columns]
#Subtype_strat_df.index.name='subtype'
#Subtype_strat_df=Subtype_strat_df.reset_index()
textable(df__,#pd.concat([highrisknamesdf1.head(30),highrisknamesdf1.tail(10)]),
         tabname='../../../../tex/overleaf3/Figures/tabdata/riskycount_v3.tex',
         FORMAT='%s',INDEX=True,
         TABFORMAT='L{.6in}|L{.9in}|L{.9in}|L{.9in}|L{.9in}',LNTERM='\\\\\n')

