import copy
import pandas as pd
import numpy as np
import Levenshtein
class SeqInfo(object):
    """Holds information regarding the sequence.
    
    """
    def __init__(self, seq, 
                 protein,
                 accession,
                 name=None,
                 subtype=None,
                 host=None, 
                 date=None, 
                 erisk=None,
                 irisk=None,
                 risk_flag=None,
                 country=None):
        self.name = name
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
    """Holds information regarding the sequences in the records.
    
    Args:
        records (list): list of records parsed from NCBI
        cov19_accessions (list): of accessions corresponding to cov19
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
                erisk=record.predicted_emergence_score,
                irisk=record.predicted_impact_score,
                risk_flag = record.predicted_emergence_score > self.risk_threshold,
                host=None,
                date=None,
                country=None)
            #print(record.predicted_emergence_score > self.risk_threshold)
            self.seq_infos[seqinfo.accession] = seqinfo
            
    
    def compute_L_diatance_matrix(self):
        highriskseq = pd.DataFrame.from_dict({key:val.seq 
                                              for (key,val) in self.seq_infos.items() 
                                              if val.risk_flag},orient='index',columns=['seq'])
        num=highriskseq.index.size
        d=np.zeros([num,num])
        for i in np.arange(num):
            for j in np.arange(num):
                if i > j:
                    d[i,j] = Levenshtein.distance(highriskseq.seq.values[i],
                                                  highriskseq.seq.values[j])
                    ds=pd.DataFrame(d)        
                    ds=(ds+ds.transpose())
                    ds.columns=highriskseq.index.values
                    self.highriskdistancematrix=ds.copy()
        return ds
    
    
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
    


#N=10000
df=pd.read_csv('./combined_results.csv',index_col=0)
#df1=df[['subtype','predicted_impact_score', 'predicted_emergence_score', 'ha', 'na']]
#df1=df1.sort_values('predicted_emergence_score',ascending=False).head(N)
#df1=df1[df1.predicted_emergence_score>6.2]
#df1.subtype.value_counts()

ALLinfoHA=MultipleSeqInfo(df.reset_index(),'ha_accession','ha',risk_threshold=6.2)
ds=ALLinfoHA.compute_L_diatance_matrix()


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


ds.to_csv('dm.csv',index=None)





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



PHYLO_DIR='./'

Phylo.convert(
    PHYLO_DIR + 'ldistanceh1n1.xml','phyloxml',
    PHYLO_DIR + 'ldistance.nhx','newick')

ltree = load_tree(
    PHYLO_DIR + 'ldistance.nhx',
    type_='newick')



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
        node.host = recordinfo.seq_infos[name].subtype
        print(node.name,node.host)
    return tree


labelled_tree=label_nodes(
    ltree, ALLinfoHA)

import copy



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

def collapse(tree, min_dist,AllrecordInfo):
    # cache the tip content of each node to reduce the number of times the tree is traversed
    
    tree = copy.deepcopy(tree)
    
    node2tips = tree.get_cached_content()
    root_distance = cache_distances(tree)

    for node in tree.get_descendants('preorder'):
        if not node.is_leaf():
            avg_distance_to_tips = mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])
            print(avg_distance_to_tips)
            if avg_distance_to_tips < min_dist:
                # do whatever, ete support node annotation, deletion, labeling, etc.

                # rename
                #                node.name += ' COLLAPSED avg_d:%g {%s}' %(avg_distance_to_tips,
                #                                                 ','.join([tip.name for tip in node2tips[node]]))
                #node.name += '{%s}' %(list(node2tips[node])[-1].name)
                #node.name = 'avg_d:%g' %(avg_distance_to_tips)
                # label
            
                closest_name = closest_node(node, node2tips, root_distance).name
                node.host = AllrecordInfo.seq_infos[closest_name].subtype
                node.name = '%s (%g)' %(closest_name,avg_distance_to_tips)
                
            
                node.add_features(collapsed=True)

                # set drawing attribute so they look collapsed when displayed with tree.show()
                node.img_style['draw_descendants'] = False

    return tree
# etc...



ltree_collapsed = collapse(
    labelled_tree, 
    min_dist=5, 
    AllrecordInfo=ALLinfoHA)

prune_nodes(ltree_collapsed)


# COLBAT='DarkRed'
# COLRAT='SteelBlue'
COLHUMAN='DarkGreen'
COLCOVID='DarkRed'
COLBAT='Red'
COLRAT='Blue'
COLCAMEL='Purple'
COLGAME='Red'
COLCATTLE='Yellow'
# COLHUMAN='Black'
FS=70
PW=20



def nodeAttribConstruct(color, node):
    N = AttrFace(
        "name", fsize=FS, 
        text_prefix=" ",penwidth=PW,ftype='Arial',
        fgcolor=color,fstyle='bold')
    faces.add_face_to_node(N, node, 1, position="branch-right")
    return N

def layout(node):
    if node.is_leaf():
        if  node.host == 'H1N1':
            N = nodeAttribConstruct(COLBAT,node)
        elif node.host == 'H3N2':
            N = nodeAttribConstruct(COLRAT,node)
        elif node.host == 'H7N9':
            N = nodeAttribConstruct(COLHUMAN,node)
        elif node.host == 'H1N2':
            N = nodeAttribConstruct(COLCATTLE,node)
        elif node.host == 'game':
            N = nodeAttribConstruct(COLGAME,node)
        elif node.host == 'camel':
            N = nodeAttribConstruct(COLCAMEL,node)
        else:
            N = nodeAttribConstruct(COLGAME,node)
            
            
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
    ts.scale=20
    ts.show_scale=True
    # ts.show_branch_length=True
    #ts.min_leaf_separation=10
    #ts.optimal_scale_level='full'
    #ts.branch_vertical_margin=0
    
    ns.hz_line_width=2
    ns.vt_line_width=1
    #ts.layout_fn = layout
    ns["vt_line_width"] = 5
    ns["hz_line_width"] = 5
    #     ns['fsize'] = 20
    for n in tree.traverse():
        n.set_style(ns)
        
    #all_accessions = all_seq_data['accessions'].values
    for n in tree:
        ts.layout_fn = layout

        
    tree.set_style(ns)
    #tree.set_style(ts)
    
    #t.show()
    tree.render(
        outfile, 
        dpi=300, 
        h=50,
        tree_style=ts)


render_tree(
    ltree_collapsed, './tmp.pdf')
#    
