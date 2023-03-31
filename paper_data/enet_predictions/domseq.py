import numpy as np
import pandas as pd
from Bio import SeqIO
from Levenshtein import distance
from quasinet.qnet import Qnet, qdistance_matrix, save_qnet, load_qnet, membership_degree
from sklearn.cluster import KMeans
from sklearn.manifold import MDS


class DomSeq(object):
    """Find and predict dominant sequences.

    Parameters
    ----------
    seq_trunc_length : int
        Length to truncate sequences in Emergenet analysis
        (Sequences used to train Emergenet and compute E-distance must be of same length)

    random_state : int
        Sets seed for random number generator
    """

    def __init__(self, seq_trunc_length, random_state=None):
        if seq_trunc_length <= 0:
            raise ValueError('Length to truncate sequences must be positive!')
        self.seq_trunc_length = seq_trunc_length

        if random_state < 0:
            raise ValueError('Seed must be between 0 and 2**32 - 1!')
        self.random_state = random_state

    def __repr__(self):
        return "emergenet.DomSeq"

    def __str__(self):
        return self.__repr__()
    
    @staticmethod
    def _count_seqs(filepath):
        """Returns number of sequences in a fasta file.

        Parameters
        ----------
        filepath : str
            File name

        Returns
        -------
        seq_df : int
            Number of sequences
        """
        with open(filepath, 'r') as f:
            fasta = SeqIO.parse(f, 'fasta')
            if not any(fasta):
                raise ValueError('The infile must be in fasta format!')
        with open(filepath, 'r') as f:
            lines = f.read()
        num_sequences = lines.count('>')
        return num_sequences
    
    @staticmethod
    def _get_acc(ID):
        """Returns accession code of a sequence.

        Parameters
        ----------
        ID : str
            Sequence metadata

        Returns
        -------
        acc : str
            Accession code
        """
        ids = ID.split('|')
        # NCBI format
        if len(ids) == 2: 
            acc = ids[0]
            return acc
        # GISAID format
        else: 
            acc = ids[4]
            return acc

    @staticmethod
    def _get_name(ID):
        """Returns name of a sequence.

        Parameters
        ----------
        ID : str
            Sequence metadata

        Returns
        -------
        name : str
            Name
        """
        ids = ID.split('|')
        # NCBI format
        if len(ids) == 2: 
            name = ''
            start = False
            for c in ids[1]:
                if c == '(' and not start:
                    start = True
                elif c == '(' and start:
                    break
                elif start == True:
                    name += c
            return name
        # GISAID format
        else: 
            name = ids[0]
            return name
    
    def _parse_fasta(self, filepath):
        """Parses a fasta file and returns DataFrame with three columns.
        `acc` contains sequence accession.
        `name` contains sequence name.
        `sequence` contains sequences truncated to `seq_trunc_length`.

        Parameters
        ----------
        filepath : str
            File name

        Returns
        -------
        seq_df : pd.DataFrame
            DataFrame of sequences
        """
        if self._count_seqs(filepath) == 0:
            raise ValueError('The file contains no sequences!')
        accs = []
        names = []
        seqs = []
        for record in SeqIO.parse(filepath, 'fasta'):
            if len(record.seq) < self.seq_trunc_length:
                continue
            accs.append(self._get_acc(str(record.description)))
            names.append(self._get_name(str(record.description)))
            seqs.append(''.join(record.seq[:self.seq_trunc_length].upper()))
        seq_df = pd.DataFrame({'acc':accs, 'name':names, 'sequence':seqs})
        return seq_df

    def _sequence_array(self, seq_df, sample_size=None):
        """Extracts array of sequence arrays from DataFrame; includes target sequence.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame containing sequences

        sample_size : int
            Number of strains to sample randomly

        Returns
        -------
        seq_lst: numpy.ndarray
            Array of sequence arrays
        """
        if 'sequence' not in seq_df.columns:
            raise ValueError('The DataFrame must store sequences in `sequence` column!')
        if sample_size is None or sample_size > len(seq_df):
            sample_size = len(seq_df)
        if sample_size < len(seq_df):
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        seqs = seq_df['sequence'].apply(lambda x: np.array(list(x[:self.seq_trunc_length])))
        seq_lst = []
        for seq in seqs:
            seq_lst.append(seq)
        seq_lst = np.array(seq_lst)
        return seq_lst

    def load_data(self, filepath, outfile=None):
        """Loads fasta file data and optionally saves to CSV.

        Parameters
        ----------
        filepath : str
            File name

        outfile : str
            File name to save to ('.csv')

        Returns
        -------
        seq_df : pd.DataFrame
            DataFrame of sequences
        """
        seq_df = self._parse_fasta(filepath)
        if outfile is not None:
            if not outfile.endswith('.csv'):
                raise ValueError('The outfile must end with `.csv`!')
            seq_df.to_csv(outfile, index=False)
        return seq_df

    def train(self, seq_df, sample_size=None, n_jobs=1):
        """Trains an Emergenet model.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        sample_size : int
            Number of strains to train Emergenet on, sampled randomly

        n_jobs : int
            Number of CPUs to use when training

        Returns
        -------
        enet : Qnet
            Trained Emergenet
        """
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        seq_arr = self._sequence_array(seq_df, sample_size)
        enet = Qnet(feature_names=['x' + str(i) for i in np.arange(self.seq_trunc_length)],
                    random_state=self.random_state, n_jobs=n_jobs)
        enet.fit(seq_arr)
        return enet
    
    def compute_domseq(self, seq_df, sample_size=None):
        """Computes the dominant sequence in the edit (Levenshtein) distance metric.
        Also returns DataFrame of all sequences and their total edit distances.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        sample_size : int
            Number of strains to compute dominant strain with, sampled randomly

        Returns
        -------
        dom_id : str
            Dominant sequence metadata

        dom_seq : str
            Dominant sequence
        
        dom_df : str
            All sequences with total edit distances
        """
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if sample_size is None or sample_size > len(seq_df):
            sample_size = len(seq_df)
        if sample_size < len(seq_df):   
            seq_df = seq_df.sample(sample_size, random_state = self.random_state)
        seqs = seq_df['sequence'].values
        edit_dists = []
        for seq in seqs:
            edit_dist = 0
            for seq1 in seqs:
                edit_dist += distance(seq, seq1)
            edit_dists.append(edit_dist)
        ind_min = np.argmin(edit_dists)
        seq_df['total_edit_dist'] = edit_dists
        dom_df = seq_df.sort_values(by='total_edit_dist').reset_index()
        return dom_df

    def predict_domseq(self, seq_df, enet, n_clusters=3, sample_size=None, save_dm=None):
        """Predicts the future dominant sequence.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        enet : Qnet
            Emergenet that sequences in seq_df belong to
            
        n_clusters : int
            Number of clusters to predict dominant strain on; default 3

        sample_size : int
            Number of strains to compute dominant strain with, sampled randomly
        
        save_dm : str
            If filepath is given, save distance matrix

        Returns
        -------
        pred_accs : list
            Accession for Emergenet recommended sequences
            
        pred_names : list
            Names for Emergenet recommended sequences

        pred_seqs : list
            Emergenet recommended sequences
            
        cluster_sizes : list
            Sizes of cluster corresponding to each prediction

        """
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if sample_size is None or sample_size > len(seq_df):
            sample_size = len(seq_df)
        if sample_size < len(seq_df):   
            seq_df = seq_df.sample(sample_size, random_state = self.random_state)
        # compute qdistance matrix
        seq_arr = self._sequence_array(seq_df, sample_size)
        dist_matrix = qdistance_matrix(seq_arr, seq_arr, enet, enet)
        # convert dist_matrix to dataframe
        columns = np.arange(0, dist_matrix.shape[1])
        index = np.arange(0, dist_matrix.shape[0])
        dm = pd.DataFrame(dist_matrix, columns=columns, index=index)
        # check for null values
        for col in dm.columns:
            # delete strains that create many null values
            # ex. if dm[i][j] == NaN, if dm[i].isna().sum() is 1 (namely dm[i][j]) but dm[j].isna().sum() is 500,
            # when we delete sequence j we also take care of the NaN from sequence i, and don't delete i
            if dm[col].isna().sum() >= sample_size/100:
                dm.drop(index=col, columns=col, inplace=True)
        for col in dm.columns:
            # delete remaining strains with null values
            if dm[col].isna().sum() > 0:
                dm.drop(index=col, columns=col, inplace=True)
        if save_dm is not None:
            dm.to_csv(save_dm)
        # convert distance matrix to embedding
        embedding = MDS(n_components=2, dissimilarity="precomputed", random_state=self.random_state)
        dm_embed = embedding.fit_transform(dm)
        # cluster the distance matrix
        clustering = KMeans(n_clusters=n_clusters, random_state=self.random_state)
        clustering_predictions = clustering.fit_predict(dm_embed)
        # find unique clusters
        unique_clusters = np.unique(clustering_predictions)
        # predictions
        pred_accs = []
        pred_names = []
        pred_seqs = []
        cluster_sizes = []
        for class_ in unique_clusters:
            # separate distance matrix into submatrices
            wanted_names = dm.columns[clustering_predictions == class_]
            sub_dist_matrix = dm.loc[wanted_names, wanted_names].values
            # find centroid
            first_term = sum(sub_dist_matrix)
            H = len(sub_dist_matrix)
            A = 0.95/(np.sqrt(8) * self.seq_trunc_length**2)
            second_term = []
            for i in wanted_names:
                second_term.append(np.log(membership_degree(seq_arr[i], enet)) * H * A)
            sums = first_term - second_term
            pred_acc = seq_df.iloc[wanted_names[np.argmin(sums)]].values[0]
            pred_name = seq_df.iloc[wanted_names[np.argmin(sums)]].values[1]
            pred_seq = seq_df.iloc[wanted_names[np.argmin(sums)]].values[2]
            # save predictions
            pred_accs.append(pred_acc)
            pred_names.append(pred_name)
            pred_seqs.append(pred_seq)
            cluster_sizes.append(len(wanted_names))
        return pred_accs, pred_names, pred_seqs, cluster_sizes

def save_model(enet, outfile, low_mem=False):
    """Saves an Emergenet model.

    Parameters
    ----------
    enet : Qnet
        An Emergenet instance

    outfile : str
        File name to save to ('.joblib')

    low_mem : bool
        If True, save the Emergenet with low memory by deleting all data attributes except the tree structure

    Returns
    -------
    None
    """
    save_qnet(enet, outfile, low_mem)

    
def load_model(filepath):
    """Loads an Emergenet model.

    Parameters
    ----------
    filepath : str
        File name

    Returns
    -------
    enet : Qnet
        An Emergenet instance
    """
    enet = load_qnet(filepath)
    return enet
