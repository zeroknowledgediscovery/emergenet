import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Levenshtein import distance
from quasinet.qnet import Qnet, qdistance_matrix, save_qnet, load_qnet, membership_degree
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.manifold import MDS
from collections import Counter


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

    def __init__(self, seq_trunc_length, random_state=42):
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
        if len(ids) == 3: 
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
        if len(ids) == 3: 
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
        
    @staticmethod
    def _get_date(ID):
        """Returns collection date of a sequence.

        Parameters
        ----------
        ID : str
            Sequence metadata

        Returns
        -------
        date : str
            Collection date
        """
        ids = ID.split('|')
        # NCBI format
        if len(ids) == 3: 
            date = ids[2]
            return date
        # GISAID format
        else:
            date = ids[3]
            return date
    
    @staticmethod
    def _dm_to_df(dm):
        """Converts a distance matrix to a DataFrame.

        Parameters
        ----------
        dm : numpy.ndarray
            Distance matrix as 2D array

        Returns
        -------
        dm : pd.Dataframe
            Distance matrix as DataFrame
        """
        columns = np.arange(0, dm.shape[1])
        index = np.arange(0, dm.shape[0])
        df = pd.DataFrame(dm, columns=columns, index=index)
        return df
    
    def _parse_fasta(self, filepath):
        """Parses a fasta file and returns DataFrame with four columns.
        FASTA metadata must be in the format:
            NCBI: Accession | GenBank Title | Collection Date
            GISAID: Isolate name | Type | Segment | Collection date | Protein Accession no. | Protein INSDC
        `acc` contains sequence accession.
        `name` contains sequence name.
        `date` contains collection date.
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
        dates = []
        seqs = []
        for record in SeqIO.parse(filepath, 'fasta'):
            date = self._get_date(str(record.description))
            if len(record.seq) < self.seq_trunc_length:
                continue
            if len(date) == 4:
                date += '-09-25'
            elif len(date) == 7:
                date += '-15'
            accs.append(self._get_acc(str(record.description)))
            names.append(self._get_name(str(record.description)))
            dates.append(date)
            seqs.append(''.join(record.seq[:self.seq_trunc_length].upper()))
        seq_df = pd.DataFrame({'acc':accs, 'name':names, 'date':dates, 'sequence':seqs})
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
        if sample_size is not None and sample_size < len(seq_df):
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
    
    def _compute_domseq_distance_matrix(self, seq_df):
        """Computes distance matrix in the Levenshtein distance.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        Returns
        -------
        dm : pd.DataFrame
            Distance Matrix
        """
        if 'sequence' not in seq_df.columns:
            raise ValueError('The DataFrame must store sequences in `sequence` column!')
        seqs = seq_df['sequence'].values
        n = len(seqs)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                dist_matrix[i, j] = distance(seqs[i][:self.seq_trunc_length], seqs[j][:self.seq_trunc_length])
        dist_matrix += dist_matrix.T
        dm = self._dm_to_df(dist_matrix)
        return dm
    
    def compute_domseq(self, seq_df, sample_size=None, save_data=None):
        """Computes distance matrix in the Levenshtein distance.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences
            
        sample_size : int
            Number of strains to sample randomly from seq_df
            
        save_data : string
            If directory is given, save the following files:
            1. seq_df.csv (post-sampling, if sample_size is given)
            2. dm.csv (the distance matrix corresponding to seq_df)

        Returns
        -------
        dom_seqs : pd.DataFrame
            Dominant sequences, with additional column for cluster sizes
        """
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if sample_size is not None and sample_size < len(seq_df):
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        # distance matrix
        dm = self._compute_domseq_distance_matrix(seq_df)
        if save_data is not None and os.path.isdir(save_data):
            seq_df.to_csv(save_data+'seq_df.csv', index=False)
            dm.to_csv(save_data+'dm.csv', index=False)
        # clustering
        embedding = MDS(n_components=2, dissimilarity="precomputed", random_state=self.random_state)
        dm_embed = embedding.fit_transform(dm)
        bandwidth = estimate_bandwidth(dm_embed, quantile=0.3, random_state=self.random_state)
        clustering = MeanShift(bandwidth=bandwidth)
        clustering_predictions = clustering.fit_predict(dm_embed)
        unique_clusters = np.unique(clustering_predictions)
        # dominant sequences
        dom_seqs = pd.DataFrame(columns=seq_df.columns)
        cluster_sizes = []
        for class_ in unique_clusters:
            wanted_names = dm.columns[clustering_predictions == class_]
            # find the centroid in this cluster
            cluster_seq_df = seq_df.iloc[wanted_names]
            cluster_dm = dm.iloc[wanted_names, wanted_names]
            argmin = np.argmin(cluster_dm.sum(axis='columns').values)
            # save sequence
            dom_seqs = dom_seqs.append(cluster_seq_df.iloc[argmin])
            cluster_sizes.append(len(wanted_names))
        dom_seqs['cluster_size'] = cluster_sizes
        return dom_seqs
    
    def _compute_risk_for_predict_domseq(self, seq_df, pred_seq_df, enet):
        """Computes risk scores for potential prediction sequences.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of current population sequences
            
        pred_seq_df : pd.DataFrame
            DataFrame of candidate sequences

        enet : Qnet
            Emergenet that sequences in seq_df belong to

        Returns
        -------
        pred_seq_df : pd.DataFrame
            Potential prediction sequences, with additional columns for risk score

        """
        # distance matrix for computing prediction
        seq_arr = self._sequence_array(seq_df)
        pred_seq_arr = self._sequence_array(pred_seq_df)
        dist_matrix = qdistance_matrix(pred_seq_arr, seq_arr, enet, enet) 
        # convert dist_matrix to dataframe
        dm = self._dm_to_df(dist_matrix)
        # find centroid data (sum across rows)
        first_term = dm.sum(axis='columns').values
        H = len(seq_df)
        A = 0.95/(np.sqrt(8) * self.seq_trunc_length**2)
        second_term = []
        for i in range(len(pred_seq_arr)):
            second_term.append(membership_degree(pred_seq_arr[i], enet) * H * A)
        sums = first_term - second_term
        pred_seq_df['first_term'] = first_term
        pred_seq_df['second_term'] = second_term
        pred_seq_df['sum'] = sums
        pred_seq_df = pred_seq_df.sort_values(by='sum',ascending=True)
        return pred_seq_df
    
    def predict_domseq(self, seq_df, pred_seq_df, enet, n_clusters=3, sample_size=None, save_data=None):
        """Predicts the future dominant sequence.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of current population sequences
            
        pred_seq_df : pd.DataFrame
            DataFrame of candidate sequences

        enet : Qnet
            Emergenet that sequences in seq_df belong to
            
        n_clusters : int
            Number of clusters to predict dominant strain on; default 3

        sample_size : int
            Number of strains to sample randomly from seq_df
            
        save_data : string
            If directory is given, save the following files:
            1. seq_df.csv (post-sampling, if sample_size is given)
            2. dm.csv (the distance matrix corresponding to seq_df)
            3. pred_seq_df_i.csv (pred_seq_df with additional columns for risk score, i = 1, 2, ... n_clusters)
            4. clusters.txt (clusters and cluster sizes)

        Returns
        -------
        pred_seqs : pd.DataFrame
            Emergenet recommended sequences, with additional column for cluster sizes

        """
        if len(seq_df) < 1 or len(pred_seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if sample_size is not None and sample_size < len(seq_df):
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        # distance matrix for clustering strains from current season
        seq_arr = self._sequence_array(seq_df)
        dist_matrix = qdistance_matrix(seq_arr, seq_arr, enet, enet)
        dm = self._dm_to_df(dist_matrix)
        if save_data is not None and os.path.isdir(save_data):
            seq_df.to_csv(save_data+'seq_df.csv', index=False)
            dm.to_csv(save_data+'dm.csv', index=False)
        # clustering
        embedding = MDS(n_components=2, dissimilarity="precomputed", random_state=self.random_state)
        dm_embed = embedding.fit_transform(dm)
        bandwidth = estimate_bandwidth(dm_embed, quantile=0.3, random_state=self.random_state)
        clustering = MeanShift(bandwidth=bandwidth)
        clustering_predictions = clustering.fit_predict(dm_embed)
        # sort unique clusters by size
        label_counts = Counter(clustering_predictions)
        unique_clusters = sorted(label_counts.items(), key=lambda x: x[1], reverse=True)
        if save_data is not None and os.path.isdir(save_data):
            with open(save_data + 'clusters.txt', 'w') as file:
                for item in unique_clusters:
                    file.write(str(item) + '\n')
        # predictions
        pred_seqs = pd.DataFrame(columns=pred_seq_df.columns)
        cluster_sizes = []
        for class_ in unique_clusters[:n_clusters]:
            wanted_names = dm.columns[clustering_predictions == class_[0]]
            # take riskiest strain using this cluster
            cluster_seq_df = seq_df.iloc[wanted_names]
            pred_seq_df_with_risk = self._compute_risk_for_predict_domseq(cluster_seq_df, pred_seq_df, enet)
            if save_data is not None and os.path.isdir(save_data):
                pred_seq_df_with_risk.to_csv(save_data+'pred_seq_df_' + str(class_) + '.csv', index=False)
            # save predictions
            pred_seqs = pred_seqs.append(pred_seq_df_with_risk.iloc[0])
            cluster_sizes.append(len(wanted_names))
        pred_seqs['cluster_size'] = cluster_sizes
        return pred_seqs

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
