import numpy as np
import pandas as pd
from Bio import SeqIO
from Levenshtein import distance
from quasinet.qnet import Qnet, qdistance, qdistance_matrix, save_qnet, load_qnet, membership_degree


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
    
    
def _parse_fasta(filepath, seq_trunc_length):
    """Parses a fasta file and returns DataFrame with two columns.
    `id` contains fasta metadata.
    `sequence` contains sequences truncated to `seq_trunc_length`.

    Parameters
    ----------
    filepath : str
        File name
        
    seq_trunc_length : int
        Length to truncate sequences
        
    Returns
    -------
    seq_df : pd.DataFrame
        DataFrame of sequences
    """
    if _count_seqs(filepath) == 0:
        raise ValueError('The file contains no sequences!')
    if seq_trunc_length <= 0:
        raise ValueError('Length to truncate sequences must be positive!')
    ids = []
    seqs = []
    for record in SeqIO.parse(filepath, 'fasta'):
        if len(record.seq) < seq_trunc_length:
            continue
        ids.append(str(record.id))
        seqs.append(''.join(record.seq[:seq_trunc_length].upper()))
    seq_df = pd.DataFrame({'id': ids, 'sequence': seqs})
    return seq_df


def _sequence_array(seq_df, seq_trunc_length, sample_size=None, random_state=42):
    """Extracts array of sequence arrays from DataFrame; includes target sequence.

    Parameters
    ----------
    seq_df : pd.DataFrame
        DataFrame containing sequences
        
    seq_trunc_length : int
        Length to truncate sequences

    sample_size : int
        Number of strains to sample randomly
        
    random_state : int
        Sets seed for random number generator

    Returns
    -------
    seq_lst: numpy.ndarray
        Array of sequence arrays
    """
    if 'sequence' not in seq_df.columns:
        raise ValueError('The DataFrame must store sequences in `sequence` column!')
    if seq_trunc_length <= 0:
        raise ValueError('Length to truncate sequences must be positive!')
    if sample_size is None or sample_size > len(seq_df):
        sample_size = len(seq_df)
    if sample_size < len(seq_df):
        seq_df = seq_df.sample(sample_size, random_state=random_state)
    seqs = seq_df['sequence'].apply(lambda x: np.array(list(x[:seq_trunc_length])))
    seq_lst = []
    for seq in seqs:
        seq_lst.append(seq)
    seq_lst = np.array(seq_lst)
    return seq_lst


def load_data(filepath, seq_trunc_length, outfile=None):
    """Loads fasta file data and optionally saves to CSV.

    Parameters
    ----------
    filepath : str
        File name
        
    seq_trunc_length : int
        Length to truncate sequences

    outfile : str
        File name to save to ('.csv')

    Returns
    -------
    seq_df : pd.DataFrame
        DataFrame of sequences
    """
    seq_df = _parse_fasta(filepath, seq_trunc_length)
    if outfile is not None:
        if not outfile.endswith('.csv'):
            raise ValueError('The outfile must end with `.csv`!')
        seq_df.to_csv(outfile, index=False)
    return seq_df


def train(seq_df, seq_trunc_length, sample_size=None, random_state=42, n_jobs=1):
    """Trains an Emergenet model.

    Parameters
    ----------
    seq_df : pd.DataFrame
        DataFrame of sequences
        
    seq_trunc_length : int
        Length to truncate sequences (number of features)

    sample_size : int
        Number of strains to train Emergenet on, sampled randomly
        
    random_state : int
        Sets seed for random number generator

    n_jobs : int
        Number of CPUs to use when training

    Returns
    -------
    enet : Qnet
        Trained Emergenet
    """
    if len(seq_df) < 1:
        raise ValueError('The DataFrame contains no sequences!')
    seq_arr = _sequence_array(seq_df, seq_trunc_length, sample_size, random_state)
    enet = Qnet(feature_names=['x' + str(i) for i in np.arange(seq_trunc_length)],
                random_state=random_state, n_jobs=n_jobs)
    enet.fit(seq_arr)
    return enet


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


def compute_domseq(seq_df, sample_size=None, random_state=42):
    """Computes the dominant sequence in the edit (Levenshtein) distance metric.

    Parameters
    ----------
    seq_df : pd.DataFrame
        DataFrame of sequences

    sample_size : int
        Number of strains to compute dominant strain with, sampled randomly

    random_state : int
        Sets seed for random number generator

    Returns
    -------
    dom_id : str
        Dominant sequence metadata

    dom_seq : str
        Dominant sequence
    """
    if len(seq_df) < 1:
        raise ValueError('The DataFrame contains no sequences!')
    if random_state < 0:
        raise ValueError('Seed must be between 0 and 2**32 - 1!')
    if sample_size is None or sample_size > len(seq_df):
        sample_size = len(seq_df)
    if sample_size < len(seq_df):   
        seq_df = seq_df.sample(min(sample_size, len(seq_df)), random_state = random_state)
    seqs = seq_df['sequence'].values
    edit_dists = []
    for seq in seqs:
        edit_dist = 0
        for seq1 in seqs:
            edit_dist += distance(seq, seq1)
        edit_dists.append(edit_dist)
    ind_min = np.argmin(edit_dists)
    dom_id = seq_df.iloc[ind_min].values[0]
    dom_seq = seq_df.iloc[ind_min].values[1]
    return dom_id, dom_seq


def predict_domseq(seq_df, enet, seq_trunc_length, sample_size=None, random_state=42):
    """Computes the dominant sequence in the edit (Levenshtein) distance metric.

    Parameters
    ----------
    seq_df : pd.DataFrame
        DataFrame of sequences
        
    enet : Qnet
        Emergenet that sequences in seq_df belong to
        
    seq_trunc_length : int
        Length to truncate sequences

    sample_size : int
        Number of strains to compute dominant strain with, sampled randomly

    random_state : int
        Sets seed for random number generator

    Returns
    -------
    pred_id : str
        Emergenet recommended sequence metadata

    pred_seq : str
        Emergenet recommended sequence
    """
    if len(seq_df) < 1:
        raise ValueError('The DataFrame contains no sequences!')
    if random_state < 0:
        raise ValueError('Seed must be between 0 and 2**32 - 1!')
    if sample_size is None or sample_size > len(seq_df):
        sample_size = len(seq_df)
    H = min(sample_size, len(seq_df))
    if sample_size < len(seq_df):   
        seq_df = seq_df.sample(H, random_state = random_state)
    seq_arr = _sequence_array(seq_df, seq_trunc_length, sample_size, random_state)
    dist_matrix = qdistance_matrix(seq_arr, seq_arr, enet, enet)
    first_term = sum(dist_matrix)
    A = 0.95/(np.sqrt(8) * seq_trunc_length**2)
    second_term = np.ones(len(seq_arr)) * H * A
    for i in range(len(seq_arr)):
        second_term[i] *= membership_degree(seq_arr[i], enet)
    sums = first_term - second_term
    pred_id = seq_df.iloc[np.argmin(sums)].values[0]
    pred_seq = seq_df.iloc[np.argmin(sums)].values[1]
    return pred_id, pred_seq
