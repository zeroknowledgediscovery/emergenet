import numpy as np
import pandas as pd
from Bio import SeqIO
from Levenshtein import distance as levdistance
from scipy.spatial.distance import pdist, squareform
from quasinet.qnet import Qnet, qdistance, qdistance_matrix, save_qnet, load_qnet


def parse_fasta(filepath, seq_trunc_length, seq_array=False):
    """Parses a fasta file and returns DataFrame with two columns.
    `id` contains fasta metadata.
    `sequence` contains sequences truncated to `seq_trunc_length` as character array.

    Parameters
    ----------
    filepath : str
        File name

    seq_trunc_length : int
        Length to truncate sequences in Qnet analysis
        (Sequences used to train Qnet and compute q-distance must be of same length)
    
    seq_array : bool
        If true, store the sequences as an array of characters
        (Required for quasinet functions)

    Returns
    -------
    seq_df : pd.DataFrame
        DataFrame of sequences
    """
    with open(filepath, 'r') as f:
        fasta = SeqIO.parse(f, 'fasta')
        if not any(fasta):
            raise ValueError('The infile must be in fasta format!')
    ids = []
    seqs = []
    for record in SeqIO.parse(filepath, 'fasta'):
        if len(record.seq) < seq_trunc_length:
            continue
        ids.append(str(record.id))
        if seq_array: 
             seqs.append(np.array(record.seq[:self.seq_trunc_length].upper()))
        else:
            seqs.append(''.join(record.seq[:seq_trunc_length].upper()))
    seq_df = pd.DataFrame({'id': ids, 'sequence': seqs})
    return seq_df
    
    
def dominant_seq(seq_df, metric):
    """Finds the dominant sequence, computed as the centroid in the edit distance metric.

    Parameters
    ----------
    seq_df : pd.DataFrame
        DataFrame containing sequences
        
    metric : str
        Either `qdistance` or `levenshtein`
    
    Returns
    -------
    dom_id : str
        Metadata of dominant sequence
    
    dom_seq : str
        Dominant sequence
    """
    if len(seq_df) < 1:
        raise ValueError('The DataFrame contains no sequences!')
    if 'sequence' not in seq_df.columns:
        raise ValueError('The DataFrame must store sequences in `sequence` column!')
    if 'id' not in seq_df.columns:
        raise ValueError('The DataFrame must metadata in `id` column!')
    
    seqs = np.array(seq_df['sequence'].values).reshape(-1,1)
    if metric == 'qdistance':
        dist_matrix = squareform(pdist(seqs, lambda x,y: qdistance(x[0],y[0])))
    elif metric == 'levenshtein':
        dist_matrix = squareform(pdist(seqs, lambda x,y: levdistance(x[0],y[0])))
    else:
        raise ValueError('The distance function must be either `qdistance` or `levenshtein`')
    ind_min = np.argmin(sum(dist_matrix))
    dom_id = seq_df.iloc[ind_min]['id']
    dom_seq = seq_df.iloc[ind_min]['sequence']
    return dom_id, dom_seq
