import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from quasinet.qnet import Qnet, qdistance, save_qnet, load_qnet


class Enet(object):
    ''' Emergenet architecture.

    Parameters
    ----------
    seq : str
        The target sequence to be analysed by Emergenet
        Either nucleotide/amino-acid sequence or fasta file path (containing '.fasta')

    seq_trunc_length : int
        Length to truncate sequences in Emergenet analysis
        (Sequences used to train Emergenet and compute E-distance must be of same length)

    seq_metadata : str
        Describes the sequence; added automatically if 'seq' is a fasta file path

    random_state : int
        Sets seed for random number generator
    ''' 

    def __init__(self, seq, seq_trunc_length, seq_metadata=None, random_state=None):
        if seq.endswith('.fasta'):
            if self._count_seqs(seq) != 1:
                raise ValueError('The file must contain exactly 1 sequence!')
            for record in SeqIO.parse(seq, 'fasta'):
                self.seq = str(record.seq.upper())
                self.seq_metadata = str(record.description)
        else:
            self.seq = seq.upper()
            self.seq_metadata = seq_metadata

        if seq_trunc_length > len(self.seq):
            raise ValueError('Length to truncate sequences must not be greater than target sequence length!')
        self.seq_trunc_length = seq_trunc_length

        if random_state is not None and random_state < 0:
            raise ValueError('Seed must be between 0 and 2**32 - 1!')
        self.random_state = random_state

        
    def __repr__(self):
        return "emergenet.Emergenet"

    
    def __str__(self):
        return self.__repr__()

    
    @staticmethod
    def _count_seqs(filepath):
        ''' Returns number of sequences in a fasta file.

        Parameters
        ----------
        filepath : str
            File name

        Returns
        -------
        num_sequences : int
            Number of sequences
        ''' 
        with open(filepath, 'r') as f:
            fasta = SeqIO.parse(f, 'fasta')
            if not any(fasta):
                raise ValueError('The infile must be in fasta format!')
        with open(filepath, 'r') as f:
            lines = f.read()
        num_sequences = lines.count('>')
        return num_sequences
    

    def _sequence_array(self, seq_df, sample_size=None, include_target=True):
        ''' Extracts array of sequence arrays from DataFrame.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame containing sequences

        sample_size : int
            Number of strains to sample randomly
            
        include_target : bool
            If true, includes target sequence

        Returns
        -------
        seq_lst: numpy.ndarray
            Array of sequence arrays
        ''' 
        if 'sequence' not in seq_df.columns:
            raise ValueError('The DataFrame must store sequences in `sequence` column!')
        if sample_size is not None:
            sample_size = min(sample_size, len(seq_df))
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        seqs = seq_df['sequence'].values
        seq_lst = []
        if include_target:
            seq_lst.append(np.array(list(self.seq[:self.seq_trunc_length])))
        for seq in seqs:
            seq_lst.append(np.array(list(seq[:self.seq_trunc_length])))
        seq_lst = np.array(seq_lst)
        return seq_lst

    
    def train(self, seq_df, sample_size=None, n_jobs=1, include_target=True):
        ''' Trains an Emergenet model.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        sample_size : int
            Number of strains to train Emergenet on, sampled randomly

        n_jobs : int
            Number of CPUs to use when training
            
        include_target : bool
            If true, includes target sequence

        Returns
        -------
        enet : Qnet
            Trained Emergenet
        ''' 
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        seq_arr = self._sequence_array(seq_df, sample_size, include_target)
        enet = Qnet(feature_names=['x' + str(i) for i in np.arange(self.seq_trunc_length)],
                    random_state=self.random_state, n_jobs=n_jobs)
        enet.fit(seq_arr)
        return enet

    
    def risk(self, seq_df, enet):
        ''' Computes risk score with qdistance.

        Parameters
        ----------
        seq_df : pd.DataFrame
            DataFrame of sequences

        enet : Qnet
            Emergenet that sequences in seq_df belong to

        Returns
        -------
        seq_df : float
            Original seq_df with extra risk column
        ''' 
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        seq_arr = self._sequence_array(seq_df, include_target=False)
        target_seq = np.array(list(self.seq[:self.seq_trunc_length]))
        risks = []
        for i in range(len(seq_arr)):
            qdist = qdistance(target_seq, seq_arr[i], enet, enet)
            if np.isnan(qdist):
                risks.append(-1)
                continue
            risks.append(qdist)
        seq_df['risk'] = risks
        return seq_df


def save_model(enet, outfile, low_mem=False):
    ''' Saves an Emergenet model.

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
    ''' 
    save_qnet(enet, outfile, low_mem)

    
def load_model(filepath):
    ''' Loads an Emergenet model.

    Parameters
    ----------
    filepath : str
        File name

    Returns
    -------
    enet : Qnet
        An Emergenet instance
    ''' 
    enet = load_qnet(filepath)
    return enet
