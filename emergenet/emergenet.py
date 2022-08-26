import numpy as np
import pandas as pd
from Bio import SeqIO
from quasinet.qnet import Qnet, qdistance, save_qnet, load_qnet


class Emergenet(object):
    """Emergenet"""

    def __init__(self, seq, seq_length, random_state=None):
        self.seq = seq
        self.seq_length = seq_length
        self.random_state = random_state

    def __repr__(self):
        return "emergenet.Emergenet"

    def __str__(self):
        return self.__repr__()

    @staticmethod
    def _is_fasta(filepath):
        with open(filepath, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)

    def _parse_fasta(self, filepath):
        if not self._is_fasta(filepath):
            raise ValueError('The infile must be in `.fasta` format')
        acc = []
        seq = []
        for record in SeqIO.parse(filepath, 'fasta'):
            if len(record.seq) < self.seq_length:
                continue
            acc.append(record.id.split('|')[0])
            seq.append(np.array(record.seq[:self.seq_length].upper()))
        seq_df = pd.DataFrame({'name': acc, 'sequence': seq})
        return seq_df

    def _sequence_array(self, seq_df, sample_size=None, strain=None):
        if sample_size is None or sample_size > len(seq_df):
            sample_size = len(seq_df)
        seqs = seq_df['sequence'].sample(sample_size, random_state=self.random_state).values
        seq_lst = []
        for seq in seqs:
            seq_lst.append(seq)
        if strain is not None:
            seq_lst.append(np.array(list(strain)))
        return np.array(seq_lst)

    def load_data(self, filepath):
        seq_df = self._parse_fasta(filepath)
        return seq_df

    def train(self, seq_df, sample_size=None):
        seq_arr = self._sequence_array(seq_df, sample_size, self.seq)
        qnet = Qnet(feature_names=['x' + str(i) for i in np.arange(self.seq_length)], n_jobs=1)
        qnet.fit(seq_arr)
        return qnet

    @staticmethod
    def save_models(qnet, filepath):
        save_qnet(qnet, filepath)

    @staticmethod
    def load_models(filepath):
        qnet = load_qnet(filepath)
        return qnet

    def emergence_risk(self, seq_df, qnet, sample_size=None):
        seq_arr = self._sequence_array(seq_df, sample_size)
        qdist_sum = 0
        for i in range(len(seq_arr)):
            qdist_sum += qdistance(self.seq, seq_arr[i], qnet, qnet)
        return qdist_sum/len(seq_arr)
