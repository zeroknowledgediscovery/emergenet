#!/usr/bin/python3.10
import numpy as np
import pandas as pd
from quasinet.qnet import qdistance, load_qnet
from emergenet.domseq import DomSeq, load_model
import os
import argparse
import glob
import shap
import multiprocessing as mulpro
from scipy import stats as st
from concurrent.futures import ProcessPoolExecutor
import itertools

SAMPLE_FRACTION=0.1

parser = argparse.ArgumentParser(description='Get SHAP theta')

parser.add_argument('--sequence_name', type=str, default='A:American wigeon:South Carolina:AH0195145:2021_ha',
                    help='Path to raw sequences')
parser.add_argument('--path_to_raw_seq', type=str, default='../../paper_data/irat_enet/raw_data/irat_sequences/',
                    help='Path to raw sequences')
parser.add_argument('--path_to_gisaid_seq', type=str, default='../../paper_data/irat_enet/raw_data/gisaid/',
                    help='Path to GISAID sequences')
parser.add_argument('--path_to_enets', type=str, default='../../paper_data/irat_enet/enet_models/irat_enets/',
                    help='Path to enets')
parser.add_argument('--TRUNCATED', type=int, default=550,
                    help='Truncated value')
parser.add_argument('--NUM_CORES', type=int, default=28,
                    help='number of cores')

args = parser.parse_args()

Q_PATH=path_to_enets+sequence_name+'.joblib.gz'
qnet__=load_qnet(Q_PATH,gz=True)
RAW_PATH=glob.glob(path_to_raw_seq+sequence_name+'*fasta')[0]

r0 = DomSeq(seq_trunc_length=TRUNCATED, random_state=42).load_data(filepath=RAW_PATH)
r0=r0.set_index('name').sequence.values
r0_=np.array(list(r0[0]))

RAW_PATH_GISAID=glob.glob(path_to_gisaid_seq+sequence_name+'*fasta')[0]

df = DomSeq(seq_trunc_length=TRUNCATED, random_state=42).load_data(filepath=RAW_PATH_GISAID)

seq=df.set_index('name').sequence.values
S=[np.array(list(x)) for x in seq]

S_consensus=r0_
s0__=S_consensus


def calculate_shap_values(S_chunk):
    explainer = shap.KernelExplainer(f, np.array([S_consensus]))
    return explainer.shap_values(S_chunk, nsamples=TRUNCATED)

def f(s_array):
    # Split s_array into chunks
    chunks = np.array_split(s_array, NUM_CORES)

    # Use a ProcessPoolExecutor to compute SHAP values in parallel
    with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
        results = executor.map(calculate_shap_values, chunks)

    # Concatenate the results
    shap_values = np.concatenate(list(results))

    return shap_values


S_=np.array([x for x in S if np.random.rand()<SAMPLE_FRACTION ])
shap_values = f(S_)

shpc=pd.DataFrame(pd.DataFrame(shap_values).abs().mean().sort_values(ascending=False),columns=['shpc'])
shpc.index.name='features'
shpc.to_csv(sequence_name+'_SHAP.csv')



