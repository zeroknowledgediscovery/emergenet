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
PARALLEL=True

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

Q_PATH=args.path_to_enets+args.sequence_name+'.joblib.gz'
qnet__=load_qnet(Q_PATH,gz=True)
RAW_PATH=glob.glob(args.path_to_raw_seq+args.sequence_name+'*fasta')[0]

domseq = DomSeq(seq_trunc_length=args.TRUNCATED, random_state=42)
r0 = domseq.load_data(filepath=RAW_PATH)
r0=r0.set_index('name').sequence.values
r0_=np.array(list(r0[0]))
S_consensus=r0_
s0__=S_consensus

RAW_PATH_GISAID=glob.glob(args.path_to_gisaid_seq+args.sequence_name+'*fasta')[0]
df = domseq.load_data(filepath=RAW_PATH_GISAID)
seq=df.set_index('name').sequence.values
S=[np.array(list(x)) for x in seq]
S_=np.array([x for x in S if np.random.rand()<SAMPLE_FRACTION ])

if PARALLEL:
    def fpar(s):
        return qdistance(s0__,s,qnet__,qnet__)

    def f(s_array):
        pool = mulpro.Pool(processes=28)
        return np.array(pool.map(fpar, s_array))

else:
    def f(s_array):
        return np.fromiter([qdistance(s0__,
                                  s,qnet__,qnet__)
                        for s in s_array],dtype=np.float,count=len(s_array))



explainer = shap.KernelExplainer(f, np.array([s0__]))
shap_values = explainer.shap_values(S_, nsamples=args.TRUNCATED)

shpc=pd.DataFrame(pd.DataFrame(shap_values).abs().mean().sort_values(ascending=False),columns=['shpc'])
shpc.index.name='features'
shpc.to_csv(args.sequence_name+'_SHAP.csv')
