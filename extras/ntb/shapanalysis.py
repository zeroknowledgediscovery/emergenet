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

PARALLEL=True

parser = argparse.ArgumentParser(description='Get SHAP theta')

parser.add_argument('--path_to_raw_seq', type=str,
                    default='../../paper_data/enet_predictions/raw_data/gisaid/',
                    help='Path to raw sequences')
parser.add_argument('--hemisphere_subtype_protein_tag',
                    type=str, default='south_h3n2_ha',
                    help='Hemisphere subtype protein tag')
parser.add_argument('--year', type=int, default='21',
                    help='Year')
parser.add_argument('--path_to_enets', type=str,
                    default='../../paper_data/enet_predictions/enet_models/',
                    help='Path to enets')
parser.add_argument('--TRUNCATED', type=int, default=565,
                    help='Truncated value')
parser.add_argument('--NUM_CORES', type=int, default=28,
                    help='number of cores')
parser.add_argument('--SAMPLE_FRACTION', type=float, default=0.05,
                    help='fractio of samples used to evaluate shap')
args = parser.parse_args()

print(args.path_to_enets+args.hemisphere_subtype_protein_tag+'/'+args.hemisphere_subtype_protein_tag+'*'+str(args.year)+'.joblib.gz')

Q_PATH=glob.glob(args.path_to_enets+args.hemisphere_subtype_protein_tag+'/'+args.hemisphere_subtype_protein_tag+'*'+str(args.year)+'.joblib.gz')[0]
qnet__=load_qnet(Q_PATH,gz=True)
RAW_PATH=glob.glob(args.path_to_raw_seq+'/'+args.hemisphere_subtype_protein_tag+\
               '/'+args.hemisphere_subtype_protein_tag+'_*'+\
                   str(args.year)+'.fasta')[0]

domseq = DomSeq(seq_trunc_length=args.TRUNCATED, random_state=42)
df = domseq.load_data(filepath=RAW_PATH)
seq=df.set_index('name').sequence.values
S=[np.array(list(x)) for x in seq]
S_consensus=st.mode(S)[0][0]
S_=np.array([x for x in S if np.random.rand()<args.SAMPLE_FRACTION ])

if PARALLEL:
    def fpar(s):
        return qdistance(S_consensus,s,qnet__,qnet__)

    def f(s_array):
        pool = mulpro.Pool(processes=args.NUM_CORES)
        return np.array(pool.map(fpar, s_array))

else:
    def f(s_array):
        return np.fromiter([qdistance(S_consensus,
                                  s,qnet__,qnet__)
                        for s in s_array],dtype=np.float,count=len(s_array))


explainer = shap.KernelExplainer(f, np.array([S_consensus]))
shap_values = explainer.shap_values(S_, nsamples=args.TRUNCATED)

shpc=pd.DataFrame(pd.DataFrame(shap_values).abs().mean().sort_values(ascending=False),columns=['shpc'])
shpc.index.name='features'
shpc.to_csv(args.hemisphere_subtype_protein_tag+'_SHAP.csv')
