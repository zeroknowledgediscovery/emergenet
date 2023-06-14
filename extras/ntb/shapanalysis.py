#!/usr/bin/python
import numpy as np
import pandas as pd
from quasinet.qnet import qdistance, load_qnet
from emergenet.domseq import DomSeq, load_model
import os
import argparse
import glob
import shap
import multiprocessing as mulpro
import itertools
from scipy import stats as st
from functools import partial


def fpar(s, s0__, qnet__):
    return qdistance(s0__, s, qnet__, qnet__)


def getSHAP_theta(path_to_raw_seq,
                  hemisphere_subtype_protein_tag,
                  year,
                  TRUNCATED=565):

    RAWPATH=glob.glob(path_to_raw_seq+'/'+hemisphere_subtype_protein_tag+\
               '/'+hemisphere_subtype_protein_tag+'_*'+str(year)+'.fasta')
    RAW_PATH=RAWPATH[0]
    
    Q_PATH='../../paper_data/enet_predictions/enet_models/'+\
        hemisphere_subtype_protein_tag+'/'+hemisphere_subtype_protein_tag+'*'+str(year)+'.joblib.gz'
    Q_PATH=glob.glob(Q_PATH)[0]
   

    qnet__=load_qnet(Q_PATH,gz=True)
    domseq = DomSeq(seq_trunc_length=TRUNCATED, random_state=42)
    df = domseq.load_data(filepath=RAW_PATH)

    seq=df.set_index('id').sequence.values
    S=[np.array(list(x))[:TRUNCATED] for x in seq]

    S_consensus=st.mode(S)[0][0]
    s0__=S_consensus


    def f(s_array):
        pool = mulpro.Pool(processes=10)
        fpar_partial = partial(fpar, s0__=s0__, qnet__=qnet__)
        return np.array(pool.map(fpar_partial, s_array))

    explainer = shap.KernelExplainer(f, np.array([S_consensus]))

    S_=np.array([x for x in S if np.random.rand()<.6 ])
    shap_values = explainer.shap_values(S_, nsamples=TRUNCATED)

    shpc=pd.DataFrame(pd.DataFrame(shap_values).abs().mean().sort_values(ascending=False),columns=['shpc'])
    shp.index.name='features'
    shpc.to_csv(TAG+'_C.csv')


    nulldata=np.array([np.array(['']*TRUNCATED)])
    explainerN = shap.KernelExplainer(f,nulldata)
    shap_valuesN = explainerN.shap_values(S_, nsamples=TRUNCATED)

    shp_N=pd.DataFrame(pd.DataFrame(shap_valuesN).abs().mean().sort_values(ascending=False),columns=['shp'])
    shp_N.index.name='features'
    shp_N.to_csv(TAG+'_N.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get SHAP theta')

    parser.add_argument('--path_to_raw_seq', type=str, default='../../paper_data/enet_predictions/raw_data/gisaid/',
                        help='Path to raw sequences')
    parser.add_argument('--hemisphere_subtype_protein_tag', type=str, default='south_h3n2_ha',
                        help='Hemisphere subtype protein tag')
    parser.add_argument('--year', type=int, default='21',
                        help='Year')
    parser.add_argument('--TRUNCATED', type=int, default=565,
                        help='Truncated value')

    args = parser.parse_args()

    getSHAP_theta(args.path_to_raw_seq, args.hemisphere_subtype_protein_tag, args.year, args.TRUNCATED)




#TRUNCATED=565
#YR='21'
#TAG_='south_h3n2_ha'
#RAWPATHPREF='../../paper_data/enet_predictions/raw_data/gisaid/'

