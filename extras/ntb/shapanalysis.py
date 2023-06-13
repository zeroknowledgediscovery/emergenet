import numpy as np
import pandas as pd
from quasinet.qnet import qdistance, load_qnet
from emergenet.domseq import DomSeq, save_model, load_model
import os
import shutil

TRUNCATED=565
TAG='south_h3n2_ha_21'
RAW_PATH='../../paper_data/enet_predictions/raw_data/gisaid/south_h3n2_ha/'
Q_PATH='../../paper_data/enet_predictions/enet_models/south_h3n2_ha/'+TAG+'.joblib.gz'
qnet__=load_qnet(Q_PATH,gz=True)

# initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNCATED, random_state=42)
df_north = domseq.load_data(filepath=RAW_PATH+TAG+'.fasta')

seq=df_north.set_index('name').sequence.values
S=[np.array(list(x))[:TRUNCATED] for x in seq]

from scipy import stats as st
S_consensus=st.mode(S)[0][0]

import multiprocessing as mulpro
import itertools

s0__=S_consensus

def fpar(s):
    return qdistance(s0__,s,qnet__,qnet__)

def f(s_array):
    pool = mulpro.Pool(processes=10)
    return np.array(pool.map(fpar, s_array))

import shap
explainer = shap.KernelExplainer(f, np.array([S_consensus]))

S_=np.array([x for x in S if np.random.rand()<.6 ])
shap_values = explainer.shap_values(S_, nsamples=TRUNCATED)

#import pylab as plt
#plt.style.use('dark_background')
#shap.summary_plot(shap_values,S_, plot_type="bar",max_display=566)

shpc=pd.DataFrame(pd.DataFrame(shap_values).abs().mean().sort_values(ascending=False),columns=['shpc'])
shp.index.name='features'
shpc.to_csv('shpc.csv')


nulldata=np.array([np.array(['']*TRUNCATED)])
explainerN = shap.KernelExplainer(f,nulldata)
shap_valuesN = explainerN.shap_values(S_, nsamples=TRUNCATED)

shp=pd.DataFrame(pd.DataFrame(shap_valuesN).abs().mean().sort_values(ascending=False),columns=['shp'])
shp.index.name='features'
shp.to_csv('shp.csv')
