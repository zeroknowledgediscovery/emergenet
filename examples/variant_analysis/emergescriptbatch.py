from emergenet.emergenet import Enet,  load_model, irat_risk
import numpy as np
import cProfile
import redis
import dill as pickle
import pandas as pd

VERBOSE=False
S=100

# Connect to Redis server
redis_client = redis.Redis(host='localhost', port=6379, db=0)

def load_model_with_cache(filepath):
    """
    Loads an Emergenet model from cache or disk.

    Parameters
    ----------
    filepath : str
        File name

    Returns
    -------
    enet : Qnet
        An Emergenet instance
    """

    # Check if the model is in the Redis cache
    cached_model = redis_client.get(filepath)
    
    if cached_model is not None:
        enet = pickle.loads(cached_model)
    else:
        if VERBOSE:
            print('model loading from disk')
        enet = load_model(filepath)
        
        redis_client.set(filepath, pickle.dumps(enet))
        
    return enet


def predict(target_seq,
            extensions,
            human_data,HA_ENET,NA_ENET,
            protein_length={'HA':550,'NA':449}):
    '''
    predict IRAT emergence and impact scores.
    '''

    enet_ha = Enet(seq=target_seq+extensions['HA'],
                   seq_trunc_length=protein_length['HA'])

    df_ha = enet_ha.load_data(filepath=DATA_DIR+'ha_sequences.fasta')   
    enet_ha_H = load_model_with_cache(filepath=HA_ENET)
    #save_model(enet=enet_ha_H, outfile=DATA_DIR+'ha_enet.joblib.gz')
    
    risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=df_ha,
                                                        enet=enet_ha_H,
                                                        sample_size=S,minimum=True)
    #print(risk_score_ha, variance_ha)
    
    enet_na = Enet(seq=target_seq+extensions['NA'],
                   seq_trunc_length=protein_length['NA'])
    df_na = enet_na.load_data(filepath=DATA_DIR+'na_sequences.fasta')
    
    enet_na_H = load_model_with_cache(filepath=NA_ENET)

    risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na,
                                                        enet=enet_na_H,
                                                        sample_size=S,minimum=True)
    #print(risk_score_na, variance_na)

    geom_mean_risk_score = np.sqrt((risk_score_ha+variance_ha) * (risk_score_na+variance_na))
    irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha,
                                                                  risk_score_na)


    return geom_mean_risk_score, irat_emergence_prediction, irat_impact_prediction


HA_ext='ha.fasta'
NA_ext='na.fasta'
HA_len=550
NA_len=449
DATA_DIR = '../example_data/emergenet/'
DATA_DIR = './human_sequences/'
EXTDICT={'HA':HA_ext,'NA':NA_ext}
PROTEIN_LENGTH={'HA':HA_len,'NA':NA_len}

HA_ENET={}
NA_ENET={}
HA_ENET['H1N1']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h1n1_ha.joblib.gz'
NA_ENET['H1N1']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h1n1_na.joblib.gz'
#HA_ENET['H1N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h1n2_ha.joblib.gz'
#NA_ENET['H1N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h1n2_na.joblib.gz'
#HA_ENET['H5N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h5n2_ha.joblib.gz'
#NA_ENET['H5N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h5n2_na.joblib.gz'
HA_ENET['H3N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h3n2_ha.joblib.gz'
NA_ENET['H3N2']='/home/ishanu/ZED/Research/emergenet/paper_data_v2/irat_enet/enet_models/current_enets/'+'h3n2_na.joblib.gz'

import argparse
import glob
import pandas as pd
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Specify the target set path.')
parser.add_argument('--target_seq', type=str, required=True, help='The path to the target sequence.')
args = parser.parse_args()


TGTS=glob.glob(args.target_seq+'*ha.fasta')
RESULT={}

for tgt in TGTS:
    result_thisseq={}
    tgtseq=tgt.replace('ha.fasta','')
    result_thisseq= {subtype: predict(tgtseq,
                                      extensions=EXTDICT,
                                      human_data=DATA_DIR,
                                      HA_ENET=HA_ENET[subtype],
                                      NA_ENET=NA_ENET[subtype])
                     for subtype in HA_ENET.keys()}
    result_thisseq=pd.DataFrame(result_thisseq).transpose()
    result_thisseq.columns=['geom_mean_risk_score',
                            'irat_emergence_prediction',
                            'irat_impact_prediction']
    result_thisseq=result_thisseq.sort_values('irat_emergence_prediction',ascending=False).head(1)
    RESULT[tgt] ={'irat_emergence_prediction':result_thisseq.irat_emergence_prediction.values[0],
                  'irat_impact_prediction':result_thisseq.irat_impact_prediction.values[0],
                  'geom_mean_risk_score':result_thisseq.geom_mean_risk_score.values[0],
                  'nearest_subtype':result_thisseq.index.values[0]}
    #print(RESULT[tgt])


Rf=pd.DataFrame(RESULT).transpose()

print(Rf.index.values[0].split('/')[-1].replace('_ha.fasta','').replace(':','/'),',IRAT Emergence score:',
      Rf.irat_emergence_prediction.values[0])


Rf.to_csv('variant_risk.csv')
