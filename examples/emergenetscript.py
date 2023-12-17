from emergenet.emergenet import Enet, save_model, load_model, irat_risk
import numpy as np

def predict(target_seq,
            extensions,
            human_data,
            protein_length={'HA':550,'NA':449}):
    '''
    predict IRAT emergence and impact scores.
    '''

    enet_ha = Enet(seq=target_seq+extensions['HA'],
                   seq_trunc_length=protein_length['HA'], random_state=22)
    df_ha = enet_ha.load_data(filepath=DATA_DIR+'ha_sequences.fasta')
    enet_ha_H = load_model(filepath=DATA_DIR+'ha_enet.joblib.gz')
    save_model(enet=enet_ha_H, outfile=DATA_DIR+'ha_enet.joblib.gz')

    risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=df_ha,
                                                        enet=enet_ha_H,
                                                        sample_size=1000)

    enet_na = Enet(seq=target_seq+extensions['NA'],
                   seq_trunc_length=protein_length['NA'], random_state=22)
    df_na = enet_na.load_data(filepath=DATA_DIR+'na_sequences.fasta')
    enet_na_H = load_model(filepath=DATA_DIR+'na_enet.joblib.gz')
    risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na,
                                                        enet=enet_na_H,
                                                        sample_size=1000)

    geom_mean_risk_score = np.sqrt(risk_score_ha * risk_score_na)
    irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha,
                                                                  risk_score_na)
    return geom_mean_risk_score, irat_emergence_prediction, irat_impact_prediction




target_seq='../extras/variants/data/variants/A:Alberta:01:2020_(H1N2)v_'
HA_ext='ha.fasta'
NA_ext='na.fasta'
HA_len=550
NA_len=449
DATA_DIR = 'example_data/emergenet/'
EXTDICT={'HA':HA_ext,'NA':NA_ext}
PROTEIN_LENGTH={'HA':HA_len,'NA':NA_len}

geom_mean_risk_score, irat_emergence_prediction, irat_impact_prediction = predict(target_seq,
                                                                                  extensions=EXTDICT,
                                                                                  human_data=DATA_DIR)


print('IRAT emergence:', irat_emergence_prediction)
print('IRAT impact:', irat_impact_prediction)
print('geometric risk:', geom_mean_risk_score)

#avg_ha, min_ha, max_ha, var_ha = enet_ha.emergence_risk_qsampling(seq_df=df_ha, enet=enet_ha1, sample_size=1000, qsamples=10, steps=10)

#avg_na, min_na, max_na, var_na = enet_na.emergence_risk_qsampling(seq_df=df_na, enet=enet_na1, sample_size=1000, qsamples=10, steps=10)


