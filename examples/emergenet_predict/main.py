from emergenet.emergenet import Enet, save_model, load_model, irat_risk
import numpy as np
from flask import jsonify, abort



def emergenet_predict(request):
    '''
    predict IRAT emergence and impact scores.
    POST request with {'HA': ha_sequence, 'NA': na_sequence}
    '''

    if request.method != "POST":
        abort(405, "Only POST requests are accepted")

    request_json = request.get_json(silent=True)
    if (
        not request_json
        or "HA" not in request_json
        or "NA" not in request_json
    ):
        return 'Missing "HA" or "NA" in the request', 400

    protein_length={'HA':550,'NA':449}
    DATA_DIR='./'

    HAseq=request_json["HA"]
    NAseq=request_json["NA"]

    enet_ha = Enet(seq=HAseq,
                   seq_trunc_length=protein_length['HA'], random_state=22)
    df_ha = enet_ha.load_data(filepath=DATA_DIR+'ha_sequences.fasta')
    enet_ha_H = load_model(filepath=DATA_DIR+'ha_enet.joblib')
    risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=df_ha,
                                                        enet=enet_ha_H,
                                                        sample_size=1000)

    enet_na = Enet(seq=NAseq,
                   seq_trunc_length=protein_length['NA'], random_state=22)
    df_na = enet_na.load_data(filepath=DATA_DIR+'na_sequences.fasta')
    enet_na_H = load_model(filepath=DATA_DIR+'na_enet.joblib.gz')
    risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na,
                                                        enet=enet_na_H,
                                                        sample_size=1000)

    geom_mean_risk_score = np.sqrt(risk_score_ha * risk_score_na)
    irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha,
                                                                  risk_score_na)

    return jsonify({'IRAT_emergence': irat_emergence_prediction,
                    'IRAT_impact': irat_impact_prediction}), 201


