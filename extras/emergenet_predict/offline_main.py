from emergenet.emergenet import Enet, save_model, load_model, irat_risk
import numpy as np
import os
import json


import argparse
from emergenet.emergenet import Enet, save_model, load_model, irat_risk
import numpy as np
import os
import json

def emergenet_predict(HAfasta,NAfasta,
                      DATA_DIR='./',
                      protein_length={'HA':550,
                                      'NA':449}):

    enet_ha = Enet(seq=HAfasta,
                   seq_trunc_length=protein_length['HA'], random_state=22)
    df_ha = enet_ha.load_data(filepath=DATA_DIR+'ha_sequences.fasta')
    enet_ha_H = load_model(filepath=DATA_DIR+'ha_enet.joblib.gz')
    risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=df_ha,
                                                        enet=enet_ha_H,
                                                        sample_size=1000)

    enet_na = Enet(seq=NAfasta,
                   seq_trunc_length=protein_length['NA'], random_state=22)
    df_na = enet_na.load_data(filepath=DATA_DIR+'na_sequences.fasta')
    enet_na_H = load_model(filepath=DATA_DIR+'na_enet.joblib.gz')
    risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na,
                                                        enet=enet_na_H,
                                                        sample_size=1000)

    geom_mean_risk_score = np.sqrt(risk_score_ha * risk_score_na)
    irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha,
                                                                  risk_score_na)

    data={'IRAT_emergence': irat_emergence_prediction,
                    'IRAT_impact': irat_impact_prediction}
#    with open(FILE+'.json', 'w') as file:
#        json.dump(data, file, indent=4)
        
    
    return data


parser = argparse.ArgumentParser(description='Emergenet Prediction Script')

# Add the arguments with flags
parser.add_argument('-ha', '--HAfasta', type=str, required=True, help='Path to HA fasta file')
parser.add_argument('-na', '--NAfasta', type=str, required=True, help='Path to NA fasta file')
#parser.add_argument('-d', '--DATA_DIR', type=str, default='./', help='Directory path for data')
parser.add_argument('-o', '--OUTFILE', type=str, default='result', help='Output file path')

# Execute the parse_args() method
args = parser.parse_args()

#print(args.DATA_DIR)

results = emergenet_predict(args.HAfasta, args.NAfasta)
#print(results)



#curl -X POST https://us-central1-pkcsaas-01.cloudfunctions.net/emergenet_predict -H "Content-Type: application/json" -d '{"HA":"MKAILLVLLHTFAATSADTICVGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGKAPLYLGKCNIAGWLLGNPECELPLTVSSWSYIVETSDSDNGTCYPGDFTNYEELREQLSSVSSFERFEMFPKESSWPNHETNKSVTAACPYAGASSFYRNLIWLVKKDDSYPMLNISYVNNKGKEVLVLWGIHHPPTEDDQKWLYKNADAYVFVGTSTYSQKFEPEIATRPRVRDQTGRMNYYWTLVKPGDKITFEATGNLVVPRYAFAMNRGSESGIIISDAPVHDCNTICQTPKGALNTSLPFQNVHPVTIGECPKYIKSTRLKMATGLRNTPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADQKSTQRAVDGITNKVNSIIERMNSQFTAVGKEFSNLERRIENLNKKVDDGFLDVWTYNAELLILLENERTLDFHDSNVKNLYERVRNQLRNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEESKLNREEIDGVKLDSTKVYQILAIYSTVASSLVVLVSLGALSFWMCSNGSL","NA":"MNPNQKIITIGSVSLIIATICFLMQIAILVTTVTLHFKQHDCNSSSNNQVMLCEPIIIERNKTEIVYLTNTTVEKEICPKPAEYRNWSKPQCNIAGFAPFSKDNSIRLSAGGDIWVTREPYVSCDLDKCYQFALGQGTTLNNRHSNDTVHDRTPYRTLLMNELGVPFHLGTRQVCVAWSSSSCHDGKAWLHVCITGDDNNATASFIYNGRLVDSVVSWSKNILRTQESECVCINGTCTVVMTDGSASGKADTRILFIVEGKIIHVSKLSGSAQHVEECSCYPRYPGVRCVCRDNWKGSNRPIVDINMKDYSIVSSYVCSGLVGDTPRKTDSLSSSNCLDPNNEEGDHGVKGWAFDDGDDVWMGRTINETLRLGYETFKVIKGWSKPNSKLQTNRQVIVKGGNRSGYSGIFSVEGKNCINRCFYVELIRGRREETRVWWTSNSIVVFCGTSGTYGTGSWPDGADINLMPI"}'
