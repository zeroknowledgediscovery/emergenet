#!/usr/bin/python

import argparse, time
from emergenet.emergenet import Enet, predict_irat_emergence


#HA=MKTIIAFSCILCLIFAQKLPGSDNSMATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGGICNSPHQILDGKNCTLIDALLGDPHCDDFQNKEWDLFVERSTAYSNCYPYYVPDYATLRSLVASSGNLEFTQESFNWTGVAQGGSSYACRRGSVNSFFSRLNWLYNLNYKYPEQNVTMPNNDKFDKLYIWGVHHPGTDKDQTNLYVQASGRVIVSTKRSQQTVIPNIGSRPWVRGVSSIISIYWTIVKPGDILLINSTGNLIAPRGYFKIQSGKSSIMRSDAHIDECNSECITPNGSIPNDKPFQNVNKITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAINQITGKLNRVIKKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAEILVALENQHTIDLTDSEMSKLFERTRRQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDIYRNEALNNRFQIKGVQLKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKGNIRCNICI

#NA=MNPNQKIITIGSVSLIIATICFLMQIAILVTTVTLHFKQHDYNSPPNNQAMLCEPTIIERNTTEIVYLTNITIEKEICPKLAEYRNWSKPQCNITGFAPFSKDNSIRLSAGGDIWVTREPYVSCDPDKCYQFALGQGTTLNNGHSNNTVHDRTPYRTLLMNELGVPFHLGTRQVCMAWSSSSCHDGKAWLHVCITGNDNNATASFIYNGRLVDSIGSWSKNILRTQESECVCINGTCTVVMTDGSASGKADTKILFVEEGKIVHISTLSGSAQHVEECSCYPRFPGVRCVCRDNWKGSNRPIVDINVKNYSIVSSYVCSGLVGDTPRKSDSVSSSYCLDPNNEKGGHGVKGWAFDDGNDVWMGRTINETLRLGYETFKVIEGWSKANSKLQTNRQVIVEKGDRSGYSGIFSVEGKSCINRCFYVELIRGRKEETKVWWTSNSIVVFCGTSGTYGTGSWPDGADINLMPI


def main():
    parser = argparse.ArgumentParser(description="Estimate IRAT Emergence Score using Enet.")
    parser.add_argument('ha_seq', type=str, help='HA sequence for the analysis')
    parser.add_argument('na_seq', type=str, help='NA sequence for the analysis')
    parser.add_argument('--analysis_date', '-a', type=str, default='PRESENT', 
                        help='Analysis date. If not "PRESENT", it will take longer as new Enets must be trained.')
    parser.add_argument('--risk_sample_size', '-r', type=int, default=None, 
                        help='Risk sample size. Takes ~30 seconds when analysis_date is PRESENT and risk_sample_size is 100.')
    parser.add_argument('--save_data_dir', '-s', type=str, default=None,
                        help='Directory to save data.')
    parser.add_argument('--random_state', '-d', type=int, default=None,
                        help='Random state for reproducibility.')
    args = parser.parse_args()
    if args.analysis_date != 'PRESENT':
        print('Training Enet models ...')
    if args.save_data_dir is not None:
        print(f'Saving data to {args.save_data_dir}')
    #start_time = time.time()

    # Initialize the Enet
    enet = Enet(analysis_date=args.analysis_date, 
                ha_seq=args.ha_seq, 
                na_seq=args.na_seq, 
                save_data=args.save_data_dir,
                random_state=args.random_state)

    # Estimate the Enet risk scores
    ha_risk, na_risk = enet.risk(risk_sample_size=args.risk_sample_size)

    # Map the Enet risk scores to the IRAT risk scale
    irat, irat_low, irat_high = predict_irat_emergence(ha_risk=ha_risk, 
                                                       na_risk=na_risk)
    
    print(f'Emergence Score: {irat:.2f}')
    print(f'Emergence Score high: {irat_high:.2f}')
    print(f'Emergence Score low: {irat_low:.2f}')
    #end_time = time.time()
    #elapsed_time = end_time - start_time
    #print(f'Time taken: {elapsed_time:.2f} seconds')

if __name__ == '__main__':
    main()
