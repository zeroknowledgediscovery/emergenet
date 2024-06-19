import re, os, json, joblib
import numpy as np
import pandas as pd
from typing import Tuple
from datetime import date, datetime
from collections import Counter
from pkg_resources import resource_filename
from quasinet.qnet import Qnet, qdistance
from .utils import filter_by_date_range, load_model, save_model


# Lengths to truncate sequences
HA_TRUNC = 560
NA_TRUNC = 460


class Enet(object):
    ''' Emergenet architecture for predicting emergence risk.
    ''' 

    def __init__(self, analysis_date:str, ha_seq:str, na_seq:str, 
                 pretrained_enet_path:str=None, save_data:str=None, random_state:int=None):
        ''' Initializes an Emergenet instance.

        Parameters
        ----------
        analysis_date - `PRESENT` or date of analysis (YYYY-MM-DD), supported from 2010-01-01 to 2024-01-01

        ha_seq - The target's HA sequence to be analysed by Emergenet
            
        na_seq - The target's NA sequence to be analysed by Emergenet

        save_data - Directory to save data to (Enet models, sequences used for training, etc.)

        random_state - Sets seed for random number generator
        
        pretrained_enet_path - Base path for all pretrained Enet models
        '''
        # Date
        if analysis_date != 'PRESENT':
            if not re.match(r'^[0-9]{4}-[0-9]{2}-[0-9]{2}$', analysis_date):
                raise ValueError('Date must be in format YYYY-MM-DD! or "PRESENT"')
            if analysis_date > '2024-01-01' or analysis_date < '2010-01-01':
                raise ValueError('Emergenet only supports sequences from 2010-01-01 to 2024-01-01!')
        self.analysis_date = analysis_date
        # Sequences
        if HA_TRUNC > len(ha_seq):
            print('HA sequence is shorter than required length (560), padding the end with Xs')
            ha_seq = ha_seq.ljust(HA_TRUNC, 'X')
        self.ha_seq = ha_seq.upper()[:HA_TRUNC]
        if NA_TRUNC > len(na_seq):
            print('NA sequence is shorter than required length (460), padding the end with Xs')
            na_seq = na_seq.ljust(NA_TRUNC, 'X')
        self.na_seq = na_seq.upper()[:NA_TRUNC]
        # Pretrained Enet path
        self.pretrained_enet_path = pretrained_enet_path
        # Save data
        self.save_data = save_data
        if save_data is not None:
            if not os.path.exists(save_data):
                os.makedirs(save_data, exist_ok=True)
            if self.analysis_date != 'PRESENT' and self.pretrained_enet_path is None:
                self.save_model = os.path.join(save_data, 'enet_models')
                self.save_sequences = os.path.join(save_data, 'data')
                self.save_results = os.path.join(save_data, 'results')
                os.makedirs(self.save_model, exist_ok=True)
                os.makedirs(self.save_sequences, exist_ok=True)
                os.makedirs(self.save_results, exist_ok=True)
        # Random state
        if random_state is not None and random_state < 0:
            raise ValueError('Seed must be between 0 and 2**32 - 1!')
        self.random_state = random_state

        
    def __repr__(self):
        return 'emergenet.Emergenet'


    def __str__(self):
        return self.__repr__()
    

    def _load_sequences(self, yearsbefore:int) -> pd.DataFrame:
        ''' Loads human sequences within yearsbefore years of the analysis date.

        Parameters
        ----------
        yearsbefore - Number of years prior to analysis_date to consider 

        Returns
        -------
        filtered - DataFrame of human sequences within one year of the analysis date
        '''
        filepath = resource_filename('emergenet', 'data/human.csv')
        human = pd.read_csv(filepath, na_filter=False)
        end = datetime.strptime(self.analysis_date, '%Y-%m-%d').date()
        start = date(end.year - yearsbefore, end.month, end.day)
        filtered = filter_by_date_range(human, 'date', str(start), str(end))
        print('Number of human sequences:', len(filtered))
        print(Counter(filtered[filtered['segment'] == 'HA']['HA']))
        print(Counter(filtered[filtered['segment'] == 'NA']['NA']))
        return filtered
    

    def _sequence_array(self, segment:str, seq_df:pd.DataFrame, 
                        sample_size:int=None, include_target:bool=True) -> np.ndarray:
        ''' Extracts array of sequence arrays from DataFrame.

        Parameters
        ----------
        segment - Either 'HA' or 'NA'
            
        seq_df - DataFrame containing sequences

        sample_size - Number of sequences to sample randomly
            
        include_target - If true, includes target sequence

        Returns
        -------
        seq_lst - Array of sequence arrays
        ''' 
        if 'sequence' not in seq_df.columns:
            raise ValueError('The DataFrame must store sequences in `sequence` column!')
        if sample_size is not None:
            sample_size = min(sample_size, len(seq_df))
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        seqs = seq_df['sequence'].values
        seq_lst = []
        if segment == 'HA':
            TRUNC = HA_TRUNC
            if include_target:
                seq_lst.append(np.array(list(self.ha_seq[:TRUNC])))
        elif segment == 'NA':
            TRUNC = NA_TRUNC
            if include_target:
                seq_lst.append(np.array(list(self.na_seq[:TRUNC])))
        for seq in seqs:
            if len(seq) < TRUNC:
                continue
            seq_lst.append(np.array(list(seq[:TRUNC])))
        seq_lst = np.array(seq_lst)
        return seq_lst

    
    def _compute_risks(self, segment:str, seq_df:pd.DataFrame, enet:Qnet) -> pd.DataFrame:
        ''' Computes risk score with qdistance.

        Parameters
        ----------
        segment - Either 'HA' or 'NA'

        seq_df - DataFrame of sequences

        enet - Emergenet that sequences in `seq_df` belong to

        Returns
        -------
        seq_df - The input `seq_df` with extra risk column
        ''' 
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        seq_arr = self._sequence_array(segment, seq_df, include_target=False)
        if segment == 'HA':
            TRUNC = HA_TRUNC
            target_seq = np.array(list(self.ha_seq[:TRUNC]))
        elif segment == 'NA':
            TRUNC = NA_TRUNC
            target_seq = np.array(list(self.na_seq[:TRUNC]))
        risks = []
        for i in range(len(seq_arr)):
            qdist = qdistance(target_seq, seq_arr[i], enet, enet)
            if np.isnan(qdist):
                risks.append(-1)
                continue
            risks.append(qdist)
        seq_df['risk'] = risks
        return seq_df
    

    def train(self, segment:str, seq_df:pd.DataFrame, sample_size:int=None, 
              include_target:bool=True, n_jobs:int=1) -> Qnet:
        ''' Trains an Emergenet model.

        Parameters
        ----------
        segment - Either 'HA' or 'NA'

        seq_df - DataFrame of sequences

        sample_size - Number of sequences to train Emergenet on, sampled randomly

        include_target - If true, includes target sequence

        n_jobs - Number of CPUs to use when training

        Returns
        -------
        enet - Trained Emergenet
        ''' 
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if segment not in ['HA', 'NA']:
            raise ValueError('Segment must be either HA or NA!')
        if segment == 'HA':
            TRUNC = HA_TRUNC
        elif segment == 'NA':
            TRUNC = NA_TRUNC
        seq_arr = self._sequence_array(segment, seq_df, sample_size, include_target)
        enet = Qnet(feature_names=['x' + str(i) for i in np.arange(TRUNC)],
                    random_state=self.random_state, n_jobs=n_jobs)
        enet.fit(seq_arr)
        return enet
    

    def risk(self, yearsbefore:int=1, enet_sample_size:int=None, risk_sample_size:int=None) -> Tuple[float, float]:
        ''' Computes risk scores for the target sequence.
        If `save_data` is not None, `analysis_date` is not 'PRESENT' and pretrained_enet_path is None, saves the following:
        1. Emergenet models: `save_data/enet_models/<subtype>.joblib`
        2. DataFrames of sequences used for training: `save_data/data/<subtype>.csv`
        3. Results (sequence DataFrames with extra `risk_score` column): `save_data/results/<subtype>.csv`
        4. Minimum risks for each HA and NA subtype: `save_data/results/<segment>_min_risks.csv`

        If `save_data` is not None, and `analysis_date` is 'PRESENT' or pretrained_enet_path is not None, saves the following:
        1. Results (sequence DataFrames with extra `risk_score` column): `save_data/<subtype>.csv`
        2. Minimum risks for each HA and NA subtype: `save_data/<segment>_min_risks.csv`

        Parameters
        ----------
        yearsbefore - Number of years prior to analysis_date to consider 
        
        enet_sample_size - Number of sequences of each subtype to train Emergenet on, sampled randomly
        
        risk_sample_size - Number of sequences of each subtype to sample for risk estimation

        Returns
        -------
        ha_risk - Risk score for HA segment
        
        na_risk - Risk score for NA segment
        '''
        if self.analysis_date == 'PRESENT':
            filepath = resource_filename('emergenet', 'data/current_subtypes.json')
            with open(filepath, 'r') as file:
                current_subtypes = json.load(file)
            for segment in ['HA', 'NA']:
                risks = pd.DataFrame()
                if segment == 'HA':
                    TRUNC = HA_TRUNC
                elif segment == 'NA':
                    TRUNC = NA_TRUNC
                for subtype in current_subtypes[segment]:
                    # Load human sequences for and pretrained Enet models for current subtype
                    human_filepath = resource_filename('emergenet', f'data/current/{subtype}.csv')
                    model_filepath = resource_filename('emergenet', f'models/{subtype}.joblib')
                    df = pd.read_csv(human_filepath, na_filter=False)
                    # Sample from human sequences if needed
                    if risk_sample_size is not None:
                        sample_size = min(risk_sample_size, len(df) // 2)
                        df = df.sample(n=sample_size, replace=False, random_state=self.random_state)
                    # Load Enet and compute risks
                    enet = load_model(model_filepath)
                    df = self._compute_risks(segment, df, enet)
                    if self.save_data is not None:
                        df.to_csv(os.path.join(self.save_data, subtype + '.csv'), index=False)
                    # Save minimum risk for current subtype
                    risks[subtype] = [np.min(df['risk'])]
                # Save overall minimum risk
                if self.save_data is not None:
                    risks.to_csv(os.path.join(self.save_data, segment + '_min_risks.csv'), index=False)
                if segment == 'HA':
                    ha_risk = risks.min(axis=1).values[0] + 1e-5
                elif segment == 'NA':
                    na_risk = risks.min(axis=1).values[0] + 1e-5
            return ha_risk, na_risk

        elif self.pretrained_enet_path is not None:
            human = self._load_sequences(yearsbefore)
            for segment in ['HA', 'NA']:
                risks = pd.DataFrame()
                if segment == 'HA':
                    TRUNC = HA_TRUNC
                elif segment == 'NA':
                    TRUNC = NA_TRUNC
                human1 = human[human['segment'] == segment]
                human1 = human1[human1['sequence'].str.len() >= TRUNC]
                subtypes = Counter(human1[segment])
                for subtype in subtypes:
                    # Skip subtypes with less than 15 sequences
                    if subtypes[subtype] < 15:
                        continue                    
                    # Use only unique sequences for inference
                    df = human1[human1[segment] == subtype].drop_duplicates(subset=['sequence'])
                    # Sample from human sequences if needed
                    if risk_sample_size is not None:
                        sample_size = min(risk_sample_size, len(df) // 2)
                        df = df.sample(n=sample_size, replace=False, random_state=self.random_state)
                    # Load Enet and compute risks
                    enet = load_model(os.path.join(self.pretrained_enet_path, subtype + '.joblib'))
                    df = self._compute_risks(segment, df, enet)
                    if self.save_data is not None:
                        df.to_csv(os.path.join(self.save_data, subtype + '.csv'), index=False)
                    # Save minimum risk for current subtype
                    risks[subtype] = [np.min(df['risk'])]
                # Save overall minimum risk
                if self.save_data is not None:
                    risks.to_csv(os.path.join(self.save_data, segment + '_min_risks.csv'), index=False)
                if segment == 'HA':
                    ha_risk = risks.min(axis=1).values[0] + 1e-5
                elif segment == 'NA':
                    na_risk = risks.min(axis=1).values[0] + 1e-5
            return ha_risk, na_risk
        
        else:
            human = self._load_sequences(yearsbefore)
            for segment in ['HA', 'NA']:
                risks = pd.DataFrame()
                if segment == 'HA':
                    TRUNC = HA_TRUNC
                elif segment == 'NA':
                    TRUNC = NA_TRUNC
                human1 = human[human['segment'] == segment]
                human1 = human1[human1['sequence'].str.len() >= TRUNC]
                subtypes = Counter(human1[segment])
                for subtype in subtypes:
                    # Skip subtypes with less than 15 sequences
                    if subtypes[subtype] < 15:
                        continue
                    # Use entire population for constructing Enet
                    df = human1[human1[segment] == subtype]
                    if self.save_data is not None:
                        df.to_csv(os.path.join(self.save_sequences, subtype + '.csv'), index=False)
                    # Train Enet
                    enet = self.train(segment, df, sample_size=enet_sample_size)
                    if self.save_data is not None:
                        save_model(enet, os.path.join(self.save_model, subtype + '.joblib'))
                    # Use only unique sequences for inference
                    df = df.drop_duplicates(subset=['sequence'])
                    # Sample from human sequences if needed
                    if risk_sample_size is not None:
                        sample_size = min(risk_sample_size, len(df) // 2)
                        df = df.sample(n=sample_size, replace=False, random_state=self.random_state)
                    # Compute risks
                    df = self._compute_risks(segment, df, enet)
                    if self.save_data is not None:
                        df.to_csv(os.path.join(self.save_results, subtype + '.csv'), index=False)
                    # Save minimum risk for current subtype
                    risks[subtype] = [np.min(df['risk'])]
                # Save overall minimum risk
                if self.save_data is not None:
                    risks.to_csv(os.path.join(self.save_results, segment + '_min_risks.csv'), index=False)
                if segment == 'HA':
                    ha_risk = risks.min(axis=1).values[0] + 1e-5
                elif segment == 'NA':
                    na_risk = risks.min(axis=1).values[0] + 1e-5
            return ha_risk, na_risk
        

def predict_irat_emergence(ha_risk:float, na_risk:float) -> Tuple[float, float, float]:
    ''' Computes IRAT emergence risk score.

    Parameters
    ----------
    ha_risk - Risk score for HA segment

    na_risk - Risk score for NA segment

    Returns
    -------
    irat_emergence - Predicted IRAT emergence risk score

    irat_emergence_low - Lower bound IRAT emergence risk score

    irat_emergence_high - Upper bound IRAT emergence risk score
    '''
    geom_mean = np.sqrt(ha_risk * na_risk)
    x = -np.log(geom_mean + 5e-4)
    emergence_model = joblib.load(resource_filename('emergenet', 'models/emergence_model.joblib'))
    emergence_low_model = joblib.load(resource_filename('emergenet', 'models/emergence_low_model.joblib'))
    emergence_high_model = joblib.load(resource_filename('emergenet', 'models/emergence_high_model.joblib'))
    irat_emergence = emergence_model['intercept'] + emergence_model['slope'] * x
    irat_emergence_low = emergence_low_model['intercept'] + emergence_low_model['slope'] * x
    irat_emergence_high = emergence_high_model['intercept'] + emergence_high_model['slope'] * x
    return irat_emergence, irat_emergence_low, irat_emergence_high
        