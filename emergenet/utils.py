import pandas as pd
from Bio import SeqIO
from quasinet.qnet import Qnet, save_qnet, load_qnet


def parse_fasta(filepath:str) -> pd.DataFrame:
    ''' Parses a GISAID fasta file into a dataframe. Metadata should be in the format: 
    
    `Isolate name|Type|Gene name|Collection date|Protein Accession no.|Isolate ID|Lineage|Clade`

    Parameters
    ----------
    filepath: File path of fasta file
    
    Returns
    -------
    df: Dataframe with columns: `[name, subtype, segment, date, accession, sequence, HA, NA]`
    '''
    name = []
    subtype = []
    segment = []
    dates = []
    accession = []
    sequence = []
    for record in SeqIO.parse(filepath, 'fasta'):
        metadata = record.id.split('|')
        if not metadata[1].startswith('A_/_') or len(metadata[1].split('_')[2]) < 4:
            continue
        name.append(metadata[0])
        subtype.append(metadata[1].split('_')[2])
        segment.append(metadata[2])
        dates.append(metadata[3])
        accession.append(metadata[4])
        sequence.append(str(record.seq.upper()))
    df = pd.DataFrame({'name':name, 
                       'subtype':subtype,
                       'segment':segment, 
                       'date':dates,
                       'accession':accession,
                       'sequence':sequence})
    df[['HA', 'NA']] = df['subtype'].str.extract(r'H(\d+)N(\d+)')
    df['HA'] = df['HA'].apply(lambda x: 'H' + str(x))
    df['NA'] = df['NA'].apply(lambda x: 'N' + str(x))
    return df


def filter_by_date_range(df:pd.DataFrame, date_column:str, 
                         start_date:str, end_date:str) -> pd.DataFrame:
    ''' Filters a DataFrame by a date range.

    Parameters
    ----------
    df - DataFrame to filter

    date_column - Name of date column in df, entries must be in format 'YYYY-MM-DD'

    start_date - Start date in format 'YYYY-MM-DD'

    end_date - End date in format 'YYYY-MM-DD'
    
    Returns
    -------
    filtered_df - Filtered DataFrame sorted by date in descending order
    '''
    df[date_column] = pd.to_datetime(df[date_column])
    filtered_df = df[(df[date_column] >= start_date) & (df[date_column] <= end_date)]
    filtered_df.sort_values(by=[date_column], inplace=True, ascending=False)
    return filtered_df


def save_model(enet:Qnet, outfile:str, low_mem:bool=False):
    ''' Saves an Emergenet model.

    Parameters
    ----------
    enet - An Emergenet instance

    outfile - File name to save to ('.joblib')

    low_mem - If true, save the Emergenet with low memory by deleting all data attributes
    ''' 
    save_qnet(enet, outfile, low_mem)

    
def load_model(filepath:str) -> Qnet:
    ''' Loads an Emergenet model.

    Parameters
    ----------
    filepath - File name

    Returns
    -------
    enet - An Emergenet instance
    ''' 
    enet = load_qnet(filepath)
    return enet
