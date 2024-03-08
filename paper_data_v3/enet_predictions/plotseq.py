import os 
import sys
import random
import numpy as np
np.random.seed(42)
import pandas as pd
from tqdm.notebook import trange
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
plt.rcParams["figure.dpi"] = 200
from itertools import cycle
from sklearn.manifold import MDS
from sklearn.cluster import MeanShift, estimate_bandwidth
import shapely
import alphashape
from distance import hamming as distance
from quasinet.qnet import qdistance_matrix, load_qnet


TRUNC = 565 # 2 less than official length of 567
AREA_MULTIPLIER = 100


# Construct seasons
NORTH_YEARS = []
for i in np.arange(2, 23):
    YEAR = ''
    if i < 10:
        YEAR += '0' + str(i)
    else:
        YEAR += (str(i))
    if i + 1 < 10:
        YEAR += '_0' + str(i + 1)
    else:
        YEAR += '_' + str(i + 1)
    NORTH_YEARS.append(YEAR)
        
SOUTH_YEARS = []
for i in np.arange(2, 23):
    if i < 10:
        SOUTH_YEARS.append('0' + str(i))
    else:
        SOUTH_YEARS.append(str(i))

            
def _dm_to_df(dm):
    ''' Converts a distance matrix to a DataFrame.

    Parameters
    ----------
    dm : numpy.ndarray
        Distance matrix as 2D array

    Returns
    -------
    df : pd.Dataframe
        Distance matrix as DataFrame
    '''
    columns = np.arange(0, dm.shape[1])
    index = np.arange(0, dm.shape[0])
    df = pd.DataFrame(dm, columns=columns, index=index)
    return df
    
    
def _compute_distance_matrix(seqs):
    ''' Computes distance matrix in the Hamming metric.

    Parameters
    ----------
    seqs : numpy.ndarray
        Array of sequences

    Returns
    -------
    dm : pd.DataFrame
        Distance Matrix
    '''
    n = len(seqs)
    dist_matrix = np.zeros((n, n))
    for i in trange(n):
        for j in range(i+1, n):
            dist_matrix[i, j] = distance(seqs[i][:TRUNC], seqs[j][:TRUNC])
    dist_matrix += dist_matrix.T
    dm = _dm_to_df(dist_matrix)
    return dm


def _compute_qdistance_matrix(seqs, enet):
    ''' Computes distance matrix in the Qdistance metric.

    Parameters
    ----------
    seqs : numpy.ndarray
        Array of sequences

    Returns
    -------
    dm : pd.DataFrame
        Distance Matrix
    '''
    seq_arr = np.array([np.array(list(seq[:TRUNC])) for seq in seqs])
    dist_matrix = qdistance_matrix(seq_arr, seq_arr, enet, enet) 
    dm = _dm_to_df(dist_matrix)
    return dm


def plot_embedding(dm_embed, clustering_predictions, n_predictions, title, alpha):
    ''' Plots the embedding and shows Enet and WHO predictions.

    Parameters
    ----------
    dm_embed : numpy.ndarray
        Embedding points, the first point is WHO, the next n_predictions points are Enet
    
    clustering_predictions : numpy.array
        Array of clustering predictions
        
    n_predictions : int
        Number of Enet predictions to use
    
    title : str
        Title of plot
    
    alpha : int
        Alpha parameter
    '''
    # Normalize points
    mean = np.mean(dm_embed, axis=0)
    std = np.std(dm_embed, axis=0)
    dm_embed_normalized = (dm_embed - mean) / std
    if sum(std) == 0:
        dm_embed_normalized = dm_embed
        
    # Add noise to offset equal points
    noise = np.random.normal(0, 0.02, size=dm_embed_normalized.shape)
    dm_embed_normalized += noise
    
    # Find unique clusters
    unique_clusters = np.unique(clustering_predictions)
    pts = []
    for _class in unique_clusters:
        cluster_idx = np.where(clustering_predictions == _class)[0]
        dm_embed_class = dm_embed_normalized[cluster_idx]
        pts.append(dm_embed_class)
        
    # Plot clusters
    plt.figure(figsize=(10, 10))
    colors = cycle(['pink','limegreen','orange','skyblue','gold'])
    markers = cycle(['o', 'v', 's', 'D'])
    for pt in pts:
        x = pt[:,0]
        y = pt[:,1]
        points = np.column_stack((x, y))
        if len(x) < 3:
            continue
        alpha_shape = alphashape.alphashape(points, alpha=alpha)
        area = alpha_shape.area * AREA_MULTIPLIER
        if isinstance(alpha_shape, shapely.geometry.multipolygon.MultiPolygon) or isinstance(alpha_shape, shapely.GeometryCollection):
            for poly in alpha_shape.geoms:
                plt.plot(*poly.exterior.xy, 'k-', linewidth=1)
        else:
            plt.plot(*alpha_shape.exterior.xy, 'k-', linewidth=1)
        plt.scatter(x, y, color = next(colors), alpha = 0.9, marker = next(markers),
                    label = f'Size: {len(pt)}, Area: {round(area, 3)}')
        
    # Enet, WHO
    colors = cycle(['blue'])
    markers = cycle(['o', 'v', 's', 'D'])
    for i in range(2, n_predictions + 1):
        plt.scatter(dm_embed_normalized[i][0], dm_embed_normalized[i][1], color = next(colors), 
                    marker = next(markers), label = f'Enet {i-2}', s = 150)
    plt.scatter(dm_embed_normalized[1][0], dm_embed_normalized[1][1], 
                color = 'purple', marker = next(markers), label = 'Enet Single', s = 150)
    plt.scatter(dm_embed_normalized[0][0], dm_embed_normalized[0][1], 
                color = 'red', marker = next(markers), label = 'WHO', s = 150)
        
    # Axis limits
    plt.axis('auto')
    plt.title(title)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1));
    
    
def plot_predictions(hemisphere, subtype, year, metric,
                     sample_size=None, previous_season=False, alpha=4):
    ''' Plots the embedding and shows Enet and WHO predictions.

    Parameters
    ----------
    hemisphere : str
        'north' or 'south'
    
    subtype : str
        'hin1' or 'h3n2'
        
    year : int
        Index corresponding to sequence in NORTH_YEARS or SOUTH_YEARS, [1, 21]
    
    metric : str
        'hamming' or 'qdistance'
    
    sample_size : int
        Number of points to sample from the strain population to plot
    
    previous_season : bool
        If true, plot previous season rather than current season strain population
    
    alpha : int
        Alpha parameter, default 4
        
    Returns
    -------
    pred : pd.DataFrame
        Dataframe of WHO prediction and WHO and Enet errors, multi-cluster
        
    pred_single : pd.DataFrame
        Dataframe of WHO prediction and WHO and Enet errors, single-cluster
    '''
    if hemisphere == 'north':
        YEARS = NORTH_YEARS
    elif hemisphere == 'south':
        YEARS = SOUTH_YEARS
    SEASON = YEARS[year - previous_season]
    DIR = hemisphere + '_' + subtype
    NAME = DIR + '_' + SEASON 

    # Load sequences and predictions
    seq_df = pd.read_csv('data/merged/' + DIR + '/' + NAME + '.csv')
    seqs = list(seq_df['sequence'].values)
    pred = pd.read_csv('results/enet_who_comparison/' + DIR + '.csv', converters={'season': str})
    pred = pred[pred['season'] == YEARS[year]]
    pred_single = pd.read_csv('results/enet_who_comparison/' + DIR + '_single_cluster.csv', converters={'season': str})
    pred_single = pred_single[pred_single['season'] == YEARS[year]]
    
    # WHO seq
    who_seq = pred['ha_seq_who'].values[0]
    
    # Enet seqs
    n_predictions = 3
    enet_seqs = [pred_single['ha_seq'].values[0], pred['ha_seq_0'].values[0], pred['ha_seq_1'].values[0]]

    # In embedding, WHO seq at index 0, Enet seq at index 1 - n_predictions
    if sample_size is not None and sample_size < len(seqs):
        seqs = random.sample(seqs, sample_size)
    seqs = np.array([who_seq] + enet_seqs + seqs)

    # DM embedding
    if metric == 'hamming':
        dm = _compute_distance_matrix(seqs)
    elif metric == 'qdistance':
        enet = load_qnet('enet_models/' + DIR + '/' + NAME + '.joblib')
        dm = _compute_qdistance_matrix(seqs, enet)
    embedding = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    dm_embed = embedding.fit_transform(dm)

    # Clustering
    bandwidth = max(estimate_bandwidth(dm_embed, quantile=0.3, random_state=42), 1e-5)
    clustering = MeanShift(bandwidth=bandwidth)
    clustering_predictions = clustering.fit_predict(dm_embed)

    # Show errors
    if not previous_season and metric == 'hamming':
        print('Multi-cluster')
        print(f'\tWHO error: {pred["ha_who_error"].values[0]:.3f}')
        print(f'\tEnet error: {pred["ha_enet_error"].values[0]:.3f}')
        print(f'\tError difference: {pred["ha_who_error"].values[0] - pred["ha_enet_error"].values[0]:.3f}')
        print('Single-cluster')
        print(f'\tWHO error: {pred_single["ha_who_error"].values[0]:.3f}')
        print(f'\tEnet error: {pred_single["ha_enet_error"].values[0]:.3f}')
        print(f'\tError difference: {pred_single["ha_who_error"].values[0] - pred_single["ha_enet_error"].values[0]:.3f}')
    
    # Plot
    title = f'Prediction: {hemisphere.upper()} {subtype.upper()} {YEARS[year]}\nMetric: {metric.upper()}\nBackground Population: {SEASON}'
    plot_embedding(dm_embed, clustering_predictions, n_predictions, title, alpha)    
    return pred, pred_single
