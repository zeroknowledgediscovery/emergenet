import os
import numpy as np
import pandas as pd
from distance import hamming
from itertools import cycle
import shapely
import alphashape
from quasinet.qnet import Qnet, qdistance_matrix, membership_degree
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.manifold import MDS
from collections import Counter
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')


# Area multiplier for alpha shapes
AREA_MULTIPLIER = 100


class DomSeq(object):
    ''' DomSeq architecture for predicting dominant sequences.
    ''' 

    def __init__(self, seq_trunc_length:int, random_state:int=42):
        ''' Initializes an DomSeq instance.

        Parameters
        ----------
        seq_trunc_length - Length to truncate sequences in Emergenet analysis
            (Sequences used to train Emergenet and compute E-distance must be of same length)

        random_state - Sets seed for random number generator
        '''
        if seq_trunc_length <= 0:
            raise ValueError('Length to truncate sequences must be positive!')
        self.seq_trunc_length = seq_trunc_length

        if random_state is not None and random_state < 0:
            raise ValueError('Seed must be between 0 and 2**32 - 1!')
        self.random_state = random_state

        
    def __repr__(self):
        return 'emergenet.DomSeq'

    
    def __str__(self):
        return self.__repr__()
    
    
    @staticmethod
    def _dm_to_df(dm:np.ndarray) -> pd.DataFrame:
        ''' Converts a distance matrix to a DataFrame.

        Parameters
        ----------
        dm - Distance matrix as 2D array

        Returns
        -------
        df - Distance matrix as DataFrame
        ''' 
        columns = np.arange(0, dm.shape[1])
        index = np.arange(0, dm.shape[0])
        df = pd.DataFrame(dm, columns=columns, index=index)
        return df

    
    def _sequence_array(self, seq_df:pd.DataFrame, sample_size:int=None) -> np.ndarray:
        ''' Extracts array of sequence arrays from DataFrame; includes target sequence.

        Parameters
        ----------
        seq_df - DataFrame containing sequences

        sample_size - Number of strains to sample randomly

        Returns
        -------
        seq_lst - Array of sequence arrays
        ''' 
        if 'sequence' not in seq_df.columns:
            raise ValueError('The DataFrame must store sequences in `sequence` column!')
        if sample_size is not None and sample_size < len(seq_df):
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        seqs = seq_df['sequence'].apply(lambda x: np.array(list(x[:self.seq_trunc_length])))
        seq_lst = []
        for seq in seqs:
            seq_lst.append(seq)
        seq_lst = np.array(seq_lst)
        return seq_lst

    
    def train(self, seq_df:pd.DataFrame, sample_size:int=None, n_jobs:int=1) -> Qnet:
        ''' Trains an Emergenet model.

        Parameters
        ----------
        seq_df - DataFrame of sequences

        sample_size - Number of strains to train Emergenet on, sampled randomly

        n_jobs - Number of CPUs to use when training

        Returns
        -------
        enet - Trained Emergenet
        ''' 
        if len(seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        seq_arr = self._sequence_array(seq_df, sample_size)
        enet = Qnet(feature_names=['x' + str(i) for i in np.arange(self.seq_trunc_length)],
                    random_state=self.random_state, n_jobs=n_jobs)
        enet.fit(seq_arr)
        return enet
    
    
    def _compute_risk_for_predict_domseq(self, seq_df:pd.DataFrame, 
                                         pred_seq_df:pd.DataFrame, enet:Qnet) -> pd.DataFrame:
        ''' Computes risk scores for potential prediction sequences.

        Parameters
        ----------
        seq_df - DataFrame of current population sequences
            
        pred_seq_df - DataFrame of candidate sequences

        enet - Emergenet that sequences in seq_df belong to

        Returns
        -------
        pred_seq_df - Potential prediction sequences, with additional columns for risk score
        ''' 
        # Distance matrix for computing prediction
        seq_arr = self._sequence_array(seq_df)
        pred_seq_arr = self._sequence_array(pred_seq_df)
        dist_matrix = qdistance_matrix(pred_seq_arr, seq_arr, enet, enet) 
        
        # Convert dist_matrix to dataframe
        dm = self._dm_to_df(dist_matrix)
        assert len(dm) == len(pred_seq_arr)
        
        # Find centroid data (sum across rows)
        first_term = dm.sum(axis='columns').values
        H = len(seq_df)
        A = 0.95/(np.sqrt(8) * self.seq_trunc_length**2)
        second_term = np.array([membership_degree(seq, enet) * H * A for seq in pred_seq_arr])
        sums = first_term - second_term
        pred_seq_df['first_term'] = first_term
        pred_seq_df['second_term'] = second_term
        pred_seq_df['sum'] = sums
        pred_seq_df = pred_seq_df.sort_values(by='sum',ascending=True)
        return pred_seq_df
    
    @staticmethod
    def _plot_embedding(dm_embed:np.ndarray, unique_clusters:Counter, clustering_predictions:np.array, 
                        n_clusters:int, save_data:str=None, alpha:int=4) -> list:
        ''' Plots the embedding and shows Enet and WHO predictions.

        Parameters
        ----------
        dm_embed - Embedding points
            
        unique_clusters - Counter (class, size) of unique clusters, sorted by size

        clustering_predictions - Array of clustering predictions

        n_clusters - Number of cluster alphashapes to draw

        save_data - Directory to save plot to
        
        alpha - Alpha parameter, default 4
            
        Returns
        -------
        cluster_areas - Areas of largest n_clusters
        '''
        # Add noise to offset equal points
        noise = np.random.normal(0, 0.02, size=dm_embed.shape)
        dm_embed += noise

        # Find unique clusters
        pts = []
        for class_ in unique_clusters:
            cluster_idx = np.where(clustering_predictions == class_[0])[0]
            pt = dm_embed[cluster_idx]
            x = pt[:,0]
            y = pt[:,1]
            points = np.column_stack((x, y))
            pts.append(points)

        # Plot clusters
        plt.figure(figsize=(10, 10))
        colors = cycle(['skyblue','red','orange','green', 'purple'])
        markers = cycle(['o', 'v'])
        
        # Plot at most 5 clusters
        cluster_areas = []
        for i in range(min(len(pts), 5)):
            pt = pts[i]
            x = pt[:,0]
            y = pt[:,1]
            points = np.column_stack((x, y))
            
            # Draw alphashape for the first n clusters if >= 3 points
            # Otherwise just scatter the points
            if i < n_clusters:
                if len(x) < 3:
                    cluster_areas.append(0)
                    plt.scatter(x, y, color = next(colors), alpha = 0.7, 
                                marker = next(markers), label = f'Size: {len(pt)}')
                else:
                    alpha_shape = alphashape.alphashape(points, alpha=alpha)
                    area = alpha_shape.area * AREA_MULTIPLIER
                    cluster_areas.append(area)
                    if isinstance(alpha_shape, shapely.geometry.multipolygon.MultiPolygon) or \
                       isinstance(alpha_shape, shapely.GeometryCollection):
                        for poly in alpha_shape.geoms:
                            plt.plot(*poly.exterior.xy, 'k-', linewidth=1)
                    else:
                        plt.plot(*alpha_shape.exterior.xy, 'k-', linewidth=1)
                    plt.scatter(x, y, color = next(colors), alpha = 0.7, marker = next(markers),
                                label = f'Size: {len(pt)}, Area: {round(area, 3)}')
            else:
                if len(x) == 1:
                    continue
                plt.scatter(x, y, color = next(colors), alpha = 0.7, 
                            marker = next(markers), label = f'Size: {len(pt)}')
        # Save figure
        plt.axis('auto')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1));
        if save_data is not None and os.path.isdir(save_data): 
            plt.savefig(os.path.join(save_data, 'embedding.png'), dpi=300, 
                        bbox_inches = 'tight', pad_inches = 0, transparent=True)
        return cluster_areas
    
    
    def predict_domseq(self, seq_df:pd.DataFrame, pred_seq_df:pd.DataFrame, enet:Qnet, 
                       n_clusters:int=3, sample_size:int=None, save_data:str=None) -> pd.DataFrame:
        ''' Predicts the future dominant sequence with multiple clusters.

        Parameters
        ----------
        seq_df - DataFrame of current population sequences
            
        pred_seq_df - DataFrame of candidate sequences

        enet - Emergenet that sequences in seq_df belong to
            
        n_clusters - Number of clusters to predict dominant strain on; default 3

        sample_size - Number of strains to sample randomly from seq_df
            
        save_data - If directory is given, save the following files:
            1. seq_df.csv (post-sampling, if sample_size is given)
            2. dm.csv (the distance matrix corresponding to seq_df)
            3. pred_seq_df_i.csv (pred_seq_df with columns for risk score, i = 1,2,...n_clusters)
            4. clusters.txt (clusters and cluster counts)
            5. embedding.png (embedding plot with alpha shapes drawn for n_clusters)

        Returns
        -------
        pred_seqs - Emergenet recommended sequences, with additional two columns for cluster counts and areas
        ''' 
        if len(seq_df) < 1 or len(pred_seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
        if sample_size is not None and sample_size < len(seq_df):
            seq_df = seq_df.sample(sample_size, random_state=self.random_state)
        if n_clusters <= 1:
            raise ValueError('For single cluster, please use predict_single_domseq!')
        if len(pred_seq_df) > 20000:
            pred_seq_df = pred_seq_df.sample(20000, random_state=self.random_state)
        
        # Parameter for alphashape
        ALPHA = 4
            
        # Distance matrix for clustering strains from current season
        seq_arr = self._sequence_array(seq_df)
        dist_matrix = qdistance_matrix(seq_arr, seq_arr, enet, enet)
        dm = self._dm_to_df(dist_matrix)
        if save_data is not None and os.path.isdir(save_data):
            seq_df.to_csv(os.path.join(save_data, 'seq_df.csv'), index=False)
            dm.to_csv(os.path.join(save_data, 'dm.csv'), index=False)
            
        # Clustering
        embedding = MDS(n_components=2, dissimilarity='precomputed', random_state=self.random_state)
        dm_embed = embedding.fit_transform(dm)
        bandwidth = estimate_bandwidth(dm_embed, quantile=0.3, random_state=self.random_state)
        bandwidth = max(bandwidth, 1e-5)
        clustering = MeanShift(bandwidth=bandwidth)
        clustering_predictions = clustering.fit_predict(dm_embed)
        
        # Normalize points
        mean = np.mean(dm_embed, axis=0)
        std = np.std(dm_embed, axis=0)
        dm_embed_normalized = (dm_embed - mean) / std
        if sum(std) == 0:
            dm_embed_normalized = dm_embed
        
        # Sort unique clusters by size
        label_counts = Counter(clustering_predictions)
        unique_clusters = sorted(label_counts.items(), key=lambda x: x[1], reverse=True)
        if save_data is not None and os.path.isdir(save_data):
            with open(os.path.join(save_data, 'clusters.txt'), 'w') as file:
                for item in unique_clusters:
                    file.write(str(item) + '\n')
        
        # Predictions
        pred_seqs = pd.DataFrame(columns=pred_seq_df.columns)
        cluster_counts = []
        for class_ in unique_clusters[:n_clusters]:
            # Take riskiest strain using this cluster
            wanted_names = dm.columns[clustering_predictions == class_[0]]
            cluster_seq_df = seq_df.iloc[wanted_names]
            pred_seq_df_with_risk = self._compute_risk_for_predict_domseq(cluster_seq_df, pred_seq_df, enet)
            if save_data is not None and os.path.isdir(save_data):
                pred_seq_df_with_risk.to_csv(os.path.join(save_data, 'pred_seq_df_'+str(class_)+'.csv'), index=False)
            
            # Save predictions
            pred_seqs = pred_seqs.append(pred_seq_df_with_risk.iloc[0])
            cluster_counts.append(len(wanted_names))
        
        # Cluster areas
        cluster_areas = self._plot_embedding(dm_embed_normalized, unique_clusters, 
                                             clustering_predictions, n_clusters, save_data, ALPHA)
            
        pred_seqs['cluster_count'] = cluster_counts
        pred_seqs['cluster_area'] = cluster_areas
        return pred_seqs


    def predict_single_domseq(self, pred_seqs:pd.DataFrame, pred_seq_df:pd.DataFrame) -> pd.DataFrame:
        ''' Predicts a single future dominant sequence with a single cluster.
        
        Parameters
        ----------
        pred_seqs - Emergenet recommended sequences, with additional column 'cluster_area'
        
        pred_seq_df - DataFrame of candidate sequences
            
        Returns
        -------
        pred_seq - Emergenet recommended sequence
        ''' 
        if len(pred_seqs) < 2:
            raise ValueError('Not enough prediction sequences!')
        if len(pred_seq_df) < 1:
            raise ValueError('The DataFrame contains no sequences!')
            
        # Sort by cluster area and get top two cluster predictions
        pred_seqs = pred_seqs.sort_values(by='cluster_count',ascending=False)
        cluster_sizes = pred_seqs['cluster_count'].values
        cluster_ratio = cluster_sizes[1]/cluster_sizes[0]
        seq1 = pred_seqs['sequence'][0][:self.seq_trunc_length]
        seq2 = pred_seqs['sequence'][1][:self.seq_trunc_length]
        dist = hamming(seq1, seq2)
        
        # If there is one dominant cluster, return the prediction from that cluster
        if seq1 == seq2 or cluster_ratio <= 1/10:
            return pred_seqs.head(1)
        ratios = []
        
        # Find the strain x closest to satisfying dist(seq1,x)/dist(seq2,x) = cluster2/cluster1
        for seq in pred_seq_df['sequence'].values:
            dist1 = hamming(seq1, seq[:self.seq_trunc_length])
            dist2 = hamming(seq2, seq[:self.seq_trunc_length])
            if dist2 == 0 or dist1 + dist2 > dist:
                ratios.append(1e5)
                continue
            ratios.append(abs(dist1/dist2 - cluster_ratio))
        pred_seq_df['ratio'] = ratios
        pred_seq_df = pred_seq_df.sort_values(by='ratio', ascending=True)
        pred_seq = pred_seq_df.head(1)
        return pred_seq
