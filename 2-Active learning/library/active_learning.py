import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from matplotlib import pyplot as plt

def sample_new_combination(parameter, nb_sample = 100000, seed=None):
    """
    Generates random combinations of parameter values by sampling from each parameter space independently.

    Parameters:
        parameter (list of lists or arrays): Each sublist contains possible values for one parameter.
        nb_sample (int): Number of combinations to sample.
        seed (int or None): Random seed for reproducibility.

    Returns:
        np.ndarray: An array of shape (nb_sample, len(parameter)) containing the sampled combinations.
    """
    samples = []
    np.random.seed(seed)
    for i in range(len(parameter)):
        choice = np.random.choice(parameter[i] , size = nb_sample, replace = True)
        samples.append(choice)
    samples = np.array(samples).transpose()
    return samples 


def sampling_without_repeat(sampling_condition, nb_sample, exited_data, seed = None):
     """
    Samples new parameter combinations while ensuring they are not already present in the existing dataset.

    Parameters:
        sampling_condition (list of lists or arrays): Parameter space for sampling.
        nb_sample (int): Number of unique combinations to sample.
        exited_data (array-like): Existing combinations to avoid.
        seed (int or None): Random seed for reproducibility.

    Returns:
        np.ndarray: Array of new, non-redundant combinations of shape (nb_sample, len(sampling_condition)).
    """
     if not isinstance(exited_data, np.ndarray):
        exited_data = np.array(exited_data)

     new_comb = sample_new_combination(sampling_condition, nb_sample, seed=seed)

     #check if new data already exits and resample
     while True:
        matches = np.all(new_comb[:, np.newaxis] == exited_data, axis=-1)
        rows_to_drop = np.any(matches, axis=1)

        if not np.any(rows_to_drop):
                # No rows to be dropped, so we are done
                break
        
        # resampling the same repeated number
        nb_repeat = sum(rows_to_drop)
        resample = sample_new_combination(sampling_condition, nb_sample = nb_repeat, seed= None)
        new_comb = new_comb[~rows_to_drop]
        new_comb = np.concatenate([new_comb, resample])

     return new_comb


def find_top_element(X, y, cluster_list, condition, n, return_ratio = False, verbose = True):
    """
    Selects the top 'n' elements based on a specified condition (e.g., acquisition function).

    Parameters:
        X (array-like): Input data points.
        y (array-like): Predicted target values for X.
        cluster_list (array-like): Cluster labels corresponding to X.
        condition (array-like): Scores (e.g., acquisition function values) used to rank X.
        n (int): Number of top elements to select.
        return_ratio (bool): If True, also return condition values of selected points.
        verbose (bool): If True, print the maximum predicted yield.

    Returns:
        choosen_X (np.ndarray): Selected top data points.
        choosen_y (np.ndarray): Predicted yields of selected points.
        ratio (list or np.ndarray): Condition values of selected points (empty if return_ratio is False).
        choosen_cluster (np.ndarray): Cluster assignments of selected points.
    """
        # Convert X and y to NumPy arrays if they are not already
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    # Sort the list in ascending order (smallest to largest)
    idx_top_n = np.argsort(-condition)[:n]
    choosen_X = X[idx_top_n,:]
    choosen_y = y[idx_top_n]
    choosen_cluster = cluster_list[idx_top_n]

    if return_ratio:
        ratio = condition[idx_top_n]
#        ratio = choosen_y/(ucb - choosen_y)
    else:
        ratio = []

    if verbose:
        print(f"Maximum yield prediction = {max(choosen_y)}")
    return choosen_X, choosen_y, ratio, choosen_cluster


def find_update_top_element(X, y, condition, n, verbose = True):
    """
    Selects the top 'n' elements based on a given condition and returns both selected and unselected sets.

    Parameters:
        X (array-like): Input data points.
        y (array-like): Predicted target values for X.
        condition (array-like): Scores used to rank X (e.g., acquisition function values).
        n (int): Number of top elements to select.
        verbose (bool): If True, print the maximum yield of selected samples.

    Returns:
        choosen_X (np.ndarray): Selected top data points.
        choosen_y (np.ndarray): Corresponding target values.
        unselected_X (np.ndarray): Remaining data points.
        unselected_y (np.ndarray): Remaining target values.
    """
        # Convert X and y to NumPy arrays if they are not already
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    if not isinstance(y, np.ndarray):
        y = np.array(y)

    # Sort the list in ascending order (smallest to largest)
    idx_top_n = np.argsort(-condition)[:n]
    choosen_X = X[idx_top_n,:]
    choosen_y = y[idx_top_n]
    
    unselected_X = np.delete(X, idx_top_n, axis=0)
    unselected_y = np.delete(y, idx_top_n)

    if verbose:
        print(f"Maximising sample has yield = {max(choosen_y)}")
        
    return choosen_X, choosen_y, unselected_X, unselected_y



def active_found(ensemble_models, sampling_condition, nb_new_data, X, scaler,verbose = True):
    """
    Identifies promising new experiments using an ensemble model and UCB-based active learning strategy.

    Parameters:
        ensemble_models (list): List of trained models for prediction.
        sampling_condition (DataFrame or array): Space of possible new experiments.
        nb_new_data (int): Number of new experiments to select.
        X (array-like): Existing experimental data to avoid repeats.
        scaler (object): Scaler used to normalize input features.
        verbose (bool): If True, print selection summaries.

    Returns:
        ucb_top (np.ndarray): Top points selected using UCB strategy.
        ratio (np.ndarray): Acquisition scores of UCB-selected points.
        exploit_top (np.ndarray): Top points based on predicted mean yield.
        explore_top (np.ndarray): Top points based on prediction uncertainty.
    """
    X_new= sampling_without_repeat(sampling_condition, nb_sample = nb_new_data*5, exited_data=X)

    y_pred = []
    X_new_normalized = scaler.transform(X_new)
    for i in range(len(ensemble_models)):
        model = ensemble_models[i]
        y_pred.append(model.predict(X_new_normalized))

    y_pred = np.array(y_pred)
    y_mean = np.average(y_pred, axis=0)
    y_stdv = np.std(y_pred, axis=0)

    ucb = 1*y_mean + 1.14*y_stdv 

    print("For UCB:")
    ucb_top, _, ratio = find_top_element(X_new, y_mean,ucb, nb_new_data, return_ratio= True, verbose = verbose)
    print("For exploitation:")
    exploit_top, _, _ = find_top_element(X_new, y_mean, y_mean, nb_new_data, verbose)
    print("For exploration:")
    explore_top, _, _ = find_top_element(X_new,y_mean, y_stdv, nb_new_data, verbose)

    return ucb_top, ratio, exploit_top, explore_top


def probability_of_improvement(mu, sigma, current_best):
    """
    Calculate Probability of Improvement (PI) for Gaussian process predictions.

    Parameters:
    - mu: Mean of the Gaussian process prediction.
    - sigma: Standard deviation of the Gaussian process prediction.
    - current_best: Current best observed value.
    Returns:
    - pi: Probability of Improvement.
    """

    # Avoid division by zero
    sigma = sigma + 1e-4

    # Calculate standard normal cumulative distribution function
    z = (mu - current_best) / sigma
    pi = norm.cdf(z)

    return pi


def expected_improvement(mu, sigma, current_best):
    """
    Calculate Expected Improvement (EI) for Gaussian process predictions.

    Parameters:
    - mu: Mean of the Gaussian process prediction.
    - sigma: Standard deviation of the Gaussian process prediction.
    - current_best: Current best observed value.
    - epsilon: Small positive constant to avoid division by zero.

    Returns:
    - ei: Expected Improvement.
    """

    # Avoid division by zero
    sigma = sigma + 1e-4

    # Calculate standard normal cumulative distribution function
    z = (mu - current_best) / sigma
    ei = (mu - current_best) * norm.cdf(z) + sigma * norm.pdf(z)

    return ei

def upper_confident_bound(mu, sigma, theta, r2):
    """
    Calculate UCB for Gaussian process predictions.

    Parameters:
    - mu: Mean of the Gaussian process prediction.
    - sigma: Standard deviation of the Gaussian process prediction.
    - epsilon: Small positive constant to avoid division by zero.

    Returns:
    - ucb: 
    """
    if r2 <= 0.8:
        ucb = mu + theta*sigma
    else:
        ucb = mu
    return ucb


def cluster(X, n_clusters, plot = False):
    """
    Performs hierarchical clustering on input data and assigns cluster labels.

    Parameters:
        X (array-like): Input feature matrix.
        n_clusters (int): Number of clusters to form.
        plot (bool): If True, display the dendrogram.

    Returns:
        clusters (np.ndarray): Cluster labels for each sample.
    """
    Z = linkage(X, 'average')
    clusters = fcluster(Z, n_clusters, criterion='maxclust')

    if plot: # Plot dendrogram       
        plt.figure(figsize=(10, 5))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index')
        plt.ylabel('distance')
        dendrogram(Z)
        plt.show()    
    
    return clusters

def pick_by_cluster(choosen_clust, n_sample, seed = 42, verbose = False):
    """
    Selects a specified number of samples by drawing one from each cluster in a round-robin fashion.

    Parameters:
        choosen_clust (array-like): Cluster labels for each sample.
        n_sample (int): Total number of samples to select.
        seed (int): Random seed for reproducibility.
        verbose (bool): If True, print cluster details.

    Returns:
        selected_indices (list): Indices of the selected samples.
    """
    np.random.seed(seed)

    # Initialize a dictionary to store cluster information
    cluster_info = {}
    # Iterate over unique cluster names
    for cluster_name in np.unique(choosen_clust):
        # Find indices of data points belonging to the current cluster
        indices = np.where(choosen_clust == cluster_name)[0]
        # Count the number of points in the cluster
        num_points = len(indices)
        # Store cluster information in the dictionary
        cluster_info[cluster_name] = {'num_points': num_points, 'indices': indices}

    sorted_clusters = sorted(cluster_info.items(), key=lambda x: x[1]['num_points'], reverse=False)

    if verbose:
        # Print cluster information
        for cluster_name, info in sorted_clusters:
            print(f"Cluster {cluster_name}: {info['num_points']} points, Indices: {info['indices']}")
            
    # Initialize a list to store the selected indices
    selected_indices = []

    # Initialize a counter for the number of selected points
    selected_count = 0

    # Iterate over sorted clusters and select one index from each cluster
    while selected_count < n_sample:
        for cluster_name, info in sorted_clusters:
            # Check if there are any points left in this cluster
            if info['num_points'] > 0:
                # Randomly select one index from this cluster
                selected_index = np.random.choice(info['indices'])
                # Add the selected index to the list
                selected_indices.append(selected_index)
                # Decrement the number of points in this cluster
                cluster_info[cluster_name]['num_points'] -= 1
                # Increment the total count of selected points
                selected_count += 1
                # Check if we have selected enough points
                if selected_count == n_sample:
                    break

    return selected_indices

def pick_by_cluster(choosen_clust, n_sample, seed=42, verbose=False):
    """
    Selects samples with highest UCB values from each cluster in a round-robin fashion.

    Parameters:
        choosen_clust (np.ndarray): 2D array where column 0 is cluster label and column 1 is UCB value.
        n_sample (int): Number of total samples to select.
        seed (int): Random seed for reproducibility.
        verbose (bool): If True, print number of points per cluster.

    Returns:
        selected_indices (list): Indices of the selected high-UCB samples.
    """
    np.random.seed(seed)

    # Initialize a dictionary to store cluster information
    cluster_info = {}
    # Iterate over unique cluster names
    for cluster_name in np.unique(choosen_clust[:, 0]):
        # Find indices of data points belonging to the current cluster
        indices = np.where(choosen_clust[:, 0] == cluster_name)[0]
        # Extract UCB values for the points in this cluster
        ucb_values = choosen_clust[indices, 1]
        # Count the number of points in the cluster
        num_points = len(indices)
        # Store cluster information in the dictionary
        cluster_info[cluster_name] = {'indices': indices, 'ucb_values': ucb_values, 'num_points': num_points}

    sorted_clusters = sorted(cluster_info.items(), key=lambda x: x[1]['num_points'], reverse=False)

    if verbose:
        # Print cluster information
        for cluster_name, info in sorted_clusters:
            print(f"Cluster {cluster_name}: {info['num_points']} points")

    # Initialize a list to store the selected indices
    selected_indices = []

    # Initialize a counter for the number of selected points
    selected_count = 0

    # Iterate over sorted clusters and select one index from each cluster with highest UCB
    while selected_count < n_sample:
        for cluster_name, info in sorted_clusters:
            # Check if there are any points left in this cluster
            if info['indices'].size > 0:
                # Select the index with the highest UCB from this cluster
                selected_index = info['indices'][np.argmax(info['ucb_values'])]
                # Add the selected index to the list
                selected_indices.append(selected_index)
                # Remove the selected index from the list of available indices in this cluster
                info['indices'] = np.delete(info['indices'], np.argmax(info['ucb_values']))
                # Decrement the total count of selected points
                selected_count += 1
                # Check if we have selected enough points
                if selected_count == n_sample:
                    break

    return selected_indices