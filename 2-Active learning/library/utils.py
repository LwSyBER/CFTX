import glob
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.preprocessing import MaxAbsScaler


def import_data(folder_for_data, verbose = True):
    """
    Imports and concatenates data from all CSV files in a specified folder.

    Parameters:
    - folder_for_data (str): Path to the folder containing the CSV files.
    - verbose (bool, optional): If True, prints the number and names of files read and displays the concatenated DataFrame. Default is True.

    Returns:
    - concatenated_data (DataFrame): A pandas DataFrame containing all data concatenated from the CSV files.
    - size_list (list of int): A list containing the number of rows in each individual file.

    The function reads all `.csv` files in the specified folder, checks for column consistency, and concatenates the data. 
    It assumes that all files have the same structure and appends them into a single DataFrame. 
    If `verbose` is enabled, it outputs a summary of the import process and displays the full dataset.
    """
    files = glob.glob(folder_for_data + "\\*.csv")

    # Initialize a DataFrame with NaN values
    concatenated_data = pd.DataFrame()
    size_list = []

    # Initialize a flag to check if columns are consistent across files
    consistent_columns = True

    # Iterate through files
    for file in files:
        df = pd.read_csv(file)

        # Check if columns are consistent across files
        if concatenated_data.columns.empty:
            concatenated_data = df
            size_list.append(len(df))
#        elif not concatenated_data.columns.equals(df.columns):
#            consistent_columns = False
#            print(f"Ignoring file {file}: Column orders are not consistent.")
#            files.remove(file)
        else: 
            concatenated_data = pd.concat([concatenated_data, df], ignore_index=True)
            size_list.append(len(df))

    if verbose:
        print("Read ", len(files), " files: ")
        for file in files:
            data_name = os.path.basename(file)
            print("- ", data_name)
        display(concatenated_data)

        if consistent_columns:
            print("All files have consistent column orders.")
        else:
            print("Some files have inconsistent column orders.")

    return concatenated_data, size_list


def select_from_iteration (data, selected_plate = (0,1)):
    """
    Extracts and concatenates data from a specified range of iterations (plates).

    Parameters:
    - data (list of arrays): A list where each element corresponds to data from one iteration or plate.
    - selected_plate (tuple of int, optional): A tuple indicating the start (inclusive) and end (exclusive)
      indices of the plates to select. Defaults to (0, 1).

    Returns:
    - current_data (ndarray): A single NumPy array containing the concatenated data from the selected plates.

    This function is useful for focusing analysis on a specific subset of iterations, rather than processing all data at once.
    """
    selected_plates = data[selected_plate[0]:selected_plate[1]]
    current_data = np.concatenate(selected_plates, axis = 0)
    return(current_data)


def import_parameter(parameter_file, nb_new_data_predict, sep = ';', verbose = False):
    """
    Imports and processes parameter data from a CSV file to generate concentration matrices 
    for combinatorial experiments.

    Parameters:
    - parameter_file (str): Path to the CSV file containing parameter information.
    - nb_new_data_predict (int): Number of new data points to be predicted; used to compute the coverage ratio.
    - sep (str, optional): Column separator in the CSV file (default is ';').
    - verbose (bool, optional): If True, prints detailed information about the imported data and displays possible combinations.

    Returns:
    - element_list (list): Names of the elements/metabolites.
    - element_max (list): Maximum concentration values for each element.
    - concentration (ndarray): Matrix of scaled concentration values across all conditions.
    
    The function reads the input file, scales relative concentration values by their respective maximum values,
    calculates the total number of possible combinations, and optionally displays a preview of the generated concentration matrix.
    """
    df = pd.read_csv(parameter_file, sep = sep)
    element_list = df.iloc[:,0].to_list()
    element_max = df.iloc[:,1].to_list()
    concentration = df.iloc[:,2:].to_numpy()

    multi_max = np.repeat([element_max], concentration.shape[1], axis = 0).T
    concentration = multi_max*concentration

    row,col = concentration.shape
    combi = col**row
    if verbose:
        print(f"Number of metabolites : {len(element_list)}")
        print(f"Number of combinations - poolsize : {combi}")
        print(f"Searching ratio : {round(nb_new_data_predict*100/combi, 2)} %")
        print(f"Possible concentrations: ")
        df = pd.DataFrame(concentration.T, columns= element_list)
        display(df)
    return element_list, element_max, concentration   


def check_column_names(data,target,element_list):
    """
    Verifies that the column names in the data match the expected element_list.

    Parameters:
    - data (DataFrame): The input DataFrame containing experimental or observational data.
    - target (str): The name of the column to exclude from the comparison (e.g., the target/output variable).
    - element_list (list of str): The list of expected column names (typically from a parameter file).

    Behavior:
    - Compares the column names in `data` (excluding `target`) with `element_list`.
    - Prints the names of columns that are present only in the data or only in the parameter list.
    - Raises and prints a ValueError if there is any mismatch; otherwise, confirms that all names match.

    This function is useful for ensuring consistency between data and parameter definitions before downstream processing.
    """
    element_data = data.columns.drop(target)

    # Find the different column names
    columns_only_in_data = set(element_data).difference(element_list)

    # Find columns in df2 but not in df1
    columns_only_in_parameter = set(element_list).difference(element_data)

    try:
        if columns_only_in_data or columns_only_in_parameter:
            print("- Columns only in data files:", columns_only_in_data)
            print("- Columns only in parameter files:", columns_only_in_parameter)
            raise ValueError("Columns names are not matched, please modify parameter column names")
        else:
            print("All column names matched!")
    except ValueError as e:
        text = "{}: {}".format(type(e).__name__, e)
        print("\033[1;91m{}\033[0m".format(text))


#def normalize_data(data):
#     """
#     Normalizes input features using Min-Max scaling across concatenated data arrays.
#
#     Parameters:
#     - data (list of arrays): A list where each element is a NumPy array containing input features 
#       and target values. Each array is expected to have the last two columns as target and possibly 
#       an additional output or metadata column.
#
#     Returns:
#     - scaled (ndarray): The Min-Max scaled feature matrix (excluding the last two columns).
#     - y (ndarray): The target values extracted from the second-to-last column.
#
#     The function concatenates all input arrays, extracts the feature matrix `X` (excluding the last two columns),
#     and the target vector `y` (assumed to be in the second-to-last column). It applies Min-Max scaling
#     to `X` to transform features to the [0, 1] range.
#     """
#    data_array = np.concatenate(data,axis = 0)
#    X = data_array[:,:-2]
#    y = data_array[:,-2]
#    scaler = MinMaxScaler()
#    scaled = scaler.fit_transform(X)
#    return scaled, y

def import_split_data(data_folder, element_list, target, type = 'first', idx = (102,204), seed = None):
    """
    Imports data and splits it into training, pool, and test sets based on the specified method.

    Parameters:
    - data_folder (str): Path to the folder containing the data files.
    - element_list (list of str): List of column names to be used as features for training and prediction.
    - target (list of str): List containing the target column(s) for prediction.
    - type (str, optional): Method for selecting the training data. Options are:
        - 'first': Selects the first data segment based on the index range `idx` for training.
        - 'random': Randomly splits the data into training and pool sets. Default is 'first'.
    - idx (tuple of int, optional): The start and end indices (inclusive) for the training data when `type='first'`. Default is (102, 204).
    - seed (int, optional): Random seed for reproducibility of the split. Default is None.

    Returns:
    - X_train (DataFrame): Features for the training set.
    - X_pool (DataFrame): Features for the pool set.
    - X_test (DataFrame): Features for the test set.
    - y_train (Series): Target values for the training set.
    - y_pool (Series): Target values for the pool set.
    - y_test (Series): Target values for the test set.

    The function imports the data from the specified folder and splits it into training, pool, and test sets.
    The training set can be selected either from the first portion of the data (`type='first'`) or randomly (`type='random'`).
    The pool and test sets are further split randomly.
    """
    data = import_data(data_folder, verbose = False)
    if type == 'first':
        train = data.loc[idx[0]:idx[1]]
        data = data.drop(range(idx[0],idx[1]))

        X_train, y_train = train[element_list], train[target[0]]
        X_pool, y_pool = data[element_list], data[target[0]]

    if type == 'random':
        X_pool, y_pool = data[element_list], data[target[0]]
        X_pool, X_train, y_pool, y_train = train_test_split(X_pool, y_pool, test_size=0.1, random_state=seed)

    X_pool, X_test, y_pool, y_test = train_test_split(X_pool, y_pool, test_size=0.11, random_state=seed)

    return X_train, X_pool, X_test, y_train, y_pool, y_test


def normalized(X_train = [], X_test = [], X_pool = []):
    """
    Normalizes the input features using MaxAbsScaler.

    Parameters:
    - X_train (array-like, optional): Feature matrix for the training set. Default is an empty list.
    - X_test (array-like, optional): Feature matrix for the test set. Default is an empty list.
    - X_pool (array-like, optional): Feature matrix for the pool set. Default is an empty list.

    Returns:
    - X_normalized (ndarray): The normalized feature matrix for the training set.
    - X_test_normalized (ndarray): The normalized feature matrix for the test set, or an empty list if not provided.
    - X_pool_normalized (ndarray): The normalized feature matrix for the pool set, or an empty list if not provided.

    This function applies MaxAbs scaling to the training, test, and pool feature matrices. The MaxAbsScaler scales each feature
    by its maximum absolute value, which is useful when the data is sparse and needs scaling without shifting the data.
    If `X_test` or `X_pool` are not provided, their normalized versions will be returned as empty lists.
    """
    scaler = MaxAbsScaler()
    scaler.fit(X_train)
    X_normalized = scaler.transform(X_train)
    X_test_normalized, X_pool_normalized = [], []
    
    if len(X_test) > 0:
        X_test_normalized = scaler.transform(X_test)

    if len(X_pool) > 0:
        X_pool_normalized = scaler.transform(X_pool)

    return X_normalized, X_test_normalized, X_pool_normalized


def find_outlier(experiment,condi):
    """
    Identifies outliers in the experiment data based on the given condition.

    Parameters:
    - experiment (array-like): The data from the experiment, such as a NumPy array or pandas DataFrame.
    - condi (array-like): A boolean array or condition that indicates which data points are outliers.

    Returns:
    - result (array-like): A list of data points from the `experiment` that are considered outliers based on the condition `condi`.
      If the sum of `condi` is greater than 1, returns an empty list. Otherwise, returns the data points where `condi` is False.

    The function checks whether the sum of the boolean condition `condi` exceeds 1. If so, it returns an empty list, assuming
    that there are no outliers. Otherwise, it returns the data points from `experiment` where the condition is False, which are
    assumed to be the outliers.
    """
    if sum(condi) > 1:  
        result= []
    else: 
        result=experiment[~condi]
    return result


def remove_outlier_strict(sample):
    """
    Removes outliers from the given sample data based on Z-scores, applying stricter thresholds for standard deviations.

    Parameters:
    - sample (array-like): A 2D array or list containing the sample data to be processed. The rows represent different data points,
      and the columns represent different features for each data point.

    Returns:
    - result (list of lists): A list where each entry corresponds to the outliers detected in each row of the sample data. The outliers 
      are identified by their Z-scores being above a dynamic threshold based on the standard deviation.

    This function calculates the Z-scores for each feature in the sample data. It compares the Z-scores to a threshold that is stricter 
    when the standard deviation of the data row is below 0.2 (threshold of 1.4) and less strict (threshold of 1) when the standard deviation
    is above 0.2. The function uses the `find_outlier` function to identify the outliers and returns them in a list.
    """
    if not isinstance(sample, np.ndarray):
        sample = np.array(sample)

    mean = np.mean(sample, axis = 1)
    mean = np.tile(mean,(sample.shape[1], 1)).T

    std = np.std(sample, axis = 1)
    std = np.tile(std,(sample.shape[1], 1)).T

    z_scores = np.abs((sample-mean)/std)
#    outliers_z = np.abs(z_scores) > threshold
    result = []
    for i in range(len(z_scores)):
        z_row = z_scores[i]
        data_row = sample[i]
        std_row = std[i][0]
        if std_row > 0.2:
            condi = z_row > 1
        else:
            condi = z_row > 1.4
        result.append(find_outlier(data_row,condi))
    return result

def remove_outlier(sample, threshold = 0.2):
    """
    Removes outliers from the given sample data based on Z-scores, applying a dynamic threshold.

    Parameters:
    - sample (array-like): A 2D array or list containing the sample data to be processed. The rows represent different data points,
      and the columns represent different features for each data point.
    - threshold (float, optional): The Z-score threshold above which values are considered outliers. Default is 0.2.

    Returns:
    - result (list of lists): A list where each entry corresponds to the outliers detected in each row of the sample data. 
      If the standard deviation for a row exceeds the threshold, outliers are identified based on a Z-score greater than 1.

    This function calculates the Z-scores for each feature in the sample data. It compares the Z-scores to a threshold that is applied
    to identify outliers. If the standard deviation of a row exceeds the threshold, outliers are determined by Z-scores greater than 1.
    The `find_outlier` function is used to isolate these outliers, and the result is a list where outliers are separated from the rest of the data.
    """
    if not isinstance(sample, np.ndarray):
        sample = np.array(sample)

    mean = np.mean(sample, axis = 1)
    mean = np.tile(mean,(sample.shape[1], 1)).T

    std = np.std(sample, axis = 1)
    std = np.tile(std,(sample.shape[1], 1)).T

    z_scores = np.abs((sample-mean)/std)
    outliers_z = np.abs(z_scores) > threshold
    result = []
    for i in range(len(z_scores)):
        z_row = z_scores[i]
        data_row = sample[i]
        std_row = std[i][0]
        if std_row > threshold:
            condi = z_row > 1
            result.append(find_outlier(data_row,condi))
        else:
            result.append(data_row)
    return result


def flatten_X_y(X, y):
    """
    Flattens the input data arrays X and y by expanding the target values and removing NaNs.

    Parameters:
    - X (array-like): A 2D array or list containing feature data. Each row corresponds to a data point, and each column represents a feature.
    - y (array-like): A 2D array or list where each row corresponds to the target values associated with the data point in X.

    Returns:
    - X_flat (array): A flattened array where each data point from X is repeated for each corresponding target value in y.
    - y_flat (array): A flattened array of the target values from y, with NaNs removed.

    This function processes the input arrays X and y by flattening them. For each row in X, the corresponding target values from y
    are expanded, and then both X and y are flattened into 1D arrays. Any NaN values in y are removed, and the corresponding rows
    in X are also discarded to maintain consistency.
    """
    X_flat = []
    for i in range(X.shape[0]):
        x_loop = X[i]
        y_loop = y[i]
        for _ in y_loop:
            X_flat.append(x_loop)
    X_flat = np.array(X_flat)
    y_flat = np.array([element for sublist in y for element in sublist])
    
    idx = np.isnan(y_flat)
    y_flat = y_flat[~idx]
    X_flat = X_flat[~idx]

    return X_flat, y_flat


def average_and_drop_na(X,y):
    """
    Averages the target values in y, removes rows with NaN values, and returns the cleaned data.

    Parameters:
    - X (array-like): A 2D array or list containing feature data. Each row represents a data point, and each column represents a feature.
    - y (array-like): A 2D array or list where each row contains the target values associated with the corresponding data point in X.

    Returns:
    - X (array): A 2D array with rows corresponding to data points from X, with NaN rows removed.
    - y (array): A 1D array of averaged target values from y, with NaN values removed.

    This function computes the mean of the target values for each row in y (ignoring NaNs), removes any rows where the target value is NaN,
    and returns the cleaned version of both X and y.
    """
    if len(y)== 0:
        return X, y
    
    y = np.nanmean(y, axis = 1)
    idx = np.isnan(y)
    y = y[~idx]
    X = X[~idx]
    return X, y


def split_and_flatten(medium, yield_array, ratio = 0.2, seed = None, flatten = True):
    """
    Splits the data into training and testing sets, flattens or averages the data, and handles missing values.

    Parameters:
    - medium (array-like): A 2D array or list containing feature data. Each row represents a data point, and each column represents a feature.
    - yield_array (array-like): A list or 2D array containing target values associated with the data points in `medium`.
    - ratio (float, optional): The proportion of data to include in the test set (default is 0.2, meaning 20% test data).
    - seed (int, optional): The random seed for splitting the data into training and testing sets (default is None).
    - flatten (bool, optional): Whether to flatten the data after splitting (default is True). If False, the function will average the target values and drop NaNs.

    Returns:
    - X_train (array): The training data features after splitting and processing.
    - X_test (array): The testing data features after splitting and processing.
    - y_train (array): The training target values after processing (flattened or averaged).
    - y_test (array): The testing target values after processing (flattened).

    This function performs the following steps:
    1. Pads the `yield_array` to ensure all arrays are of equal length.
    2. Splits the data (`medium` and `yield_array`) into training and test sets based on the specified ratio.
    3. Optionally flattens or averages the data.
    4. Removes NaN values from the target array.

    If `flatten` is set to `True`, the function flattens the training data and corresponding targets. If set to `False`, it averages the target values and removes NaN entries.
    """
    medium = np.array(medium)
    yield_array = np.array([np.pad(arr, (0, len(max(yield_array, key=len)) - len(arr)), 'constant', constant_values=np.nan) for arr in yield_array])
    _, unique_indices = np.unique(medium, axis=0, return_index=True)
    if len(unique_indices) != len(medium):
        print('BE CAREFUL! medium has repeatition')

    # Split 
    if ratio == 0:
        X_train, X_test = medium, np.array([])
        y_train, y_test = yield_array, np.array([])
    else:
        train_indices, test_indices = train_test_split(unique_indices, test_size=ratio, random_state=seed)
        X_train, X_test = medium[train_indices],  medium[test_indices]
        y_train = [yield_array[i] for i in train_indices]   # yield_array is list sometime
        y_test = [yield_array[i] for i in test_indices]

    # Flatten 
    if flatten:
       X_train, y_train = flatten_X_y(X_train, y_train) 
    else: 
        X_train, y_train = average_and_drop_na(X_train, y_train)
    
    X_test, y_test = flatten_X_y(X_test, y_test)
    return X_train, X_test, y_train, y_test

def euclidean_distance(point1, point2):
    """
    Computes the Euclidean distance between two points in a multi-dimensional space.

    Parameters:
    - point1 (array-like): A 1D or 2D array representing the first point(s) in the space. 
      If a 2D array is provided, each row corresponds to a separate point.
    - point2 (array-like): A 1D or 2D array representing the second point(s) in the space.
      It should have the same shape as `point1`. 

    Returns:
    - distance (array): The Euclidean distance(s) between `point1` and `point2`.
      If both `point1` and `point2` are arrays of points, the function returns an array of distances for each pair of corresponding points.

    The function calculates the Euclidean distance using the formula:
    `sqrt(sum((point1 - point2)**2))`, where `point1` and `point2` can be vectors or sets of points.
    """
    return np.sqrt(np.sum((point1 - point2)**2, axis = 1))


#########################################################################################
def plot_loss(model_list, highlight=False):
    """
    Plots the loss curve for multiple models, comparing their training losses over iterations.

    Parameters:
    - model_list (list): A list of models that have a `loss_curve_` attribute (e.g., models trained using sklearn).
                          Each model in the list should have a `loss_curve_` attribute representing the training loss at each iteration.
    - highlight (bool, optional): If set to True, the model with the minimum final loss will be highlighted in red, 
                                  and others will be shown in blue. Default is False (all models in blue).

    Returns:
    - None: Displays the plot of training loss curves for all models.

    The function creates subplots for each model, where the x-axis represents the iteration number, 
    and the y-axis represents the corresponding loss value. If `highlight` is True, the model with the 
    minimum loss at the end of the training will be shown in red.
    """
    if not isinstance(model_list, list):
        model_list = [model_list]

    num_models = len(model_list)
    num_cols = round(np.sqrt(num_models))
    num_rows = (num_models + num_cols - 1) // num_cols

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 8))
    for i in range(num_cols*num_rows): # delete empty subplots
        if i >= num_models:
            fig.delaxes(axes.flatten()[i])
    
    min_loss_index = np.argmin([model.loss_ for model in model_list])

    for i, (model, ax) in enumerate(zip(model_list, np.ravel(axes))):
        losses = model.loss_curve_
        iterations = range(1, len(losses) + 1)

        if highlight:
            line_color = 'red' if i == min_loss_index else 'blue'
        else:
            line_color = 'blue'

        ax.plot(iterations, losses, color=line_color, label='Training loss')
        ax.set_xlabel('Iterations')
        ax.set_ylabel('Loss')
        ax.set_title(f'Loss of Model {i + 1} : {round(model.loss_, 5)}')

    plt.tight_layout()
    plt.show()  


def plot_feature(data, label):
    """
    Plots histograms of the features in the provided data.

    Parameters:
    - data (pandas.DataFrame): The DataFrame containing the features to be visualized. 
      Each column of the DataFrame represents a feature.
    - label (str): A label for the plot title, used to specify the feature set being visualized.

    Returns:
    - None: Displays the histograms of the features in the DataFrame.
    """
    data.hist(bins=4, color='blue', edgecolor='black', figsize=(8, 6))

    # Add labels and title
    plt.suptitle(f'Histograms for {label} features')
    plt.tight_layout()
    plt.show()


def plot_r2_curve(y_true, y_pred):
    """
    Plots a scatter plot comparing the experimental values (`y_true`) with the predicted values (`y_pred`), 
    and includes a line representing the ideal prediction (where experimental values equal predicted values). 
    The R-squared (R²) value is also displayed on the plot.

    Parameters:
    - y_true (array-like): The ground truth or experimental values.
    - y_pred (array-like): The predicted values from the model.

    Returns:
    - None: Displays the plot with the scatter plot and R² value.
    """
    r2 = r2_score(y_true, y_pred)
    TRUE, PRED = y_true, y_pred
    sns.set(
        font="arial",
        palette="colorblind",
        style="whitegrid",
        font_scale=1.5,
        rc={"figure.figsize": (5, 5), "axes.grid": False},
    )
    sns.regplot(
        x=TRUE,
        y=PRED,
        fit_reg=0,
        marker="+",
        color="black",
        scatter_kws={"s": 40, "linewidths": 0.7},
    )
    plt.plot([min(TRUE), max(TRUE)], [min(TRUE), max(TRUE)], 
             linestyle='--', 
             color='blue',
             linewidth=1)
    plt.xlabel("Experiment ground truth ")
    plt.ylabel("Model prediction")
    plt.title(f'R2: {r2:.2f}', fontsize=14)
    plt.xlim(min(TRUE) - 0.2, max(TRUE) + 0.5)
    plt.ylim(min(PRED) - 0.2, max(PRED) + 0.5)
    plt.show()


def plot_hist_yield_std(mean, training_mean, std, training_std):
    """
    Plots side-by-side histograms comparing the prediction and training yield and standard deviation.

    Parameters:
    - mean (array-like): Predicted yield values.
    - training_mean (array-like): Training yield values.
    - std (array-like): Predicted standard deviation values.
    - training_std (array-like): Training standard deviation values.

    Returns:
    - None: Displays histograms comparing prediction and training yield and standard deviation.
    """
    # Plotting histograms side by side
    plt.figure(figsize=(10, 4))

    plt.subplot(2, 2, 1) 
    hist_range = [min(mean), max(mean)]
    plt.hist(mean, bins= 10, color='red', alpha=0.7, label='New points')
    plt.legend(prop={'size': 10})
    plt.title('Histogram of prediction yield')

    plt.subplot(2, 2, 3) 
    plt.hist(training_mean, bins=10, range=hist_range, color='black', alpha=0.5, label='Training points')
    plt.legend(prop={'size': 10})
    plt.title('Histogram of training yield')

    plt.subplot(2, 2, 2)  
    hist_range = [min(std), max(std)]
    plt.hist(std, bins=10, color='green', alpha=0.7, label='New points')
    plt.legend(prop={'size': 10})
    plt.title('Histogram of prediction std')

    plt.subplot(2, 2, 4)  
    plt.hist(training_std, bins=10, range=hist_range, color='black', alpha=0.5, label='Training points')
    plt.legend(prop={'size': 10})
    plt.title('Histogram of training std')
    plt.tight_layout()  # Adjust layout for better spacing
    plt.show()


def plot_hist_selected_yield(mean, selected_mean, title = 'Add title'):
    """
    Plots a histogram comparing the predicted yield and selected points' yield.

    Parameters:
    - mean (array-like): Predicted yield values.
    - selected_mean (array-like): Selected yield values to compare.
    - title (str, optional): Title for the plot. Default is 'Add title'.

    Returns:
    - None: Displays the histogram comparing predicted and selected yield.
    """
    hist_range = [min(mean), max(mean)]
    plt.figure(figsize=(6, 3))
    plt.hist(mean, bins=20, color='green', alpha=0.7, label='Prediction')
    plt.hist(selected_mean, bins=20, range=hist_range, color='red', alpha=0.5, label='Selected points')
    
    plt.title(title)
    plt.ylabel('Frequency')
    plt.legend(loc='upper left')
    plt.tight_layout()  
    plt.show()


def plot_distance(ax, X_new_norm,ucb_top_norm, ucb, ratio_ucb, color, title):
    """
    Plots the distance between selected and unselected points in a scatter plot.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes to plot on.
    - X_new_norm (array-like): The normalized new points.
    - ucb_top_norm (array-like): The normalized UCB top points.
    - ucb (array-like): UCB values corresponding to the points.
    - ratio_ucb (float): Ratio value for the selected UCB points.
    - color (str): Color for the selected points.
    - title (str): Title for the plot.

    Returns:
    - None: Displays the scatter plot of distances.
    """
    top_ucb = ucb_top_norm[0]
    distance_ucb = euclidean_distance(top_ucb, ucb_top_norm)
    distance_all_ucb = euclidean_distance(top_ucb,X_new_norm)

    ax.scatter(distance_all_ucb, ucb, color= 'grey', alpha = 0.5, label='Unselected points')
    ax.scatter(distance_ucb, ratio_ucb, color= color, alpha = 0.5, label='Selected points')
    ax.set_title(title)
    ax.legend(loc='upper right',  prop={'size': 10})


def plot_selected_point(ax, y_pred, std_pred, condition, title):
    """
    Plots selected and unselected points based on predicted values and standard deviations.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes to plot on.
    - y_pred (array-like): Predicted yield values.
    - std_pred (array-like): Predicted standard deviation values.
    - condition (array-like): The condition to select specific points.
    - title (str): Title for the plot.

    Returns:
    - None: Displays a scatter plot with selected and unselected points.
    """
    # Specify the X positions where you want different colors
    position = np.where(np.isin(y_pred, condition))[0]
    selected_std = std_pred[position]
    selected_y = y_pred[position]

    # not selected points
    y_not = np.delete(y_pred,position)
    std_not = np.delete(std_pred,position)

    # Create a scatter plot with different colors at special X positions
    ax.scatter(y_not, std_not, c='grey', label='Unselected points', alpha = 0.5)
    ax.scatter(selected_y, selected_std, c='red', alpha = 0.5, label='Selected_point')

    # Customize the plot
    ax.set_title(title)
    ax.set_xlabel('Predicted yield')
    ax.set_ylabel('Predicted std')
    ax.legend(loc='upper left',  prop={'size': 10})


def plot_heatmap(ax, X_new_norm, y_pred, element_list, title):
    """
    Plots a heatmap of normalized data, sorted by predicted values.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes to plot on.
    - X_new_norm (array-like): Normalized input data.
    - y_pred (array-like): Predicted values used for sorting.
    - element_list (list): List of element labels for the y-axis.
    - title (str): Title for the plot.

    Returns:
    - None: Displays the heatmap with sorted rows based on predicted values.
    """
    # Get the indices that would sort the array based on the values list
    n = len(y_pred)
    sorted_indices = np.argsort(y_pred)
    # Use the sorted indices to rearrange the rows of the array
    sorted_X = X_new_norm[sorted_indices]
    sorted_X = sorted_X.T

    # Create a heatmap
    heatmap=ax.imshow(sorted_X, aspect='auto', cmap='viridis')

    ax.set_xticks([])
    ax.set_title(f"{title}, sample = {n}", size = 12)
    ax.set_yticks(np.arange(len(element_list)), element_list, size = 12)
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.set_label('ratio with max concentrations', size = 12)

def plot_rmse(model):
    """
    Plots a histogram of the RMSE (Root Mean Square Error) values from the model's cross-validation scores.

    Parameters:
    - model (sklearn model): The model object with cross-validation results.

    Returns:
    - None: Displays a histogram of RMSE values during cross-validation.

    Note: RMSE (Root Mean Square Error) is a metric that measures the average magnitude 
    of prediction errors in a regression model. It is calculated by taking the square root 
    of the average of squared differences between predicted and actual values. 
    Lower RMSE values indicate better model performance, as they mean predictions are closer to actual values.
    """
    scores = model.cv_score
    score_list = scores[scores.rank_test_score == 1].iloc[:, 6:-3]
    score_list = np.array(score_list)*-1
    plt.hist(score_list[0], bins = 20, color='orange')
    plt.title(f'Histogram of RMSE during cross validation, mean = {round(np.mean(score_list),2)}', size = 12)


def plot_each_round(y,size_list, predict = False):
    """
    Visualizes the yield evolution through each active learning query using box plots.
    
    Parameters:
        y (array-like): The yield data for all rounds.
        size_list (list): A list of sizes representing the number of data points in each round.
        predict (bool, optional): If True, highlights the last round's box in silver. Default is False.
    
    Returns:
        None: Displays a boxplot showing the yield evolution.
    """
    # Calculate the cumulative sum of the size_list
    cumulative_sizes = np.cumsum(size_list)

    # Split the array into subarrays based on the cumulative sizes
    subarrays = np.split(y, cumulative_sizes)

    # Flatten each subarray
    flattened_arrays = [subarray.flatten() for subarray in subarrays]

    # Create a DataFrame 
    y_by_file = {}
    for i in range(len(size_list)):
        name = 'round ' + str(i)
        y_by_file[name] = flattened_arrays[i]

    y_by_file = pd.DataFrame.from_dict(y_by_file, orient='index').transpose()

    boxprops = dict(linewidth=0)
    medianprops = dict(linewidth=1, color='red', linestyle='dashed')
    plt.figure(figsize=(10, 5))
    ax = sns.boxplot(y_by_file, color = 'yellow', width=0.3, boxprops=boxprops, medianprops= medianprops)

    # Add markers for the maximum values
    #max_values = y_by_file.max()
    #for i, value in enumerate(max_values):
    #    ax.annotate(f'{value:.2f}', (i, value), textcoords="offset points", xytext=(0,5),
    #                ha='center', fontsize=8, color='black')

    # Add markers for the median values
    median_values = y_by_file.median()
    for i, value in enumerate(median_values):
        ax.annotate(f'{value:.2f}', (i, value), textcoords="offset points", xytext=(0,3),
                    ha='center', fontsize=8, color='black')

    if predict:
        # Get the last box patch
        last_box = ax.patches[-1]

        # Change the color of the last box
        last_box.set_facecolor('silver')

    plt.ylabel('Yield')
    plt.title('Yield evolution through each active learning query')
    # Show the plot
    plt.show()


def plot_train_test(train, test, element_list):
    """
    Plots histograms comparing training and test data distributions for each element in `element_list`.

    Parameters:
        train (array-like): The training data.
        test (array-like): The test data.
        element_list (list): List of column names to plot histograms for.

    Returns:
        None: Displays histograms for each element comparing training and test data.
    """
    test = pd.DataFrame(test, columns= element_list)
    train = pd.DataFrame(train, columns= element_list)

    # Plot histograms for each column in both DataFrames on the same figure
    no_element = len(element_list)
    nrows = int(np.sqrt(no_element))
    ncols = no_element//nrows
    _, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 6))

    for i, column in enumerate(element_list):
        row = i // ncols
        col = i % ncols
        hist_range = [min(train[column]), max(train[column])]
        axes[row, col].hist(train[column], alpha=0.3, label='Previous data', bins=10)
        axes[row, col].hist(test[column], range=hist_range, alpha=1, label='New data', bins=10)
        axes[row, col].set_title(column)
        axes[row, col].legend()
plt.tight_layout()
plt.show()