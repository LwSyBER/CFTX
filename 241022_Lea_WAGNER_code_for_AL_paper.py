#%% [protein+MG] IMPORTS 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string
import seaborn as sns
import re
from matplotlib.patches import Rectangle
from sklearn.feature_selection import mutual_info_regression
from sklearn.decomposition import PCA

plt.rcParams['font.size'] = plt.rcParamsDefault['font.size']

#%% FUNCTIONS

def get_plate_variables(Name_plate):
    """
    This function retrieves specific experimental data and metadata for a given plate 
    based on its name. It extracts information such as the corresponding Excel file,
    fluorescence data for GFP and malachite green aptamer, time step, number of replicates, 
    and control conditions.
    
    Input:
    - Name_plate (str): The name of the plate, which must be one of the valid names 
      ('plate 0' to 'plate 10').
    
    Output:
    - plate_info (dict): A dictionary containing the following keys:
        - 'Excel_file': The name of the Excel file associated with the plate.
        - 'DF_data_GFP': DataFrame of GFP fluorescence data for each well.
        - 'DF_data_MG': DataFrame of malachite green aptamer fluorescence data for each well.
        - 'Time_step': The time step between fluorescence measurements.
        - 'Name_plate': The name of the plate.
        - 'nbr_replicate': The number of replicates on the plate.
        - 'nbr_controls': The number of control wells on the plate.
        - 'name_ref_plus_MG': The reference condition for malachite green.
        - 'name_ref_plus_GFP': The reference condition for GFP.
    
    If the provided Name_plate is not valid, an error message is printed, and the function returns None.
    """
    plate_info = {}
    
    # Check if Name_plate is valid
    valid_plate_names = ['plate 0', 'plate 1', 'plate 2', 'plate 3', 'plate 4', 'plate 5', 'plate 6', 'plate 7', 'plate 8', 'plate 9', 'plate 10']
    if Name_plate not in valid_plate_names:
        print("Error: a proper Name_plate input is required.")
        return
    
    if Name_plate == 'plate 0':
        plate_info['Excel_file'] = '231107_ECHO_DNA_plate0_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 0'
        plate_info['nbr_replicate'] = 3
        plate_info['nbr_controls'] = 3
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"
        
    elif Name_plate == 'plate 1':
        plate_info['Excel_file'] = '231221_ECHO_DNA_plate1_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 1'
        plate_info['nbr_replicate'] = 3
        plate_info['nbr_controls'] = 10
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 2':
        plate_info['Excel_file'] = '240115_ECHO_DNA_plate2_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 2'
        plate_info['nbr_replicate'] = 3
        plate_info['nbr_controls'] = 10
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 3':
        plate_info['Excel_file'] = '240131_ECHO_DNA_plate3_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 3'
        plate_info['nbr_replicate'] = 3
        plate_info['nbr_controls'] = 10
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 4':
        plate_info['Excel_file'] = '240212_ECHO_DNA_plate4_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 4'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 11
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 5':
        plate_info['Excel_file'] = '240314_ECHO_DNA_plate5_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 5'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 11
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 6':
        plate_info['Excel_file'] = '240319_ECHO_DNA_plate6_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 6'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 11
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 7':
        plate_info['Excel_file'] = '240327_ECHO_DNA_plate7_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 7'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 6
        plate_info['name_ref_plus_MG'] = "P0_Jove+"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 8':
        plate_info['Excel_file'] = '240403_ECHO_DNA_plate8_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 8'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 12
        plate_info['name_ref_plus_MG'] = "P5_Buffer 21"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 9':
        plate_info['Excel_file'] = '240410_ECHO_DNA_plate9_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 9'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 12
        plate_info['name_ref_plus_MG'] = "P5_Buffer 21"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"

    elif Name_plate == 'plate 10':
        plate_info['Excel_file'] = '240417_ECHO_DNA_plate10_fluoMax_data.xlsm'
        plate_info['DF_data_GFP'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_GFP', header=None, index_col=0)
        plate_info['DF_data_MG'] = pd.read_excel(plate_info['Excel_file'], sheet_name='Data_MGapt', header=None, index_col=0)
        plate_info['Time_step'] = 6.116667
        plate_info['Name_plate'] = 'plate 10'
        plate_info['nbr_replicate'] = 6
        plate_info['nbr_controls'] = 12
        plate_info['name_ref_plus_MG'] = "P5_Buffer 21"
        plate_info['name_ref_plus_GFP'] = "P0_Jove+"
        
    return plate_info

def legend(Excel_file):
    """
    Generate a legend mapping each well to a strain name.

    Parameters:
        Excel_file (str): Path to the Excel file containing the layout and strains information, comes from get_plate_variables
    Returns:
        - dico_strains_wells: Dictionary mapping strain names to lists of corresponding wells.
    """
    
    # Legend: associate each well to a strain name
    layoutP = pd.read_excel(Excel_file, sheet_name='layoutW')

    # List of every well for each strain
    list_strains = list(pd.read_excel(Excel_file, sheet_name='strains')['Name'])
    dico_strains_wells = {strain: [] for strain in list_strains}

    # Find the wells corresponding to each strain
    alphabet = string.ascii_uppercase
    alphabet = alphabet[1:15]
    for j in range(len(dico_strains_wells)):
        for cells in layoutP:
            if layoutP[cells].to_string(index=False) == list(dico_strains_wells.keys())[j]:
                dico_strains_wells[list(dico_strains_wells.keys())[j]].append(cells)

    # Iterate through the dictionary and filter out 'nan' values from each list
    for key, value in dico_strains_wells.items():
        dico_strains_wells[key] = [item for item in value if not pd.isna(item)]
        
    return list_strains, dico_strains_wells

def get_initial_conditions(name_plate):
    """
    Get initial conditions for each buffer in a specific plate.

    Parameters:
        name_plate (str): Name of the plate, defined at the beginning of the code

    Returns:
        DataFrame: DataFrame containing initial conditions for each buffer in the plate.
        
    Raises:
        ValueError: If an invalid plate name is provided.
    """
    if name_plate == "plate 0":
        DF_sampling = pd.read_excel('231107_ECHO_DNA_plate0_fluoMax_data.xlsm', sheet_name='sampling')
        DF_sampling.set_index('Buffer', inplace=True)
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 1":
        DF_sampling = pd.read_excel('231221_ECHO_DNA_plate1_fluoMax_data.xlsm', sheet_name='sampling')        
        DF_sampling.set_index('Buffer', inplace=True)
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 2":
        DF_sampling = pd.read_excel('240115_ECHO_DNA_plate2_fluoMax_data.xlsm', sheet_name='sampling')
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 3":
        DF_sampling = pd.read_excel('240131_ECHO_DNA_plate3_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 4":
        DF_sampling = pd.read_excel('240212_ECHO_DNA_plate4_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 5":
        DF_sampling = pd.read_excel('240314_ECHO_DNA_plate5_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 6":
        DF_sampling = pd.read_excel('240319_ECHO_DNA_plate6_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 7":
        DF_sampling = pd.read_excel('240327_ECHO_DNA_plate7_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 8":
        DF_sampling = pd.read_excel('240403_ECHO_DNA_plate8_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 9":
        DF_sampling = pd.read_excel('240410_ECHO_DNA_plate9_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    elif name_plate == "plate 10":
        DF_sampling = pd.read_excel('240417_ECHO_DNA_plate10_fluoMax_data.xlsm', sheet_name='sampling')        
        nbr_bff_wo_triplicates = len(DF_sampling) // nbr_replicate
        DF_sampling = DF_sampling.head(nbr_bff_wo_triplicates)
        DF_sampling = DF_sampling.set_index(pd.Index(list_strains))
        DF_sampling = DF_sampling.iloc[:, :8]
        return DF_sampling
    else:
        raise ValueError("Invalid plate name. Please provide a valid plate name.")

def calc_median_std_per_buffer(DF_data_corr, dico_strains_wells):
    """
    Calculate the median and standard deviation for a buffer using the relevant fluoresence value of the wells corresponding to the considered buffer.
    
    Parameters:
        DF_data_corr (pd.Series): Series containing the relevant fluorescence point for each well (maximum for MG, final value for GFP).
        dico_strains_wells (dict): Dictionary containing buffer names as keys and lists of corresponding well names as values.
        
    Returns:
        pd.DataFrame: DataFrame containing the median and standard deviation of fluorescence for each buffer.
    """
    # Create an empty DataFrame to store the results
    DF_results_fluo = pd.DataFrame()

    # Calculate the median and standard deviation for each buffer
    medians = []
    std_devs = []
    for buffer, well_list in dico_strains_wells.items():
        buffer_data = DF_data_corr.loc[well_list]
        median_value = buffer_data.median()
        std_dev_value = buffer_data.std()
        medians.append(median_value)
        std_devs.append(std_dev_value)

    # Create a DataFrame with the medians and standard deviations
    DF_results_fluo.index = dico_strains_wells.keys()
    DF_results_fluo['Median'] = medians
    DF_results_fluo['Std Dev'] = std_devs

    return DF_results_fluo
    

def barplot_one_plate(DF_results_yield, fluo_or_yield, MG_or_GFP='MG'):
    """
    Plot bar charts of median values with standard deviation for a given plate.
    
    Parameters:
        DF_results_yield (DataFrame): DataFrame containing the median and standard deviation of yields or fluorescence.
        fluo_or_yield (str): Indicates whether to plot fluorescence ('fluo') or yield ('yield').
        MG_or_GFP (str, optional): Indicates whether the data is for MG or GFP. Default is 'MG'.
    
    Returns:
        None
    """
    # Initialize colors
    colors = []
    colors_sorted = []

    # Plot by buffer name
    plt.figure(figsize=(21, 6))
    if MG_or_GFP == 'MG':
        colors = ['red' if 'Jove' in buffer_name else 'chocolate' if name_ref_plus_MG in buffer_name else 'coral' for buffer_name in DF_results_yield.index]
    elif MG_or_GFP == 'GFP':
        colors = ['green' if 'Jove' in buffer_name else 'olivedrab' if name_ref_plus_GFP in buffer_name else 'yellowgreen' for buffer_name in DF_results_yield.index]
        
    plt.bar(DF_results_yield.index, DF_results_yield['Median'], yerr=DF_results_yield['Std Dev'], capsize=5, color=colors)
    plt.xlabel('Buffer')
    plt.ylabel('Median')

    # Set the plot title based on MG_or_GFP and fluo_or_yield
    if MG_or_GFP == 'MG':
        title_part = 'mRNA'
    elif MG_or_GFP == 'GFP':
        title_part = 'protein'
        
    plt.title(f'Maximum of {title_part} production {fluo_or_yield} per buffer for {Name_plate}')
    plt.xticks(rotation=90, ha='center')
    plt.show()
    
    # Plot ranked values
    DF_results_yield_sorted = DF_results_yield.sort_values(by='Median')
    if MG_or_GFP == 'MG':
        colors_sorted = ['red' if 'Jove' in buffer_name else 'chocolate' if name_ref_plus_MG in buffer_name else 'coral' for buffer_name in DF_results_yield_sorted.index]
    elif MG_or_GFP == 'GFP':
        colors_sorted = ['green' if 'Jove' in buffer_name else 'olivedrab' if name_ref_plus_GFP in buffer_name else 'yellowgreen' for buffer_name in DF_results_yield_sorted.index]
    
    plt.figure(figsize=(21, 6))
    plt.bar(DF_results_yield_sorted.index, DF_results_yield_sorted['Median'], yerr=DF_results_yield_sorted['Std Dev'], capsize=5, color=colors_sorted)
    plt.xlabel('Buffer')
    plt.ylabel('Median')
    
    plt.title(f'Maximum {title_part} production {fluo_or_yield} per buffer for {Name_plate}')
    plt.xticks(rotation=90, ha='center')
    plt.show()
  
def list_strains_controls(Name_plate):
    """
    Generate a list of strains or controls to trace for a given plate.
    
    Parameters:
        Name_plate (str): The name of the plate (e.g., 'plate 5').
    
    Returns:
        list: List of strains or controls to trace.
    """
    if Name_plate == 'plate 5':
        list_strains_to_trace = ['P0_Buffer 6','P0_Jove-','P0_Buffer 68','P0_H20','P0_Buffer 18','P0_Jove+','P1_Buffer 44','P1_Buffer 25','P2_Buffer 38','P3_Buffer 9']
    else :
        list_strains_to_trace = list(dico_strains_wells.keys())[-nbr_controls:]
    return list_strains_to_trace


def plot_control_buffers(dico_strains_wells, DF_results_fluo, nbr_controls, name_ref_plus, list_strains_to_trace, Name_plate, MG_or_GFP='MG'):
    """
    Plot bar charts of median fluorescence for control buffers.
    
    Parameters:
        dico_strains_wells (dict): Dictionary containing the mapping of strains to wells.
        DF_results_fluo (DataFrame): DataFrame containing the median and standard deviation of fluorescence results.
        nbr_controls (int): Number of control buffers, comes from get_plate_variables.
        name_ref_plus (str): Name of the positive control buffer, comes from get_plate_variables.
        list_strains_to_trace (list): List of strains or buffers to trace, comes from list_strains_controls.
        Name_plate (str): The name of the plate (e.g., 'plate 3', 'plate 2'), defined at the beginning of the code.
        MG_or_GFP (str, optional): Indicates whether the data is for MG or GFP. Default is 'MG'.
    
    Returns:
        None
    """
    # Copy and reindex DataFrame based on the strains to trace
    DF_results_fluo_MG_ctrl = DF_results_fluo.copy()
    DF_results_fluo_MG_ctrl = DF_results_fluo_MG_ctrl.reindex(list_strains_to_trace)

    # Specific handling for 'plate 3'
    if Name_plate == 'plate 3':
        DF_results_fluo_MG_ctrl = DF_results_fluo_MG_ctrl.drop('P1_Buffer 38', axis=0)  # SPECIAL P3 (control poorly chosen)

    # Define title and colors based on MG_or_GFP
    if MG_or_GFP == 'MG':
        mid_title = 'mRNA'
        colors = ['red' if 'Jove' in buffer_name else 'chocolate' if name_ref_plus in buffer_name else 'coral' for buffer_name in DF_results_fluo_MG_ctrl.index]
    elif MG_or_GFP == 'GFP':
        mid_title = 'protein'
        colors = ['green' if 'Jove' in buffer_name else 'olivedrab' if name_ref_plus in buffer_name else 'yellowgreen' for buffer_name in DF_results_fluo_MG_ctrl.index]

    # Plot each control buffer
    for index_pos in [i for i, val in enumerate(DF_results_fluo_MG_ctrl.index) if val in list_strains_to_trace]:
        buffer_name = DF_results_fluo_MG_ctrl.index[index_pos]
        plt.bar(buffer_name, DF_results_fluo_MG_ctrl['Median'][index_pos], yerr=DF_results_fluo_MG_ctrl['Std Dev'][index_pos], capsize=5, color=colors[index_pos])

    plt.xlabel('Control buffers')
    plt.ylabel('Median mRNA production fluo (n=3)')
    
    # Set plot title
    plt.title(f'Maximum {mid_title} production per control buffer (n=3) for {Name_plate}')
    plt.xticks(rotation=90, ha='center')
    plt.show()

 
def calc_yield_NaN(dico_strains_wells, DF_data_corr, DF_results_fluo, name_ref_plus):
    """
    Calculate the yield for each buffer by normalizing the fluorescence data and handling NaN values.
    
    Parameters:
        dico_strains_wells (dict): Dictionary containing the mapping of strains to wells.
        DF_data_corr (DataFrame): DataFrame containing corrected data.
        DF_results_fluo (DataFrame): DataFrame containing the median and standard deviation of fluorescence results.
        name_ref_plus (str): Reference name to identify specific buffers for normalization.

    Returns:
        Tuple[DataFrame, DataFrame]: A tuple containing:
            - DF_results_yield: DataFrame with median and standard deviation of yields.
            - DF_yield_WperW: DataFrame with yield values per well.
    """
    # Create empty DataFrames to store the results for raw and adjusted data
    DF_results_yield = pd.DataFrame()
    DF_yield_WperW = pd.DataFrame(index=dico_strains_wells.keys(), columns=[f'yield_{i}' for i in range(1, 7)])

    # Loop over each buffer and calculate yield for raw data
    for buffer, well_list in dico_strains_wells.items():
        buffer_data = DF_data_corr.loc[well_list]

        # Calculations for raw data
        Ref_minus_median = DF_results_fluo.loc["P0_Jove-", "Median"]
        Ref_plus_median = DF_results_fluo.loc[name_ref_plus, "Median"]

        yields = (buffer_data - Ref_minus_median) / (Ref_plus_median - Ref_minus_median)
        median_yield = np.nanmedian(yields)
        std_dev_yield = np.nanstd(yields)

        nan_array = np.full(DF_yield_WperW.shape[1], np.nan)
        nan_array[:len(yields)] = yields

        DF_yield_WperW.loc[buffer] = nan_array
        DF_yield_WperW = DF_yield_WperW.apply(pd.to_numeric, errors='coerce')

        # Store results in DF_results_yield
        DF_results_yield.loc[buffer, 'Median'] = median_yield
        DF_results_yield.loc[buffer, 'Std Dev'] = std_dev_yield
    
    return DF_results_yield, DF_yield_WperW


def correct_gain_change(name_ref_plus, name_plate, DF_results_yield, DF_yield_WperW, MG_or_GFP='MG'):
    """
    Correct the yield values for a given plate by adjusting for gain changes.
    
    Parameters:
        name_ref_plus (str): Reference name for the buffer used for normalization.
        name_plate (str): The name of the plate being processed.
        DF_results_yield (DataFrame): DataFrame containing the median and standard deviation of yield results.
        DF_yield_WperW (DataFrame): DataFrame with yield values per well.
        MG_or_GFP (str, optional): Type of data ('MG' or 'GFP'). Default is 'MG'.
    
    Returns:
        Tuple[DataFrame, DataFrame, float]: A tuple containing:
            - DF_results_yield: Corrected DataFrame with median and standard deviation of yields.
            - DF_yield_WperW: Corrected DataFrame with yield values per well.
            - Ratio_Refplus_Joveplus: Calculated ratio used for correction.
    """
    # Adapt name of Name_plate_compil for the compilated files
    Last_plate = name_plate
    all_Name_plate_compil = 'P0_P1_P2_P3_P4_P5_P6_P7_P8_P9_P10'
    last_part = Last_plate.split()[-1]
    number = int(last_part.split()[-1]) - 1
    last_part = str(number)
    
    if number < 10:
        Name_plate_compil = all_Name_plate_compil[:all_Name_plate_compil.index(last_part) + len(last_part)]
    elif number > 10 and last_part in all_Name_plate_compil:
        Name_plate_compil = all_Name_plate_compil[:all_Name_plate_compil.index(last_part) + len(last_part)]
    else:
        Name_plate_compil = all_Name_plate_compil[:all_Name_plate_compil.index(last_part) + len(last_part)]
    
    # Calculation of the ratio
    DF_gain100_all_data = pd.read_csv(Name_plate_compil + "_AL_corr_everything_MG_GFP.csv")
    DF_gain_change = DF_gain100_all_data[(DF_gain100_all_data['Buffer'] == name_ref_plus) | (DF_gain100_all_data['Buffer'] == 'P0_Jove+')]
    DF_gain_change_piv = pd.pivot_table(
        DF_gain_change,
        values=[f'{MG_or_GFP}_yield_median'],
        index=['name_plate'],
        columns=['Buffer'],
        aggfunc='first'
    )
    DF_gain_change_piv['Ref_plus/P0_Jove+'] = DF_gain_change_piv[f'{MG_or_GFP}_yield_median', name_ref_plus] / DF_gain_change_piv[f'{MG_or_GFP}_yield_median', 'P0_Jove+']
    Ratio_Refplus_Joveplus = DF_gain_change_piv.loc['plate 7', (f'{MG_or_GFP}_yield_median', name_ref_plus)]

    # Open file for writing output
    with open(f'correct_gain_change_params_{MG_or_GFP}_{Name_plate}.txt', 'w') as file:
        file.write(DF_gain_change_piv.to_string())
        file.write(f"\nRatio_Refplus_Joveplus = {Ratio_Refplus_Joveplus}\n")

    # Multiply yield (median and std) of the plate with lower gain by this ratio
    DF_results_yield = DF_results_yield * Ratio_Refplus_Joveplus
    Ratio_Refplus_Joveplus_newP = DF_results_yield.loc[name_ref_plus, 'Median'] / DF_results_yield.loc['P0_Jove+', 'Median']
    
    # Append new results to the file
    with open(f'correct_gain_change_params_{MG_or_GFP}_{Name_plate}.txt', 'a') as file:
        file.write(f"Ratio_Refplus_Joveplus_newP = {Ratio_Refplus_Joveplus_newP} in {Name_plate}\n")

    # Multiply yield (_1,_2,_3,_4,_5,_6) of the plate with lower gain by this ratio
    DF_yield_WperW = DF_yield_WperW * Ratio_Refplus_Joveplus

    return DF_results_yield, DF_yield_WperW, Ratio_Refplus_Joveplus



def process_fluo_DF_data_corr(dico_strains_wells, DF_data, MG_or_GFP):
    """
    Process fluorescence data and create DataFrames with structured values and well information.
    
    Parameters:
        dico_strains_wells (dict): Dictionary mapping buffer names to lists of wells.
        DF_data (DataFrame): DataFrame containing the fluorescence data for each well.
        MG_or_GFP (str): Prefix to use for column names in the resulting DataFrames, either 'MG' or 'GFP'.
    
    Returns:
        Tuple[DataFrame, DataFrame]: A tuple containing:
            - DF_fluo_values: DataFrame with structured fluorescence values for each buffer.
            - DF_wells: DataFrame with well identifiers for each buffer.
    """
    list_dico_fluo_values = []
    list_dico_wells = []
    
    max_wells = max(len(wells) for wells in dico_strains_wells.values())
    
    for buffer, well_list in dico_strains_wells.items():
        values = DF_data.loc[well_list].tolist()

        # Pad values with NaNs if they are shorter than max_wells
        values += [np.nan] * (max_wells - len(values))
        dico_values = {f'{MG_or_GFP}_fluo_{i+1}': values[i] for i in range(max_wells)}
        dico_wells = {f'well_{i+1}': well_list[i] for i in range(len(well_list))}
        list_dico_fluo_values.append(dico_values)
        list_dico_wells.append(dico_wells)

    DF_fluo_values = pd.DataFrame(list_dico_fluo_values)
    DF_wells = pd.DataFrame(list_dico_wells)
    
    DF_fluo_values.index = dico_strains_wells.keys()
    DF_wells.index = dico_strains_wells.keys()
    
    return DF_fluo_values, DF_wells


def process_fluo_DF_results_fluo(DF_results, MG_or_GFP):
    """
    Rename columns in the fluorescence results DataFrame for consistency with a given prefix 'fluo',
    in order to be compiled within plate_AL_X_corr_everything_MG (DataFrame).
    
    Parameters:
        DF_results (DataFrame): DataFrame containing the fluorescence results with 'Median' and 'Std Dev' columns.
        MG_or_GFP (str): Prefix to use for column names, either 'MG' or 'GFP'.
    
    Returns:
        DataFrame: DataFrame with renamed columns for median and standard deviation.
    """
    return DF_results.rename(columns={'Median': f'{MG_or_GFP}_fluo_median', 'Std Dev': f'{MG_or_GFP}_fluo_std'})

def process_yield_DF_WperW(DF_yield, MG_or_GFP):
    """
    Rename columns in the yield WperW DataFrame for consistency with a given prefix 'yield',
    in order to be compiled within plate_AL_X_corr_everything_MG (DataFrame).
    
    Parameters:
        DF_yield (DataFrame): DataFrame containing yield data for each well.
        MG_or_GFP (str): Prefix to use for column names, either 'MG' or 'GFP'.
    
    Returns:
        DataFrame: DataFrame with renamed columns for yield values.
    """
    return DF_yield.rename(columns={col: f'{MG_or_GFP}_yield_{col.split("_")[-1]}' for col in DF_yield.columns})

def process_yield_DF_results_yield(DF_results, MG_or_GFP):
    """
    Rename columns in the yield results DataFrame for consistency with a given prefix 'yield',
    in order to be compiled within plate_AL_X_corr_everything_MG (DataFrame).
    
    Parameters:
        DF_results (DataFrame): DataFrame containing the yield results with 'Median' and 'Std Dev' columns.
        MG_or_GFP (str): Prefix to use for column names, either 'MG' or 'GFP'.
    
    Returns:
        DataFrame: DataFrame with renamed columns for median and standard deviation.
    """
    return DF_results.rename(columns={'Median': f'{MG_or_GFP}_yield_median', 'Std Dev': f'{MG_or_GFP}_yield_std'})


def subplot_controls_across_plates(DF_all_data_MG_GFP, Name_plate, MG_or_GFP):
    
    if Name_plate == 'plate 0':
        print('No control comparison for plate 0')
        
    else:
        # Filter out row where 'Buffer' is not duplicated to remove all those that have been tested only once
        filtered_df = DF_all_data_MG_GFP[DF_all_data_MG_GFP['Buffer'].duplicated(keep=False)]
        
        # Plot
        unique_buffers = filtered_df['Buffer'].unique()
        
        if MG_or_GFP == 'MG':
            columns_to_plot = ['MG_fluo_median', 'MG_yield_median'] # List of columns to plot
            my_col = 'coral'
            tit_mid = 'mRNA'
        elif MG_or_GFP == 'GFP':
            columns_to_plot = ['GFP_fluo_median', 'GFP_yield_median'] # List of columns to plot
            my_col = 'yellowgreen'
            tit_mid = 'protein'
            
        # Iterate over each column to create separate figures
        for col in columns_to_plot:
            # Calculate the number of rows and columns for the subplots
            num_rows = int(np.ceil(len(unique_buffers) / 4))  # Adjusted to accommodate four columns
            num_cols = min(len(unique_buffers), 4)            # Adjusted to ensure enough space
        
            # Calculate the figure size to make it square
            figsize = (15, 3 * num_rows)
        
            # Create subplots
            fig, axs = plt.subplots(num_rows, num_cols, figsize=figsize, sharex=False, sharey=True)
        
            # Flatten axs if it's not already flat
            axs = np.array(axs).flatten()
        
            # Iterate over each buffer
            for j, buffer_name in enumerate(unique_buffers):
                # Filter the DataFrame for rows 'Buffer' is the current buffer
                filtered_data = filtered_df[filtered_df['Buffer'] == buffer_name]
        
                # Define width of each bar
                bar_width = 0.35
        
                # Set x coordinates for bars
                x_mg = np.arange(len(filtered_data['name_plate']))
        
                # Plotting with error bars for the current column
                axs[j].bar(x_mg, filtered_data[col], width=bar_width,
                           yerr=filtered_data[col.replace('median', 'std')],
                           capsize=5, color=my_col, label=col)
        
                # Set labels and title
                axs[j].set_title(f'{buffer_name} - {col}')
                axs[j].set_xticks(np.arange(len(filtered_data['name_plate'])))
                axs[j].set_xticklabels(filtered_data['name_plate'])
                axs[j].tick_params(axis='x', rotation=45)
                #axs[j].set_ylim(0,4)
        
                # Add grid
                # axs[j].grid(True, zorder=0)
        
            # Adjust layout
            plt.tight_layout()
            fig.suptitle(f'{col} - Yield and fluorescence of {tit_mid} production across plates per buffer', fontsize=18, y=1.02)
        
            plt.show()
            
        return(unique_buffers)

def get_label_color(label_text):
    """
    Determine the color for a label based on its text,
    in order to make the legend of the heatmap colored by name_plate.

    Parameters:
        label_text (str): Text of the label.

    Returns:
        str: Color for the label.
    """
    # Default color for other labels
    label_color = 'black'
    tab10_palette = sns.color_palette('tab10')
    # Iterate through plate numbers and set label color accordingly
    for i, ending in enumerate(['in plate 0', 'in plate 1', 'in plate 2', 'in plate 3', 'in plate 4', 'in plate 5', 'in plate 6', 'in plate 7', 'in plate 8', 'in plate 9', 'in plate 10']):
        if label_text.endswith(ending):
            # Set label color using tab10_palette with modulo operation
            label_color = tab10_palette[(i) % len(tab10_palette)]
            break  # Exit loop once color is set
    return label_color

def plot_barplot_heatmap_compil(DF_all_data, name_ref_plus, MG_or_GFP='MG', figsize_x=150, figsize_y=15):
    """
    Plot a bar plot and a heatmap side by side.
    
    Parameters:
        MG_or_GFP (str, optional): Prefix to use for colors and title, either 'MG' or 'GFP'. Default is 'MG'.
        figsize_x (int, optional): Width of the figure in inches. Default is 150.
        figsize_y (int, optional): Height of the figure in inches. Default is 15.
    
    Returns:
        None
    """
    # Create subplots with 2 rows and 1 column
    fig, axs = plt.subplots(2, 1, figsize=(figsize_x, figsize_y))
    
    # Set a suptitle for the whole figure
    if MG_or_GFP == 'MG':
        mid_title = 'mRNA'
    elif MG_or_GFP == 'GFP':
        mid_title = 'protein'
    else:
        raise ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")
    
    plt.suptitle(f'Maximum {mid_title} production yield per buffer for {Name_plate_compil}', fontsize=20)
    
    # First subplot for the bar plot
    ax_barplot = axs[0]

    # Data for barplot
    DF_barplot = DF_all_data.copy()
    DF_barplot_sorted = DF_barplot.sort_values(by=f'{MG_or_GFP}_yield_median')

    # Define colors based on buffer names
    pattern = re.compile(re.escape(name_ref_plus) + r'\s')
    colors = ['red' if 'Jove' in buffer_name else 'chocolate' if pattern.search(buffer_name) else 'coral' for buffer_name in DF_barplot_sorted.index] if MG_or_GFP == 'MG' else ['green' if 'Jove' in buffer_name else 'olivedrab' if pattern.search(buffer_name) else 'yellowgreen' for buffer_name in DF_barplot_sorted.index]

    # Create the bar plot
    ax_barplot.bar(DF_barplot_sorted.index, DF_barplot_sorted[f'{MG_or_GFP}_yield_median'], yerr=DF_barplot_sorted[f'{MG_or_GFP}_yield_std'], capsize=5, color=colors)
    ax_barplot.set_ylabel('Median')
    ax_barplot.set_xlim(-0.5, len(DF_barplot_sorted.index) - 0.5)
    ax_barplot.set_xticks([])
    ax_barplot.grid(True, color='white')  # Set grid color to white
    ax_barplot.yaxis.set_label_position("right")
    ax_barplot.yaxis.tick_right()

    # Second subplot for the heatmap
    ax_heatmap = axs[1]

    # Data for heatmap
    Excel_file_param_max_across_plates = 'Params_concentration_range.csv'  # The final one after plate 10
    DF_data_param = pd.read_csv(Excel_file_param_max_across_plates)
    list_max_conc = DF_data_param['maxValue'].tolist()
    list_igdt = DF_data_param['Parameter'].tolist()

    DF_heatmap = DF_all_data.copy()
    # Normalize each column by its corresponding maximum concentration
    for col, max_concentration in zip(list_igdt, list_max_conc):
        DF_heatmap[col] = (DF_heatmap[col] / max_concentration) * 100

    DF_heatmap = DF_heatmap.T
    DF_heatmap_sorted = DF_heatmap.sort_values(by=f'{MG_or_GFP}_yield_median', axis=1)
    DF_heatmap_sorted_ = DF_heatmap_sorted.iloc[:8].apply(pd.to_numeric, errors='coerce')

    # Create the heatmap
    heatmap_plot = sns.heatmap(DF_heatmap_sorted_, cmap='viridis', linewidths=0.5, linecolor='grey', annot=False, cbar=None, robust=True, ax=ax_heatmap)

    # x-axis labels
    xtick_labels = heatmap_plot.get_xticklabels()
    for label in xtick_labels:
        label.set_color(get_label_color(label.get_text()))
    heatmap_plot.set_xticklabels(xtick_labels, ha='center')

    # y-axis ticks
    ax_heatmap.yaxis.tick_right()
    # y-axis labels
    ytick_labels = heatmap_plot.get_yticklabels()
    for label in ytick_labels:
        label.set_color(get_label_color(label.get_text()))
        label.set_rotation(0)  # Set rotation to 0 degrees for horizontal orientation
    heatmap_plot.set_yticklabels(ytick_labels, ha='left')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the plot
    plt.show()
    
def filter_remove_controls(DF):
    # Extract numbers from Buffer and name_plate
    DF['buffer_number'] = DF['Buffer'].str.extract(r'P(\d+)_')[0].astype(int)  # Extracting number between "P" and "_"
    DF['plate_number'] = DF['name_plate'].str.extract(r'(\d+)')[0].astype(int)  # Extracting number from name_plate
    
    # Keep only the rows where buffer_number matches plate_number
    filtered_df = DF[DF['buffer_number'] == DF['plate_number']]
    
    # Drop the temporary columns used for filtering
    filtered_df = filtered_df.drop(columns=['buffer_number', 'plate_number'])
    
    return filtered_df

    
def plot_barplot_heatmap_compil_plate_per_plate(DF_all_data, name_ref_plus, MG_or_GFP='MG', figsize_x=40, figsize_y=15):
    """
    Plot a bar plot and a heatmap side by side.

    Parameters:
        MG_or_GFP (str, optional): Prefix to use for colors and title, either 'MG' or 'GFP'. Default is 'MG'.
        figsize_x (int, optional): Width of the figure in inches. Default is 40.
        figsize_y (int, optional): Height of the figure in inches. Default is 15.
    
    Returns:
        None
    """
    # Create subplots with 2 rows and 1 column
    fig, axs = plt.subplots(2, 1, figsize=(figsize_x, figsize_y))
    
    # Set a suptitle for the whole figure
    if MG_or_GFP == 'MG':
        mid_title = 'mRNA'
    elif MG_or_GFP == 'GFP':
        mid_title = 'protein'
    else:
        raise ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")

    plt.suptitle(f'Maximum {mid_title} production yield per buffer for {Name_plate_compil}', fontsize=20)

    # First subplot for the bar plot
    ax_barplot = axs[0]

    # Data for barplot
    plates = [f'plate {i}' for i in range(11)]  # List of plates
    DF_complete_sorted_data = pd.DataFrame()  # Initialize an empty DataFrame
    DF_all_data = filter_remove_controls(DF_all_data) # Remove controls

    for plate in plates:
        # Filter data for the current plate
        plate_data = DF_all_data[DF_all_data['name_plate'] == plate]
        # Sort the data by 'yield_median'
        #plate_data_sorted = plate_data.sort_values(by=f'{MG_or_GFP}_yield_median')
        plate_data_sorted = plate_data.sort_values(by='MG_yield_median')

        # Concatenate the sorted data
        DF_complete_sorted_data = pd.concat([DF_complete_sorted_data, plate_data_sorted])
 
    DF_barplot_sorted = DF_complete_sorted_data.copy()

    # Define colors based on buffer names
    colors = ['coral' for buffer_name in DF_barplot_sorted.index] if MG_or_GFP == 'MG' else ['#3cb371' for buffer_name in DF_barplot_sorted.index]

    # Create the bar plot
    ax_barplot.bar(DF_barplot_sorted.index, DF_barplot_sorted[f'{MG_or_GFP}_yield_median'], yerr=DF_barplot_sorted[f'{MG_or_GFP}_yield_std'], capsize=5, color=colors)
    ax_barplot.set_ylabel('Median')
    ax_barplot.set_xlim(-0.5, len(DF_barplot_sorted.index) - 0.5)
    ax_barplot.set_xticks([])
    ax_barplot.grid(True, color='white')
    ax_barplot.yaxis.set_label_position("right")
    ax_barplot.yaxis.tick_right()

    # Second subplot for the heatmap
    ax_heatmap = axs[1]

    # Data for heatmap
    Excel_file_param_max_across_plates = 'Params_concentration_range.csv'  # Final file after plate 10
    DF_data_param = pd.read_csv(Excel_file_param_max_across_plates)
    list_max_conc = DF_data_param['maxValue'].tolist()
    list_igdt = DF_data_param['Parameter'].tolist()

    DF_heatmap = DF_complete_sorted_data.copy()
    # Normalize each column by its corresponding maximum concentration
    for col, max_concentration in zip(list_igdt, list_max_conc):
        DF_heatmap[col] = (DF_heatmap[col] / max_concentration) * 100

    DF_heatmap = DF_heatmap.T
    DF_heatmap_sorted_ = DF_heatmap.iloc[:8].apply(pd.to_numeric, errors='coerce')

    # Create the heatmap
    heatmap_plot = sns.heatmap(DF_heatmap_sorted_, cmap='viridis', linewidths=0.5, linecolor='grey', annot=False, cbar=None, robust=True, ax=ax_heatmap)

    # Hide x-ticks on the heatmap
    ax_heatmap.set_xticks([])

    # y-axis ticks
    ax_heatmap.yaxis.tick_right()
    # y-axis labels
    ytick_labels = heatmap_plot.get_yticklabels()
    for label in ytick_labels:
        label.set_color(get_label_color(label.get_text()))
        label.set_rotation(0)
    heatmap_plot.set_yticklabels(ytick_labels, ha='left')

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Show the plot
    plt.show()

    
def boxplot_yield_across_plates(MG_or_GFP='MG'):
    """
    Plot a boxplot of yield across plates with swarmplot overlay.

    Parameters:
        MG_or_GFP (str, optional): Prefix to use for the title, either 'MG' or 'GFP'. Default is 'MG'.

    Returns:
        None
    """
    
    # Create a color dictionary for unique buffers
    color_dict = {buffer_name: color for buffer_name, color in zip(unique_buffers, sns.color_palette('Paired', len(unique_buffers)))}
    
    # Create a boxplot using seaborn with tab10 colors
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=DF_all_data_MG_GFP, ax=ax, palette='tab10')
    
    # Plotting the median points on top of the boxes in black
    sns.swarmplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=DF_all_data_MG_GFP, color='black', size=3, ax=ax)
    
    # Plotting the points associated with unique_buffers with unique colors
    for buffer_name in unique_buffers:
        buffer_data = DF_all_data_MG_GFP[DF_all_data_MG_GFP['Buffer'] == buffer_name]
        sns.swarmplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=buffer_data, color=color_dict[buffer_name], size=3, ax=ax)
    
    # Set labels and title
    ax.set_xlabel('Plate')
    ax.set_ylabel('Average yield')
    
    mid_title = 'mRNA' if MG_or_GFP == 'MG' else 'protein' if MG_or_GFP == 'GFP' else ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")
    ax.set_title(f'Yield for maximum {mid_title} production across plates')
    
    # Add custom legend annotations
    legend_height = 0.95
    for buffer_name, color in color_dict.items():
        rect = Rectangle((1.01, legend_height - 0.02), 0.05, 0.03, color=color, transform=ax.transAxes, clip_on=False)
        ax.add_patch(rect)
        ax.text(1.08, legend_height, buffer_name, color='black', fontsize=8, verticalalignment='center', transform=ax.transAxes)
        legend_height -= 0.05
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_xlabel('')  # Delete the legend axis name 'Plate'

    # Adjust the layout to ensure everything fits and is centered
    plt.tight_layout()

    # Render the plot
    plt.draw()
    
    # Show the plot
    plt.show()
    
    
def boxplot_yield_across_plates_sleek(MG_or_GFP='MG'):
    """
    Plot a boxplot of yield across plates with swarmplot overlay.

    Parameters:
        MG_or_GFP (str, optional): Prefix to use for the title, either 'MG' or 'GFP'. Default is 'MG'.

    Returns:
        None
    """
    
    # Define color and title based on input parameters
    my_col, tit_mid = ('coral', 'mRNA') if MG_or_GFP == 'MG' else ('yellowgreen', 'protein') if MG_or_GFP == 'GFP' else (None, None)
    if my_col is None:
        raise ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")

    # Create a boxplot
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=DF_all_data_MG_GFP_sleek, ax=ax, color=my_col)
    
    # Plotting the median points on top of the boxes in black
    sns.swarmplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=DF_all_data_MG_GFP_sleek, color='black', size=3, ax=ax)
    
    # Plotting the points associated with unique_buffers with unique colors
    for buffer_name in unique_buffers:
        buffer_data = DF_all_data_MG_GFP_sleek[DF_all_data_MG_GFP_sleek['Buffer'] == buffer_name]
        sns.swarmplot(x='name_plate', y=f'{MG_or_GFP}_yield_median', data=buffer_data, color='cyan', size=3, ax=ax)
    
    # Set labels and title
    ax.set_xlabel('Plate')
    ax.set_ylabel('Average yield')
    ax.set_title(f'Yield for maximum {tit_mid} production across plates')

    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='center')
    ax.set_xlabel('')  # Delete the legend axis name 'Plate'
    
    # Adjust the layout to ensure everything fits and is centered
    plt.tight_layout()
    plt.draw()
    plt.show()

    
def plot_barplot_all_data_plate_by_plate(DF_all_data_MG_GFP, MG_or_GFP):
    """
    Plot bar plots of median yield for each plate, grouped by buffer name.

    Parameters:
        DF_all_data_MG_GFP (DataFrame): DataFrame containing all data.
        MG_or_GFP (str, optional): Prefix to use for the title, either 'MG' or 'GFP'.

    Returns:
        None
    """
    # Define titles based on input parameters
    tit_mid = 'mRNA' if MG_or_GFP == 'MG' else 'protein' if MG_or_GFP == 'GFP' else ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")
    # List of unique plate names
    plates = [f'plate {i}' for i in range(11)]
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(40, 15))
    
    # Position tracker
    pos = 0
    # Plotting data for each plate
    for plate in plates:
        # Filter data for the current plate
        plate_data = DF_all_data_MG_GFP[DF_all_data_MG_GFP['name_plate'] == plate]
        
        # Sort the data by 'yield_median'
        plate_data_sorted = plate_data.sort_values(by=f'{MG_or_GFP}_yield_median')
        
        if MG_or_GFP == 'MG':
            colors = ['red' if 'Jove' in buffer_name else 'chocolate' if name_ref_plus_MG in buffer_name else 'coral' for buffer_name in plate_data_sorted.index]
        elif MG_or_GFP == 'GFP':
            colors = ['green' if 'Jove' in buffer_name else 'olivedrab' if name_ref_plus_GFP in buffer_name else 'yellowgreen' for buffer_name in plate_data_sorted.index]
        
        # Plotting the bar for the current plate
        bars = ax.bar([i + pos for i in range(len(plate_data_sorted))], plate_data_sorted[f'{MG_or_GFP}_yield_median'], yerr=plate_data_sorted[f'{MG_or_GFP}_yield_std'], capsize=5, color=colors)
        
        # Print plate name below the middle bar
        if plate_data_sorted.shape[0] > 0:
            middle_bar_index = len(bars) // 2  # Find the index of the middle bar
            ax.text(pos + middle_bar_index, -0.6, plate, ha='center', va='top', fontsize=30, color='black')
    
        # Update position for the next plate
        pos += len(plate_data_sorted) + 2  # Adding some space between plates
    
    # Remove x-axis ticks
    ax.set_xticks([])
    # Set y-axis label
    ax.set_ylabel('Median')
    # Set plot title
    ax.set_title(f'Yield for maximum {tit_mid} production for all plates')
    # Hide legend
    ax.legend().remove()
    # Move y-axis label to the right
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    # Adjust x-axis limits to reduce blank spaces on the sides
    plt.xlim(-0.5, pos - 0.5)
    
    plt.show()

def plot_correlation_by_plate_labels(DF_all_data_MG_GFP, MG_or_GFP_x, MG_or_GFP_y, ort='', selected_buffers=None):
    """
    Plot correlation scatter plot across plates.
    
    Parameters:
    - DF_all_data_MG_GFP (DataFrame): DataFrame containing all data.
    - MG_or_GFP_x (str): Prefix for x-axis, either 'MG' or 'GFP'.
    - MG_or_GFP_y (str): Prefix for y-axis, either 'MG' or 'GFP'.
    - ort (str, optional): Whether to enforce an orthogonal plot with same scale on both axes. Default is ''.
    
    Returns:
    None
    """
    # Define titles based on input parameters
    tit_mid_x = 'mRNA' if MG_or_GFP_x == 'MG' else 'protein'
    tit_mid_y = 'mRNA' if MG_or_GFP_y == 'MG' else 'protein'
    
    deb_tit_x = 'maximum'
    deb_tit_y = 'maximum'
    
    # Number and list of unique plates
    list_name_plates = DF_all_data_MG_GFP['name_plate'].unique()
    nbr_names_plates = len(list_name_plates)
    
    # Create a color dictionary mapping each plate name to a color from the viridis palette
    large_palette = sns.color_palette("viridis", n_colors=nbr_names_plates)
    color_dict = dict(zip(list_name_plates, large_palette))
    
    plt.figure(figsize=(8, 6))
    
    # Plot data for each plate with transparency
    for name_plate in list_name_plates:
        plate_data_x = DF_all_data_MG_GFP[DF_all_data_MG_GFP['name_plate'] == name_plate][f'{MG_or_GFP_x}_yield_median']
        plate_data_y = DF_all_data_MG_GFP[DF_all_data_MG_GFP['name_plate'] == name_plate][f'{MG_or_GFP_y}_yield_median']
        
        plt.scatter(plate_data_x, plate_data_y, 
                    color=color_dict.get(name_plate, 'black'), 
                    s=50, alpha=0.3)  # Default points
    
    # Highlight selected buffers with custom colors, transparency, and no legend entry
    if selected_buffers is not None:
        highlight_colors = {
            "P0_Jove-": 'grey',
            "P0_Jove+": 'black',
            "P10_Buffer 9": 'red',
            "P8_Buffer 9": 'purple'
        }
        
        highlighted_data = DF_all_data_MG_GFP[DF_all_data_MG_GFP['Buffer'].isin(selected_buffers)]
        for idx, row in highlighted_data.iterrows():
            buffer_name = row['Buffer']
            plate_data_x = row[f'{MG_or_GFP_x}_yield_median']
            plate_data_y = row[f'{MG_or_GFP_y}_yield_median']
            color = highlight_colors.get(buffer_name, 'black')
            
            plt.scatter(plate_data_x, plate_data_y, 
                        color=color, s=80, alpha=0.9, label='_nolegend_')  # Suppress legend
            
            # Add label next to highlighted points
            plt.text(plate_data_x + 0.05, plate_data_y, buffer_name, fontsize=9, ha='left', alpha=0.9)
    
    # Set up the plot appearance
    plt.xlabel(f'Yield for {deb_tit_x} {tit_mid_x} production')
    plt.ylabel(f'Yield for {deb_tit_y} {tit_mid_y} production')
    plt.title('Correlation plot across plates')
    plt.grid(visible=True)
    plt.gca().set_axisbelow(True)  # Set the axis below the grid
    
    if ort == 'ort':
        # Calculate the maximum and minimum values in x and y
        max_x = DF_all_data_MG_GFP[f'{MG_or_GFP_x}_yield_median'].max() + 0.5
        min_x = DF_all_data_MG_GFP[f'{MG_or_GFP_x}_yield_median'].min() - 0.5
        max_y = DF_all_data_MG_GFP[f'{MG_or_GFP_y}_yield_median'].max() + 0.5
        min_y = DF_all_data_MG_GFP[f'{MG_or_GFP_y}_yield_median'].min() - 0.5
        
        # Set the same scale on the x and y axis
        plt.gca().set_xlim([min(min_x, min_y), max(max_x, max_y)])
        plt.gca().set_ylim([min(min_x, min_y), max(max_x, max_y)])
        
        # Make the plot orthogonal
        plt.gca().set_aspect('equal', adjustable='box')
    
    plt.show()

    
def plot_correlation_by_plate(DF_all_data_MG_GFP, MG_or_GFP_x, MG_or_GFP_y, ort='', selected_buffers=None, show_labels=True):
    """
    Plot correlation scatter plot across plates.

    Parameters:
    - DF_all_data_MG_GFP (DataFrame): DataFrame containing all data.
    - MG_or_GFP_x (str): Prefix for x-axis, either 'MG' or 'GFP'.
    - MG_or_GFP_y (str): Prefix for y-axis, either 'MG' or 'GFP'.
    - ort (str, optional): Whether to enforce an orthogonal plot with same scale on both axes. Default is ''.
    - selected_buffers (list, optional): List of buffers to highlight. Default is None.
    - show_labels (bool, optional): Whether to show text labels for highlighted points. Default is True.
    
    Returns:
    None
    """
    # Define titles based on input parameters
    tit_mid_x = 'mRNA' if MG_or_GFP_x == 'MG' else 'protein'
    tit_mid_y = 'mRNA' if MG_or_GFP_y == 'MG' else 'protein'
    
    deb_tit_x = 'maximum'
    deb_tit_y = 'maximum'

    # Number and list of unique plates
    list_name_plates = DF_all_data_MG_GFP['name_plate'].unique()
    nbr_names_plates = len(list_name_plates)

    # Create a color dictionary mapping each plate name to a color from the viridis palette
    large_palette = sns.color_palette("viridis", n_colors=nbr_names_plates)
    color_dict = dict(zip(list_name_plates, large_palette))

    plt.figure(figsize=(8, 6))

    # Plot data for each plate with transparency
    for name_plate in list_name_plates:
        plate_data_x = DF_all_data_MG_GFP[DF_all_data_MG_GFP['name_plate'] == name_plate][f'{MG_or_GFP_x}_yield_median']
        plate_data_y = DF_all_data_MG_GFP[DF_all_data_MG_GFP['name_plate'] == name_plate][f'{MG_or_GFP_y}_yield_median']
        
        plt.scatter(plate_data_x, plate_data_y, 
                    color=color_dict.get(name_plate, 'black'), 
                    s=50, alpha=0.3, label=name_plate)  # Add label for legend

    # Highlight selected buffers with custom colors and no legend entry
    if selected_buffers is not None:
        highlighted_data = DF_all_data_MG_GFP[DF_all_data_MG_GFP['Buffer'].isin(selected_buffers)]

        for idx, row in highlighted_data.iterrows():
            plate_data_x = row[f'{MG_or_GFP_x}_yield_median']
            plate_data_y = row[f'{MG_or_GFP_y}_yield_median']

            # All highlighted points in black
            plt.scatter(plate_data_x, plate_data_y, 
                        color='black', s=50, alpha=0.9, label='_nolegend_')  # Suppress legend

            # Add label next to highlighted points only for plate 'P10' if show_labels is True
            if show_labels and row['name_plate'] == 'P10':
                plt.text(plate_data_x + 0.05, plate_data_y, row['Buffer'], fontsize=9, ha='left', alpha=0.9)
    
    # Set up the plot appearance
    plt.xlabel(f'Yield for {deb_tit_x} {tit_mid_x} production')
    plt.ylabel(f'Yield for {deb_tit_y} {tit_mid_y} production')
    plt.title('Correlation plot across plates')
    plt.grid(visible=True)
    plt.gca().set_axisbelow(True)  # Set the axis below the grid

    # Show the legend excluding highlighted points
    plt.legend(markerscale=1)
    
    if ort == 'ort':
        # Calculate the maximum and minimum values in x and y
        max_x = DF_all_data_MG_GFP[f'{MG_or_GFP_x}_yield_median'].max() + 0.5
        min_x = DF_all_data_MG_GFP[f'{MG_or_GFP_x}_yield_median'].min() - 0.5
        max_y = DF_all_data_MG_GFP[f'{MG_or_GFP_y}_yield_median'].max() + 0.5
        min_y = DF_all_data_MG_GFP[f'{MG_or_GFP_y}_yield_median'].min() - 0.5
        
        # Set the same scale on the x and y axis
        plt.gca().set_xlim([min(min_x, min_y), max(max_x, max_y)])
        plt.gca().set_ylim([min(min_x, min_y), max(max_x, max_y)])
        
        # Make the plot orthogonal
        plt.gca().set_aspect('equal', adjustable='box')
    
    plt.show()

    
def plot_correlation_by_ingredient(DF_all_data_MG_GFP, MG_or_GFP_x, MG_or_GFP_y, ort=''):
    """
    Plot correlation scatter plots for each ingredient.

    Parameters:
    - DF_all_data_MG_GFP (DataFrame): DataFrame containing all data.
    - MG_or_GFP_x (str): Prefix for x-axis title, either 'MG' or 'GFP'.
    - MG_or_GFP_y (str): Prefix for y-axis title, either 'MG' or 'GFP'.
    - ort (str, optional): Whether to enforce orthogonality. Default is ''.

    Returns:
    None
    """
    list_ingredients = ['Mg-glutamate', 'K-glutamate', 'Amino acid', 'Spermidine', '3-PGA', 'NTPs', 'PEG-8000', 'DNA']
    
    # Define titles based on input parameters
    tit_mid_x = 'mRNA' if MG_or_GFP_x == 'MG' else 'protein'
    tit_mid_y = 'mRNA' if MG_or_GFP_y == 'MG' else 'protein'

    deb_tit_x = 'maximum'
    deb_tit_y = 'maximum'

    # Set number of rows and columns for subplots
    n_rows = 2
    n_cols = 4

    # Create subplots with shared x-axis
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, ingredient in enumerate(list_ingredients):
        # Sort data by the ingredient values
        DF_sorted = DF_all_data_MG_GFP.sort_values(by=ingredient)
        
        # Create color palette for the ingredient values
        unique_values = DF_sorted[ingredient].unique()
        palette = sns.color_palette("viridis", n_colors=len(unique_values))
        color_dict = dict(zip(unique_values, palette))
        
        ax = axes[idx]
        
        # Plot data for each ingredient value
        for value in unique_values:
            data_x = DF_sorted[DF_sorted[ingredient] == value][f'{MG_or_GFP_x}_yield_median']
            data_y = DF_sorted[DF_sorted[ingredient] == value][f'{MG_or_GFP_y}_yield_median']
            
            ax.scatter(data_x, data_y, color=color_dict.get(value, 'black'), label=value, s=5)

        ax.set_title(ingredient)
        ax.grid(True)
        ax.set_axisbelow(True)

        if idx % n_cols == 0:
            ax.set_ylabel(f'Yield for {deb_tit_y} {tit_mid_y} production')
        if idx >= (n_rows - 1) * n_cols:
            ax.set_xlabel(f'Yield for {deb_tit_x} {tit_mid_x} production')

        # Add legend outside the plot
        handles, labels = ax.get_legend_handles_labels()
        
        if ort == 'ort':
            # Calculate the maximum and minimum values in x and y with margins
            max_x = DF_sorted[f'{MG_or_GFP_x}_yield_median'].max() + 0.5
            min_x = DF_sorted[f'{MG_or_GFP_x}_yield_median'].min() - 0.5
            max_y = DF_sorted[f'{MG_or_GFP_y}_yield_median'].max() + 0.5
            min_y = DF_sorted[f'{MG_or_GFP_y}_yield_median'].min() - 0.5
            
            # Set the same scale on the x and y axis
            ax.set_xlim([min(min_x, min_y), max(max_x, max_y)])
            ax.set_ylim([min(min_x, min_y), max(max_x, max_y)])
            
            # Make the plot orthogonal
            ax.set_aspect('equal', adjustable='box')

    fig.suptitle('Correlation plot by ingredient', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.98])  # [left, bottom, right, top]
    plt.show()


def mutual_info_score(DF_all_data, MG_or_GFP):
    """
    Plot the mutual information scores between features and target variable.
    
    Parameters:
    - DF_all_data (DataFrame): DataFrame containing all data.
    - MG_or_GFP (str): Prefix for target variable, either 'MG' for mRNA or 'GFP' for protein.
    
    Returns:
    None
    """
    # Define titles and colors based on the MG_or_GFP argument
    if MG_or_GFP == 'MG':
        mid_title = 'mRNA'
        color = 'coral'
    elif MG_or_GFP == 'GFP':
        mid_title = 'protein'
        color = 'yellowgreen'
    else:
        raise ValueError("A correct MG_or_GFP argument is required: 'MG' for mRNA, 'GFP' for protein data")
        
    deb_title = 'maximum' 

    # Extract features and target variable
    features = DF_all_data[['Mg-glutamate', 'K-glutamate', 'Amino acid', 'Spermidine', '3-PGA', 'NTPs', 'PEG-8000', 'DNA']]
    target = DF_all_data[f'{MG_or_GFP}_yield_median']

    # Calculate mutual information scores
    mi_scores = mutual_info_regression(features, target)
    
    # Plot the mutual information scores
    plt.figure(figsize=(10, 6))
    plt.bar(features.columns, mi_scores, color=color)
    plt.xlabel('Features')
    plt.ylabel('Mutual information score')
    plt.title(f'Mutual information score for {deb_title} {mid_title}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

#%% [protein+MG] INPUTS
      
Name_plate = 'plate 10'

# Call the function to fill the dictionary with the correct values depending on Name_plate
plate_info = get_plate_variables(Name_plate)  # Example usage

# Unwrap the dictionary if plate_info is not empty
if plate_info:
    Excel_file = plate_info.get('Excel_file')
    DF_data_GFP = plate_info.get('DF_data_GFP')
    DF_data_MG = plate_info.get('DF_data_MG')
    Time_step = plate_info.get('Time_step')
    Name_plate = plate_info.get('Name_plate')
    nbr_replicate = plate_info.get('nbr_replicate')
    nbr_controls = plate_info.get('nbr_controls')
    name_ref_plus_MG = plate_info.get('name_ref_plus_MG')
    name_ref_plus_GFP = plate_info.get('name_ref_plus_GFP')
    
#%% [protein+MG] LEGEND

# Call the function for processing data MG
list_strains, dico_strains_wells = legend(Excel_file)

#%% INITIAL CONDITIONS : concentrations of each ingredients of the buffer per well for GFP

# Call the functiun importing initial condition 
DF_sampling = get_initial_conditions(Name_plate)

#%% [MG] CALCULATION OF FLUO - MAXIMUM FOR MG FROM RAW DATA
'''This part calculates the fluo for each buffer, from raw MG data.
The fluo is defined as the max within the two firsts hours of our timelines. 
INPUTS:  - DF_data_MG_corr (all wells over time, preincub subtracted)
OUTPUTS: - DF_results_fluo_MG (buffer as index, one column median and one column stdev)
'''
# Create a Dataframe with only the relevant timepoint fluorescence for each well (for MG: max within 10h)
DF_data_MG_corr_maximum = DF_data_MG.squeeze()

# Call the function with your input data
DF_results_fluo_MG = calc_median_std_per_buffer(DF_data_MG_corr_maximum, dico_strains_wells)

# Call the function for FLUO values of MAXIMUM mRNA production
barplot_one_plate(DF_results_fluo_MG, 'fluo', 'MG')

#%% [GFP] CALCULATION OF FLUO - MAXIMUM FOR GFP FROM RAW DATA
'''This part calculates the fluo for each buffer, from raw GFP data.
The fluo is defined as the vale of the last point. 
INPUTS:  - DF_data_GFP_corr (all wells over time, preincub subtracted)
OUTPUTS: - DF_results_fluo_GFP (buffer as index, one column median and one column stdev)
'''
# Create a Dataframe with only the relevant timepoint fluorescence for each well (for GFP: mean of the final points)
DF_data_GFP_corr_maximum = DF_data_GFP.squeeze()

# Call the function with your input data
DF_results_fluo_GFP = calc_median_std_per_buffer(DF_data_GFP_corr_maximum, dico_strains_wells)

# Call the function for FLUO values of MAXIMUM mRNA production
barplot_one_plate(DF_results_fluo_GFP, 'fluo', 'GFP')

#%% [MG] [FLUO] CONTROLS ANALYSIS IN PLATE n+1

### For MG
plot_control_buffers(dico_strains_wells, DF_results_fluo_MG, nbr_controls, name_ref_plus_MG, list_strains_controls(Name_plate), Name_plate, 'MG')
### For GFP
plot_control_buffers(dico_strains_wells, DF_results_fluo_GFP, nbr_controls, name_ref_plus_GFP, list_strains_controls(Name_plate), Name_plate, 'GFP')

#%% CALCULATION OF YIELD WITH NaN - calcul

### For MG
DF_results_yield_MG, DF_yield_MG_WperW = calc_yield_NaN(dico_strains_wells, DF_data_MG_corr_maximum, DF_results_fluo_MG, name_ref_plus_MG)
### For GFP
DF_results_yield_GFP, DF_yield_GFP_WperW = calc_yield_NaN(dico_strains_wells, DF_data_GFP_corr_maximum, DF_results_fluo_GFP, name_ref_plus_GFP)

#%% CORRECTION OF GAIN CHANGE

# Call the function if necessary (gain has been changes in plate 8 and following)
if re.search(r'\d+', Name_plate) and int(re.search(r'\d+', Name_plate).group()) > 7:
    # Call the function correct_gain_change for correcting the yield based on maximum value of mRNA production
    DF_results_yield_MG, DF_yield_MG_WperW, Ratio_Refplus_Joveplus_MG = correct_gain_change(name_ref_plus_MG, Name_plate, DF_results_yield_MG, DF_yield_MG_WperW, 'MG')
else:
    print("Plate number is not greater than 7. Skipping the call of the function correct_gain_change.")
    
#%% CALCULATION OF YIELD WITH NaN - barplots

### For MG
barplot_one_plate(DF_results_yield_MG, 'yield', 'MG')
### For GFP
barplot_one_plate(DF_results_yield_GFP, 'yield', 'GFP')

#%%[MG] COMPILATION OF DATA FOR THE CURRENT PLATE

# Iterate through the dictionary and filter out 'nan' values from each list
for key, value in dico_strains_wells.items():
    dico_strains_wells[key] = [item for item in value if item is not np.nan]

# PROCESS DATA USING MAXIMUM FOR MG
DF_fluo_values_max_MG, DF_wells_max_MG = process_fluo_DF_data_corr(dico_strains_wells, DF_data_MG_corr_maximum, 'MG')
DF_fluo_median_std_max_MG = process_fluo_DF_results_fluo(DF_results_fluo_MG, 'MG')
DF_yield_values_max_MG = process_yield_DF_WperW(DF_yield_MG_WperW, 'MG')
DF_yield_median_std_max_MG = process_yield_DF_results_yield(DF_results_yield_MG, 'MG')

# PROCESS DATA USING MAXIMUM FOR GFP
DF_fluo_values_max_GFP, DF_wells_max_GFP = process_fluo_DF_data_corr(dico_strains_wells, DF_data_GFP_corr_maximum, 'GFP')
DF_fluo_median_std_max_GFP = process_fluo_DF_results_fluo(DF_results_fluo_GFP, 'GFP')
DF_yield_values_max_GFP = process_yield_DF_WperW(DF_yield_GFP_WperW, 'GFP')
DF_yield_median_std_max_GFP = process_yield_DF_results_yield(DF_results_yield_GFP, 'GFP')

# CONCATENATE
plate_AL_X_corr_everything_MG_GFP = pd.concat([DF_sampling, 
                                           DF_fluo_values_max_MG, DF_fluo_median_std_max_MG, 
                                           DF_yield_values_max_MG, DF_yield_median_std_max_MG,
                                           DF_fluo_values_max_GFP, DF_fluo_median_std_max_GFP, 
                                           DF_yield_values_max_GFP, DF_yield_median_std_max_GFP], axis=1)
plate_AL_X_corr_everything_MG_GFP['name_plate'] = Name_plate
plate_AL_X_corr_everything_MG_GFP['Buffer'] = list_strains
plate_AL_X_corr_everything_MG_GFP = pd.concat([plate_AL_X_corr_everything_MG_GFP, DF_wells_max_MG], axis=1)

#Converting a dataframe to a csv Excel file
plate_AL_X_corr_everything_MG_GFP.to_csv(Name_plate+'_AL_corr_everything_MG_GFP.csv', index=False)

#%% [protein/MG] COMPILATION OF COMPLETE DATA FILES ACROSS PLATES

Last_plate = Name_plate

### Adapt name of Name_plate_compil for the compilated files
all_Name_plate_compil = 'P0_P1_P2_P3_P4_P5_P6_P7_P8_P9_P10'
# Extract the last part from Last_plate
last_part = Last_plate.split()[-1]
# Find the corresponding substring in all_Name_plate_compil
Name_plate_compil = all_Name_plate_compil[:all_Name_plate_compil.index(last_part) + len(last_part)]

### Adapt the source files to use depending on Last_plate
# Define the list of plate names
list_Name_plate = ['plate 0', 'plate 1', 'plate 2', 'plate 3', 'plate 4', 'plate 5', 'plate 6', 'plate 7','plate 8','plate 9','plate 10']
# Generate list_Excel_files containing plates up to Name_plate
list_Excel_files = [plate + '_AL_corr_everything_MG_GFP' for plate in list_Name_plate[:list_Name_plate.index(Last_plate)+1]]

list_index = ['Mg-glutamate', 'K-glutamate', 'Amino acid', 'Spermidine', '3-PGA', 'NTPs', 'PEG-8000', 'DNA', 
       'MG_fluo_1', 'MG_fluo_2', 'MG_fluo_3', 'MG_fluo_4', 'MG_fluo_5', 'MG_fluo_6','MG_fluo_median', 'MG_fluo_std',
       'MG_yield_1', 'MG_yield_2','MG_yield_3', 'MG_yield_4', 'MG_yield_5', 'MG_yield_6','MG_yield_median', 'MG_yield_std',
       'GFP_fluo_1','GFP_fluo_2', 'GFP_fluo_3', 'GFP_fluo_4','GFP_fluo_5', 'GFP_fluo_6','GFP_fluo_median', 'GFP_fluo_std',
       'GFP_yield_1', 'GFP_yield_2', 'GFP_yield_3', 'GFP_yield_4','GFP_yield_5', 'GFP_yield_6', 'GFP_yield_median', 'GFP_yield_std',
       'name_plate', 'Buffer', 'well_1', 'well_2', 'well_3','well_4', 'well_5', 'well_6']

# Read the first file to initialize the DataFrame
DF_all_data_MG_GFP = pd.read_csv(list_Excel_files[0] + '.csv')

# Iterate over the rest of the files and concatenate them
for file_name in list_Excel_files[1:]:
    df = pd.read_csv(file_name + '.csv')
    DF_all_data_MG_GFP = pd.concat([DF_all_data_MG_GFP, df], ignore_index=True)
DF_all_data_MG_GFP = DF_all_data_MG_GFP.reindex(columns=list_index)

# Export
DF_all_data_MG_GFP.to_csv(Name_plate_compil + '_AL_corr_everything_MG_GFP.csv', index=False)

### For next analysis
## Informative index
DF_all_data_MG_GFP['Concatenated'] = DF_all_data_MG_GFP['Buffer']+' in '+DF_all_data_MG_GFP['name_plate']
DF_all_data_MG_GFP.set_index('Concatenated', inplace=True)

#%% CONTROLS COMPARISON ACROSS PLATES

### For MG
unique_buffers = subplot_controls_across_plates(DF_all_data_MG_GFP, Name_plate, 'MG')
### For GFP
unique_buffers = subplot_controls_across_plates(DF_all_data_MG_GFP, Name_plate, 'GFP')

#%% FIGURE OF COMPILATION WITH ALIGNED HEATMAP AND BARPLOT

### For MG
plot_barplot_heatmap_compil(DF_all_data_MG_GFP, name_ref_plus_MG, 'MG')
### For GFP
plot_barplot_heatmap_compil(DF_all_data_MG_GFP,name_ref_plus_GFP, 'GFP')


#%% FIGURE OF COMPILATION WITH ALIGNED HEATMAP AND BARPLOT - light legend

# Change text size for all plots
plt.rcParams.update({'font.size': 28})  # Adjust the value as needed

### For MG
plot_barplot_heatmap_compil(DF_all_data_MG_GFP, name_ref_plus_MG, 'MG',40,15)
### For GFP
plot_barplot_heatmap_compil(DF_all_data_MG_GFP, name_ref_plus_GFP, 'GFP',40,15)

# Back to the default parameter
plt.rcParams['font.size'] = plt.rcParamsDefault['font.size']


#%% EVOLUTION OF YIELD PLATE PER PLATE SLEEK

if Name_plate != 'plate 0':
    # Change text size for all plots
    plt.rcParams.update({'font.size': 25})  # Adjust the value as needed
    
    ### For MG
    plot_barplot_all_data_plate_by_plate(DF_all_data_MG_GFP, 'MG')
    ### For GFP
    plot_barplot_all_data_plate_by_plate(DF_all_data_MG_GFP, 'GFP')
    
    # Back to the default parameter
    plt.rcParams['font.size'] = plt.rcParamsDefault['font.size']
    
#%% EVOLUTION OF YIELD PLATE PER PLATE SLEEK WITH HEATMAP

if Name_plate != 'plate 0':
    # Change text size for all plots
    plt.rcParams.update({'font.size': 25})  # Adjust the value as needed
    
    ### For MG
    plot_barplot_heatmap_compil_plate_per_plate(DF_all_data_MG_GFP, name_ref_plus_MG, 'MG')
    ### For GFP
    plot_barplot_heatmap_compil_plate_per_plate(DF_all_data_MG_GFP, name_ref_plus_GFP, 'GFP')
    
    # Back to the default parameter
    plt.rcParams['font.size'] = plt.rcParamsDefault['font.size']
    
#%% BOXPLOT ACROSS PLATES

if Name_plate != 'plate 0':
    ### For MG
    boxplot_yield_across_plates('MG')
    ### For GFP
    boxplot_yield_across_plates('GFP')

#%% BOXPLOT ACROSS PLATES - SLEEK VERSION
if Name_plate != 'plate 0':
    DF_all_data_MG_GFP_sleek = DF_all_data_MG_GFP.copy()
    
    ### For MG
    boxplot_yield_across_plates_sleek('MG')
    ### For GFP
    # Call the function for a detailed figure based on MAXIMUM protein
    boxplot_yield_across_plates_sleek('GFP')

#%% PLOT CORRELATION COLORED BY PLATE

# FOR MG_MAX VERSUS GFP_MAX
selection = ["P0_Jove-", "P0_Jove+", "P10_Buffer 9", "P8_Buffer 9"]
plot_correlation_by_plate_labels(DF_all_data_MG_GFP, 'MG', 'GFP', selected_buffers=selection)

#%% PLOT CORRELATION COLORED BY INGREDIENT

# FOR MG_MAX VERSUS GFP_MAX
plot_correlation_by_ingredient(DF_all_data_MG_GFP, 'MG', 'GFP')

#%% PCA (Principal Component Analysis)

Name_data = 'mRNA data'

# Feature definition
feature_columns = ['Mg-glutamate', 'K-glutamate', 'Amino acid', 'Spermidine', '3-PGA', 'NTPs', 'PEG-8000', 'DNA','MG_yield_median']
# Target definition
target_column = 'MG_yield_median'
# Extract features (X) and target variable (y)
X = DF_all_data_MG_GFP[feature_columns]
y = DF_all_data_MG_GFP[target_column]

# Standardize the features (important for PCA)
X_standardized = (X - X.median()) / X.std()

# Perform PCA
pca = PCA()
X_pca = pca.fit_transform(X_standardized)

# Plot the explained variance ratio using a barplot
explained_var_ratio = pca.explained_variance_ratio_
plt.bar(range(1, len(explained_var_ratio) + 1), explained_var_ratio, color='coral')
plt.xlabel('Principal components')
plt.ylabel('Explained variance ratio')
plt.title('Explained variance ratio for each principal component for '+Name_data+' - '+Name_plate)
plt.show()

# Choose the number of components based on the explained variance ratio
n_components = np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.95) + 1
print(f"Number of components to explain 95% variance: {n_components}")

# Plot PCA projection with heatmap
pc_df = pd.DataFrame(data=X_pca[:, :2], columns=['PC1', 'PC2'])
pc_df[target_column] = y.values

plt.figure(figsize=(12, 8))
sns.scatterplot(x='PC1', y='PC2', hue=pc_df[target_column], data=pc_df, cmap='rocket')
plt.xlabel('Component 1 ({:.2%})'.format(pca.explained_variance_ratio_[0]))
plt.ylabel('Component 2 ({:.2%})'.format(pca.explained_variance_ratio_[1]))
plt.title('PCA projection for '+Name_data+' - '+Name_plate)

# Plot the correlation circle
coefficients = np.transpose(pca.components_[:2, :])
labels = X.columns

fig, ax = plt.subplots(figsize=(10, 10))
circle = plt.Circle((0, 0), 1, facecolor='none', edgecolor='b')
ax.add_artist(circle)

for i, (label, coeff) in enumerate(zip(labels, coefficients)):
    ax.arrow(0, 0, coeff[0], coeff[1], head_width=0.05, head_length=0.1, fc='r', ec='r')
    ax.text(coeff[0] * 1.15, coeff[1] * 1.15, label, color='k', ha='center', va='center')

plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Correlation circle for '+Name_data+' - '+Name_plate)
plt.grid()
plt.show()

#%% MUTUAL INFORMATION SCORE

### For MG
mutual_info_score(DF_all_data_MG_GFP,'MG')
### For GFP
mutual_info_score(DF_all_data_MG_GFP,'GFP')
