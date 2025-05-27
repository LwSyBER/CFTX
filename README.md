# CFTX (Cell-free transcriptomics)
This repository contains the code and data associated to the publication XXX. It contains two folders:
1-From biotek data to active learning: contains the code for making biotek export into a list of yield associated to rach buffer, to feed the machine learning algorithm. 
2-Active learning: using the output of the previous code, contains the code that trains and runs the machine learning model to suggest new compositions to test in next plate.
To go from this list of compositions to Echo instructions, see github https://github.com/brsynth/icfree-ml

# 1-From biotek data to active learning
## General description of the data
Each cell-free reaction is performed in a 384-well plate in 10.5 µL and consists of :
- a fixed quantity of *E. coli* BL21 DE3 cell lysate (3.5 µ)
- a variable buffer
- a reporter plasmid of varying concentration (corresponding to the "DNA" ingredient in the buffer).

This dataset measures over time :
- the red fluorescence (590/635 nm) representing the quantity of messenger RNA thanks to the malachite green aptamer encoded on the reporter plasmid
- the green fluorescence (485/528 nm) representing the quantity of protein thanks to the deGFP encoded on the reporter plasmid

The aim of this experiment was to optimise the buffer for cell-free mRNA production. To achieve this, an active learning algorithm was used to suggest at each iteration (1 plate = 1 iteration) the next most interesting buffer compositions to be tested. The initialization was done by a Latin hypercube sampling algorithm (plate 0), and the active learning corresponds to plates 1 to 10. This result was obtained (20-fold improvement in the mRNA quantity) and in the process we collected a large dataset that could be exploited.

A total of 650 different buffer compositions were tested overall:
![image](https://github.com/user-attachments/assets/d80b2178-7cbc-4462-8b69-a56bf60c7db6)

## Inputs
### Raw data
The dataset is provided as .xlsm Excel files, with one file per plate. Each file consists of several sheets:
- **Data_GFP**: Contains the maximum fluorescence values for GFP in each well, used to track translation over time.
- **Data_MGap**t: Contains the maximum fluorescence values for the malachite green aptamer in each well, used to track transcription over time.
- **layout**: Displays the plate layout, linking each well to its corresponding buffer condition.
- **layoutW**: A version of the layout formatted for easier use in Python.
- **strains**: Lists the conditions that will be considered and analyzed when executing the code.
- **sampling**: Details the composition of each buffer in terms of Mg-glutamate (mM), K-glutamate (mM), amino acids (mM), spermidine (mM), 3-PGA (mM), NTPs (mM), PEG-8000 (%), DNA (nM), HEPES (mM), and malachite green (mM).

*Nota bene:* Each buffer composition is referenced by a unique identifier that includes the plate number and buffer number. 
For example:
**P1_Buffer 1** refers to the first buffer tested on plate 1.
**P2_Buffer 1** refers to the first buffer tested on plate 2.
The “P0_Jove-” and “P0_Jove+” buffers are the negative and positive controls used to calculate the yield.

### Additional file
**Params_concentration_range.csv**: This file is used to generate heatmaps representing the buffer compositions.

## Outputs
The code provides:
- **plate X_AL_corr_everything_MG_GFP.csv**: analysed data for the plate X.
  - Each row corresponds to a buffer that is identified by the columns **Buffer** (*e.g.* P5_Buffer 21) and **name_plate** (*e.g.* plate 8).
  - The wells were this buffer has been tested are indicated in the columns **well_1,	well_2,	well_3,	well_4,	well_5,	well_6** according to the number of replicates.
  - The composition of the buffer is given by the first eight columns: **Mg-glutamate,	K-glutamate,	Amino acid,	Spermidine,	3-PGA,	NTPs,	PEG-8000,	DNA**. All concentrations are in mM, except for PEG-8000 (%) and DNA (nM). The concentration of HEPES, malachite green and lysate are not indicated here because they are fixed in all experiments: 50 mM HEPES, 0.010 mM malachite green and 3.5µL of lysate among the 10.5µL final.
  - _For the malachite green aptamer data:_
  - The maximum fluorescence value for each well (corresponding to the raw data) are given in the columns **MG_fluo_1,	MG_fluo_2,	MG_fluo_3,	MG_fluo_4,	MG_fluo_5,	MG_fluo_6** according to the number of replicates.
  - The median and the standard deviation of these fluorescence values are calculated per buffer in the columns **MG_fluo_median** and	**MG_fluo_std**.
  - The yield (normalised fluorescence values) is calculated for each well and given in the columns **MG_yield_1,	MG_yield_2,	MG_yield_3,	MG_yield_4,	MG_yield_5,	MG_yield_6** according to the number of replicates.
  - The median and the standard deviation of these yields values are calculated per buffer in the columns **MG_yield_median** and	**MG_yield_std**.
  - _For the GFP data:_
  - The maximum fluorescence value for each well (corresponding to the raw data) are given in the columns **GFP_fluo_1,	GFP_fluo_2,	GFP_fluo_3,	GFP_fluo_4,	GFP_fluo_5,	GFP_fluo_6** according to the number of replicates.
  - The median and the standard deviation of these fluorescence values are calculated per buffer in the columns **GFP_fluo_median** and	**GFP_fluo_std**.
  - The yield (normalised fluorescence values) is calculated for each well and given in the columns **GFP_yield_1,	GFP_yield_2,	GFP_yield_3,	GFP_yield_4,	GFP_yield_5,	GFP_yield_6** according to the number of replicates.
  - The median and the standard deviation of these yields values are calculated per buffer in the columns **GFP_yield_median** and	**GFP_yield_std**.

- **P0_P1_P2_P3_PX_AL_corr_everything_MG_GFP.csv**: analysed data for the plate X compiled with all previous data.
  - This file contains the same columns and information as described above, but for all plates together. It allows to plot graphs representing the evolution of MG_yield and GFP_yield between plates.
 
- **correct_gain_change_params_MG_plate X.txt** for plates 8, 9 and 10 for which the yield is calculated differently to correct for the change of gain of the plate reader.

- **Numerous figures .svg** among which the figures of the paper.

## A few subtelties
### Experimental biaises and the way there are corrected
- between plate 0 and plates 1-10, the plate sealer has been changed, reducing the value of the red background. Consequently, it is normal for the raw fluorescence of plate 0 to be higher than that of subsequent plates. This is corrected by the calculation of the yield, that normalized all fluorescence values to make them comparable across plates.
- from plate 7 onwards, the red fluorescence rose so much that we had to change a parameter in the plate reader called the gain. This forced us to calculate the yield a little differently. It is therefore normal for the raw fluorescence of plate 7 to be higher than that of plates 8, 9 and 10.
![image](https://github.com/user-attachments/assets/c336d41c-cb84-4617-9e0d-ef04a4175a1f)

- for plates 9 and 10 only, some wells are removed as outliers by the remove_outliers function. These are due to experimental error (I know that in some wells/columns I failed to add the lysate correctly, and this is visible when looking at the data). This should not hinder the analysis as there are always at least 5 replicates per buffer left.

### Controls
- in each plate, a certain number of controls (= buffers already tested in previous plates) were tested. They are used to:
1) check the comparability of data from different plates
2) correct the data of plates 8, 9 and 10 after changing the gain.
Thus, the buffer "P5_Buffer 21" was created in plate 5, and was tested as control in the following plates (plate 6, 8, 9, 10).

These controls validate the use of the yield to correct the experimental biaises listed above:
![image](https://github.com/user-attachments/assets/911c0b2f-444f-4430-981a-45f345c5e367)

## Authors
This code was written by Léa WAGNER. 

# 2-Active learning
The script “active_learning_bayesian_style_vf.py” generates a list of buffer compositions to test in the next plate according to active learning. Yield optimization was performed using a Bayesian optimization active learning approach, with a Gaussian Process (GP) model to probabilistically estimate the relationship between component concentrations and yield.

## Authors
This code was written by An HOANG. 

