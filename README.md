# CFTX (Cell-free transcriptomics)
This repository contains the code and data associated to the publication XXX.
It allows to reproduce the figures of the paper and plot that data of each plate. 
As an active learning code, this code was designed to be run plate per plate, to feed the machine learning algorithm. You must run the different plates from 0 to 10 to be able to generate a compiled file of all data (corresponding to supplementary table X).

# Inputs
## Raw data
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

## Additional file
**Params_concentration_range.csv**: This file is used to generate heatmaps representing the buffer compositions.

# Outputs
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





  The first 8 columns contain the composition

# Authors
This code was written by Léa WAGNER (PhD student). 
