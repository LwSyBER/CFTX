# CFTX (Cell-free transcriptomics)
This repository contains the code and data associated to the publication XXX.

# Inputs
## Raw data
The dataset is provided as .xlsm Excel files, with one file per plate. Each file consists of several sheets:
- **Data_GFP**: Contains the maximum fluorescence values for GFP in each well, used to track translation over time.
- **Data_MGap**t: Contains the maximum fluorescence values for the malachite green aptamer in each well, used to track transcription over time.
- **layout**: Displays the plate layout, linking each well to its corresponding buffer condition.
- **layoutW**: A version of the layout formatted for easier use in Python.
- **strains**: Lists the conditions that will be considered and analyzed when executing the code.
- **sampling**: Details the composition of each buffer in terms of Mg-glutamate (mM), K-glutamate (mM), amino acids (mM), spermidine (mM), 3-PGA (mM), NTPs (mM), PEG-8000 (%), DNA (nM), HEPES (mM), and malachite green (mM).

*Nota bene:* Each buffer composition is referenced by a unique identifier that includes the plate number and buffer number. For example:
P1_Buffer 1 refers to the first buffer tested on plate 1.
P2_Buffer 1 refers to the first buffer tested on plate 2.

## Additional file
Params_concentration_range.csv: This file is used to generate heatmaps representing the buffer compositions.

