# PKD_eDNA
This repository contains scripts used for the analysis of qPCR and ddPCR data to evaluate methods for the detection of _T. bryosalmonae_ eDNA, for reproducing results presented in Stelzer _et al_ (2023).

## LOD_LOQ_qPCR_ddPCR_JO130123.R
Script used to calculate limits of detection and quantification for qPCR and ddPCR assays, based on methods used previously by Klymus _et al_ (2020) (https://doi.org/10.1002/edn3.29)
### Input files:
- data/LOD_LOQ/ddPCR_for_LOD_LOQ_small_assaysrenamed.csv
- data/LOD_LOQ/qPCR_for_LOD_LOQ_corrected_assaysrenamed.csv

## data_processing_env_samples_200123.R
The script estimates copy numbers from qPCR CQ values using a linear model, then compares these to copy numbers derived by ddPCR.
It also generates the detections matrix used for the occupancy modelling (data/environmental/detections_qpcr_ddpcr_filtered_200123.txt).
### Input files:
- data/environmental/ddpcr_counts_final.csv
- data/environmental/qPCR_data_merged_200123.csv

## internal_control_CQ_deviation_250823.R
Calculates the difference in CQ values of internal spike-in controls (IC) between field samples and standards/negatives to test the possibility of qPCR inhibition in field samples.
The IC-CQ of a technical replicate is subtracted from the average standard/negative IC-CQ for the same run, and the average difference calculated for each sample.
### Input files
- data/environmental/qPCR_data_with_IC_CQ.xlsx

## msocc_models_200123.R
Performs multilevel occupancy modeling to estimate sample-level capture probabilities and replicate-level detection probabilities as functions of sample volume, filter type, and assay type (qPCR or ddPCR).
It used models implemented in the msocc package (https://github.com/StrattonCh/msocc)
### Input files
- data/environmental/detections_qpcr_ddpcr_filtered_200123.txt

## sample_size_sim_030223.R
Based on estimates of sample-level capture probability and replicate-level detection proability from a multilevel model fit with msocc, calculates the cumulative probsbility of PKD eDNA detection for different numbers of collected samples.
### Input files
- detections_qpcr_ddpcr_filtered_200123.txt

Note: data/environmental/qPCR_data_merged_200123.csv and data/environmental/qPCR_data_with_IC_CQ.xlsx contion data from the same samples. The former was modified so that sample IDs of the Sterivex samples would match those used for the same samples in the ddPCR data.
