# PET/fMRI Integration

This academic project investigates the relationship between dynamic Positron Emission Tomography (PET) and resting-state functional MRI (fMRI) data using MATLAB. The analysis combines tracer kinetic modeling with functional connectivity metrics to assess how glucose metabolism relates to brain network organization.

## Project Description
This project was completed as part of the Imaging for Neuroscience course. The goal was to explore possible associations between PET-derived metabolic activity and fMRI-derived connectivity in a healthy subject, using data acquired from a hybrid PET/MRI system.

## Aims of the Study
- Quantify tracer uptake (Ki) from PET data using Patlak graphical analysis at both ROI and voxel levels.
- Extract functional connectivity (FC) matrices and compute clustering coefficients from resting-state fMRI data.
- Investigate potential correlations between metabolic activity (Ki) and functional network segregation (clustering coefficient).

## Technologies & Tools
- Language: MATLAB
- Toolboxes:
  - Brain Connectivity Toolbox
  - SPM12-compatible filters and tools
  - Custom utility functions for NIfTI I/O and data processing

The dataset used is too heavy to be uploaded here, but it's available [here](https://drive.google.com/drive/folders/1LaK51OAIB6LonDEf_wV95aNQOutJtqCT?usp=share_link).

## Folder Structure
```bash
PET_fMRI_integration/
│
├── main_script.m                    # MATLAB script for PET/fMRI analysis
├── PET_data.mat                     # Ki ROI and grey mask matrices computed during the first PET analysis, used in the final integration with fMRI
├── PET_results_HW1_group2.mat       # Results of the PET analysis
├── fMRI_results_HW1_group2.mat      # Results of the fMRI analysis
├── solution_report.pdf              # Final PDF presentation
├── assignment.pdf                   # Original assignment instructions
└── README.md
```

## Main Steps of the Analysis

### PET Analysis
- Thresholded GM mask creation from segmentation
- Extraction of TACs at ROI level
- Patlak graphical modeling for Ki estimation
- 3D Ki maps and ROI-level summaries
- Comparison of voxel vs ROI average estimates

### fMRI Analysis
- WM/CSF signal extraction and PCA denoising
- ROI time-series extraction and motion regression
- High-pass filtering (cutoff = 1/128 Hz, Butterworth)
- FC matrix computation (Pearson correlation + Fisher’s Z)
- Clustering coefficient extraction from binarized FC matrix

### Integration
Scatterplot and correlation analysis between PET Ki and fMRI clustering coefficient at ROI level

## Key Findings
- Strong agreement (r = 0.99) between voxel-averaged and ROI-based Ki estimates.
- No significant correlation found between Ki values and clustering coefficients, suggesting they represent different physiological phenomena (local metabolism vs network integration).
