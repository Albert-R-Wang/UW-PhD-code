# Image Analysis

This folder contains MATLAB scripts for microscopy image analysis, including tools for quantification of various cellular features. The goal is to provide flexible and reproducible pipelines for extracting meaningful data from fluorescence imaging experiments.

## Contents

- `fociCount.m`: A semi-automated pipeline for quantifying γH2AX foci at the single-cell level, with support for dynamic and universal thresholding, per-nucleus metric extraction, and annotated outputs.
  - `fociCount_workflow.pdf`: A workflow overview document that describes the foci counting pipeline and guides users through image processing.
- `AR_translocation_analysis.m`: This script quantifies androgen receptor (AR) nuclear translocation by calculating nucleus-to-cytoplasm signal ratios from TIFF images, with user-defined thresholding, region segmentation, and export of annotated overlays.
  - `AR_translocation_analysis.pdf`
- `ScratchWoundAnalysis.m`: This script analyzes scratch wound assay images by detecting manually drawn boundary lines, filling the enclosed wound region, and quantifying both wound area and length for each image.
  - `ScratchWoundAnalysis.pdf`

## Notes

- The code is provided as-is for reference and documentation purposes.
- Project-specific data and raw input files are not included due to privacy and data sharing restrictions.
