# Project：Preliminary Ecological Analysis

## Overview
This repository contains the preliminary ecological analysis for **Preliminary Ecological Analysis**. 
We investigated the potential co-selection between Colistin resistance and alternative antibiotics (TMP-SMX, Tetracycline) in bovine *E. coli* using surveillance data from the RESAPATH network (France).

## Hypothesis
We hypothesized that the use of alternative antibiotics drives the persistence of Colistin resistance through chromosomal integration of mobile genetic elements.

## Data Source
- **RESAPATH Network:** Clinical bovine *E. coli* isolates.
- **Years Analyzed:** 2024 (Current) and 2019 (Historical).
- **Level:** Departmental level (French administrative regions).

## Methodology
- Data cleaning and extraction using custom R scripts ("Surgical Extraction" method).
- Pearson correlation analysis between resistance rates of:
  - Colistin vs. TMP-SMX
  - Colistin vs. Tetracycline

## Results (Summary)
Contrary to the co-selection hypothesis, **no significant correlation** was found at the departmental level in either 2019 or 2024.

| Year | Antibiotic Pair | Pearson r | p-value | Result |
| :--- | :--- | :--- | :--- | :--- |
| **2024** | Colistin vs TMP-SMX | -0.12 | 0.52 | No Correlation |
| **2024** | Colistin vs Tetracycline | -0.26 | 0.17 | No Correlation |
| **2019** | Colistin vs TMP-SMX | 0.09 | 0.61 | No Correlation |
| **2019** | Colistin vs Tetracycline | 0.28 | 0.13 | No Correlation |

TMP-SMX vs Colistin: Pearson r=−0.12,p=0.52
![Correlation Plot](plots/correlation_plot_tmp_2019.jpg)
![Correlation Plot](plots/correlation_plot_tmp_2024.jpg)

Tetracycline vs Colistin: Pearson r=−0.26,p=0.17
![Correlation Plot](plots/correlation_plot_terra_2019.jpg)
![Correlation Plot](plots/correlation_plot_terra_2024.jpg)

## Discussion 
The consistent lack of correlation in aggregated regional data empahasizes the **"Ecological Challenge"** of macroscopic surveillance. It suggests that mechanisms of chromosomal integration of mobile genetic elements are invisible to current monitoring strategies. Also, the lack of correlation with antibiotic usage indicates that resistance is no longer driven by co-selection pressure, but has likely stabilized via chromosomal integration. These support the need for the genomic approach to be adopted to the current surveillance system as proposed in my PhD proposal. 

## Files
- `data/`: Cleaned CSV datasets.
- `scripts/`: R analysis scripts.
- `plots/`: Generated correlation scatter plots.

---
*Author: Makiko Fujita*
