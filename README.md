# calculate_effeX
Calculate effect size precursors to effect size explorer

## master.R includes all steps
To summarize...
1. Clean the data with clean_data.R
2. Perform QC on the cleaned data with qc.R (produces an html report to visually inspect matrices)
3. Calculate effect sizes with calc_d.R
4. Estimate simultaneous confidence intervals with estimate_simci.R
