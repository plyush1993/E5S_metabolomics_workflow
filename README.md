# E5S_metabolomics_workflow.R
E5S metabolomics workflow.R was used for fully processing peak table in .csv format. All .csv files were placed in a working directory, each first column - sample name, second - "Label".
To simplify the calculations, name the data tables according to:
- raw data (ds_raw.csv)
- ds_raw after EigenMS (dsr.csv)
- dsr after UVF+MVI (ds.csv)
- ds after ML+SFE (ds_d.csv)
- ds_d after RFE (ds_rfe.csv)
- raw data after UVF+MVI (ds_raw_ds.csv)
