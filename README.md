# E5S_metabolomics_workflow.R
E5S metabolomics workflow.R was used for fully processing peak table in .csv format. All .csv files were placed in a working directory, each first column - sample name, second - "Label".
To simplify the calculations, name the data tables according to:
- raw data (ds_raw)
- ds_raw after EigenMS (dsr)
- dsr after UVF+MVI (ds)
- ds after ML+SFE (ds_d)
- ds_d after RFE (ds_rfe)
