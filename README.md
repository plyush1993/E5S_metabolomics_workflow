# E5S_metabolomics_workflow.R
Easy 5 Steps metabolomics workflow

E5S metabolomics workflow.R was used for fully processing peak table in .csv format. All .csv files were placed in a working directory, each first column - sample name, second - "Label", other - numeric values or NA. Row and column names should be presented.

raw data after integration.csv is data table of own dataset after integration.

To simplify the calculations, name the data tables according to:
- raw data (ds_raw.csv)
- ds_raw after EigenMS (dsr.csv)
- dsr after UVF+MVI (ds.csv)
- ds after ML+SFE (ds_d.csv)
- ds_d after RFE (ds_rfe.csv)
- raw data after UVF+MVI (ds_raw_ds.csv)
