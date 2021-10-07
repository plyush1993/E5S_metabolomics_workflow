# E5S_metabolomics_workflow.R
Easy 5 Steps metabolomics workflow

E5S metabolomics workflow.R was used for fully processing peak table in .csv format. All .csv files were placed in a working directory, each first column - sample name, second - "Label", other - numeric values or NA. Row and column names should be presented.

ds_raw.csv is data table of own dataset after integration and HM MVI.

To simplify the calculations, name the data tables according to:
- raw data (ds_raw.csv)
- ds_raw after EigenMS (dsr.csv)
- dsr after UVF+MVI (ds.csv)
- ds after ML+SFE (ds_d.csv)
- ds_d after RFE (ds_rfe.csv)
- raw data after UVF+MVI (ds_raw_ds.csv)

## Citation:
[Plyushchenko Ivan, et al. "An approach for feature selection with data modelling in LC-MS metabolomics." Analytical Methods 12.28 (2020): 3582-3591.](https://pubs.rsc.org/en/content/articlelanding/2020/ay/d0ay00204f#!divAbstract)

## Dataset Ref:
Metabolomics Workbench [project PR000857, study ST001271](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000857)

## Contact:
Please send any comment, suggestion or question you may have to the author (Mr. Ivan Plyushchenko), email: plyushchenko.ivan@gmail.com.
