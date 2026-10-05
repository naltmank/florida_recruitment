# Florida adult, recruit, and juvenile coral dynamics
Understanding the drivers of coral recruits and bottlenecks for succesful transitions to the juvenile stage.

Repo must be cloned as a .Rproj in order for pathing to work correctly. 

# Data
Data are organized into three folders: raw_data, clean_data, and summary_tables.

raw_data contain the original recruit, juvenile, adult density, and adult living tissue area data prior to any manipulations. They are cleaned and organized into long-form data in the script 01_data_cleaning. The resultant csvs are in the clean_data folder, which are the only ones used for analyses. One could theoretically skip running any of the 01_ script and jump straight to 02_ and fully reproduce all figures and analyses. 

summary_tables contain statistical outputs for creating clean tables (in excel) or as reference for constructing the SEM path diagrams.

# Code
Code is organized in order of the analyses that are being conducted in the manuscript, starting with multivariate analyses, then taxa-level recruit glmms, then taxa-level juvenile glmms, etc. 

All code uses data from the clean_data folder. glmm and sem scripts output csvs to the summary_tables folder. 
