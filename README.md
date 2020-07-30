# Rapid Online Assessment of Reading ability (ROAR)


## Measuring reading ability in the web-browser with a lexical decision task
Jason D. Yeatman, Kenny An Tang, Patrick M. Donnelly, Maya Yablonski, Maha Ramamurthy, Iliana I. Karipidis, Sendy Caffarra,  Megumi E. Takada, Michal Ben-Shachar, Benjamin W. Domingue

## Summary 
This analysis code agregates data from three versions of the online experiment (Stanford Child, UW Child, Stanford Adults) and generates figures and analyses for the manuscript.

## Repo Organization
* **data_UW and data_Stanford** .csv files for each subject that ran the experiment
* **data_allsubs** contains data organized in .csv files that has been concatenated across subjects and organized as [tidy data](https://r4ds.had.co.nz/tidy-data.html) tables.
  - **LDT_alldata_long.csv** All the data in long format. Each row is a trial of the LDT task. Raction time (rt), Accuracy (acc), stimulus (word), word length (worrdLength), whether it's a real/pseudoword (realpseuudo), and subject ID (subj) are coded for each trial.
  - **LDT_alldata_wide.csv** Accuracy data in wide format. Each row is a subject. Each column is a stimulus (word). Each entry is correct/incorrect for that stimulus/subject. Subject ID (subj) is coded in the first column. This data is also displayed in figs/WordMatrix.eps
  - **LDT_alldata_wide_sorted.csv** Same data as above but the rows are order based on each subject's percent correct and the columns are ordered based on each columns percent correct. See figs/WordMatrix.eps
  - **metadata_all.csv** Subject metadata including ages and reading scores. This needs to be filtered based on subject ID to match rows to the LDT data.
  - **wordStatistics.csv** Statistics on each stimulus such as word frequency, bigram and trigram statistics, etc. For definitions of word stats see: https://www.neuro.mcw.edu/mcword/definitions.html#UN1
* **analysis** Directory of analysis code 
  - **LDT_Run_All.sh** This shell script runs the whole analysis pipeline. In a shell type "source LDT_Run_All.sh" and it will run each step of the analysis pipeline and make figures.
  - **LDT_AnalyzeData.m** This MATLAB script concatenates all of the raw data files into a table and performs some initial analyses. Must run this first to add new data to the tables that every other function uses
  - **LDT_concatmetadata.m** MATLAB script to combine metadata downloaded from redcap into tidy tables concatenating UW and Stanford subjects and making sure subject id is coded properly
  - **analysis/LDT_IRT_Models.R** Script that does main analysis of item response data using varios item response theory models. Makes Figure 3 and outputs various csv files with subject ability estimates and word statistics


