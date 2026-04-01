# Seasonality of anti-rotavirus titers in maternal sera in a cohort in Bamako, Mali

Administered by Marina Antillon and Meagan C. Fitzpatrick  
School of Medicine  
University of Maryland  

------

# Project Objective 

The project seeks to understand how anti-rotavirus titers in mothers at the time of delivery may fluctuate through the year in during the pre-vaccine era in Mali. The sera of 329 mothers who gave birth in 2012 was assayed for anti-rotavirus IgA and IgG titers by ELISA at the Laboratory for Specialized Clinical Studies, Cincinnati Children’s Hospital Medical Center, Cincinnati, Ohio. Concurrently, the Global Enteric Multicenter Study-1a (GEMS-1a) collected data that allowed us to calculate the incidence study of diarrheal disease attributable to rotavirus. We therefore undertook a descriptive statistical analysis to detect seasonality in the maternal titers and we tested the cross-correlation with various time-lags between the fluctuations in titers and incidence. 

Maternal flu vaccine in 2012: https://clinicaltrials.gov/study/NCT01430689; https://doi.org/10.1016/S1473-3099(16)30054-8.
Global Enteric Multicenter Study-1a in 2012: https://doi.org/10.1016/S2214-109X(19)30076-2; https://doi.org/10.1016/S2214-109X(19)30541-8.

------

# Computational considerations

## Software considerations
The tools for the economic model are coded in R. While we would highly recommend using the code within RStudio environment (in part because of it's features to manage the project with .RProj and renv) this is not strictly necessary and the benefits of renv are available from a classic R interface or a shell command line.

For detailed information, see 'Installation to-do list' below.

## Hardware considerations

For reference, a single run of the code for all strategies in one health zone takes about 2-3 minutes in a MacBook Pro (2020 model) with an Apple M1 processor and 16 GB of RAM. 

## Installation to-do list (all free)

### R
- [R](https://www.r-project.org) (required)
- [RStudio](https://www.rstudio.com/) (highly recommended though not required; it integrates well with the other management tools used in the project. The free version is more than enough.)

### R Packages
The packages to run the whole analysis are available in a private repository constructed via renv.

In order to install the packages from the private library, open R and change the working directory to the project directory final_code.

Then type in renv::restore().

------

# Overview of analysis 

The analysis is broadly defined in three parts, each of these parts is executed in various R code files. Because the titer data cannot be posted in this repository, a fourth part is included at the beginning to create the synthetic data that is provided with this repository.
I. Creating synthetic data
II. Mandolo data re-analysis
III. Validation of cross-correlation functions
IV. Analysis of data

## Files

`./data/synthetic_rotavirus_titers.csv` - the data that would stand in the place of the real data. 
`./data/epi_desc_1a_lsd.csv` - the incidence data from the GEMS-1a study for less severe diarrhea.
`./data/epi_desc_1a_msd.csv` - the incidence data from the GEMS-1a study for moderate-to-severe diarrhea.
`./data/Mandolo_data_new.xlsx` - data from (CITE) Figure 2a digitized by automeris.io. We used these data to calculate the odds-ratio of vaccine sero-conversion based on a doubling of IgG titers in the child before vaccination.

`./functions/` - various useful functions for the analysis, primarily the portion included in the Rmarkdown.

`./00a_IgA and IgG synthetic data.R` - the code that makes `./data/synthetic_rotavirus_titers.csv`. The real data would be available for intrested researchers by the University of Maryland by emailing Prof. Fitzpatrick at meagan.fitzpatrick@som.umaryland.edu.
`./00b_Mandolo_baby_data.R` - the code that re-analyses the data digitized from [Mandolo et al 2025](https://doi.org/10.1371/journal.pmed.1004734). The results are reported in the supplement.
`./01a_cross_corr_interpretation.R` - Basic figure showing interpretation for Fig X. Not real data or analysis.
`./01b_cross_corr_comp_validation.R` - Fig 1 in the paper.
`./02_cross_corr_catalytics_analysis.R` - analysis that gives Fig 2c in the paper. This portion MUST be run before `03_titers_moms_2batches_word.Rmd`. 
`./03_titers_moms_2batches_word.Rmd` - Rmd analysis pulling together the whole analysis for all tables and figures in the manuscript. The Rmarkdown document is meant to be used with Rstudio and it outputs a Word document. When you click the **Knit** button at the top left, a Word document will be generated that includes both content as well as the output of any embedded R code chunks within the document. One may run individual chunks by hitting the green right-facing triangle on the top right of each chunk. By default, the R code is excluded from the Rmarkdown, as given by the first line of code of the following chunk. Note that the `echo = FALSE` or `echo = TRUE` parameter can be added to each chunk specifically to prevent or dis/enable printing of the R code that generated the table or plot.

## Troubleshooting issues with the renv repository
If there is an issue with the repository, try typing `renv:: rinit()` and when prompted type the number 1, for Restored the project from lockfile.
Adding a package to the renv library: just use the command `install.packages()`. Note that doing this won't install the package for use with other projects. Then update the lockfile by typing `renv::snapshot()`.
