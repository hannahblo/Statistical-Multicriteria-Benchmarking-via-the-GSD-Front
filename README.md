# Robust Statistical Comparison of Classification Algorithms

## Introduction
This anonymous repository contains R-code and data sets corresponding to the "Robust Statistical Comparison of Classification Algorithms"

The structure of the repository is as follows:
- File TODOOOO
- Folder R/ contains the needed functions (permutation test, sampling, definition of the constraint matrix etc) for the analysis.
- TODOOO Folder data/ 
- File _setup_session.R installs all (except for gurobi) needed R-packages.

The code was tested with
- R version 4.2.1
- R version 4.2.2

on

- Linux Ubuntu 20.04.5
- Windows 10 

## Setup
First, please install all necessary R-packages:
- For the computation of the linear programs, we used the R interface of gurobi optimizer, see [here](https://www.gurobi.com/) (accessed: 08.02.2023). This is a commercial
solver that offers a free academic licenses which can be found [here](https://www.gurobi.com/features/academic-named-user-license/) (accessed: 08.02.2023). To install this package, please follow the instructions there. A documentation can be found [here](https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf) (page 643ff) (accessed: 08.02.2023).
- Afterwards, please install all dependencies by sourcing the file _setup_session.R and source the files 

Then download the following files and save them in a folder named 'R':
- constraints_r1_r2.R
- sample_permutation_test.R

Afterwards, please download the data sets and save them in a folder named 'data':
TODOOO


In order to reproduce the papers' key results (and visualizations thereof) further download these scripts and save in respective folder:
- TODOOO (estimated runtime: TODOO days)

Running these three files sources the files in 'R/' and reproduces our result. Note that the folder 'R/', 'data/' and the three xxx_permutation_tests.R files must be in the same folder structure. Also make sure that the working directory of the R session links in the same folder.

## References to the data sets:
TODOOO

