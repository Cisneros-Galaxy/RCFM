# Rotation Curve Fitting Model

This repository is the working code for the Rotation Curve Fitting Model (RCFM), a method for interpreting galactic rotation curves without the need for dark matter or modifying gravity. This includes several sets of measured galactic rotation curves used to develop and test the model, and in the primary publications laying out the model.

### Organization

The main code for the model is in `Neros.py`. The files `DataAid.py` and `DataImporter.py` contain utilities related to reading the rotation curve data files. `Model.ipynb` is an example of using the model, this will eventually be simplified to require less "wrapper code" to read files and create plots.

The `data` directory contains the rotation curve data for multiple Milky Way models (`McGaugh` and `XueSofue`) and several collections of galaxies, including Sparc and Little Things. 

The `dev` directory is a collection of files used in developing and testing the model. It is a storage space for no longer used files.

The `documents` directory contains papers the authors referenced or found useful. The `documents/notes` subdirectory has python notebooks containing explanations of some parts of the model or techniques that were used.

The `fit-analysis` directory has the code (a notebook) used to determine some goodness of fit measures as part of the process. The output graphs appear in the paper.

The `graphs` directory is a placeholder. The `Model.ipynb` saves the graphs it generates here when run, but these aren't included in the repo. But you can reproduce them by running that notebook.

The `imported-data` is also a placeholder, the `Model.ipynb` saves a file of fit parameters there. We should rename this and make it an actual placeholder.

The `utils` directory contains some utility files. These aren't currently used by the main notebook or py files.

### References

This is the code powering the analysis in the papers:
- [2024] https://arxiv.org/pdf/2310.04372
- [2015] https://arxiv.org/abs/1506.04587 
- [2014] https://arxiv.org/abs/1407.7583



## To download this project

Go to where you'd like to clone the respoditory (a good place is documents/github) and type:
`git clone https://github.com/Cisneros-Galaxy/RCFM.git`

## To download latest changes, updates, new branches (etc.) onto your local machine

`git pull`

## Installing Dependancies

Before running anything, you need to install all of the dependencies. There are two ways to do this (do only one or the other):

### (1) Manually

[Install Python3](https://www.python.org/downloads/)  
 [Install Jupyter Notebook](http://jupyter.org/install)  
 [Install MatPlotLib](https://matplotlib.org/faq/installing_faq.html)

### (2) With Anaconda - RECOMMENDED

Anaconda is a platform that makes it easy to work with Python and Notebook. It comes with all the dependencies already installed! Just open it up, click "Launch" under "Jupyter Notebook" then navigate to the project's folder.  
 [Install Anaconda Python 3.6 Version](https://matplotlib.org/faq/installing_faq.html)

## Running the code

1. Open Anaconda
2. Launch Jupyter Notebook
3. Navigate to where the code is located
4. Click any of the .ipynb files

## Contributing

1. Request to be added as a collaborator to the RCFM project.


# Working with the model

BLOCK 1: Line 12
choose version of Neros file (composes vNeros function)
Block 3: Line 5
choose version of Neros file

##  choosing MW

Block 3:   line 5 & 6
 
## choosing sample

 Options are in Block 2
 
 input choice to Block 4, line 2.

## Editing Neros file (functions)

### (1) Manually

write your own function and get somebody to check it, then add it to the list of existing functions

### (2) use functions already established

hack the existing code

## Running the code

1. Open Anaconda
2. Launch Jupyter Notebook
3. Navigate to where the code is located
4. Click on the model.ipynb file
5. hit "run-all" or clear the cashe and run. 
6. record resulting table values
7. plot alpha correlation function

## Contributing

1. Got something really good? [criteria: functions which have good chisquare (<4), coefficients a & b (0.6-0.8), correlation fcn]
