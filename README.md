# SCOPA (Software for COrrelated Phenotype Analysis): BGEN Implementation

SCOPA (Software for COrrelated Phenotype Analysis) is a software that implements the reverse regression model to perform GWAS analysis of multiple correlated phenotypes. The mathematical model and its example usage is documented in the paper SCOPA and META-SCOPA: software for the analysis and aggregation of genome-wide association studies of multiple correlated phenotypes (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1437-3).

The installation instructions and documentation can be found on: https://genomics.ut.ee/en/tools

This project aims to adapt SCOPA to support the BGEN format (https://www.well.ox.ac.uk/~gav/bgen_format/) to make use of the data from the UK Biobank (BGEN v1.2).

## Instructions
### Download
The code can be downloaded with `git clone https://github.com/jiaq8877/SCOPA-BGEN-extension.git`
Alternatively, the zipped version can be downloaded and unzipped.
### Installation
To compile SCOPA program, use command: 
`make` 

in the folder where files have been unpacked. The program can be run by typing: 
`./SCOPA
`
### Input files
SCOPA requires specification of input files - a genotype file (either BGEN or GEN) and a phenotype file (SAMPLE).

## To Be Updated
### Command line options
`./SCOPA  [--debug] [--print_covariance] [--print_complex] [--betas]

            [--print_all] [--remove_missing] --pheno_name <string> ... 

            [--imp_threshold <double>] [--missing_phenotype <string>] [-e

            <string>] -o <string> -g <string> [--chr <int>] -s <string>

            [--] [--version] [-h]`
Where: 
`   --debug
`        Debug mode on (default OFF)
        
`   --print_covariance`        
Print covariance matrix data for the model with all phenotypes. This is necessary for METASCOPA and can only be used with "--print_complex" option
        (default OFF)

### SCOPA output columns
  
