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
Print covariance matrix data for the model with all phenotypes. This is necessary for METASCOPA and can only be used with "`--print_complex`" option (default OFF)

`   --print_complex
`
Print only the model with all phenotypes. These ful models can be meta-analysed with METASCOPA (default OFF)

`   --betas
`
Print each phenotype's effect size and stderr info of all selected models into separate output file (default OFF)

`   --print_all
`
Print out all models (default OFF)

`   --remove_missing
`
Remove sample if any of the phenotype values is missing. This is necessary if you want to compare models based on BIC scores (default OFF)

`   --pheno_name <string>  (accepted multiple times)
`
**(required)**  Name of phenotype to use (use this command multiple times i.e. --pheno_name BMI --pheno_name HEIGHT etc.)

`   --imp_threshold <double>
`
Imputation quality threshold (default 0)

`   --missing_phenotype <string>
`
This specifies missing data value (default NA)

`   -e <string>,  --exclusion <string>
`
This specifies marker exclusion list

`   -o <string>,  --out <string>
`
**(required)**  This specifies output root

`   -g <string>,  --gen <string>
`
**(required)**  This specifies genotype file.

`   --chr <int>
`
This specifies chromosome to be printed into chromosome column

`   -s <string>,  --sample <string>
`
**(required)** This specifies sample file

`   --,  --ignore_rest
`
Ignores the rest of the labeled arguments following this flag

`   --version
`
Displays version information and exits

`   -h,  --help
`
Displays usage information and exits
            
### SCOPA output columns
  
