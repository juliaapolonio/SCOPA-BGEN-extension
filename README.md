## SCOPA (Software for COrrelated Phenotype Analysis)
### Overview
SCOPA (Software for COrrelated Phenotype Analysis) is a software that implements the reverse regression model to perform GWAS analysis of multiple correlated phenotypes. The mathematical model and its example usage is documented in the paper [SCOPA and META-SCOPA: software for the analysis and aggregation of genome-wide association studies of multiple correlated phenotypes](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1437-3). Further installation instructions and documentation can be found [here](https://genomics.ut.ee/en/tools). 

### BGEN Implementation
This project aims to adapt SCOPA to support the [BGEN format](https://www.well.ox.ac.uk/~gav/bgen_format/) to make use of the data from the UK Biobank (BGEN v1.2). It is primarily written in C++ and utilises files from the ALGLIB, TCLAP and BGEN libraries. 

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
SCOPA requires specification of input files - a genotype file (either BGEN or GEN) and a phenotype file (SAMPLE). Genotype file can be gzipped, if it has *.gz extension.

### Command line options
            ./SCOPA  [--debug] [--print_covariance] [--print_complex] [--betas]
            
            [--print_all] [--remove_missing] --pheno_name <string> ... 

            [--imp_threshold <double>] [--missing_phenotype <string>] [-e

            <string>] -o <string> -g <string> [--chr <int>] -s <string>

            [--] [--version] [-h]
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
1     Chromosome - chromosome of variant if set with --chr option. Otherwise 0

2     Position - position of variant

3     MarkerName - variant name

4     EffectAllele - effect allele (necessary for meta-analysis)

5     OtherAllele - non-effect allele (necessary for meta-analysis)

6     InfoScore -  Imputation quality measurement calculated similarly to IMPUTE2

7     HWE - p-value for HWE

8     MAF - minor allele frequency

9     N - samplesize

10    AA - genotype counts from imputed data

11    AB - genotype counts from imputed data

12    BB - genotype counts from imputed data

13    PhenotypeCount - number of phenotypes in model

14    Mask - binary mask showing the phenotypes used in current model (1-usd, 0-unused)

15    LogLikelihood - model likelihood

16    nullLogLikelihood - null model likelihood

17    LikelihoodRatio - likelyhood ratio

18    P-value - model p-value

19    BIC - Bayesian information score

20    BICnull - Bayesan iformation score for null model

21    Model - phenotypes in the order they were used in model (important for selecting covariance matrix for meta-analysis)

22    sortedModel - phenotypes in model in alphabetical order

23    beta_1 - effect size for phenotype 1

24    se_1 - stderr of effect for phenotype 1

25    beta_2 - effect size for phenotype 2

26    se_2 - stderr of effect for phenotype 2

27    beta_3 - effect size for phenotype 3

28    se_3 - stderr of effect for phenotype 3

29    cov_1_1 - inverted covariance matrix values

30    cov_1_2 - inverted covariance matrix values

31    cov_1_3 - inverted covariance matrix values

32    cov_2_2 - inverted covariance matrix values

33    cov_2_3 - inverted covariance matrix values

34    cov_3_3 - inverted covariance matrix values
