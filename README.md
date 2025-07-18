# SCOPA (Software for COrrelated Phenotype Analysis) – BGEN Extension

SCOPA is a GWAS tool designed for the joint analysis of multiple correlated phenotypes using a reverse regression model. This repository extends SCOPA with support for the (BGEN format)[https://www.well.ox.ac.uk/~gav/bgen_format/], enabling the analysis of large-scale biobank data such as the **UK Biobank**.

Developed by Mägi *et al.*, the original method and application are detailed in the publication:

> [Mägi et al. (2017). SCOPA and META-SCOPA: software for the analysis and aggregation of genome-wide association studies of multiple correlated phenotypes.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1437-3)

---

## Features

- Reverse regression model for correlated phenotypes.
- Efficient parsing of **BGEN v1.2 genotype files**.
- Compatible with **UK Biobank** data.
- Support for meta-analysis via **META-SCOPA**.
- Outputs full likelihood model details and inverted covariance matrix.

---

## Installation

Clone this repository and compile the code:

```bash
git clone https://github.com/juliaapolonio/SCOPA-BGEN-extension.git
cd SCOPA-BGEN-extension
make
```

---

## Dependencies

The project uses:

- [**BGEN**](https://github.com/biobank-uk/BGEN) library for genotype parsing.
- [**TCLAP**](http://tclap.sourceforge.net/) for command-line argument parsing.
- [**ALGLIB**](https://www.alglib.net/) for matrix operations and statistics.

These dependencies are already in the cloned repository. Make sure headers and compiled libraries of the above are accessible during compilation (`make`).

---

## Input Files

- `--gen` / `-g` : Genotype file in **BGEN v1.2** format
- `--sample` / `-s` : Phenotype/sample file in plain text format

---

## Basic Usage

```bash
./SCOPA \
  --pheno_name <phenotype_1> \
  --pheno_name <phenotype_2> \
  -g data/genotypes.bgen \
  -s data/phenotypes.sample \
  -o results/output_prefix
```

---

### Common Options

| Option | Description |
|--------|-------------|
| `--pheno_name` | **Required** String. Phenotype(s) to include. Use multiple times for multivariate models |
| `-o, --out` | **Required** String. Specifies output directory |
| `-g, --gen` | **Required** String. Specifies genotype file |
| `-s, --sample` | **Required** String. Specifies sample file |
| `-e, --exclusion` | String. Specifies marker exclusion list |
| `--chr` | Integer. Specifies chromosome to be printed into chromosome column |
| `--missing_phenotype` | String. Specifies missing data value (default is NA) |
| `--remove_missing` | Binary. Removes samples with missing phenotype(s) (default is off) |
| `--imp_threshold` | Double. Minimum INFO score for imputed genotypes |
| `--betas` | Binary. Outputs per-phenotype betas and SEs (default is off) |
| `--print_complex` | Binary. Outputs full model with all phenotypes for META-SCOPA (default is off) |
| `--print_covariance` | Binary. Outputs covariance matrix (requires `--print_complex`; default is off) |
| `--print_all` | Binary. Print out all modes (default is off) |
| `--, --ignore_rest` | Binary. Ignores the arguments following this flag |
| `--debug` | Binary, enable debug mode (default is off) |
| `--version` | Binary. Displays version information and exits |
| `-h, --help` | binary. Displays usage information and exits | 

Run `./SCOPA --help` for full list of options.

---

## Output

The main output is a tab-delimited file with the following columns:

| Column | Description |
|--------|-------------|
| Chromosome | Chromosome of variant |
| Position | Position of variant |
| MarkerName | Variant name |
| EffectAllele | Effect allele (for meta-analysis)  |
| OtherAllele | Other allele (for meta-analysis) |
| InfoScore | Imputation quality measurement |
| HWE | p-value for HWE |
| MAF | minor allele frequency |
| N | Sample size |
| AA | Genotype counts from imputed data |
| AB | Genotype counts from imputed data |
| BB | Genotype counts from imputed data |
| PhenotypeCount | Number of phenotypes in model |
| Mask | Binary mask showing the phenotypes used in current model (1-used, 0-unused) |
| LogLikelihood | Model likelihood |
| nullLogLikelihood | null model likelihood |
| LikelihoodRatio | Likelyhood ratio |
| P-value | GWAS association p-value |
| BIC | Bayesian information score |
| BICnull | Bayesan iformation score for null model |
| Model | Phenotypes in the order they were used in model |
| sortedModel | phenotypes in model in alphabetical order |
| beta_* | Effect sizes |
| se_* | Standard errors |
| cov_\*_* | Inverted covariance matrix values |

You can check for sample outputs in `SAMPLE_SCOPA_OUTPUT`

## Example Dataset

Example command with dummy files:

```bash
./SCOPA \
--pheno_name pheno1 \
--pheno_name pheno2 \
-g SAMPLE_SCOPA_INPUT_FILES/cohort1_0X.bgen \
-s SAMPLE_SCOPA_INPUT_FILES/cohort1.sample \
-o $PWD/results
```

---

## Reference

If you use SCOPA in your research, please cite:

- **Mägi, R. et al. (2017)** SCOPA and META-SCOPA: software for the analysis and aggregation of genome-wide association studies of multiple correlated phenotypes. *Bioinformatics*, 33(15), 2314–2316. https://doi.org/10.1093/bioinformatics/btx153

- **Lagou, V. et al. (2025)** SCOPA software for analysis of correlated phenotypes identifies shared genetic loci between type 2 diabetes and colorectal cancer in the UK Biobank. *Under review*

---

## License

This software is distributed for academic use only. For other licensing, please contact the repository maintainer.

---

## Contact

For issues or contributions, please open an [issue](https://github.com/jiaq8877/SCOPA-BGEN-extension/issues) or reach out to the original developers.

