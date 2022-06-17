# UK Biobank PRS Evaluation Tool

This repository contains the source code for the UK Biobank PRS Evaluation Tool developed by Genomics plc as part of
the Thompson et al publication https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1.

This tool will be made available via the [UK Biobank Research Analysis Platform](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform) and may be used to compare the performance of Polygenic Scores against the models published in Thompson et al.

You may only use this tool under the terms of the [licence](LICENCE).

## Introduction

`ukb-pret` is a software package that will evaluate and compare an input PRS against the scores and phenotypes 
published in the UK Biobank PRS Release.
In order to enable like-for-like performance evaluation for different PRSs for the same disease or trait, all evaluation metrics are calculated in the multi-ancestry UK Biobank PRS Release Testing subgroup.


## Usage

`ukb-pret` can be used by any registered user of the 
[UK Biobank Research Analysis Platform](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform) 
with access to the UK Biobank data and the UK Biobank PRS Release data. Simply load a project with access to these datasets, 
select the `ukb-pret` workflow from the "Tool Library" dropdown and follow the instructions to begin your analysis.

The workflow outputs a comprehensive PDF report that containing useful evaluation metrics and
plots describing the data. A CSV report of the metrics is also produced. 

### Running the tool locally

`ukb-pret` can be run and installed locally to compare any two PRS for UKB individuals, evaluated in the UK Biobank PRS Release Testing subgroup. 

#### Installation 

_Installation instructions to follow via pip_

#### Usage

Instructions on how to use `ukb-pret`'s `evaluate-prs` CLI tool can be viewed by typing `evaluate-prs --help` 
into your terminal: 

```
usage: evaluate-prs [-h] --prs-files PRS_FILES [PRS_FILES ...] --pheno-file
                    PHENO_FILE [--pcs-file PCS_FILE] [--output-dir OUTPUT_DIR]
                    [--list-phenotypes]

A CLI tool for evaluating a set of PRS against a set of phenotypes.

optional arguments:
  -h, --help            show this help message and exit
  --prs-files PRS_FILES [PRS_FILES ...]
                        Paths to two files, each containing Polygenic Risk
                        Score (PRS) and participant eIDs in CSV format.
                        Headers should be [eid,<data_tag>], where <data_tag>
                        is the field used to identify the PRS in the output
                        (REQUIRED)
  --pheno-file PHENO_FILE
                        Paths to a file containing phenotype data and
                        participant eIDs in CSV format. Headers should contain
                        at least [eid,<trait_code>,sex,in_ukb_wbu_testing],
                        where <trait_code> matches a Genomics plc phenotype
                        definition (type "evaluate-prs --list-phenotypes" for
                        supported phenotypes). Additionally, this file can
                        include the following columns to enable survival
                        analysis in binary traits: ['age_at_first_assessment',
                        'date_of_diagnosis', 'date_of_first_assessment',
                        'date_of_death', 'incident', 'sex',
                        'in_ukb_wbu_testing'] (REQUIRED)
  --pcs-file PCS_FILE   Path to a file containing the first 4 genetically
                        inferred principal components. Headers should be
                        [eid,pc1,pc2,pc3,pc4] (OPTIONAL) (when omitted,
                        evaluation is carried out across all ancestries & the
                        report does not contain a quality control section)
  --output-dir OUTPUT_DIR
                        Output directory for evaluation report and CSV
                        containing metrics (OPTIONAL) (default is current
                        working directory)
  --list-phenotypes     List supported Genomics plc phenotype definitions
```

### Contact

Please feel free to contact us at research+prs@genomicsplc.com
