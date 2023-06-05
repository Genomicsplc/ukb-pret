# UK Biobank PRS Evaluation Tool

This repository contains the source code for the UK Biobank PRS Evaluation Tool developed by Genomics plc as part of
the Thompson et al publication https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1.

This tool will be made available via the [UK Biobank Research Analysis Platform](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform)
and may be used to compare the performance of Polygenic Scores against the models published in
 [Thompson et al](https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1).

You may only use this tool under the terms of the [licence](LICENCE).

## Introduction

`ukb-pret` is a software package which compares the performance of two PRS in predicting an outcome phenotype.
This tool has been designed to evaluate an input set of PRS against a PRS selected from the UKB PRS release when
used as predictors for a set of input phenotype values for a relevant trait.


## Usage

`ukb-pret` can be used by any registered user of the 
[UK Biobank Research Analysis Platform](https://www.ukbiobank.ac.uk/enable-your-research/research-analysis-platform) 
with access to the UK Biobank data and the UK Biobank PRS Release data. Simply load a project with access to the 
UK Biobank data and the UK Biobank PRS Release data, and follow the steps below.

### Running the tool on the Research Analysis Platform (RAP)

#### Overview
In order to compare your own PRS against one from the UKB PRS release, you will need the following:
- A RAP project with dispensed data including the PRS release fields
- `dx` CLI tool installed
- docker to build the ukb-pret code for use in the RAP (with docker engine >v17.05)
- Your PRS scores in a csv file
- Your phenotype values in a csv file

Specific instructions on how to install/prepare each of these prerequisites can be found below

#### Installing the `dx` client

If you don't already have this installed, type `pip3 install dxpy` into the command 
line.

`dxpy` can also be installed 
[from a tarball](https://documentation.dnanexus.com/downloads#installation-of-java-r-and-c++-sdk-from-tarball) 
on Ubuntu, MacOS, Windows, CentOS, Red Hat Enterprise Linux and Scientific Linux 5/6/7 machines.

Details can be found on the 
[DNANexus website](https://documentation.dnanexus.com/downloads), 
and full documentation for `dxpy` can be found  [here](http://autodoc.dnanexus.com/bindings/python/current/dxpy.html).


#### Loading the UKB PRS Release Data

The first step in running `ukb-pret` on the RAP is to export the PRS of interest and other required data to a CSV file. 

This is done using the DNANexus tool `table-exporter`. A helper script is provided to run this tool and generate the
input CSV: 

1. Choose the PRS you are interested in comparing against from the list 
   [here](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=302) and note the "trait code" identifier.  For example, 
   the trait code for `"PRS for enhanced age at menopause (AAM)"` is `AAM`
2. Log into DNANexus by typing `dx login` to the command line, and select the project containing the UKB PRS dataset 
3. Execute the [`run_table_exporter.sh`](./run_table_exporter.sh) script. The script takes two arguments to identify
   your desired comparison PRS from step 1 and your Spark dataset ID. For example:

    ```shell script
    ./run_table_exporter.sh AAM app1234_20221008123456
    ```
    The pipeline should now run automatically, creating a CSV file named `table_exporter_output_<trait code>.csv` in your 
    project directory containing some data that are required to run `ukb-pret`


#### Generating the ukb-pret Docker image

The second step is to generate a Docker image containing an installation of the `ukb-pret` tool*.

To do this, enter the top-level parent [`ukb-pret`]('.') directory and run 

    docker build -t ukb-pret-docker -f docker/Dockerfile .
    
This will generate a Docker image named `ukb-pret-docker`, which should be converted to `.tar.gz` format and uploaded to 
the RAP for use with the `ukb-pret` DNANexus app.

Convert to a compressed tarball:

    docker save ukb-pret-docker | gzip > ukb-pret-docker.tar.gz

and upload to DNANexus platform 

    dx upload ukb-pret-docker.tar.gz

*_An installation of docker engine >v17.05 is required_

#### Prepare your input PRS

The third step is create a CSV file containing the PRS you wish to compare against the UKB PRS release, 
and upload it to your DNANexus project. 

This file should contain two columns `[eid,<data_tag>]`, where `<data_tag>` is the field used to 
identify the PRS in the output.

Next upload the file to DNANexus with `dx upload <my_prs_file>.csv`, replacing `<my_prs_file>` with your own CSV path.


#### Building & running the ukb-pret app

The final step is to build and run the `ukb-pret` app on the platform.

To build the app, go to the [`ukb-pret/app`](./app) directory and type

    dx build ukb-pret

The app should now be built and available to view on the DNANexus platform. 

Go to the online GUI, click on the app, and select your inputs through the GUI. These inputs should be: 

- The PRS file generated by the table-exporter app
- The PRS you wish to compare to the UKB PRS release (as per instructions in the 
[previous section](#prepare-your-input-prs))
- The phenotype file you wish to evaluate your PRS performance in (the headers should be be `[eid,<trait_code>]`, 
where `<trait_code>` can either refer to a UKB PRS release phenotype definition or to a user-defined trait. 
Additionally, this file can include the following columns to enable survival analysis in binary traits: 
`[age_at_first_assessment,date_of_diagnosis,date_of_first_assessment,date_of_death,incident]`)* 
- The docker tarball you uploaded to the platform

Click the "Run" button to commence the analysis. Once the job is finished, you will find a PDF named 
`prs_evaluation_report.pdf` containing useful plots, metrics and summary statistics that describe and compare the two 
input PRS. A CSV named `evaluation_metrics.csv` is also produced which can be used to create your own bespoke plots 
and figures from the raw data.

*Notebooks for generating some of the Enhanced PRS phenotypes can be found in the 
[pheconstructors directory](./pheconstructors) (see the [README](./pheconstructors/README.md) for instructions on how
to run these in the RAP)


### Running the tool locally

`ukb-pret` can be run and installed locally to compare the performance of two PRS in predicting a 
phenotype, provided the inputs conform to the formats defined below.

Please note that survival analysis will only be performed if the fields 
`[age_at_first_assessment,date_of_diagnosis,date_of_first_assessment,date_of_death,incident]` are provided in the 
phenotype file.

Ancestry-stratified analysis will only be performed if a valid principal components file is provided.

#### Installation 

_Installation instructions to follow via pip_

#### `ukb-pret` usage outside the RAP

The `ukb-pret` tool can also be used more generally to compare the performance of any two sets of PRS scores,
for prediction of a phenotype. In this feature, all inputs come directly from the user and are called using the `evaluate-prs` entrypoint.

Instructions on how to use `ukb-pret`'s `evaluate-prs` CLI tool can be viewed by typing `evaluate-prs --help` 
into your terminal: 

```
usage: evaluate-prs [-h] --prs-files PRS_FILES PRS_FILES --pheno-file
                    PHENO_FILE [--pcs-file PCS_FILE] [--output-dir OUTPUT_DIR]

A CLI tool for evaluating a set of PRS against a phenotype.

optional arguments:
  -h, --help            show this help message and exit
  --prs-files PRS_FILES PRS_FILES
                        Paths to two files, each containing Polygenic Risk
                        Score (PRS) and participant eIDs in CSV format.
                        Headers should be [eid,<data_tag>], where <data_tag>
                        is a field without spaces or special characters that
                        is used to identify the PRS in the output (REQUIRED)
  --pheno-file PHENO_FILE
                        Paths to a file containing phenotype data and
                        participant eIDs in CSV format. Headers should contain
                        at least [eid,<trait_code>], where <trait_code> is a
                        field without spaces or special characters that can
                        either correspond to an existing Gplc phenotype
                        definition or be defined by the user. [sex] can also
                        be included as a header for stratified analysis using
                        coding {0: female, 1: male}. Additionally, this file
                        can include the following columns to enable survival
                        analysis in binary traits: ['age_at_first_assessment',
                        'date_of_diagnosis', 'date_of_first_assessment',
                        'date_of_death', 'incident'] (REQUIRED)
  --pcs-file PCS_FILE   Path to a file containing the first 4 genetically
                        inferred principal components. Headers should be
                        [eid,pc1,pc2,pc3,pc4] (OPTIONAL) (when omitted,
                        evaluation is carried out across all ancestries & the
                        report does not contain a quality control section)
  --output-dir OUTPUT_DIR
                        Output directory for evaluation report and CSV
                        containing metrics (default is current working
                        directory) (OPTIONAL)
```

### Contact

Please feel free to contact us at research+prs@genomicsplc.com
