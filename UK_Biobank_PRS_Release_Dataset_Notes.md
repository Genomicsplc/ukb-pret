# UK Biobank PRS Release

This document contains notes regarding the UK Biobank PRS Release data

## Dataset release notes, issues and clarifications

- PC centering issues
- Missing enhanced PRS for CD and UC

## Original Release Notes

This Return-of-Results package relates to the following paper:
Selzam et al. (2022). UK Biobank release and systematic evaluation of polygenic risk scores for 53 diseases and
quantitative traits. Submitted.

The results data were generated under UK Biobank project number 9659. This project originally focused on methods for
studying pleiotropy (a phenomenon where genetic variants influence multiple human traits), but the project was
subsequently granted an extension to include more clinical data and aims, including “providing information relevant
to improving the efficacy and safety of therapeutic interventions that are needed to improve human health and healthcare”.

### UK BIOBANK PRS RELEASE BRIEF DESCRIPTION
Our paper describes the UK Biobank PRS Release: a set of polygenic scores for 53 disease and quantitative traits,
together with allied phenotypic and genetic information.

A separate tool (UKB-PRET) is provided to facilitate evaluation in a testing subset of the UK Biobank. See
https://github.com/Genomicsplc/ukb-pret

PRS scores were generated by the application of a modified version of LDpred
(Vilhjálmsson et al 2015, Am. J. Hum. Genet. 97, 576–592) to meta-analysed summary statistic GWAS data,
obtained either entirely from external GWAS data (the Standard PRS set) or from a combination of external and internal
UKB data (the Enhanced PRS set). A subsequent PC-based ancestry centering step
(Khera et al 2019, Circulation 139, 1593–1602) was applied to approximately centre the score distributions on
zero across all ancestries, although it should be noted that residual deviations of mean(PRS) from zero were still
observed for some traits in some ancestries. Score distributions were also standardised to have approximately unit
variance within ancestry groups, as determined by a geometric inference in PC space. For further details, please refer
to the paper.

### UK BIOBANK PRS RELEASE DATA FILES
The Release is packaged as a series of data files, with fields as described below. For mapping of phenotype codes to
phenotypes, please refer to the Binary_trait_mapping.csv and Quantitative_trait_mapping.csv.

### DATA FILES
Common data file (PCs, ancestry): pc_ancestry_info_ukb.csv  
Phenotype files (one for each trait):  
e.g. phenotype_info_ukb_BC.csv  
PRS files (one for each trait and PRS type (Enhanced (also known as UKB-WBU) and Standard (also known as UKB-Free))):  
e.g. prs_values_UKB_WBU.csv, prs_values_UKB_Free.csv

### COMMON DATA FILE CONTENT DESCRIPTION
eid: participant ID as provided to Genomics plc under project number 9659