# UK Biobank PRS Release

This document contains notes regarding the UK Biobank PRS Release data

## Dataset release notes, issues and clarifications

- All PRS in the release are meant to be centered. We note however that the prostate cancer PRS (Standard version) has suboptimal centering properties. This issue will be fixed in the next release.
- Both Crohn's disease and ulcerative colitis are primarily based on a relatively sparse ImmunoChip array. This sparsity complicated the meta-analysis process, leading to less effective Enhanced PRS. For this reason, the Enhanced version of the PRS is not provided for either of these two diseases. Please use the Standard versions in lieu of the Enhanced versions.
- Motivated by the evaluation of the enhanced PRS set, which uses the
  White British Unrelated group as training, we defined our testing
  set to exclude individuals related to the White British Unrelated
  group. This process resulted in an evaluation set of 104,621
  individuals. However, we noted that this process was not perfectly
  applied and out of the 104,621 individuals in testing, 9K show a
  relatedness level above 0.2 (and 28K with greater than third degree
  relatedness). We compared PRS performance in individuals with and
  without relatedness to training and this did not measurably impact
  the evaluation (see
  ![](img/excluded_vs_only_related_28k_binary.png?raw=true | width=250)). 
  This will be fixed in the next release. Meanwhile, users keen to investigate further can identify separately individuals related to the WBU set (using the relatedness data field provuded by Bycroft et al).


## Original Release Notes

This Return-of-Results package relates to the following paper:
Thompson et al. 2022. UK Biobank release and systematic evaluation of optimised polygenic risk scores for 53 diseases and
quantitative traits. https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1

The results data were generated under UK Biobank project number 9659 (under its extension to therapeutic applications).

### UK BIOBANK PRS RELEASE BRIEF DESCRIPTION
Our paper describes the UK Biobank PRS Release: a set of polygenic scores for 53 disease and quantitative traits,
together with associated phenotypic and genetic information (see below for list of traits).

A separate tool (UKB-PRET) is provided to facilitate evaluation in a testing subset of the UK Biobank. See
https://github.com/Genomicsplc/ukb-pret

The polygenic risk scores (PRSs) are described in Thompson et al. 2022 (https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1).
Data supporting these scores were obtained either entirely from external GWAS data (the Standard PRS set) or from a combination of external and internal UK Biobank data (the Enhanced PRS set). This work was conducted by Genomics PLC under UK Biobank project 9659."

### UK BIOBANK PRS RELEASE DATA FILES
The Release is packaged as a series of data files, with fields as described below.

#### DATA
PRS genetic principal components – field 26201  
In UKB WBU PRS testing set – field 26200  
PRS fields:  
One for each trait and PRS type (Category 302 - Enhanced (also known as UKB-WBU) and Category 301 - Standard (also known as UKB-Free)))

#### COMMON DATA CONTENT DESCRIPTION
PRS genetic principal components: the first four PC axes obtained from the common genotype data in the 1000 genomes
reference dataset using standard methods (Price et al, 2006) are provided in the four array indices in field 26201.

In UKB WBU PRS testing set: indicator of participants in the UK Biobank PRS Release testing set

#### PRS FILES (Enhanced (also known as UKB-WBU)) CONTENT DESCRIPTION
prs: polygenic score for the trait for every participant in the UK Biobank PRS Release Testing subgroup, as described in the paper

#### PRS FILES (Standard (also known as UKB-Free)) CONTENT DESCRIPTION
prs: polygenic score for the trait for all UK Biobank participants, as described in the paper

#### PHENOTYPE CODES
The following is a list of traits and their phenotype codes (as used in file naming).
See also [src/ukb_pret/data/traits.yaml](./src/ukb_pret/data/traits.yaml)

DISEASE TRAITS  
Age-related macular degeneration	AMD  
Alzheimer's disease	AD  
Asthma	AST  
Atrial fibrillation	AF  
Bipolar disorder	BD  
Bowel cancer	CRC  
Breast cancer	BC  
Cardiovascular disease	CVD  
Coeliac disease	CED  
Coronary artery disease	CAD  
Crohn's disease	CD  
Epithelial ovarian cancer	EOC  
Hypertension	HT  
Ischaemic stroke	ISS  
Melanoma	MEL  
Multiple sclerosis	MS  
Osteoporosis	OP  
Prostate cancer	PC  
Parkinson's disease	PD  
Primary open angle glaucoma	POAG  
Psoriasis	PSO  
Rheumatoid arthritis	RA  
Schizophrenia	SCZ  
Systemic lupus erythematosus	SLE  
Type 1 diabetes	T1D  
Type 2 diabetes	T2D  
Ulcerative colitis	UC  
Venous thromboembolic disease	VTE

QUANTITATIVE TRAITS  
Age at menopause	AAM  
Apolipoprotein A1	APOEA  
Apolipoprotein B	APOEB  
Body mass index	BMI  
Calcium	CAL  
Docosahexaenoic acid	DOA  
Estimated bone mineral density T-score	EBMDT  
Estimated glomerular filtration rate (creatinine based)	EGCR  
Estimated glomerular filtration rate (cystatin based)	EGCY  
Glycated haemoglobin	HBA1C_DF  
High density lipoprotein cholesterol	HDL  
Height	HEIGHT  
Intraocular pressure	IOP
