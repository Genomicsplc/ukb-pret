# Generating Phenotypes in the RAP

This directory contains iPython notebooks to reproduce some phenotypes used to train and evaluate the PRS release models 
used to generate the PRS in the UK Biobank PRS Release.*

Details on these phenotypes can be found in the Thompson et al publication 
https://www.medrxiv.org/content/10.1101/2022.06.16.22276246v1.

These notebooks can be run from within a DNANexus project to generate phenotypes that are resolved to your 
unique participant identifiers.

*Please note that these notebooks many not recapitulate the exact phenotypes used in the UKB PRS Release if the backing 
data changes (for example if UKB update data in a new release, or if samples are withdrawn)

### Constructing your phenotype

Before attempting to run a pheconstructor file, please ensure that you have a DNANexus project prepared with the 
UKB data dispensed to a spark dataset.

1. Select the file corresponding to the trait of interest and convert this file to a Jupyter notebook using `jupytext`. 
This can be done via the command line with the command 
```jupytext --to ipynb <trait_code>_pheconstructor.py```
2. Upload the notebook to your DNANexus project with `dx upload <trait_code>_pheconstructor.ipynb`
3. Start a JupyterLab instance from the DNANexus web interface, ensuring that it is running on a Spark cluster
4. Open the notebook and replace the `'<dataset_id>.dataset'` string in the second cell with the 
ID of your UKB dataset
5. Run the notebook. This should output the phenotype data to a CSV file compatible with the `ukb-pret` tool*

*_For instruction on running your phenotype with `ukb-pret` in the DNANexus Research Analysis Platform, please refer to 
[project README](../README.md)_