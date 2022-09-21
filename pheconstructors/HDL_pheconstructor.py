# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # HDL Pheconstructor
#
# This notebook contains code to reproduce the HDL phenotype used in the UK Biobank PRS release.
#
# To reproduce the phenotype:
#
# 1. Convert this file to an Jupyter notebook using jupytext
# 2. Upload the notebook to your DNANexus project with `dx upload HDL_pheconstructor.ipynb`
# 3. Start a JupyterLab instance from the DNANexus web interface, ensuring that it is running on a Spark cluster
# 4. Open the pheconstructor notebook and replace the `<dataset_id>.dataset` string in the second cell with the ID of your UKB dataset
# 5. Run the notebook - this should output the phenotype data to a CSV file compatible with the UKB-PRET tool

# +
import pyspark
import dxpy
import dxdata
from datetime import datetime

from typing import List
import pandas
import numpy
# -

# Load the database into pyspark and dxdata
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
dataset = dxdata.load_dataset('<dataset_id>.dataset')

# +
# Load the UKB data fields corresponding to HDL

hdl_ukb_code = '30760'
[participant_data] = [e for e in dataset.entities if e.name == 'participant']
codes = [f.name for f in participant_data.fields if f.name.startswith('p' + hdl_ukb_code)]
names = ['eid'] + codes

hdl_data = participant_data.retrieve_fields(names=names, engine=dxdata.connect())
hdl_df = hdl_data.toPandas()
hdl_df = hdl_df.set_index('eid')

# -

# Remove any instances where all instances for an individual are NaN
hdl_df = hdl_df.dropna(how='all')

# HDL has two instances
# Take the value of the first instance where available, and the second instance otherwise
hdl_df['HDL'] = hdl_df['p30760_i0'].fillna(hdl_df['p30760_i1'])

# Select the HDL column only and write the phenotype to CSV
out_df = hdl_df['HDL'].to_frame()
out_df.to_csv('HDL_phenotype.csv')
