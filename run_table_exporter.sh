#!/usr/bin/env bash

set -euf -o pipefail

prs_trait_code="${1?'Not enough args... Usage: ./run_table_exporter.sh PRS_TRAIT_CODE SPARK_DATASET_ID'}"
spark_dataset_id="${2?'Not enough args... Usage: ./run_table_exporter.sh PRS_TRAIT_CODE SPARK_DATASET_ID'}"

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
trait_code_map_file="${script_dir}/src/ukb_pret/data/ukb_enhanced_prs_codes.tsv"

if [[ ! -f "${trait_code_map_file}" ]]
then
    >&2 echo "cannot locate required trait code mapping file ${trait_code_map_file}"
    exit 1
fi

if ! dx ls "${spark_dataset_id%.dataset}.dataset" &> /dev/null
then
    >&2 echo "cannot find dataset for input SPARK DATASET_ID: ${spark_dataset_id%.dataset}"
    exit 1
fi

get_trait_code() {
    awk -v FS=$'\t' -v TARGET=$1 '{ if (TARGET == $1) { print $3 } }' "${trait_code_map_file}"
}

prs_field_title=$(get_trait_code "${prs_trait_code}")
if [[ -z "${prs_field_title}" ]]
then
    >&2 echo "could not resolve valid PRS field from input PRS_TRAIT_CODE: ${prs_trait_code}"
    exit 1
fi

dx run table-exporter -y \
  -idataset_or_cohort_or_dashboard="${spark_dataset_id%.dataset}.dataset" \
  -ioutput="table_exporter_output_${prs_trait_code}" \
  -ientity=participant \
  -ifield_titles="Participant ID" \
  -ifield_titles="Sex" \
  -ifield_titles="${prs_field_title}" \
  -ifield_titles="PRS genetic principle components | Array 0" \
  -ifield_titles="PRS genetic principle components | Array 1" \
  -ifield_titles="PRS genetic principle components | Array 2" \
  -ifield_titles="PRS genetic principle components | Array 3"
