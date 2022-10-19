#!/bin/bash
set -euf -o pipefail

main() {

    # Download the inputs
    dx download "$ukb_release_prs_file"
    dx download "$user_prs"
    dx download "$pheno_file"
    dx download "$ukb_pret_docker"

    # Replace the name of the fixed headers in the table-exporter output file
    sed -i "s/p26201_a0/pc1/;s/p26201_a1/pc2/;s/p26201_a2/pc3/;s/p26201_a3/pc4/;s/p31/sex/;s/Male/1/g;s/Female/0/g" $ukb_release_prs_file_name

    # Load the docker image
    docker load < $ukb_pret_docker_name

    # Run the docker image with /home/dnanexus mounted
    docker run -v /home/dnanexus:/data/dnanexus/ ukb-pret-docker evaluate-prs-rap \
     --ukb-release-prs-file /data/dnanexus/$ukb_release_prs_file_name \
     --user-prs-file /data/dnanexus/$user_prs_name \
     --pheno-file /data/dnanexus/$pheno_file_name \
     --output-dir /data/dnanexus

    # Upload the results to the Project directory
    pdf_report=$(dx upload prs_evaluation_report.pdf --brief)
    csv_output=$(dx upload evaluation_metrics.csv --brief)
    dx-jobutil-add-output pdf_report "$pdf_report" --class=file
    dx-jobutil-add-output csv_output "$csv_output" --class=file
}
