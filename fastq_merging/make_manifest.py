import os
import json
import re
import argparse
import pandas as pd

def create_json(input_folder, output_json, sample, ena_submission_file):
    # Template
    template = {
        "study": "PRJEB97253",
        "sample": "",
        "name": "",
        "platform": "ILLUMINA",
        "library-source": "SYNTHETIC",
        "library_selection": "other",
        "libraryStrategy": "WGS",
        "fastq": []
    }
    
    # Get all .fastq.gz files
    fastq_files = sorted([
        f for f in os.listdir(input_folder)
        if f.startswith(sample)
    ])

    if not fastq_files:
        raise ValueError("No .fastq.gz files found in input folder.")

    # Infer sample name from first file (everything before .R1/.R2/.R*)
    match = re.match(r"(.+)\.R[12]\.fastq\.gz", fastq_files[0])
    if match:
        sample_name = match.group(1)
    else:
        sample_name = os.path.splitext(fastq_files[0])[0]

    # Fill in name and sample
    template["name"] = sample


    ena_runs=pd.read_csv(ena_submission_file)

    template["sample"] = str(ena_runs.loc[ena_runs["title"] == sample, "id"].iloc[0])
    # Add all FASTQ files
    for f in fastq_files:
        template["fastq"].append({
            "value": f,
            "attributes": {
                "read_type": "paired"
            }
        })

    # Write JSON
    with open(output_json, "w") as outfile:
        json.dump(template, outfile, indent=2)

    print(f"✅ JSON created: {output_json}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate JSON from FASTQ files.")
    parser.add_argument("-i", "--input", required=True, help="Input folder containing FASTQ files.")
    parser.add_argument("-o", "--output", default="sample_metadata.json", help="Output JSON file name.")
    parser.add_argument("-s", "--sample", required=True, help="Sample name, fastq prefix"),
    parser.add_argument("-ef", "--ena_submission_file", required=True, help="Runs report on ENA")
    args = parser.parse_args()

    create_json(args.input, args.output, args.sample, args.ena_submission_file)
