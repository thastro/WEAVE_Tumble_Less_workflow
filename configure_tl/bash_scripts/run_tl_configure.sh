#!/bin/bash

input_path_base="./xmls/plate"
output_path_base="./xmls/configure"


epoch=2024.0
plate="A"
hemisphere=""
fibre_status_file=""

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --epoch) epoch="$2"; shift ;;
        --plate) plate="$2"; shift ;;
        --hemisphere) hemisphere="$2"; shift ;;
        --fibre_status) fibre_status_file="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate inputs
if [[ -z "$fibre_status_file" ]]; then
    echo "Error: --fibre_status must be specified."
    exit 1
fi

if [[ ! -f "$fibre_status_file" ]]; then
    echo "Error: Fibre status file '$fibre_status_file' not found."
    exit 1
fi

if [[ -n "$hemisphere" && "$hemisphere" != "N" && "$hemisphere" != "S" ]]; then
    echo "Invalid hemisphere value: $hemisphere. Expected 'N' or 'S'."
    exit 1
fi

if [[ "$plate" == "A" ]]; then
    input_path="${input_path_base}_A"
    output_path="${output_path_base}_A"
elif [[ "$plate" == "B" ]]; then
    input_path="${input_path_base}_B"
    output_path="${output_path_base}_B"
else
    echo "Unknown plate value passed: $plate. Expected 'A' or 'B'."
    exit 1
fi

# Add hemisphere suffix if set
if [[ -n "$hemisphere" ]]; then
    input_path="$input_path/$hemisphere"
    output_path="$output_path"
fi

# Run fibre status command once before loop
echo "Running twilight_configure --fibre_status with: $fibre_status_file"


# Loop through each .xml file in the input directory
count=1
for file in "$input_path"/*.xml; do
    if [[ ! -f "$file" ]]; then
        echo "No XML files found in $input_path"
        exit 1
    fi

    filename=$(basename -- "$file")
    output_file="${output_path}/${filename%.xml}_${epoch}_conf.xml"

    echo "Processing $file -> $output_file"
    /home/thajnik/Software/configure_tl/./twilight_configure --fibre_status "$fibre_status_file" --epoch "$epoch" -f "$file" -o "$output_file"
    ((count++))
done
