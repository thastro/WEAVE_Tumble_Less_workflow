#!/bin/bash

# Default values for input and output paths
input_path_base="/home/thajnik/Projects/WEAVE_workflows/tl/output/xmls"
output_path_base="/home/thajnik/Projects/WEAVE_workflows/tl/output/configure_output"

# Default epoch value
epoch=2024.0

# Default plate value
plate="A"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --epoch) epoch="$2"; shift ;;
        --plate) plate="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Update input and output paths based on the plate parameter
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

# Loop through each .xml file in the input directory
count=1
for file in "$input_path"/*.xml; do
    # Create the output filename by appending 'conf' before the .xml extension in the input file name
    filename=$(basename -- "$file")
    output_file="${output_path}/${filename%.xml}_${epoch}_conf.xml"
    #echo $count
    # Run the twilight_configure command for each file
    /home/thajnik/Software/configure_tl/configure/./twilight_configure --epoch $epoch -f "$file" -o "$output_file"
    ((count++))
done
