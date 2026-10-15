#!/bin/bash

#INPUT_DIR="/home/thajnik/Software/configure_tl/xmls/configure_output" # The first argument is the input directory containing XML files.
#OUTPUT_DIR="/home/thajnik/Software/configure_tl/xmls/splits_output" # The second argument is the output directory for split XML files.

INPUT_DIR="/home/thajnik/Projects/WEAVE_workflows/tl/output/configure_output" # The first argument is the input directory containing XML files.
OUTPUT_DIR="/home/thajnik/Projects/WEAVE_workflows/tl/output/xmls_split" # The second argument is the output directory for split XML files.



if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for input_file in "$INPUT_DIR"/*.xml; do
    filename=$(basename -- "$input_file")
    file_id="${filename%.*}"
    specific_output_dir="$OUTPUT_DIR"  # Adjusted to place all splits directly in the output directory

    # The Python script is updated to handle naming, so no changes needed here
    python xml_splitter.py -i "$input_file" -o "$specific_output_dir"
done

echo "Processing completed."
