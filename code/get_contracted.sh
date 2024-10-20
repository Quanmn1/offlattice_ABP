#!/bin/bash

# Usage: ./script.sh inputfile outputfile
name_all="pfaps_qsaps_test16_harmonic_metastability"
rf=0.07
rf_target=0.06

name="$name_all"_"$rf"
new_name="$name_all"_"$rf_target"
dir="$new_name"_video
if [ ! -d "$dir" ]; then
    mkdir "$dir"
fi

input_file="$name"_video/"$name"_last_state
output_file="$dir"/"$new_name"_last_state

# AWK script
awk 'NR==1 {print $0}  # Print the header (first line)
     NR>1 {            # For other rows
        $1 = $1 * 0.9  # Multiply first column by 0.9
        $2 = $2 * 0.9  # Multiply second column by 0.9
        print $0       # Print the modified line
     }' "$input_file" > "$output_file"
