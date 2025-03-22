#!/bin/bash

# Usage: ./script.sh inputfile outputfile
name_all="pfaps_qsaps_test17_solids"
rf=0.05
rf_target=0.055

name="$name_all"_"0.05"
new_name=pfaps_qsaps_test17_solids_"$rf_target"
dir="$new_name"_video
if [ ! -d "$dir" ]; then
    mkdir "$dir"
fi

input_file="$name"_video/"$name"_last_state
output_file="$dir"/"$new_name"_last_state

# AWK script
awk 'NR==1 {print $0}  # Print the header (first line)
     NR>1 {            # For other rows
        $1 = $1 * 1.1
        $2 = $2 * 1.1
        print $0       # Print the modified line
     }' "$input_file" > "$output_file"
