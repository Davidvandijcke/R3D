#!/bin/bash

# Define the output file
output_file="combined_R_code.R"

# Empty or create the output file
> "$output_file"

# Find all .R files recursively and concatenate their content
find . -type f -name '*.R' -exec cat {} \; >> "$output_file"

echo "All R code combined into $output_file"
