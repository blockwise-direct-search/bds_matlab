mkdir summary
find . -type f -name 'summary.pdf' -exec cp {} summary/ \;
cd summary
# Define keywords to search for PDF files
keywords=(
    "plain"
    "noisy_1_no_rotation"
    "noisy_2_no_rotation"
    "noisy_3_no_rotation"
    "noisy_4_no_rotation"
    "rotation_noisy_1"
    "rotation_noisy_2"
    "rotation_noisy_3"
    "rotation_noisy_4"
    "permuted"
    "linearly_transformed"
    "quantized"
    "perturbed_x0"
    "perturbed_x0_01"
    "perturbed_x0_10"
    "random_nan_5"
    "random_nan_10"
    "random_nan_20"
    "truncated_1"
    "truncated_2"
    "truncated_3"
    "truncated_4"
)

output_file="merged.pdf"
declare -a pdf_files  # Use an array to store PDF files

# Obtain all PDF files named 'summary.pdf' in the current directory
all_pdf_files=(summary.pdf)

# Print the array content for debugging
echo "Found these PDF files:"
for file in "${all_pdf_files[@]}"; do
    echo "  $file"
done

# Create an associative array to store keyword to file mapping
declare -A keyword_to_file

# Search for PDF files with keywords
for keyword in "${keywords[@]}"; do
    echo "Searching for keyword: $keyword"
    for file in "${all_pdf_files[@]}"; do
        if [[ $file == *"$keyword"* ]]; then
            echo "  Found file for keyword '$keyword': $file"
            keyword_to_file[$keyword]=$file
            break  # Every keyword should have only one file
        fi
    done
done

# Clear array to store PDF files in order of keywords
pdf_files=()

# Add PDF files to the array in order of keywords
for keyword in "${keywords[@]}"; do
    if [[ -n "${keyword_to_file[$keyword]}" ]]; then
        pdf_files+=("${keyword_to_file[$keyword]}")
    fi
done

# Print the array content for debugging
echo -e "\nFiles in order of keywords:"
printf '%s\n' "${pdf_files[@]}"

# Print total number of files found
echo -e "\nTotal files found: ${#pdf_files[@]}"

# Merge PDF files
if [[ ${#pdf_files[@]} -gt 0 ]]; then
    pdfunite "${pdf_files[@]}" "$output_file"
    echo "Merge successfully: $output_file"
else
    echo "There are no PDF files to merge."
    echo -e "\nAll PDF files in current directory:"
    ls summary*.pdf
fi