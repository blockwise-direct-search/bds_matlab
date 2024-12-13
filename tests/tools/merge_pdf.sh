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
    "random_nan"
    "truncated"
)

output_file="merged.pdf"
declare -a pdf_files  # Use an array to store PDF files

# Print all PDF files for debugging
echo "Found these PDF files:"
find . -maxdepth 1 -name "summary*.pdf" -type f | while read -r file; do
    echo "  $file"
done

# Create an associative array to store keyword to files mapping
declare -A keyword_to_files

# Function to add files to array in natural sort order
add_files_sorted() {
    local key=$1
    local files=("${@:2}")
    if [ ${#files[@]} -gt 0 ]; then
        # Convert array to newline-separated string, sort, and store
        keyword_to_files[$key]="$(printf "%s\n" "${files[@]}" | sort -V)"
    fi
}

# Search for PDF files with keywords and sort them
for keyword in "${keywords[@]}"; do
    echo "Searching for keyword: $keyword"
    case "$keyword" in
        "perturbed_x0")
            # Match perturbed_x0.pdf, perturbed_x0_01.pdf, perturbed_x0_10.pdf
            files=($(find . -maxdepth 1 -type f -name "summary*perturbed_x0*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        "random_nan")
            # Match random_nan_05.pdf, random_nan_10.pdf, random_nan_20.pdf
            files=($(find . -maxdepth 1 -type f -name "summary*random_nan*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        "truncated")
            # Match truncated_1.pdf through truncated_4.pdf
            files=($(find . -maxdepth 1 -type f -name "summary*truncated*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
        *)
            files=($(find . -maxdepth 1 -type f -name "summary*${keyword}.pdf" -o -name "summary*${keyword}_*.pdf"))
            add_files_sorted "$keyword" "${files[@]}"
            ;;
    esac
done

# Clear array to store PDF files in order of keywords
pdf_files=()

# Add PDF files to the array in order of keywords
for keyword in "${keywords[@]}"; do
    if [[ -n "${keyword_to_files[$keyword]}" ]]; then
        while IFS= read -r file; do
            if [[ -f "$file" ]]; then
                pdf_files+=("$file")
            fi
        done <<< "${keyword_to_files[$keyword]}"
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
    find . -maxdepth 1 -name "summary*.pdf" -type f
fi