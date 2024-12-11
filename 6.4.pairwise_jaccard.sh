#!/bin/bash

for i in .; do

sort -t, -nk 7 $i
echo This is the value 


# Check if the input folder is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <bedfiles_folder> <output_csv>"
    exit 1
fi

BED_FOLDER="$1"
OUTPUT_CSV="$2"
DEBUG_FILE="jaccard_debug.txt"

# Check if the folder exists
if [[ ! -d "$BED_FOLDER" ]]; then
    echo "Error: Folder '$BED_FOLDER' does not exist."
    exit 1
fi

# Initialize the CSV file with a header
echo "File1,File2,Intersection,Union,JaccardIndex,NumIntersections" > "$OUTPUT_CSV"

# Initialize the debug text file
echo "Jaccard Index Raw Results" > "$DEBUG_FILE"

# Get the list of BED files
FILES=("$BED_FOLDER"/*.bed)

# Loop through all possible pairwise combinations of files
for ((i=0; i<${#FILES[@]}; i++)); do
    for ((j=i+1; j<${#FILES[@]}; j++)); do
        FILE1="${FILES[$i]}"
        FILE2="${FILES[$j]}"
        
        # Get file names without path for logging
        FILE1_NAME=$(basename "$FILE1")
        FILE2_NAME=$(basename "$FILE2")
        
        # Calculate Jaccard Index and log the raw output
        JACCARD_RESULT=$(bedtools jaccard -a "$FILE1" -b "$FILE2" 2>/dev/null)
        echo -e "Comparing: $FILE1_NAME vs $FILE2_NAME\n$JACCARD_RESULT\n" >> "$DEBUG_FILE"
        
        # Parse Jaccard output
        INTERSECTION=$(echo "$JACCARD_RESULT" | awk 'NR==2 {print $1}')
        UNION=$(echo "$JACCARD_RESULT" | awk 'NR==2 {print $2}')
        JACCARD=$(echo "$JACCARD_RESULT" | awk 'NR==2 {print $3}')
        NUM_INTERSECTIONS=$(echo "$JACCARD_RESULT" | awk 'NR==2 {print $4}')

        # Handle missing or invalid values
        if [[ -z "$INTERSECTION" || -z "$UNION" || -z "$JACCARD" || -z "$NUM_INTERSECTIONS" ]]; then
            INTERSECTION=0
            UNION=0
            JACCARD=0
            NUM_INTERSECTIONS=0
            echo "Parsing error: $FILE1_NAME vs $FILE2_NAME" >> "$DEBUG_FILE"
        fi

        # Write results to the CSV
        echo "$FILE1_NAME,$FILE2_NAME,$INTERSECTION,$UNION,$JACCARD,$NUM_INTERSECTIONS" >> "$OUTPUT_CSV"
        echo "Processed: $FILE1_NAME vs $FILE2_NAME"
    done
done

echo "Pairwise Jaccard Indexes saved to '$OUTPUT_CSV'"
echo "Debug log saved to '$DEBUG_FILE'"

