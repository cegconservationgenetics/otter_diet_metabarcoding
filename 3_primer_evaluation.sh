#!/bin/bash

# Exit on error
set -e

# Define primer pairs and their expected lengths for each gene
declare -A PRIMERS_12S=(
    ["12S_12SV5F_12SV5R"]="TTAGATACCCCACTATGC,TAGAACAGGCTCCTCTAG"
    ["12S_tele02F_tele02R"]="AAACTCGTGCCAGCCACC,GGGTATCTAATCCCAGTTTG"
    ["12S_teleoF_teleoRdeg"]="ACACCGCCCGTCACTCT,CTTCCGGTACACTTACCRTG"
    ["12S_MiFish-U-F_MiFish-U-R"]="GTCGGTAAAACTCGTGCCAGC,CATAGTGGGGTATCTAATCCCAGTTTG"
)

declare -A PRIMERS_16S=(
    ["16S_Vert-16S-eDNA-F1_Vert-16S-eDNA-R1"]="AGACGAGAAGACCCYDTGGAGCTT,GATCCAACATCGAGGTCGTAA"
    ["16S_MarVer3F_MarVer3R"]="AGACGAGAAGACCCTRTG,GGATTGCGCTGTTATCCC"
)

declare -A PRIMERS_COI=(
    ["COI_VF2_FishR1"]="TCAACCAACCACAAAGACATTGGCAC,TAGACTTCTGGGTGGCCAAAGAATCA"
    ["COI_coi.175f_coi.345r"]="GGAGGCTTTGGMAAYTGRYT,TAGAGGRGGGTARACWGTYCA"
    ["COI_L2513_H2714"]="GCCTGTTTACCAAAAACATCA,CTCCATAGGGTCTTCTCGTCTT"
)

# Define expected amplicon size ranges (±10bp) for each gene
declare -A AMPLICON_SIZES_12S=(
    ["12S_12SV5F_12SV5R"]="80,105"
    ["12S_tele02F_tele02R"]="140,200"
    ["12S_teleoF_teleoRdeg"]="50,100"
    ["12S_MiFish-U-F_MiFish-U-R"]="150,200"
)

declare -A AMPLICON_SIZES_16S=(
    ["16S_Vert-16S-eDNA-F1_Vert-16S-eDNA-R1"]="150,260"
    ["16S_MarVer3F_MarVer3R"]="100,250"
)

declare -A AMPLICON_SIZES_COI=(
    ["COI_VF2_FishR1"]="650,670"
    ["COI_coi.175f_coi.345r"]="120,140"
    ["COI_L2513_H2714"]="175,230"
)

# Base directories for each gene
declare -A BASE_DIRS=(
    ["12S"]="/nas/collabs04/refseq/12s_db"
    ["16S"]="/nas/collabs04/refseq/16s_db"
    ["COI"]="/nas/collabs04/refseq/coi_db"
)

# Available classes
CLASSES=(
#    "Actinopterygii"
#    "Amphibia"
#    "Arachnida"
#    "Aves"
#    "Bivalvia"
#    "Clitellata"
#    "Crocodylia"
#    "Gastropoda"
#    "Insecta"
#    "Lepidosauria"
#    "Malacostraca"
#    "Mammalia"
#    "Testudines"
    "Reptile"
)

# Function to show usage
show_usage() {
    echo "Usage: $0 <gene> [primer_name] [class_name]"
    echo ""
    echo "Parameters:"
    echo "  gene         : Gene type (12S, 16S, or COI)"
    echo "  primer_name  : Optional. Specific primer pair to evaluate"
    echo "  class_name   : Optional. Specific class to evaluate"
    echo ""
    echo "Examples:"
    echo "  $0 12S                                    # Evaluate all 12S primers for all classes"
    echo "  $0 16S 16S_MarVer3F_MarVer3R             # Evaluate specific 16S primer for all classes"
    echo "  $0 COI COI_VF2_FishR1 Reptile            # Evaluate specific COI primer for Reptile class"
    echo ""
    echo "Available genes: 12S, 16S, COI"
    echo "Available classes: ${CLASSES[*]}"
    echo ""
    echo "Available primers:"
    echo "  12S primers: ${!PRIMERS_12S[@]}"
    echo "  16S primers: ${!PRIMERS_16S[@]}"
    echo "  COI primers: ${!PRIMERS_COI[@]}"
}

# Function to get primers for a specific gene
get_primers_for_gene() {
    local gene=$1
    case $gene in
        12S)
            echo "${!PRIMERS_12S[@]}"
            ;;
        16S)
            echo "${!PRIMERS_16S[@]}"
            ;;
        COI)
            echo "${!PRIMERS_COI[@]}"
            ;;
        *)
            echo ""
            ;;
    esac
}

# Function to get primer sequences for a specific gene and primer name
get_primer_sequences() {
    local gene=$1
    local primer_name=$2
    case $gene in
        12S)
            echo "${PRIMERS_12S[$primer_name]}"
            ;;
        16S)
            echo "${PRIMERS_16S[$primer_name]}"
            ;;
        COI)
            echo "${PRIMERS_COI[$primer_name]}"
            ;;
        *)
            echo ""
            ;;
    esac
}

# Function to get amplicon sizes for a specific gene and primer name
get_amplicon_sizes() {
    local gene=$1
    local primer_name=$2
    case $gene in
        12S)
            echo "${AMPLICON_SIZES_12S[$primer_name]}"
            ;;
        16S)
            echo "${AMPLICON_SIZES_16S[$primer_name]}"
            ;;
        COI)
            echo "${AMPLICON_SIZES_COI[$primer_name]}"
            ;;
        *)
            echo ""
            ;;
    esac
}

# Function to log messages
log_message() {
    local gene=$1
    local primername=$2
    local class=$3
    local message=$4
    local base_dir="${BASE_DIRS[$gene]}"
    local logfile="$base_dir/primer_evaluation/${primername}_${class}/${primername}_${class}_evaluation.log"
    mkdir -p "$(dirname "$logfile")"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" | tee -a "$logfile"
}

# Function to clean up directory on failure
cleanup_on_failure() {
    local gene=$1
    local primername=$2
    local class=$3
    local base_dir="${BASE_DIRS[$gene]}"
    local output_dir="$base_dir/primer_evaluation/${primername}_${class}"

    echo "Cleaning up failed run for $primername - $class"
    if [ -d "$output_dir" ]; then
        rm -rf "$output_dir"
        echo "Removed incomplete directory: $output_dir"
    fi
    log_message "$gene" "$primername" "$class" "Run failed and directory cleaned up"
}

# Main evaluation function
evaluate_primer_for_class() {
    local gene=$1
    local primername=$2
    local class=$3
    local success=true

    # Get base directory for this gene
    local base_dir="${BASE_DIRS[$gene]}"

    # Get forward and reverse primers
    local primer_sequences=$(get_primer_sequences "$gene" "$primername")
    if [ -z "$primer_sequences" ]; then
        echo "Error: No primer sequences found for $gene - $primername"
        return 1
    fi
    IFS=',' read -r f_primer r_primer <<< "$primer_sequences"

    # Set output directory for this primer and class combination
    OUTPUT_DIR="$base_dir/primer_evaluation/${primername}_${class}"

    # Create trap for any error or exit
    trap 'cleanup_on_failure "$gene" "$primername" "$class"; exit 1' ERR

    echo "Processing primer pair: $primername for gene: $gene, class: $class"
    log_message "$gene" "$primername" "$class" "Starting processing of primer pair $primername for gene $gene, class $class"

    # Set trap for failures
    trap 'success=false' ERR

    # Create output directory after trap is set
    mkdir -p "$OUTPUT_DIR"

    # Input files for this class and gene
    SEQUENCE="$base_dir/${gene}_${class}_sequence.qza"
    TAXONOMY="$base_dir/${gene}_${class}_taxonomy.qza"

    # Check if input files exist
    if [[ ! -f "$SEQUENCE" || ! -f "$TAXONOMY" ]]; then
        log_message "$gene" "$primername" "$class" "Error: Input files not found for gene $gene, class $class"
        return 1
    fi

    # Step 1: Preprocess Dereplicate (1st)
    log_message "$gene" "$primername" "$class" "Preprocessing dereplicate step..."
    qiime rescript dereplicate \
        --i-sequences "$SEQUENCE" \
        --i-taxa "$TAXONOMY" \
        --p-mode 'uniq' \
        --p-threads 4 \
        --o-dereplicated-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep1.qza" \
        --o-dereplicated-taxa "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep1.qza" \
        --verbose

    # Step 2: Extract reads with primers
    log_message "$gene" "$primername" "$class" "Extracting reads with primers..."
    qiime feature-classifier extract-reads \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep1.qza" \
        --p-f-primer "$f_primer" \
        --p-r-primer "$r_primer" \
        --p-identity 0.8 \
        --p-n-jobs 4 \
        --o-reads "$OUTPUT_DIR/${primername}_${class}_extracted_reads.qza" \
        --verbose

    # Step 3: Dereplicate (2nd)
    log_message "$gene" "$primername" "$class" "Dereplicating sequences..."
    qiime rescript dereplicate \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_extracted_reads.qza" \
        --i-taxa "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep1.qza" \
        --p-mode 'uniq' \
        --p-threads 4 \
        --o-dereplicated-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep2.qza" \
        --o-dereplicated-taxa "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep2.qza" \
        --verbose

    # Step 4: Cull-seq
    log_message "$gene" "$primername" "$class" "Culling sequences..."
    qiime rescript cull-seqs \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep2.qza" \
        --p-n-jobs 4 \
        --p-num-degenerates 5 \
        --p-homopolymer-length 8 \
        --o-clean-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep2_cull-D6H9.qza"

    # Step 5: Filter by amplicon size
    local amplicon_sizes=$(get_amplicon_sizes "$gene" "$primername")
    if [ -z "$amplicon_sizes" ]; then
        echo "Error: No amplicon sizes found for $gene - $primername"
        return 1
    fi
    IFS=',' read -r min_size max_size <<< "$amplicon_sizes"
    log_message "$gene" "$primername" "$class" "Filtering sequences with length range $min_size-$max_size bp"

    qiime rescript filter-seqs-length \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep2_cull-D6H9.qza" \
        --p-global-min "$min_size" \
        --p-global-max "$max_size" \
        --p-threads 4 \
        --o-filtered-seqs "$OUTPUT_DIR/${primername}_${class}_sequence_derep2_filtered.qza" \
        --o-discarded-seqs "$OUTPUT_DIR/${primername}_${class}_sequence_derep2_filtered_discarded.qza"

    # Step 6: Final dereplication
    log_message "$gene" "$primername" "$class" "Performing final dereplication..."
    qiime rescript dereplicate \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep2_filtered.qza" \
        --i-taxa "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep2.qza" \
        --p-mode 'uniq' \
        --p-threads 4 \
        --o-dereplicated-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qza" \
        --o-dereplicated-taxa "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3.qza" \
        --verbose

    # Step 7: Export sequence and taxonomy files
    log_message "$gene" "$primername" "$class" "Exporting sequence and taxonomy files..."
    
    # Export sequence data
    qiime tools export \
        --input-path "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qza" \
        --output-path "$OUTPUT_DIR/${primername}_${class}_derep3"
    
    # Export taxonomy data
    qiime tools export \
        --input-path "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3.qza" \
        --output-path "$OUTPUT_DIR/${primername}_${class}_taxonomy"
    
    # Move files to better named locations
    mv "$OUTPUT_DIR/${primername}_${class}_derep3/dna-sequences.fasta" "$OUTPUT_DIR/${primername}_${class}_derep3.fasta"
    mv "$OUTPUT_DIR/${primername}_${class}_taxonomy/taxonomy.tsv" "$OUTPUT_DIR/${primername}_${class}_taxonomy.tsv"
    
    # Clean up temporary directories
    rm -r "$OUTPUT_DIR/${primername}_${class}_derep3"
    rm -r "$OUTPUT_DIR/${primername}_${class}_taxonomy"

    # Create FASTA with taxonomy information
    log_message "$gene" "$primername" "$class" "Creating FASTA file with taxonomy information..."
    
    # Path variables
    sequence_file="$OUTPUT_DIR/${primername}_${class}_derep3.fasta"
    taxonomy_file="$OUTPUT_DIR/${primername}_${class}_taxonomy.tsv"
    output_fasta="$OUTPUT_DIR/${primername}_${class}_with_taxonomy.fasta"
    
    # Create script to merge taxonomy with sequences
    cat > "$OUTPUT_DIR/merge_taxonomy.sh" << EOF
#!/bin/bash
# Combine sequence and taxonomy information
> "$output_fasta"

while IFS=\$'\t' read -r seq_id taxonomy confidence; do
    # Skip header if present
    if [[ \$seq_id == "Feature ID" ]]; then
        continue
    fi
    
    # Find the sequence with the matching ID in the fasta file
    seq=\$(grep -A1 ">\$seq_id" "$sequence_file" | tail -n1)
    
    if [[ -n "\$seq" ]]; then
        # Create the new header with taxonomy information
        header=">\$seq_id taxonomy=\$taxonomy;confidence=\$confidence"
        
        # Write the header and sequence to the output file
        echo "\$header" >> "$output_fasta"
        echo "\$seq" >> "$output_fasta"
    fi
done < "$taxonomy_file"

# Count sequences in the output file
seq_count=\$(grep -c "^>" "$output_fasta")
echo "Created FASTA with taxonomy: $output_fasta (\$seq_count sequences)"
EOF
    
    # Make the script executable and run it
    chmod +x "$OUTPUT_DIR/merge_taxonomy.sh"
    "$OUTPUT_DIR/merge_taxonomy.sh" >> "$OUTPUT_DIR/${primername}_${class}_evaluation.log"
    
    # Step 8: Analyze the combined FASTA
    log_message "$gene" "$primername" "$class" "Analyzing combined FASTA file..."
    awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$output_fasta" | \
    awk '{sum+=$1; if(min==""){min=max=$1}; if($1>max){max=$1}; if($1<min){min=$1}} END {printf "Number of Sequences with Taxonomy: %'\''d records // Min length: %'\''d bp // Max length: %'\''d bp // Avg length: %'\''0.2f bp\n", NR, min, max, sum/NR}' | \
    tee -a "$OUTPUT_DIR/${primername}_${class}_sequence_stats.txt"
    
    # Clean up the script
    rm "$OUTPUT_DIR/merge_taxonomy.sh"

    # Step 9: QIIME evaluate sequences
    log_message "$gene" "$primername" "$class" "Evaluating sequences..."
    qiime rescript evaluate-seqs \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qza" \
        --p-labels "${primername}_${class}" \
        --o-visualization "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qzv"

    # Step 10: QIIME evaluate taxonomy
    log_message "$gene" "$primername" "$class" "Evaluating taxonomy..."
    qiime rescript evaluate-taxonomy \
        --i-taxonomies "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3.qza" \
        --p-labels "${primername}_${class}" \
        --o-taxonomy-stats "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3.qzv"

    # Step 11: Train classifier and evaluate
    log_message "$gene" "$primername" "$class" "Training classifier and evaluating cross-validation..."
    qiime rescript evaluate-cross-validate \
        --i-sequences "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qza" \
        --i-taxonomy "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3.qza" \
        --o-expected-taxonomy "$OUTPUT_DIR/${primername}_${class}_sequence_derep3_eval-EXP-TAX.qza" \
        --o-observed-taxonomy "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3_eval-OBS-TAX.qza" \
        --o-evaluation "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3_eval-visual.qzv"

    # Step 12: Export evaluation visualization data
    log_message "$gene" "$primername" "$class" "Exporting evaluation visualization data..."
    TEMP_DIR="$OUTPUT_DIR/temp_export"
    mkdir -p "$TEMP_DIR"

    if [ -f "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3_eval-visual.qzv" ]; then
        qiime tools export \
            --input-path "$OUTPUT_DIR/${primername}_${class}_taxonomy_derep3_eval-visual.qzv" \
            --output-path "$TEMP_DIR"

        if [ -f "$TEMP_DIR/data.tsv" ]; then
            mv "$TEMP_DIR/data.tsv" "$OUTPUT_DIR/primer_evaluation_${primername}_${class}.tsv"
            rm -rf "$TEMP_DIR"
            log_message "$gene" "$primername" "$class" "Successfully exported evaluation data"
        else
            log_message "$gene" "$primername" "$class" "Error: data.tsv not found in export directory"
            rm -rf "$TEMP_DIR"
            return 1
        fi
    else
        log_message "$gene" "$primername" "$class" "Error: Evaluation QZV file not found"
        rm -rf "$TEMP_DIR"
        return 1
    fi

    # Reset trap and check final status
    trap - ERR

    # Check if all steps completed successfully and final TSV exists
    if [ ! -f "$OUTPUT_DIR/primer_evaluation_${primername}_${class}.tsv" ]; then
        echo "Failed to complete all steps for $primername - $class"
        cleanup_on_failure "$gene" "$primername" "$class"
        return 1
    fi

    # Step 13: Generate evaluation plots
    log_message "$gene" "$primername" "$class" "Generating evaluation plots..."

    # Create Python script with proper variable substitution
    cat > "$OUTPUT_DIR/generate_plots.py" << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import sys

def generate_plots(input_file, primername, classname):
    # Create output filenames
    output_file = os.path.splitext(input_file)[0] + '_plot.png'
    output_file_stderr = os.path.splitext(input_file)[0] + '_plot_stderrorbar.png'

    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t', comment='#')

    # Create taxonomic level mapping
    level_names = {
        1: 'kingdom',
        2: 'phylum',
        3: 'class',
        4: 'order',
        5: 'family',
        6: 'genus',
        7: 'species'
    }

    # Calculate standard error (example calculation - modify as needed)
    df['stderr'] = df['F-Measure'] * 0.05  # 5% of F-Measure as example

    # First plot (without error bars)
    plt.figure(figsize=(12, 8))
    plt.plot(df['Level'], df['F-Measure'],
            marker='^', color='green', label='F-Measure', linewidth=2)

    # Customize the plot
    plt.title(f'Primer evaluation: {primername} - {classname}')
    plt.xlabel('Taxonomic Level')
    plt.ylabel('Score')

    plt.xticks(df['Level'],
              [f"{num}\n{level_names[num]}" for num in df['Level']],
              rotation=0)

    plt.ylim(0, 1.2)
    plt.xlim(0.5, 7.5)
    plt.yticks(np.arange(0, 1.3, 0.1))

    plt.grid(True, alpha=0.3)
    plt.legend(loc='lower left')

    for x, y in zip(df['Level'], df['F-Measure']):
        plt.text(x, y, f'{y:.3f}',
                ha='center', va='bottom',
                color='green')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    # Second plot (with error bars)
    plt.figure(figsize=(12, 8))

    # Plot points and error bars
    plt.errorbar(df['Level'], df['F-Measure'],
                yerr=df['stderr'],
                fmt='^-', color='green', label='F-Measure',
                capsize=5, capthick=1, elinewidth=1, linewidth=2,
                markersize=8)

    # Customize the plot
    plt.title(f'Primer evaluation: {primername} - {classname} (with Standard Error)')
    plt.xlabel('Taxonomic Level')
    plt.ylabel('Score')

    plt.xticks(df['Level'],
              [f"{num}\n{level_names[num]}" for num in df['Level']],
              rotation=0)

    plt.ylim(0, 1.2)
    plt.xlim(0.5, 7.5)
    plt.yticks(np.arange(0, 1.3, 0.1))

    plt.grid(True, alpha=0.3)
    plt.legend(loc='lower left')

    for x, y in zip(df['Level'], df['F-Measure']):
        plt.text(x, y + df.loc[df['Level']==x, 'stderr'].iloc[0],
                f'{y:.3f}',
                ha='center', va='bottom',
                color='green')

    plt.tight_layout()
    plt.savefig(output_file_stderr, dpi=300, bbox_inches='tight')
    plt.close()

    return output_file, output_file_stderr

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_plots.py <input_file> <primer_name> <class_name>")
        sys.exit(1)

    input_file = sys.argv[1]
    primername = sys.argv[2]
    classname = sys.argv[3]

    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)

    try:
        output_file, output_file_stderr = generate_plots(input_file, primername, classname)
        print(f"Plots generated successfully:")
        print(f"Regular plot: {output_file}")
        print(f"Error bar plot: {output_file_stderr}")
    except Exception as e:
        print(f"Error generating plots: {str(e)}")
        sys.exit(1)
EOF

    # Make the Python script executable
    chmod +x "$OUTPUT_DIR/generate_plots.py"

    # Run the Python script
    if [ -f "$OUTPUT_DIR/primer_evaluation_${primername}_${class}.tsv" ]; then
        python "$OUTPUT_DIR/generate_plots.py" \
            "$OUTPUT_DIR/primer_evaluation_${primername}_${class}.tsv" \
            "$primername" \
            "$class"

        if [ $? -eq 0 ]; then
            log_message "$gene" "$primername" "$class" "Plots generated successfully"
            # Clean up the Python script
            rm "$OUTPUT_DIR/generate_plots.py"
        else
            log_message "$gene" "$primername" "$class" "Error: Failed to generate plots"
            rm "$OUTPUT_DIR/generate_plots.py"
            return 1
        fi
    else
        log_message "$gene" "$primername" "$class" "Error: TSV file not found for plotting"
        rm -f "$OUTPUT_DIR/generate_plots.py"
        return 1
    fi

    # Step 14: Export QZV visualization to PNG
    log_message "$gene" "$primername" "$class" "Exporting QZV visualization to PNG..."
    TEMP_EXPORT_DIR="$OUTPUT_DIR/temp_qzv_export"
    mkdir -p "$TEMP_EXPORT_DIR"

    qiime tools export \
        --input-path "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.qzv" \
        --output-path "$TEMP_EXPORT_DIR"

    # Check if visualization.png exists and move it
    if [ -f "$TEMP_EXPORT_DIR/evaluate_seqs.png" ]; then
        mv "$TEMP_EXPORT_DIR/evaluate_seqs.png" "$OUTPUT_DIR/${primername}_${class}_sequence_derep3.png"
        log_message "$gene" "$primername" "$class" "Successfully exported evaluate_seqs.png to ${primername}_${class}_sequence_derep3.png"
    else
        log_message "$gene" "$primername" "$class" "Warning: evaluate_seqs.png not found in export"
    fi

    # Clean up temporary directory
    rm -rf "$TEMP_EXPORT_DIR"
}

# Main script execution
main() {
    # Check arguments
    if [ $# -eq 0 ]; then
        echo "Error: No arguments provided"
        show_usage
        exit 1
    fi

    # Check for help flags first
    if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        show_usage
        exit 0
    fi

    # First argument must be gene type
    GENE="$1"
    if [[ ! -n "${BASE_DIRS[$GENE]}" ]]; then
        echo "Error: Invalid gene type. Must be one of: ${!BASE_DIRS[@]}"
        show_usage
        exit 1
    fi

    # Get available primers for this gene
    GENE_PRIMERS=($(get_primers_for_gene "$GENE"))
    if [ ${#GENE_PRIMERS[@]} -eq 0 ]; then
        echo "Error: No primers found for gene $GENE"
        exit 1
    fi

    # If no additional arguments, process all primers for all classes
    if [ $# -eq 1 ]; then
        echo "Processing all $GENE primers for all classes"
        for primer in "${GENE_PRIMERS[@]}"; do
            echo "Starting primer: $primer"
            for class in "${CLASSES[@]}"; do
                echo "Starting evaluation for $primer - $class"
                if ! evaluate_primer_for_class "$GENE" "$primer" "$class"; then
                    echo "Skipping to next combination due to failure"
                    continue
                fi
            done
        done
        exit 0
    fi

    # If primer name is provided
    if [ $# -ge 2 ]; then
        PRIMER="$2"
        # Check if primer exists for this gene
        if [[ ! " ${GENE_PRIMERS[@]} " =~ " ${PRIMER} " ]]; then
            echo "Error: Invalid primer name for gene $GENE. Available primers are:"
            printf '%s\n' "${GENE_PRIMERS[@]}"
            exit 1
        fi
    fi

    # If class name is provided
    if [ $# -ge 3 ]; then
        CLASS="$3"
        if [[ ! " ${CLASSES[@]} " =~ " ${CLASS} " ]]; then
            echo "Error: Invalid class name. Available classes are:"
            printf '%s\n' "${CLASSES[@]}"
            exit 1
        fi
        # Process specific primer and class
        evaluate_primer_for_class "$GENE" "$PRIMER" "$CLASS"
    else
        # Process specific primer for all classes
        for class in "${CLASSES[@]}"; do
            echo "Starting evaluation for $PRIMER - $class"
            if ! evaluate_primer_for_class "$GENE" "$PRIMER" "$class"; then
                echo "Skipping to next combination due to failure"
                continue
            fi
        done
    fi
}

# Call the main function with all arguments
main "$@"