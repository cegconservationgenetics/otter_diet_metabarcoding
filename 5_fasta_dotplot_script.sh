#!/bin/bash

# FASTA to Excel Dot Plot Converter
# Usage: ./fasta_to_excel_dotplot.sh input.fasta START_POSITION END_POSITION output_name.xlsx

# Function to display usage
usage() {
    echo "Usage: $0 <input.fasta> <start_position> <end_position> <output_name.xlsx>"
    echo ""
    echo "Parameters:"
    echo "  input.fasta     - Aligned FASTA file (first sequence should be blocking primer)"
    echo "  start_position  - Starting position for alignment window (1-based)"
    echo "  end_position    - Ending position for alignment window (1-based)"
    echo "  output_name.xlsx- Output Excel file name"
    echo ""
    echo "Example:"
    echo "  $0 aligned_sequences.fas 1998 2027 primer_analysis.xlsx"
    echo ""
    exit 1
}

# Check if correct number of arguments provided
if [ $# -ne 4 ]; then
    echo "Error: Incorrect number of arguments!"
    echo ""
    usage
fi

# Assign command line arguments to variables
INPUT_FASTA="$1"
START_POSITION="$2"
END_POSITION="$3"
OUTPUT_EXCEL="$4"

# Validate input parameters
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input file '$INPUT_FASTA' not found!"
    exit 1
fi

# Check if positions are valid numbers
if ! [[ "$START_POSITION" =~ ^[0-9]+$ ]] || ! [[ "$END_POSITION" =~ ^[0-9]+$ ]]; then
    echo "Error: Start and end positions must be positive integers!"
    exit 1
fi

if [ "$START_POSITION" -ge "$END_POSITION" ]; then
    echo "Error: Start position must be less than end position!"
    exit 1
fi

# Check if output file has .xlsx extension
if [[ ! "$OUTPUT_EXCEL" =~ \.xlsx$ ]]; then
    echo "Warning: Output file should have .xlsx extension. Adding it..."
    OUTPUT_EXCEL="${OUTPUT_EXCEL}.xlsx"
fi

# Calculate window size
WINDOW_SIZE=$((END_POSITION - START_POSITION + 1))

echo "========================================="
echo "FASTA to Excel Dot Plot Converter"
echo "========================================="
echo "Input file: $INPUT_FASTA"
echo "Output file: $OUTPUT_EXCEL"
echo "Alignment window: positions $START_POSITION-$END_POSITION ($WINDOW_SIZE bp)"
echo ""

# Check if required Python packages are available
python3 -c "import pandas, openpyxl" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: Required Python packages not found!"
    echo "Please install with: pip install pandas openpyxl"
    exit 1
fi

echo "Processing aligned FASTA file..."

# Use Python to process the aligned FASTA file
python3 - <<EOF
import pandas as pd
import re
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils.dataframe import dataframe_to_rows

def parse_fasta_simple(filename):
    """Simple FASTA parser for aligned sequences"""
    sequences = []
    current_header = ""
    current_seq = ""
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences.append({
                        'header': current_header[1:],  # Remove '>' character
                        'sequence': current_seq.upper()
                    })
                current_header = line
                current_seq = ""
            else:
                current_seq += line.replace(' ', '')  # Remove spaces but keep gaps
        
        if current_header and current_seq:
            sequences.append({
                'header': current_header[1:],  # Remove '>' character
                'sequence': current_seq.upper()
            })
    
    return sequences

def parse_taxonomy_from_header(header):
    """Extract taxonomy information from header"""
    # Initialize taxonomy
    taxonomy = {
        'Accession': '',
        'Class': '',
        'Order': '',
        'Family': '',
        'Scientific_name': ''
    }
    
    # Try to extract accession (usually first part before space or semicolon)
    accession_match = re.match(r'^([A-Z0-9_.]+)', header)
    if accession_match:
        taxonomy['Accession'] = accession_match.group(1)
    else:
        taxonomy['Accession'] = header.split()[0] if header.split() else header[:20]
    
    # Look for taxonomic information in the header using various patterns
    if ';' in header:
        parts = [part.strip() for part in header.split(';')]
        for part in parts:
            if part.startswith('c__'):
                taxonomy['Class'] = part[3:]
            elif part.startswith('o__'):
                taxonomy['Order'] = part[3:]
            elif part.startswith('f__'):
                taxonomy['Family'] = part[3:]
            elif part.startswith('g__') and part.startswith('s__'):
                # Handle genus and species
                continue
            elif 'class' in part.lower():
                taxonomy['Class'] = part.split('_')[-1] if '_' in part else part
            elif 'order' in part.lower():
                taxonomy['Order'] = part.split('_')[-1] if '_' in part else part
            elif 'family' in part.lower():
                taxonomy['Family'] = part.split('_')[-1] if '_' in part else part
    
    # Enhanced scientific name extraction
    scientific_name = ""
    
    # Method 1: Look for g__ and s__ patterns
    if ';' in header:
        parts = [part.strip() for part in header.split(';')]
        genus = ""
        species = ""
        for part in parts:
            if part.startswith('g__'):
                genus = part[3:]
            elif part.startswith('s__'):
                species = part[3:]
        
        if genus and species and species != "sp":
            scientific_name = f"{genus} {species}"
        elif genus:
            scientific_name = f"{genus} sp."
    
    # Method 2: Look for patterns in the header text
    if not scientific_name:
        # Common patterns in sequence headers
        patterns = [
            r'([A-Z][a-z]+)\s+([a-z]+)',  # Genus species
            r'([A-Z][a-z]+)_([a-z]+)',    # Genus_species
            r'>([A-Z][a-z]+)\s([a-z]+)',  # >Genus species
        ]
        
        for pattern in patterns:
            match = re.search(pattern, header)
            if match:
                genus, species = match.groups()
                if len(genus) > 2 and len(species) > 2 and species != "sp":
                    scientific_name = f"{genus} {species}"
                    break
                elif len(genus) > 2:
                    scientific_name = f"{genus} sp."
                    break
    
    # Method 3: Extract from common header formats
    if not scientific_name:
        # Look for words that could be genus/species
        words = re.findall(r'[A-Z][a-z]+', header)
        if len(words) >= 2:
            # Skip common non-taxonomic words
            skip_words = ['GenBank', 'NCBI', 'DNA', 'RNA', 'PCR', 'Clone', 'Isolate', 'Strain']
            filtered_words = [w for w in words if w not in skip_words]
            
            if len(filtered_words) >= 2:
                scientific_name = f"{filtered_words[0]} {filtered_words[1]}"
            elif len(filtered_words) == 1:
                scientific_name = f"{filtered_words[0]} sp."
    
    # Fallback: use accession if nothing found
    if not scientific_name:
        scientific_name = taxonomy['Accession']
    
    taxonomy['Scientific_name'] = scientific_name
    
    return taxonomy

def create_alignment_excel(input_file, output_file, start_pos, end_pos):
    """Create Excel file with alignment dot plot view"""
    
    print("Loading aligned sequences...")
    sequences = parse_fasta_simple(input_file)
    print(f"Loaded {len(sequences)} sequences")
    
    if not sequences:
        print("No sequences found!")
        return
    
    # Check alignment length
    alignment_length = len(sequences[0]['sequence'])
    print(f"Alignment length: {alignment_length} bp")
    
    if end_pos > alignment_length:
        print(f"Warning: End position {end_pos} exceeds alignment length {alignment_length}")
        end_pos = alignment_length
    
    # Prepare data for DataFrame
    data_rows = []
    window_sequences = []
    
    print("Processing sequences and extracting alignment window...")
    
    # Get reference sequence (first sequence - should be blocking primer)
    reference_record = sequences[0]
    reference_taxonomy = parse_taxonomy_from_header(reference_record['header'])
    reference_window = reference_record['sequence'][start_pos-1:end_pos]
    
    print(f"\nUsing reference sequence (blocking primer): {reference_taxonomy['Scientific_name']}")
    print(f"Reference window: {reference_window}")
    
    for i, record in enumerate(sequences):
        # Parse taxonomy
        taxonomy = parse_taxonomy_from_header(record['header'])
        
        # Show first few examples of taxonomy parsing
        if i < 5:
            print(f"  Example {i+1}: {record['header'][:80]}...")
            print(f"    -> Scientific name: {taxonomy['Scientific_name']}")
        
        # Extract sequence window (convert to 0-based indexing)
        seq_window = record['sequence'][start_pos-1:end_pos]
        window_sequences.append(seq_window)
        
        # Create base row with taxonomy info
        row = {
            'Accession': taxonomy['Accession'],
            'Class': taxonomy['Class'],
            'Order': taxonomy['Order'],
            'Family': taxonomy['Family'],
            'Scientific_name': taxonomy['Scientific_name']
        }
        
        # Add individual base positions with dot notation and gap handling
        mismatch_count = 0
        gap_count = 0
        for j, base in enumerate(seq_window):
            if i == 0:  # First sequence (reference/blocking primer) - show actual bases
                row[str(j + 1)] = base
            else:  # Comparison sequences - show dots for matches
                ref_base = reference_window[j]
                # Count gaps for n/a logic
                if base == '-':
                    gap_count += 1
                    row[str(j + 1)] = '-'  # Keep gaps as gaps (don't convert to dots)
                    # Don't increment mismatch_count for gaps
                elif base == ref_base:
                    row[str(j + 1)] = '.'  # Same base = match
                    # Don't increment mismatch_count
                else:
                    row[str(j + 1)] = base  # Different base = show actual base
                    mismatch_count += 1
        
        # Add mismatch count - show "n/a" if all positions are gaps
        if i == 0:
            row['Mismatches'] = 0  # Reference has 0 mismatches
        elif gap_count == len(seq_window):
            row['Mismatches'] = 'n/a'  # All gaps = n/a
        else:
            row['Mismatches'] = mismatch_count
        
        data_rows.append(row)
    
    # Create DataFrame
    print("Creating DataFrame...")
    
    # Define column order
    base_columns = [str(i) for i in range(1, len(seq_window) + 1)]
    column_order = ['Accession', 'Class', 'Order', 'Family', 'Scientific_name'] + base_columns + ['Mismatches']
    
    df = pd.DataFrame(data_rows)
    df = df.reindex(columns=column_order, fill_value='')
    
    print("Creating Excel file with formatting...")
    
    # Create workbook and worksheet
    wb = Workbook()
    ws = wb.active
    ws.title = "Alignment_DotPlot"
    
    # Add DataFrame to worksheet
    for r in dataframe_to_rows(df, index=False, header=True):
        ws.append(r)
    
    # Define colors for bases and dots
    base_colors = {
        'A': PatternFill(start_color="FF9999", end_color="FF9999", fill_type="solid"),  # Light red
        'T': PatternFill(start_color="99FF99", end_color="99FF99", fill_type="solid"),  # Light green
        'G': PatternFill(start_color="9999FF", end_color="9999FF", fill_type="solid"),  # Light blue
        'C': PatternFill(start_color="FFFF99", end_color="FFFF99", fill_type="solid"),  # Light yellow
        '-': PatternFill(start_color="CCCCCC", end_color="CCCCCC", fill_type="solid"),  # Light gray for gaps
        'N': PatternFill(start_color="FFFFFF", end_color="FFFFFF", fill_type="solid"),  # White for N
        '.': PatternFill(start_color="F0F0F0", end_color="F0F0F0", fill_type="solid")   # Very light gray for dots
    }
    
    # Format header row
    header_fill = PatternFill(start_color="E6E6FA", end_color="E6E6FA", fill_type="solid")  # Lavender
    header_font = Font(bold=True)
    
    for col in range(1, len(column_order) + 1):
        cell = ws.cell(row=1, column=col)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")
    
    # Format base columns with colors
    base_col_start = 6  # After taxonomy columns
    mismatch_col = base_col_start + len(base_columns)  # Mismatch column position
    
    for row in range(2, len(sequences) + 2):  # Skip header row
        # Color base columns
        for col in range(base_col_start, base_col_start + len(base_columns)):
            cell = ws.cell(row=row, column=col)
            base = cell.value
            
            if base in base_colors:
                cell.fill = base_colors[base]
            
            cell.alignment = Alignment(horizontal="center")
            cell.font = Font(name="Courier New", size=10)  # Monospace font
        
        # Format mismatch column
        mismatch_cell = ws.cell(row=row, column=mismatch_col)
        mismatch_value = mismatch_cell.value
        
        # Color-code mismatches: Green for 0, Yellow for 1-3, Red for >3, Gray for n/a
        if mismatch_value == 'n/a':
            mismatch_cell.fill = PatternFill(start_color="D3D3D3", end_color="D3D3D3", fill_type="solid")  # Light gray for n/a
        elif mismatch_value == 0:
            mismatch_cell.fill = PatternFill(start_color="90EE90", end_color="90EE90", fill_type="solid")  # Light green
        elif isinstance(mismatch_value, int) and mismatch_value <= 3:
            mismatch_cell.fill = PatternFill(start_color="FFFF99", end_color="FFFF99", fill_type="solid")  # Light yellow
        elif isinstance(mismatch_value, int) and mismatch_value > 3:
            mismatch_cell.fill = PatternFill(start_color="FFB6C1", end_color="FFB6C1", fill_type="solid")  # Light pink
        
        mismatch_cell.alignment = Alignment(horizontal="center")
        mismatch_cell.font = Font(bold=True)
    
    # Auto-adjust column widths
    for column in ws.columns:
        max_length = 0
        column_letter = column[0].column_letter
        
        for cell in column:
            try:
                if len(str(cell.value)) > max_length:
                    max_length = len(str(cell.value))
            except:
                pass
        
        # Set appropriate width
        if column_letter <= 'E':  # Taxonomy columns
            adjusted_width = min(max_length + 2, 25)
        elif column_letter == chr(ord('F') + len(base_columns)):  # Mismatch column
            adjusted_width = 12
        else:  # Base columns
            adjusted_width = 3  # Narrow columns for bases
        
        ws.column_dimensions[column_letter].width = adjusted_width
    
    # Add a summary sheet with legend
    summary_ws = wb.create_sheet("Legend")
    
    # Add legend information
    legend_data = [
        ["PCR Blocking Primer Analysis - BioEdit-Style Dot Plot"],
        ["Symbol", "Color", "Meaning"],
        ["A", "Light Red", "Adenine"],
        ["T", "Light Green", "Thymine"],
        ["G", "Light Blue", "Guanine"],
        ["C", "Light Yellow", "Cytosine"],
        ["-", "Light Gray", "Gap"],
        ["N", "White", "Ambiguous"],
        [".", "Very Light Gray", "Same as blocking primer"],
        [""],
        ["Gap Handling Logic"],
        ["Gaps (-) never count as mismatches"],
        ["Gaps are shown as (-) not dots"],
        ["Gap vs any base = no mismatch penalty"],
        ["Only actual base differences count as mismatches"],
        ["All gaps = n/a in mismatch column"],
        [""],
        ["Mismatch Color Coding"],
        ["0 mismatches", "Light Green", "Perfect match (will be blocked)"],
        ["1-3 mismatches", "Light Yellow", "Partial blocking possible"],
        [">3 mismatches", "Light Pink", "No blocking (will amplify)"],
        ["n/a", "Light Gray", "All gaps (no sequence data)"],
        [""],
        ["Blocking Primer Analysis"],
        ["First sequence = blocking primer (reference)"],
        ["Identical bases shown as dots (.)"],
        ["Different bases shown as actual nucleotides"],
        ["Lower mismatch count = better blocking efficiency"],
        [""],
        ["Analysis Information"],
        ["Input File", input_file],
        ["Position Range", f"{start_pos}-{end_pos}"],
        ["Total Sequences", len(sequences)],
        ["Window Size", f"{len(base_columns)} bp"],
        ["Blocking Primer", f"{data_rows[0]['Scientific_name']} ({data_rows[0]['Accession']})"]
    ]
    
    for row_data in legend_data:
        summary_ws.append(row_data)
    
    # Format legend
    for row in range(1, len(legend_data) + 1):
        for col in range(1, 4):
            cell = summary_ws.cell(row=row, column=col)
            if row == 1 or row == 11 or row == 18 or row == 24:  # Headers
                cell.font = Font(bold=True, size=12)
            elif row == 2:  # Subheader
                cell.font = Font(bold=True)
    
    # Save workbook
    wb.save(output_file)
    
    print(f"Excel file created: {output_file}")
    
    # Print summary with mismatch analysis
    print(f"\nSummary:")
    print(f"Total sequences: {len(sequences)}")
    print(f"Alignment positions: {start_pos}-{end_pos} ({len(base_columns)} bp)")
    
    # Analyze mismatches (with gap logic)
    mismatch_counts = [row['Mismatches'] for row in data_rows[1:]]  # Skip reference sequence
    numeric_mismatches = [x for x in mismatch_counts if isinstance(x, int)]  # Only numeric values
    na_count = sum(1 for x in mismatch_counts if x == 'n/a')  # Count n/a values
    
    if mismatch_counts:
        print(f"\nBlocking Efficiency Analysis:")
        print(f"  Perfect matches (0 mismatches): {mismatch_counts.count(0)} sequences (WILL BE BLOCKED)")
        print(f"  Low mismatches (1-3): {sum(1 for x in numeric_mismatches if 1 <= x <= 3)} sequences (PARTIAL BLOCKING)")
        print(f"  High mismatches (>3): {sum(1 for x in numeric_mismatches if x > 3)} sequences (WILL AMPLIFY)")
        print(f"  All gaps (n/a): {na_count} sequences (NO DATA)")
        if numeric_mismatches:
            print(f"  Average mismatches: {sum(numeric_mismatches)/len(numeric_mismatches):.2f}")
            print(f"  Max mismatches: {max(numeric_mismatches)}")
        
        # Show species with perfect matches (will be blocked)
        perfect_matches = [row for row in data_rows[1:] if row['Mismatches'] == 0]
        if perfect_matches:
            print(f"\n  Species that WILL BE BLOCKED (0 mismatches):")
            for species in perfect_matches[:10]:  # Show first 10
                print(f"    - {species['Scientific_name']}")
            if len(perfect_matches) > 10:
                print(f"    ... and {len(perfect_matches) - 10} more")
        
        # Show species that will amplify (>3 mismatches)
        high_mismatch = [row for row in data_rows[1:] if isinstance(row['Mismatches'], int) and row['Mismatches'] > 3]
        if high_mismatch:
            print(f"\n  Species that WILL AMPLIFY (>3 mismatches):")
            for species in high_mismatch[:10]:  # Show first 10
                print(f"    - {species['Scientific_name']} ({species['Mismatches']} mismatches)")
            if len(high_mismatch) > 10:
                print(f"    ... and {len(high_mismatch) - 10} more")
    
    print(f"\nBase composition in window:")
    
    # Count bases across all sequences in the window
    base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'N': 0, 'Other': 0}
    total_bases = 0
    
    for seq_window in window_sequences:
        for base in seq_window:
            if base in base_counts:
                base_counts[base] += 1
            else:
                base_counts['Other'] += 1
            total_bases += 1
    
    for base, count in base_counts.items():
        if count > 0:
            percentage = (count / total_bases) * 100
            print(f"  {base}: {count} ({percentage:.1f}%)")
    
    return df

# Run the conversion
try:
    result_df = create_alignment_excel('$INPUT_FASTA', '$OUTPUT_EXCEL', $START_POSITION, $END_POSITION)
    print("\n========================================")
    print("Conversion completed successfully!")
    print("========================================")
    print("The Excel file contains:")
    print("  - Main sheet: Alignment with color-coded bases")
    print("  - Legend sheet: Color coding and summary information")
    print("  - Blocking analysis: Shows which species will be blocked/amplified")
    print("")
    print("INTERPRETATION:")
    print("  🟢 Green (0 mismatches): Species will be BLOCKED by primer")
    print("  🟡 Yellow (1-3 mismatches): Partial blocking possible")
    print("  🔴 Pink (>3 mismatches): Species will AMPLIFY despite primer")
    print("  ⚪ Gray (n/a): No sequence data available")
except ImportError as e:
    if 'pandas' in str(e) or 'openpyxl' in str(e):
        print("Error: Required libraries missing. Install with: pip install pandas openpyxl")
    else:
        print(f"Error: Missing required library - {e}")
except FileNotFoundError:
    print("Error: Input file '$INPUT_FASTA' not found!")
except Exception as e:
    print(f"Error: {e}")
EOF

echo ""
echo "========================================="
echo "Script completed!"
echo "========================================="
echo "Output file: $OUTPUT_EXCEL"
echo ""
echo "Next steps:"
echo "1. Open the Excel file to view the blocking analysis"
echo "2. Check the 'Legend' sheet for detailed interpretation"
echo "3. Review species in green (will be blocked) vs pink (will amplify)"
