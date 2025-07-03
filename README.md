# CEG - Otter Diet Metabarcoding Project
DNA Metabarcoding-based diet diversity and niche differences of otters in Thailand

## 1. Retrieve Reference Sequences from NCBI Using QIIME 2
This step retrieves curated reference sequences for key mitochondrial genes — 12S, 16S, and COI — from NCBI, using flexible query strings that include multiple gene name synonyms. This ensures broad and accurate capture of target sequences across diverse taxonomic groups for metabarcoding studies.

<pre><code>gene_queries["12S"] = '(12S[Title] OR 12s[Title] OR "s-rRNA"[Title] OR rrn12[Title] OR rrnS[Title] OR rrns[Title] OR srRNA[Title] OR "small subunit ribosomal RNA"[Title] OR "small ribosomal RNA subunit"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnS"[Title])'
gene_queries["16S"] = '(16S[Title] OR 16s[Title] OR rrn16[Title] OR rrnA[Title] OR rrnB[Title] OR lrRNA[Title] OR "l-rRNA"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-rRNA"[Title] OR "mt-rrnL"[Title])'
gene_queries["COI"] = '(COI[Title] OR CO1[Title] OR "cytochrome oxidase I"[Title] OR "cytochrome c oxidase subunit I"[Title] OR "cytochrome c oxidase I"[Title] OR cox1[Title] OR coxi[Title] OR "COX subunit 1"[Title] OR mtDNA[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR "mt-co1"[Title] OR "mt-coI"[Title] OR "mt-rRNA"[Title])'</code></pre>

**Target Taxonomic Groups**<br>
We focus on 11 eukaryotic classes relevant to otter and carnivore diet analysis.

**Vertebrate**
<pre><code>"Actinopterygii"  # Ray-finned fishes
"Amphibia"        # Amphibians
"Mammalia"        # Mammals
"Aves"            # Birds
"Crocodylia"      # Crocodilians (Reptilia)
"Lepidosauria"    # Lizards and snakes (Reptilia)
"Testudines"      # Turtles (Reptilia)</code></pre> 
**Invertebrate**
<pre><code>"Arachnida"       # Scorpions and spiders
"Insecta"         # Insects
"Bivalvia"        # Clams, oysters
"Gastropoda"      # Snails
"Clitellata"      # Earthworms and leeches
"Malacostraca"    # Crabs, shrimps</code></pre>

⚠️ Note: In NCBI, reptiles are often split into subgroups like Crocodylia, Lepidosauria, and Testudines. After retrieving sequences, we merge these into Reptilia for downstream analysis.

**QIIME 2 Command for Retrieval**<br>
For each group–gene combination, sequences are downloaded using qiime rescript get-ncbi-data. Here's the core command pattern:<br>
<pre><code>qiime rescript get-ncbi-data \
  --p-query "$group[Organism] AND 80[SLEN]:30000[SLEN] AND ($gene_query) NOT \"environmental sample\"[Title] NOT \"environmental samples\"[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] NOT PREDICTED[Title]" \
  --p-ranks kingdom phylum class order family genus species \
  --p-rank-propagation \
  --p-n-jobs 5 \
  --o-sequences "${gene}_${group}_sequences.qza" \
  --o-taxonomy "${gene}_${group}_taxonomy.qza" \
  --verbose > "${gene}_${group}_qiime_download.log" 2>&1</code></pre>

This query searches NCBI for sequences from the specified organism group ($group), with sequence lengths between 80 and 30,000 bases, matching the gene-specific keywords ($gene_query). It excludes records labeled as environmental samples, uncultured, unclassified, unidentified, unverified, or predicted.

**Example Output Structure**
<pre><code>├── COI_Actinopterygii_sequences.qza
├── COI_Actinopterygii_taxonomy.qza
├── 12S_Mammalia_sequences.qza
├── 12S_Mammalia_taxonomy.qza
├── ...
├── 16S_Aves_qiime_download.log</code></pre>

## 2. Primer Evaluation

**Overview**<br>
This primer evaluation pipeline assesses the performance of DNA metabarcoding primers for taxonomic classification using QIIME 2 and RESCRIPt.
It estimates:
- **Binding capacity**: how many reference sequences are amplified by each primer set.
- **Taxonomic resolution**: how effectively each primer recovers taxonomic labels (e.g., kingdom to species levels).

The goal is to identify effective DNA metabarcoding primers by evaluating their binding coverage and taxonomic resolution, supporting accurate and efficient detection of target taxa in ecological studies such as diet analysis and biodiversity monitoring.<br>

**Example Primer Evaluation Output**<br>
![image](https://github.com/user-attachments/assets/a2b9f102-2ad3-4f2d-9dcf-6b27d7447665)


**Primer candidate**
|Gene| Target | Primer F | Fw Sequence | Primer R | Re Sequence | Amplicon Size | Reference |
|-----|--------|----------|-------------|----------|-------------|---------------|-----------|
|12s| Vertebrate + Invertebrate | 12Sv5Fw | TTAGATACCCCACTATGC | 12Sv5Re | TAGAACAGGCTCCTCTAG | 100 - 137 (inc. pm) | Riaz et al. 2011 (Ta=52), adopted by Kumari et al. 2019 (Ta=50) |
|12s| Fish - Actinopterygian | MiFish-U-F | GTCGGTAAAACTCGTGCCAGC | MiFish-U-R | CATAGTGGGGTATCTAATCCCAGTTTG | 170bp (exc. pm) 209 - 232 (inc. pm) | Miya et al. 2015 (Ta=59), 2020 |
|12s| Fish   | tele02F  | AAACTCGTGCCAGCCACC | tele02R  | GGGTATCTAATCCCAGTTTG | 219 (inc. pm) | Taberlet et al. 2018 (Ta=52) |
|12s| Fish   | teleoF   | ACACCGCCCGTCACTCT | teleoRdeg | CTTCCGGTACACTTACCRTG | 60 (exc. pm) 100 - 120 (inc.pm) | Valentini et al. 2016 |
|16s| Crustacean + Invertebrate | 16SMAVF | CCAACATCGAGGTCRYAA | 16SMAVR | ARTTACYNTAGGGATAACAG | 74 (inc.pm) | De Barba et al. 2014, adopted by Kumari et al. 2019 (Ta=55) |
|16s| Freshwater vertebrate (focusing on fish + amphibian) | Vert-16S-eDNA-F1 | AGACGAGAAGACCCYDTGGAGCTT | Vert-16S-eDNA-R1 | GATCCAACATCGAGGTCGTAA | 250 - 335 (inc. pm) | Vences et al. 2016 |
|16s| Fish + Marine Mammal | MarVer3F | AGACGAGAAGACCCTRTG | MarVer3R | GGATTGCGCTGTTATCCC | 245 - 285 (inc.pm) | Valsecchi et al. 2020 |
|COI| Fish   | coi.175f | GGAGGCTTTGGMAAYTGRYT | coi.345r | TAGAGGRGGGTARACWGTYCA | 130 - 170 (inc.pm) | Collins et al. 2019 (Ta=53) |
|COI| Fish   | L2513    | GCCTGTTTACCAAAAACATCA | H2714 | CTCCATAGGGTCTTCTCGTCTT | 250 (inc.pm) | Kitano et al. 2007 (Ta=55) |

## Usage

<pre><code>
# Basic syntax
./primer_evaluation.sh <gene> [primer_name] [class_name]

# Examples
./primer_evaluation.sh 12S                           # All 12S primers, all classes
./primer_evaluation.sh 16S 16S_MarVer3F_MarVer3R     # Specific 16S primer, all classes
./primer_evaluation.sh COI COI_VF2_FishR1 Reptile    # Specific primer and class</code></pre>

## Evaluation Pipeline Steps

### 1. **Initial Dereplication**
- Removes duplicate sequences from input database
- Maintains taxonomic assignments for unique sequences
- Uses `qiime rescript dereplicate` with 'uniq' mode

### 2. **Primer-Based Read Extraction**
- Extracts sequences that match forward and reverse primer sequences
- Uses 80% identity threshold for primer matching
- Simulates PCR amplification *in silico*

### 3. **Post-Extraction Dereplication**
- Second round of deduplication after primer extraction
- Ensures only unique amplicon sequences remain

### 4. **Sequence Quality Control**
- **Cull sequences**: Removes sequences with excessive degenerates (>5) and long homopolymers (>8bp)
- **Length filtering**: Retains only sequences within expected amplicon size ranges
- Filters out sequences likely to cause classification errors

### 5. **Final Dereplication**
- Third and final deduplication step
- Produces clean, non-redundant reference database

### 6. **Data Export and Processing**
- Exports sequences and taxonomy from QIIME 2 artifacts (.qza) to standard formats
- Creates FASTA files with integrated taxonomic information
- Generates sequence statistics (count, length distribution)

### 7. **Performance Evaluation**

#### Sequence Assessment
- Evaluates sequence quality and composition
- Generates summary statistics and visualizations

#### Taxonomic Coverage Analysis  
- Assesses taxonomic representation across different levels
- Identifies gaps in taxonomic coverage

#### Cross-Validation Testing
- Trains naive Bayes classifier using processed sequences
- Performs k-fold cross-validation to test classification accuracy
- Generates F-measure scores across taxonomic levels (kingdom → species)

### 8. **Results Visualization**
- **Performance plots**: F-measure scores across taxonomic levels
- **Error bar plots**: Include confidence intervals/standard errors
- **PNG exports**: Publication-ready visualizations of sequence quality metrics

## Output Files

For each primer-class combination, the pipeline generates:

| File Type | Description |
|-----------|-------------|
| `*_derep3.fasta` | Final processed sequences |
| `*_taxonomy.tsv` | Taxonomic assignments |
| `*_with_taxonomy.fasta` | FASTA with taxonomy in headers |
| `*_sequence_stats.txt` | Sequence count and length statistics |
| `primer_evaluation_*.tsv` | Cross-validation results |
| `*_plot.png` | Performance visualization |
| `*_plot_stderrorbar.png` | Performance plot with error bars |
| `*_evaluation.log` | Detailed processing log |

## Supported Primers

### 12S rRNA
- `12S_12SV5F_12SV5R` (80-105 bp)
- `12S_tele02F_tele02R` (140-200 bp)  
- `12S_teleoF_teleoRdeg` (50-100 bp)
- `12S_MiFish-U-F_MiFish-U-R` (150-200 bp)

### 16S rRNA
- `16S_Vert-16S-eDNA-F1_Vert-16S-eDNA-R1` (150-260 bp)
- `16S_MarVer3F_MarVer3R` (100-250 bp)

### COI
- `COI_VF2_FishR1` (650-670 bp)
- `COI_coi.175f_coi.345r` (120-140 bp)
- `COI_L2513_H2714` (175-230 bp)

## Requirements

- QIIME 2 (with RESCRIPt plugin)
- Python 3 (with pandas, matplotlib, numpy)
- Input databases in QIIME 2 format (.qza files)

## Key Features

- **Automated processing**: Complete pipeline from raw sequences to evaluation metrics
- **Multi-gene support**: Handles different gene regions with appropriate parameters
- **Robust error handling**: Cleanup procedures for failed runs
- **Comprehensive logging**: Detailed logs for troubleshooting
- **Publication-ready outputs**: High-quality plots and summary statistics

## Interpretation

**F-measure scores** range from 0-1, where:
- **> 0.9**: Excellent classification performance
- **0.7-0.9**: Good performance  
- **0.5-0.7**: Moderate performance
- **< 0.5**: Poor performance, primer may not be suitable

Performance typically decreases from kingdom level (highest) to species level (lowest), which is expected due to increased taxonomic resolution requirements.
