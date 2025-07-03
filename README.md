# CEG - Otter Diet Metabarcoding Project
DNA Metabarcoding-based diet diversity and niche differences of otters in Thailand

## 1. Retrieve Reference Sequences from NCBI Using QIIME 2 (1_download_reference_sequence.sh)
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

## 2. CRABS Insilico PCR (2_crabs_insilico_pcr.py)
**Overview**
We use CRABS (Computational Rapid Amplicon Binding Simulator) to perform in silico PCR and evaluate the expected amplicon size distribution for each candidate primer set across different taxonomic groups. This step helps confirm that primers amplify within the expected size range, identify size variability across taxa, and screen for non-target amplification or unintended products. Results are summarized using dot plots that show amplicon length distributions for each primer.

**Usage**
<pre><code># Basic syntax
python3 crabs_insilico_pcr.py <gene> [options]

# Examples
python3 crabs_insilico_pcr.py 12s                           # All 12S primers
python3 crabs_insilico_pcr.py 16s --primer 16S_MarVer3      # Specific 16S primer  
python3 crabs_insilico_pcr.py coi --mismatch 15.0           # COI with custom mismatch tolerance
python3 crabs_insilico_pcr.py 12s --plot-only               # Generate plots only</code></pre>

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

**Amplicon Size Distribution Workflow Steps**

**1. Input Preparation**
- Reads taxonomic class-specific FASTA files (*_sequences.fasta)
- Processes each class separately for comprehensive coverage
- Supports multiple gene regions with appropriate primer sets

**2. Mismatch Calculation**
- Calculates allowed mismatches based on primer length and tolerance percentage
- Default: 10% of average primer pair length
- Accounts for primer degeneracies and biological sequence variation

**3. In Silico PCR Amplification**
- Uses CRABS insilico_pcr command to simulate PCR amplification
- Searches for forward and reverse primer binding sites
- Extracts sequences between primer pairs (amplicons)
- Applied per primer-class combination

**4. Quality Control**
- Identifies and removes empty output files (no successful amplifications)
- Counts sequences successfully amplified per primer-class combination
- Logs detailed processing information for troubleshooting

**5. Amplicon Size Analysis**
- Extracts sequence lengths from all amplified products
- Generates comprehensive statistics (min, max, mean lengths)
- Categorizes sequences by size thresholds for visualization

**6. Dot Plot Visualization**
- Creates scatter plots showing amplicon size distribution
- X-axis: Sequence length (bp, capped at 1000 bp for readability)
- Y-axis: Sequence number (ordered by input)
- Color coding: Blue (≤1000 bp), Red (>1000 bp)

**7. Statistical Summary**
- Total sequence count
- Size distribution breakdown
- Length statistics (min/max/mean)
- Visual quality assessment

**Example Amplicon Size Distribution Output**<br>
![coi 175f 345r_actinopteri_dot_plot](https://github.com/user-attachments/assets/bf3a5f7e-79a2-4f5b-b497-929f7857ac54)


## 3. Primer Evaluation (3_primer_evaluation.sh)

**Overview**<br>
This primer evaluation pipeline assesses the performance of DNA metabarcoding primers for taxonomic classification using QIIME 2 and RESCRIPt.
It estimates:
- **Binding capacity**: how many reference sequences are amplified by each primer set.
- **Taxonomic resolution**: how effectively each primer recovers taxonomic labels (e.g., kingdom to species levels).

The goal is to identify effective DNA metabarcoding primers by evaluating their binding coverage and taxonomic resolution, supporting accurate and efficient detection of target taxa in ecological studies such as diet analysis and biodiversity monitoring.<br>

**Requirements**
- QIIME 2 (with RESCRIPt plugin)
- Python 3 (with pandas, matplotlib, numpy)
- Input databases in QIIME 2 format (.qza files)

**Usage**
<pre><code># Basic syntax
./primer_evaluation.sh <gene> [primer_name] [class_name]

# Examples
./primer_evaluation.sh 12S                           # All 12S primers, all classes
./primer_evaluation.sh 16S 16S_MarVer3F_MarVer3R     # Specific 16S primer, all classes
./primer_evaluation.sh COI COI_VF2_FishR1 Reptile    # Specific primer and class</code></pre>

**Primer Evaluation Workflow Steps**<br>
**1) Initial Dereplication (1st derep)**
- Removes duplicate sequences from input database
- Maintains taxonomic assignments for unique sequences
- Uses `qiime rescript dereplicate` with 'uniq' mode

**2. Primer-Based Read Extraction**
- Extracts sequences that match forward and reverse primer sequences
- Uses 80% identity threshold for primer matching
- Simulates PCR amplification *in silico*

**3. Post-Extraction Dereplication (2nd derep)**
- Second round of deduplication after primer extraction
- Ensures only unique amplicon sequences remain

**4. Sequence Quality Control**
- **Cull sequences**: Removes sequences with excessive degenerates (>5) and long homopolymers (>8bp)
- **Length filtering**: Retains only sequences within expected amplicon size ranges
- Filters out sequences likely to cause classification errors

**5. Final Dereplication (3rd derep)**
- Third and final deduplication step
- Produces clean, non-redundant reference database

**6. Data Export and Processing**
- Exports sequences and taxonomy from QIIME 2 artifacts (.qza) to standard formats
- Creates FASTA files with integrated taxonomic information **(for phylogenetic tree construction and inspection)**
- Generates sequence statistics (count, length distribution)

**7. Performance Evaluation**

**7.1) Sequence Assessment**
- Evaluates sequence quality and composition
- Generates summary statistics and visualizations

**7.2) Taxonomic Coverage Analysis**
- Assesses taxonomic representation across different levels
- Identifies gaps in taxonomic coverage

**7.3) Cross-Validation Testing**
- Trains naive Bayes classifier using processed sequences
- Performs k-fold cross-validation to test classification accuracy
- Generates F-measure scores across taxonomic levels (kingdom → species)

**8. Results Visualization**
- **Performance plots**: F-measure scores across taxonomic levels
- **Error bar plots**: Include confidence intervals/standard errors
- **PNG exports**: Publication-ready visualizations of sequence quality metrics

**Output Files**

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

**Example Primer Evaluation Output**<br>
![image](https://github.com/user-attachments/assets/a2b9f102-2ad3-4f2d-9dcf-6b27d7447665)

**F-measure scores** range from 0-1, where:
- **> 0.9**: Excellent classification performance
- **0.7-0.9**: Good performance  
- **0.5-0.7**: Moderate performance
- **< 0.5**: Poor performance, primer may not be suitable

Performance typically decreases from kingdom level (highest) to species level (lowest), which is expected due to increased taxonomic resolution requirements.

## 4. Blocking Primer Design and Evaluation


