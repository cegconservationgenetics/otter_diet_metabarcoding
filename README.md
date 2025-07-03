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
