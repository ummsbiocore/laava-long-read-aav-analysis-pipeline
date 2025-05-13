LAAVA: Long-read AAV Analysis pipeline is designed for automated analysis of long-read sequencing data from adeno-associated virus (AAV) products.

The sequencing data should be from a PacBio sequencer (Sequel II, Sequel IIe or Revio) run with on-instrument AAV mode enabled, or equivalent circular consensus sequencing (CCS) reads (Travers et al., 2010)

In this analysis, reads are aligned to the given AAV, packaging, and host reference sequences using Minimap2 (Li, 2018). The reference sequences for each primary alignment and its orientation are counted and summarized to assign read type classifications, including vector, non-vector, and chimeric reads. For reads assigned to the AAV vector, the primary alignment coordinates are compared to the annotated vector region in the reference sequence, which comprises the left and right ITRs and the genomic region between them, to assign each read to a subtype classification. Sequence variants relative to the vector reference sequence are determined directly from each read's alignment, specifically the CIGAR string indicating insertions, deletions, mismatches, and gaps.

Finally, a report is generated with relevant quality metrics and analysis results in both HTML and PDF formats.


### Overview

LAAVA: Long-read AAV Analysis pipeline is designed for the analysis of recombinant adeno-associated virus (rAAV) products by PacBio long-read sequencing. It provides standardized nomenclature, rigorous QC, and comprehensive reporting to ensure comparability across production runs and inform vector design and quality control .

### Key Use cases

* **Vector Genome Integrity Assessment:** Quantify proportions of full‑length, partial, and truncated ssAAV/scAAV genomes to evaluate functional payload yield.
* **Contaminant Detection & Quantification:** Detect and report host‑cell DNA, RepCap/helper plasmid carry‑over, backbone fragments, and chimeric reads.
* **Flip/Flop Configuration Analysis:** Determine ITR orientation (“flip” vs. “flop”) distributions via local alignment (Parasail) to assess ITR integrity and QC.

### Features

* **Multi‑Reference Alignment & Classification:** Maps HiFi reads with minimap2 to assign each read a “type” (ssAAV, scAAV, host, RepCap, helper, chimeric, etc.) and “subtype” (full, left‑partial, right‑partial, vector+backbone, snapback) based on CIGAR patterns.
* **Flip/Flop & Structural Variant Detection:** Reports size distributions and variant hotspots (insertions, deletions) across the vector.
* **Scalability:** Can process [define throughput, e.g., thousands of samples] in parallel using [e.g., cloud-based execution, HPC, Nextflow].
* **Reproducibility:** Fully containerized via Docker, ensuring consistent results across environments.
* **Comprehensive QC & Integrity Checks:** Automated MultiQC–style summaries integrated into final report.
* **Automated Reporting:** Generates HTML/PDF reports with interactive plots. 

### Inputs
#### File inputs
**inputs**
* Description: PacBio AAV sequencing read set, as HiFi/CCS reads in FASTQ or unaligned BAM format. The PacBio instrument should be run in AAV mode..
* Format: .fastq.gz / .bam
* Example File Path: /samples/ss.subsample005.bam

**vector_fasta**
* Description: Vector plasmid, as a single-record FASTA.
* Format: .fasta
* Example File Path: /samples/ss.construct.fasta

**vector_bed**
* Description: Annotated vector construct region coordinates in 4-column UCSC BED format.
* This file must indicate the transgene/payload region via either two labeled Inverted Terminal Repeat (ITR) regions (see itr_label_1 and itr_label_2 below) or, as a legacy mode, one region with the label 'vector', spanning both ITRs (inclusive).
* May also include additional labeled regions, e.g. for promoter and CDS regions; these will be ignored and will not affect the output.
* Format: .bed
* Example File Path: /samples/ss.construct.bed


#### Optional file inputs
**packaging_fa**
* Description: Packaging sequences -- helper and rep/cap plasmids, and other sequences of interest (e.g. Lambda), as a multi-record FASTA.
* Format: .fa / .fasta

**host_fa**
* Description: Host genome (recommended), as a multi-record FASTA. Best to include only the canonical chromosomes and not the alternative contigs.
* Format: .fa / .fasta

**flipflop_fa**
* Description: Flip/flop ITR sequences.
* AAV2 sequences are built in and available by default; provide custom sequences here to use another serotype (and use a different value for the "AAV Serotype" option).
* Format: .fa / .fasta


#### Label inputs
**ITR labels used in the vector annotation BED file: itr_label_1, itr_label_2, mitr_label (Optional)**
* These are case-sensitive and must match exactly.
* The order does not matter: LAAVA will check for the presence of all the given labels in the annotation BED file, and treat the first match as the left (5') ITR and the second as the right (3') ITR. If the annotation uses the same label for both, e.g. 'ITR', you only need to specify itr_label_1 as 'ITR', and that label will correctly match both ITR regions in the annotation BED.
* For scAAV vector constructs, the mutant ITR (mITR) should be specified with mitr_label instead of itr_label_2. The location of this ITR relative to the payload (left or right, 5' or 3') will be used in classification to distinguish "itr-partial" (originating in the wild-type ITR, likely packaged) from "unclassified (originating in the mITR, likely an artifact).
* If these fields are left blank, LAAVA will fall back to the legacy mode of looking for a region labeled "vector" instead.

**Sequence IDs used in the packaging FASTA file: repcap_name, helper_name, lambda_name (Optional)**
* These are case-sensitive and must match exactly.
* If not specified, the reference sequence IDs in the packaging FASTA file will be counted and reported with their original names. If specified but not found in the packaging FASTA file, LAAVA will raise an error.
* If repcap_name is specified and found in the packaging FASTA, certain plots related to RepCap read alignments will be included in the output report.

**ITR AAV Serotype name (flipflop_name)**
* Selects a set of built-in ITR sequences for flip/flop analysis. Currently, the only built-in set of ITR sequences is for the AAV2 serotype.
* Specifying 'AAV2' here is equivalent to providing the same sequences as a FASTA file via the "flip/flop ITR sequences" input field above.

### Clair3 inputs
**Platform**
* Select the sequencing platform of the input. 
* Possible options: {ont, hifi, ilmn}.

**model_name**
* Selects a pre-trained model to use in variant calling.

|           Model name           |  Platform   | Option (`-p/--platform`) |                       Training samples                       |  Date   |
| :----------------------------: | :---------: | :----------------------------------------------------------: | -------------------------------- |:------: |
|      r941_prom_sup_g5014       |     ONT r9.4.1     |     `ont`     |                    HG002,4,5 (Guppy5_sup)                    | 20220112 |
|    r941_prom_hac_g360+g422     |     ONT r9.4.1    | `ont`    |                         HG001,2,4,5                          | 20210517 |
|       r941_prom_hac_g238       |     ONT r9.4.1    | `ont`    |                         HG001,2,3,4                          | 20210627 |
|              hifi_revio              | PacBio HiFi Revio | `hifi` |                         HG002,4                         | 20230522 |
|             hifi_sequel2             | PacBio HiFi Sequel II | `hifi` |                         HG001,2,4,5                          | 20210517 |
| ilmn | Illumina | `ilmn` | HG001,2,4,5 | 20210517 |

**opt_parameters**
* Optional parameters to use in Clair3. For all possible options: [Clair3](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#options)

### Outputs
#### Published outputs
**html_report**
* Description: HTML report with various plots such as read-type proportions, read-length histograms, etc.
* Format: .html
* Example File Path: /make_report/ss_report.html
* Visualization App: HTML
* Location: make_reports

**pdf_report**
* Description: PDF report with various plots such as read-type proportions, read-length histograms, etc.
* Format: .pdf
* Example File Path: /make_report/ss_report.pdf
* Visualization App: PDF reader
* Location: make_reports

**variants**
* Description: Variants in VCF version 4.2 format.
* Format: .vcf
* Example File Path: /variants/ss.vcf.gz
* Location: Clair3

## Input Files and Parameters for a Demo Run:
### Case 1: Run the example single-stranded AAV (ssAAV) sample.
* reads: "ss.subsample005.bam"
* vector_type: "ss"
* vector_bed: "ss.annotation.bed"
* vector_fa: "ss.construct.fasta"
* itr_label_1: "ITR"
* itr_label_2: "ITR"

Example bed file(ss.annotation.bed) content:
```
pAV_CMV_GFP     0       145     ITR
pAV_CMV_GFP     0       2434    vector
pAV_CMV_GFP     493     1068    CMV_promoter
pAV_CMV_GFP     1191    1909    GFP_CDS
pAV_CMV_GFP     2033    2220    polyA
pAV_CMV_GFP     2289    2434    ITR
```


### Case 2: Run the example self-complementary AAV (scAAV) sample.
* reads: "sc.subsample005.bam"
* vector_type: "sc"
* vector_bed: "sc.annotation.bed"
* vector_fa: "sc.construct.fasta"
* itr_label_1: "wtITR"
* mitr_label: "mITR"

Example bed file(sc.annotation.bed) content:
```
dsCB-GFP        661     806     wtITR
dsCB-GFP        661     2739    vector
dsCB-GFP        1604    2321    GFP_CDS
dsCB-GFP        2594    2739    mITR
```

*For both cases, default values should be selected for other parameters*
