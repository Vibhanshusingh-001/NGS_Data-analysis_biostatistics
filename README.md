# NGS_Data-analysis_biostatistics
Analyzing complex datasets
### Task 1: Data Handling and Statistical Analysis

#### Objective: Assess candidates’ ability to handle complex data and apply statistical methods effectively

#### Background: CpG methylation is an epigenetic marker that varies across tissue types. However, the methylation status of a single CpG site is unreliable as a biomarker due to
#### errors introduced by bisulfite sequencing, sampling techniques, and biological variability.

#### Definition: Phased Methylation Pattern (PMP) is a unique set of coordinates that includes the DNA strand (‘f’ for forward (+) or ‘r’ for reverse (-)), the relative positions of three CpG
#### sites on the same strand (e.g., x:y:z), and their methylation status (e.g., ‘000’ for all unmethylated or ‘111’ for all methylated). It represents a combined epigenetic signature 

#### Hypothesis: Phased methylation patterns (PMPs) can act as reliable biomarkers to differentiate tissue types, providing higher specificity compared to individual CpG sites.

#### Dataset: The dataset (Link to Data) summarizes phased methylation patterns from NGS results across two tissues. Key columns include:

#### Strand: Indicates the DNA strand (‘f’ or ‘r’).
#### CpG Coordinates: Relative positions of three CpG sites (x:y:z).
#### Methylation Status: Eight possible patterns (‘000’ to ‘111’).
#### Sample ID: Unique identifier for each sample.
#### Replicate: Indicates technical replicates.
#### Tissue: Tissue type (Tissue #1 or Tissue #2).
