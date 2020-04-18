# Summary

In this repository I am conducting a DEG analysis of *Solanum lycopersicum* cv. Ailsa Craig, *Nicotiana obtusifolia*, and *Arabidopsis thaliana*. I will supplement this DEG analysis with publicly available ChIP-chip and ChIP-seq data from *S. lycopersicum* and *A. thaliana* as well as publicly available RNA-seq data from *FRUITFULL* mutants of both speices.

# Methods

## DEG Analysis

For RNA-seq, this will consist of STAR as the read aligner and DEseq2 as the differential expression testing software. For the microarray analyses, as a preprocessing step, the probe sequences were mapped back to transcriptomic targets, and nonmapping or multimapping probes were excluded from the analyes. Differential expression testing will be done with limma/voom. In all caases, the data are fit to a spline regression model with one fewer degree of freedom than there are timepoints in the specific study. (See supplement of [Sander et al, 2017](https://www.ncbi.nlm.nih.gov/pubmed/27797772) for an implementation of this with DESeq2.)

### In-House Data
 - I am working on uploading these to NCBI so that the whole thing will be with publicly available data

Our lab has generated several transcriptome libraries that will be used for this analysis. They are symlinked to in the `SlycRNA`, `SpimpRNA`, and `NobtRNA` for tomato and tobacco, respecitively

### Public Data

#### Gene Expression Data

The RNAseq and Microarray files were all downloaded with the `1_GetData.sh` and placed in the `ExternalData` folder.

| Species | NCBI BioProject ID | Description |
| ------- | ------------------ | ----------- |
| Tomato | [PRJNA213528](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA213528)  | RNA-seq from a developmental series of fruits using wild-type and *ful1*/*ful2* knockdown lines. From Fujusawa et al, 2014 |
| Tomato | [PRJNA177429](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA177429) | Microarray from a developmental series of fruits using wild-type and *ful1*/*ful2* knockdown lines. From Bemer et al, 2010 |
| Arabidopsis | [PRJEB25745](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB25745) | RNA-seq of wild-type fruit valve tissue. From Mizzotti et al, 2018 | 
| Arabidopsis | [PRJNA316153](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA316153) | Unpublished(?) microarray of *ful-1* FUL-GR siliques under both mock and DEX treatment. Attached to Bemer et al, 2017, but not reported. |

#### Genomes

I am currently using the [SL4.0 genome and the ITAG4.0 annotation for tomato](https://solgenomics.net/organism/Solanum_lycopersicum/genome). For tobacco (*N. obtusifolia*) I am using the publicly available [genome and annotation](http://nadh.ice.mpg.de/NaDH/download/overview) but with a slight modification to add in the the MBP20 gene, which was missed in the original annotation. For *Arabidopsis* I'll be using the TAIR10 genome and annotation. In the case of tomato and tobacco these files are symlinked to copies we already have using the `SlyDNA` and `NobtDNA` directories, respecitively. 

For Arabidopsis, the data is downloaded with the `1_GetData.sh` to `ExternalData/TAIR10`, but I have excluded those files from the git to save space.

# ChIP Analysis

The ChIP-seq files are downloaded to `ExternalData/ChIPseq` and the microarray file is downloaded to `ExternalData/Microarray` both with `1_GetData.sh`. As a preprocessing step for the microarray study, array probe sequences were mapped to the SL4.0 genome, and probe locations were updated. Nonmapping or multimapping probes were excluded from the analysis.

### Public Data

| Species | NCBI BioProject ID | Description |
| ------- | ------------------ | ----------- |
| Tomato | [PRJNA213106](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA213106) | ChIP-chip of wild-type fruits using anti-FUL1 or anti-FUL2 antibodies. From Fujisawa et al, 2014 |
| Arabidopsis | [PRJNA316152](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA316152) | ChIP-seq of pistil/siliques of *ful-1 pFUL*:*FUL-GFP* using an anti-GFP antibody. From Bemer et al, 2017 |
| Arabidopsis | [PRJNA427320](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA427320) | ChIP-seq of pistil/siliques of *ful-1 pFUL*:*FUL-GFP* using an anti-GFP antibody. From Balanza et al 2018 |



