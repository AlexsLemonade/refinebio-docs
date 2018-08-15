# Source Data

## Types of Data

The current version of [refine.bio](https://refine.bio) is designed to process gene expression data.
This includes both microarray data and RNA-seq data.
We normalize data to Ensembl gene identifiers and provide abundance estimates.

More precisely, we support microarray platforms based on their GEO or ArrayExpress accessions.
We currently support Affymetrix microarrays and Illumina BeadArrays, and we are continuing to evaluate and add support for more platforms.
This [table](https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_microarray_platforms.csv) contains the microarray platforms that we support.
We process a subset of platforms using the [BrainArray Custom CDFs](http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp), which are denoted by a `y` in the `is_brainarray` column.
We also support RNA-seq experiments performed on [these](https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt) short-read platforms.
For more information on how data are processed, see the [Processing Information](#processing-information) section of this document.
If there is a platform that you would like to see processed, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues).
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

## Sources

We download gene expression information and metadata from [EBI's ArrayExpress](https://www.ebi.ac.uk/arrayexpress/), [NCBI's Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), and [NCBI's Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra).
NCBI's SRA also contains experiments from EBI's ENA ([example](https://www.ncbi.nlm.nih.gov/sra/?term=ERP000447)) and the DDBJ ([example](https://www.ncbi.nlm.nih.gov/sra/?term=DRP000017)).

## Metadata

We provide metadata that we obtain from the source repositories.
We also attempt to, where possible, perform some harmonization of metadata to improve the discoverability of useful data.

### Submitter Supplied Metadata

TODO: How this will be done & delivered is up in the air. I notice that it's on the whiteboard to tackle before release. Thus I'm just tagging this with a TODO and [linking the issue](https://github.com/AlexsLemonade/refinebio/issues/515) in question.

Protocol information, which generally contains the type of sample preparation, is handled differently by different databases.
We push this information down to the sample level and provide it that way.
In the case of data imported from ArrayExpress, the protocol _may_ be the information from the full experiment and not just the sample in question.

### refine.bio-harmonized Metadata

Scientists who upload results don't always use the same names for related values.
This makes it challenging to search across datasets.
We have put some processes in place to smooth out some of these issues.

To produce lightly harmonized metadata, we combine certain fields based on similar keys.
For example, `treatment`, `treatment group`, `treatment protocol`, `drug treatment`, and `clinical treatment` fields get collapsed down to `treatment`.
The fields that we currently collapse to includes `specimen part`, `genetic information`, `disease`, `disease stage`, `treatment`, `race`, `subject`, `development stage`, `compound`, and `time`.
We do this for convenience and to aid in searches.
We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata.
If you find that the harmonized metadata does not accurately reflect the metadata supplied by the submitter, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues) so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

# Downloadable Files

## Gene Expression Matrix

## Sample Metadata

## Experiment Metadata


# Processing Information

## refine.bio processed

### Microarray pipelines

#### Affymetrix

#### Illumina BeadArrays

### RNA-seq pipelines

We use [Salmon](https://combine-lab.github.io/salmon/) and [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) to process all RNA-seq data in refine.bio.
We obtain fastq files run on our [supported short-read platforms](https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt) from Sequence Read Archive. 
We use the library strategy and library source metadata fields to identify RNA-seq experiments.
It's possible that experiments that are inappropriate for use with Salmon will still appear in refine.bio (e.g., long-read platforms that are labeled incorrectly in the source repository).
If you find an experiment that you believe is inappropriate for use with Salmon, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues) so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

#### Salmon 

Salmon is an alignment-free method for estimating transcript abundances from RNA-seq data. 
We use it in [quasi-mapping mode](http://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-quasi-index-and-fmd-index-based-modes), which is significantly faster than alignment-based approaches and requires us to build a Salmon transcriptome index. 

##### Transcriptome index

We build a custom reference transcriptome (using [RSEM](https://github.com/deweylab/RSEM) `rsem-prepare-reference`) by filtering the Ensembl genomic DNA assembly to remove _pseudogenes_, which we expect could negatively impact the quantification of protein-coding genes. 
This means we're obtaining abundance estimates for coding as well as non-coding transcripts. 

Building a transcriptome index with `salmon index` requires us to specify a value for the parameter `-k` that determines the size of the k-mers used for the index.
The length of a read determines what k-mer size is appropriate. 
Consistent with the recommendations of the authors of Salmon, we use an index build with _k_ = 31 when quantifying samples with reads with length > 75bp.
We use _k_ = 23 for all other read lengths.

TODO: Note about how to obtain Salmon indices. Do we also make gene to transcript mapping passed to `rsem-prepare-reference` available?

##### Quantification with Salmon

When quantifying transcripts with `salmon quant`, we take advantage of options that allow Salmon to learn and attempt to correct for certain biases in sequencing data. 
We include the flags `--seqBias` to correct for random hexamer priming and, if this is a **paired-end** experiment, `--gcBias` to correct for GC content when running salmon quant. 
We set the library type parameter such that Salmon will [infer the sequencing library type automatically](http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype) for the reads it is quantifying (`-l A`).

#### tximport

Salmon quantification is at the _transcript-level_. 
To better integrate with the microarray data contained in refine.bio, we summarize the transcript-level information to the _gene-level_ with `tximport` ([Soneson, Love, and Robinson. _F1000 Research._ 2015.](http://dx.doi.org/10.12688/f1000research.7563.1)).

Our tximport implementation generates ["lengthScaledTPM"](https://www.rdocumentation.org/packages/tximport/versions/1.0.3/topics/tximport), which are gene-level count-scale values that are generated by scaling TPM using the average transcript length across samples and to the library size. 
Note that tximport is applied at the _experiment-level_ rather than to single samples. 
For additional information, see the [tximport Bioconductor page](http://bioconductor.org/packages/release/bioc/html/tximport.html), the [tximport tutorial _Importing transcript abundance datasets with tximport_](http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html), and [Soneson, Love, and Robinson. _F1000Research._ 2015.](http://dx.doi.org/10.12688/f1000research.7563.1).

TODO: If we make the gene to tx mapping available (see above under **Transcriptome index**), we might put more information about that here in this section. 
 
## Submitter processed

## Aggregations

## Transformations

### Quantile normalization

### Gene transformations


# CCDL Developed Packages

## Illumina SCAN

## Identifier Refinery


# Species compendia


# Use Cases for Downstream Analysis
