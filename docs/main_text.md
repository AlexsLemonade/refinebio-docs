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

### RNA-seq pipelines

### Microarray pipelines

#### Affymetrix

#### Illumina BeadArrays

## Submitter Processed

## Aggregations

## Transformations

### Quantile normalization

### Gene transformations


# CCDL Developed Packages

## Illumina SCAN

## Identifier Refinery


# Species compendia


# Use Cases for Downstream Analysis
