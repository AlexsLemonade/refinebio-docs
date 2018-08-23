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

![sources](https://user-images.githubusercontent.com/15315514/44533218-08b95f80-a6c3-11e8-9eb8-827086fd2d99.png)

## Metadata

We provide metadata that we obtain from the source repositories.
We also attempt to, where possible, perform some harmonization of metadata to improve the discoverability of useful data.
Note that we do not yet obtain sample metadata from the [BioSample](https://www.ncbi.nlm.nih.gov/biosample/) database, so the metadata available for RNA-seq samples is limited.

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
We do this for convenience and to aid in searches.
For example, `treatment`, `treatment group`, `treatment protocol`, `drug treatment`, and `clinical treatment` fields get collapsed down to `treatment`.
The fields that we currently collapse to includes `specimen part`, `genetic information`, `disease`, `disease stage`, `treatment`, `race`, `subject`, `development stage`, `compound`, and `time`.
See the table below for a complete set of mappings between the keys from source data and the harmonized keys.
Values are stripped of white space and forced to lowercase.

| Harmonized key | Keys from data sources |
|:----------------:|-------------------------|
| `specimen part` | `organism part`, `cell type`, `tissue`, `tissue type`, `tissue source`, `tissue origin`, `source tissue`, `tissue subtype`, `tissue/cell type`, `tissue region`,  `tissue compartment`,  `tissues`, `tissue of origin`, `tissue-type`,  `tissue harvested`, `cell/tissue type`, `tissue subregion`, `organ`, `characteristic [organism part]`, `characteristics [organism part]`, `cell_type`, `organismpart`, `isolation source`, `tissue sampled`, `cell description`
| `genetic information` | `strain/background`, `strain`,  `strain or line`, `background strain`, `genotype`, `genetic background`, `genotype/variation`, `ecotype`, `cultivar`, `strain/genotype`|
| `disease` |  `disease `, `disease state `, `disease status `, `diagnosis `, `disease `, `infection with `, `sample type ` |
| `disease stage` | `disease state `, `disease staging `, `disease stage `, `grade `, `tumor grade `,  `who grade `, `histological grade `, `tumor grading `, `disease outcome `, `subject status ` |
| `treatment` | `treatment`, `treatment group`, `treatment protocol`,  `drug treatment`, `clinical treatment` |
| `race` | `race`, `ethnicity`, `race/ethnicity`|
| `subject` |  `subject `, `subject id `, `subject/sample source id `, `subject identifier `, `human subject anonymized id `, `individual `, `individual identifier `,  `individual id `, `patient `, `patient id `, `patient identifier `,  `patient number `, `patient no `,  `donor id `, `donor `, `sample_source_name `|
| `development stage` | `developmental stage`,  `development stage`, `development stages` |
| `compound` | `compound`, `compound1`, `compound2`, `compound name`, `drug`, `drugs`, `immunosuppressive drugs` |
| `time` | `initial time point`, `start time`, `stop time`, `time point`, `sampling time point`, `sampling time`, `time post infection` |
| `age` | `age`, `patient age`, `age of patient`, `age (years)`, `age at diagnosis`, `age at diagnosis years`, `characteristic [age]`, `characteristics [age]` |

We type-cast age values to doubles.
Sex is a special case; we map to `female` and `male` values if the values are one of the following:

| Harmonized `sex` value | Values |
|:-----------------------:|-------|
| `female` | `f`, `female`, `woman`|
| `male` | `m`, `male`, `man` |

We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata. If you find that the harmonized metadata does not accurately reflect the metadata supplied by the submitter, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues) so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

# Processing Information

## refine.bio processed

Because refine.bio is designed to be consistently updated, we use processing and normalization methods that operate on a single sample wherever possible.
Processing and normalization methods that require multiple samples (e.g., Robust Multi-array Average or RMA) generally have outputs that are influenced by whatever samples are included and rerunning these methods whenever a new sample is added to the system would be impractical.

### Microarray pipelines

#### Affymetrix

SCAN (Single Channel Array Normalization) is a normalization method for develop for single channel Affymetrix microarrays that allows us to process individual samples. 
SCAN models and corrects for the effect of technical bias, such as GC content, using a mixture-modeling approach. 
For more information about this approach, see the primary publication ([Piccolo, et al. _Genomics._ 2012.](http://dx.doi.org/10.1016/j.ygeno.2012.08.003)) and the [SCAN.UPC Bioconductor package](https://www.bioconductor.org/packages/release/bioc/html/SCAN.UPC.html) documentation. 
We specifically use the `SCANfast` implementation of SCAN and the Brainarray packages as probe-summary packages when available.

##### Platform detection

We have encountered instances where the platform label from the source repository and the metadata included in the sample's raw data file (`.CEL` file) itself do not match.
In these cases, we take the platform information included in the raw data (`.CEL`) file header to be the true platform label.

#### Illumina BeadArrays

Dr. Stephen Piccolo, the developer of SCAN, has adapted the algorithm for use with Illumina BeadArrays for refine.bio. 
Because this Illumina SCAN methodology is not yet incorporated into the SCAN.UPC package, we briefly summarize the methods below.

We require that non-normalized or raw expression values and detection p-values to be present in Illumina non-normalized data.
If we infer that background correction has not occurred in the non-normalized data (e.g., there are no negative expression values), the data are background corrected using the [`limma::nec`](http://web.mit.edu/%7Er/current/arch/i386_linux26/lib/R/library/limma/html/nec.html) function ([Shi, Oshlack, and Smyth. _Nucleic Acids Research._ 2010.](https://doi.org/10.1093/nar/gkq871)). 
Following background correction -- either upstream presumably in the Illumina BeadStudio software or in our processor, arrays are normalized with SCAN.
SCAN requires probe sequence information obtained from the [Illumina BeadArray Bioconductor annotation packages](https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip) (e.g., [`illuminaHumanv1.db`](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html)).
We only retain probes that have a "Good" or "Perfect" rating in these packages; this quality rating is in reference to how well a probe is likely to measure its target transcript.

##### Platform detection

We infer the Illumina BeadArray platform that a sample is likely to be run on by comparing the probe identifiers in the unprocessed file to probes for each of the Illumina expression arrays for a given organism. 
We again use the Illumina Bioconductor annotation packages for this step.
For instance, the overlap between the probe identifiers in a human sample and the probe identifiers in each human platform ([`v1`](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html), [`v2`](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv2.db.html), [`v3`](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv3.db.html), and [`v4`](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html)) is calculated.
The platform with the highest overlap (provided it is >75%) is inferred to be the true platform.
Some analyses around this platform detection procedure can be found in [this repository](https://github.com/jaclyn-taroni/beadarray-platform-detection).

### RNA-seq pipelines

TODO: add pipeline illustrations from processing information modals?

We use [Salmon](https://combine-lab.github.io/salmon/) and [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) to process all RNA-seq data in refine.bio.
We obtain fastq files run on our [supported short-read platforms](https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt) from Sequence Read Archive. 
We use the library strategy and library source metadata fields to identify RNA-seq experiments.
It's possible that experiments that are inappropriate for use with Salmon will still appear in refine.bio (e.g., long-read platforms that are labeled incorrectly in the source repository).
If you find an experiment that you believe is inappropriate for use with Salmon, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues) so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

#### Salmon 

Salmon is an alignment-free method for estimating transcript abundances from RNA-seq data ([Patro, et al. _Nature Methods_. 2017.](http://dx.doi.org/10.1038/nmeth.4197)). 
We use it in [quasi-mapping mode](http://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-quasi-index-and-fmd-index-based-modes), which is significantly faster than alignment-based approaches and requires us to build a Salmon transcriptome index. 

##### Transcriptome index

We build a custom reference transcriptome (using [RSEM](https://github.com/deweylab/RSEM) `rsem-prepare-reference`) by filtering the Ensembl genomic DNA assembly to remove _pseudogenes_, which we expect could negatively impact the quantification of protein-coding genes. 
This means we're obtaining abundance estimates for coding as well as non-coding transcripts. 

Building a transcriptome index with `salmon index` requires us to specify a value for the parameter `-k` that determines the size of the k-mers used for the index.
The length of a read determines what k-mer size is appropriate. 
Consistent with the recommendations of the authors of Salmon, we use an index build with _k_ = 31 when quantifying samples with reads with length > 75bp.
We use _k_ = 23 for shorter read lengths.

The refine.bio processed Salmon indices are available for download.
You can make use of our API like so:

```
https://api.staging.refine.bio/transcriptome_indices/?organism=<ORGANISM>&length=<LENGTH>
``` 
Where `ORGANISM` is the scientific name of the species in all caps separated by underscores and `LENGTH` is either `SHORT` or `LONG`. 

To obtain the zebrafish (_Danio rerio_) index used for >75bp reads, use:

```
https://api.staging.refine.bio/transcriptome_indices/?organism=DANIO_RERIO&length=LONG
```

The `s3_url` field will allow you to download the index.

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
 
TODO: Caveats for samples that are part of multiple experiments?

## Submitter processed

Sometimes raw data for a sample is either unavailable at the source repository or exists in a form that we can not process.
For microarray platforms that we support, we obtain the submitter processed expression data and use these values in refine.bio with some modification (e.g., log2-transformation where we detect it has not been performed).

As noted above, we use Ensembl gene identifiers throughout refine.bio. 
Submitter processed data may use other gene (or probe) identifiers that we must convert to Ensembl gene identifiers. 
We describe the processes for Affymetrix and Illumina data below. 
Note in the case of one-to-many mappings when going from the ID used by the submitter to the Ensembl gene ID, expression values are duplicated: 
if a probe maps to two Ensembl gene ids, those two Ensembl gene ids will both have the probe's expression value following conversion.

### Affymetrix

We have created custom gene mapping files for most of the Affymetrix platforms we support. 
Briefly, for Brainarray supported platforms, we use the Brainarray (e.g., `hgu133plus2hsensgprobe`) and the platform-specific annotation package from Bioconductor (e.g., `hgu133plus2.db`) to generate a platform-specific mapping file that includes probe IDs, Ensembl gene IDs, gene symbols, Entrez IDs, RefSeq and Unigene identifiers. 
The rationale for only using probes or IDs that are accounted for in the Brainarray package is two-fold: 1) Brainarray packages are updated as we learn more about the genome and 2) it allows for these submitter processed data to be more consistent with refine.bio processed data. 
We support identifier conversion for a limited number of platforms that either do not have a Brainarray or Bioconductory annotation packages.

The code for deriving these mappings and more details are available at https://github.com/AlexsLemonade/identifier-refinery.
If you find an issue with these mappings, please [file an issue on GitHub](https://github.com/AlexsLemonade/identifier-refinery/issues) so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

### Illumina

We support conversion from Illumina BeadArray probe IDs to Ensembl gene IDs using 
[Bioconductor Illumina Beadarray expression packages](https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip), 
allowing for one-to-many mappings.

TODO: Filtering for probe quality?

## Aggregations

refine.bio allows users to aggregate their selected samples in two ways: by experiment or by species.
We use the term aggregate or aggregation to refer to the process of combining _individual samples_ to form a _multi-sample_ gene expression matrix (see also: [Downloadable Files](#downloadable-files)).

* **By experiment:** Samples that belong to the same experiment will become a single gene expression matrix. 
If you have selected all samples from two experiments with 10 and 15 samples, respectively, and have chosen the `by experiment` option, you will receive two gene expression matrices with 10 and 15 samples, respectively.

* **By species:** All samples assaying the same species will be aggregated into a single gene expression matrix.
If you have selected three experiments each from human and mouse and the `by species` option, you receive two gene expression matrices that contain all human and all mouse samples, respectively.

For either aggregation method, we summarize duplicate Ensembl gene IDs to the mean expression value and only include genes (rows) that are represented in **all** samples being aggregated.
This is also known as an inner join and is illustrated below.
![inner join](https://user-images.githubusercontent.com/15315514/44534751-7a46dd00-a6c6-11e8-9760-e8daa91a500f.png)
Note that some early generation microarrays measure fewer genes than their more recent counterparts, so their inclusion when aggregating `by species` may result in a small number of genes being returned. 


The aggregation methodology for species compendia is different; see [Species compendia](#species-compendia) for more information.

TODO: illustration of structure of download zip file

## Transformations

### Quantile normalization

TODO: Linking to QN [issue](https://github.com/AlexsLemonade/refinebio/issues/488) and [PR](https://github.com/AlexsLemonade/refinebio/pull/519) for now.

### Gene transformations

In some cases, it may be useful to row-normalize or transform the gene expression values in a matrix (e.g., following aggregation and quantile normalization). 
We offer the following options for transformations:

* **None:** No row-wise transformation is performed.

* **Z-score:** Row values are [z-scored](https://en.wikipedia.org/wiki/Standard_score) using the [`StandardScaler`](http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html) from [`scikit-learn`](http://scikit-learn.org/stable/index.html). 
This transformation is useful for examining samples' gene expression values relative to the rest of the samples in the  expression matrix (either all selected samples from that _species_ when aggregating by species or all selected samples in an _experiment_ when aggregating by experiment). 
If a sample has a positive value for a gene, that gene is more highly expressed in that sample compared to the mean of all samples; if that value is negative, that gene is less expressed compared to the population.
It assumes that the data are normally distributed.

* **Zero to one:** Rows are scaled to values `[0,1]` using the [`MinMaxScaler`](http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html) from [`scikit-learn`](http://scikit-learn.org/stable/index.html).
We expect this transformation to be most useful for certain machine learning applications (e.g., those using cross-entropy as a loss function).

In the plot below, we demonstrate the effect of different scaling options on gene expression values (using a randomly selected human dataset, microarray platform, and gene):

<img src="https://user-images.githubusercontent.com/19534205/44432215-1a89ee00-a56f-11e8-9327-9b5cca438e39.png" width="480">

Note that the distributions retain the same general _shape_, but the range of values and the density are altered by the transformations.

# Downloadable Files

Users can download gene expression data and associated sample and experiment metadata from refine.bio.
Below we describe the files included in the delivered zip file.

## Gene Expression Matrix

Gene expression matrices are delived in [tab-separated value](https://en.wikipedia.org/wiki/Tab-separated_values) (TSV) format.
In these matrices, rows are _genes_ or _features_ and columns are _samples_.
Note that this format is consistent with the input expected by many programs specifically designed for working with gene expression matrices, but some machine learning libraries will expect this to be transposed.
The column names or header will contain values corresponding to sample titles. 
You can use these values in the header to map between a sample's gene expression data and its metadata (e.g., disease label or age). See also [Use Cases for Downstream Analysis](#use-cases-for-downstream-analysis).

TODO: I would like to write the metadata sections after working with some data once https://github.com/AlexsLemonade/refinebio/pull/526 lands.

## Sample Metadata

## Experiment Metadata


# Species compendia

TODO: Note about release schedule and/or presence on Zenodo

refine.bio periodically releases compendia comprised of all the samples from a species that we were able to process.
We refer to these as **species compendia**. 
We process these compendia in a manner that is different from the options that are available via the web user interface. 
Instead of selecting only genes available in all samples, we take the union of all genes, filling in any missing values with `NA` (e.g., perform a full outer join as illustrated below).

![outer join](https://user-images.githubusercontent.com/15315514/44534241-4dde9100-a6c5-11e8-8a9c-aa147e294e81.png)


We drop any genes that have missing values in greater than 30% of samples.
We impute the remaining missing values with KNN impute.
We then quantile normalize all samples as described above.

TODO: More information about QN


# Use Cases for Downstream Analysis
