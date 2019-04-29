# Source Data

## Types of Data

The current version of <a href ="https://www.refine.bio" target = "blank">refine.bio </a>is designed to process gene expression data.
This includes both microarray data and RNA-seq data.
We normalize data to Ensembl gene identifiers and provide abundance estimates.

More precisely, we support microarray platforms based on their GEO or ArrayExpress accessions.
We currently support Affymetrix microarrays and Illumina BeadArrays, and we are continuing to evaluate and add support for more platforms.
This <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_microarray_platforms.csv" target = "blank">table </a> contains the microarray platforms that we support.
We process a subset of platforms using the <a href = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp" target = "blank">BrainArray Custom CDFs</a>, which are denoted by a `y` in the `is_brainarray` column.
We also support RNA-seq experiments performed on <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">these</a>  short-read platforms.
For more information on how data are processed, see the [Processing Information](#processing-information) section of this document.
If there is a platform that you would like to see processed, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a>.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

## Sources

We download gene expression information and metadata from <a href = "https://www.ebi.ac.uk/arrayexpress/" target = "blank">EBI's ArrayExpress</a>, <a href = "https://www.ncbi.nlm.nih.gov/geo/" target = "blank">NCBI's Gene Expression Omnibus (GEO)</a>, and <a href = "https://www.ncbi.nlm.nih.gov/sra" target = "blank">NCBI's Sequence Read Archive (SRA)</a>.
NCBI's SRA also contains experiments from EBI's ENA (<a href = "https://www.ncbi.nlm.nih.gov/sra/?term=ERP000447" target = "blank">example</a>) and the DDBJ (<a href = "https://www.ncbi.nlm.nih.gov/sra/?term=DRP000017" target = "blank">example</a>).

![sources](https://user-images.githubusercontent.com/15315514/44533218-08b95f80-a6c3-11e8-9eb8-827086fd2d99.png)

## Metadata

We provide metadata that we obtain from the source repositories.
We also attempt to, where possible, perform some harmonization of metadata to improve the discoverability of useful data.
Note that we do not yet obtain sample metadata from the <a href = "https://www.ncbi.nlm.nih.gov/biosample/" target = "blank">BioSample</a> database, so the metadata available for RNA-seq samples is limited.

### refine.bio-harmonized Metadata

Scientists who upload results don't always use the same names for related values.
This makes it challenging to search across datasets.
We have put some processes in place to smooth out some of these issues.

![harmonized-metadata](https://user-images.githubusercontent.com/15315514/44549202-5eefc800-a6ee-11e8-8a7b-57826f0153f2.png)

To produce lightly harmonized metadata, we combine certain fields based on similar keys.
We do this for convenience and to aid in searches.
For example, `treatment`, `treatment group`, `treatment protocol`, `drug treatment`, and `clinical treatment` fields get collapsed down to `treatment`.
The fields that we currently collapse to includes `specimen part`, `genetic information`, `disease`, `disease stage`, `treatment`, `race`, `subject`, `development stage`, `compound`, `cell_line`, and `time`.

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
| `cell_line` | `cell line`, `sample strain` |

We type-cast age values to doubles.
Sex is a special case; we map to `female` and `male` values if the values are one of the following:

| Harmonized `sex` value | Values |
|:-----------------------:|-------|
| `female` | `f`, `female`, `woman`|
| `male` | `m`, `male`, `man` |

Only harmonized values are displayed in the sample table on the web interface.
When downloading refine.bio data, these harmonized metadata are denoted with the `refinebio_` prefix.

We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata. If you find that the harmonized metadata does not accurately reflect the metadata supplied by the submitter, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

### Submitter Supplied Metadata

We also capture the metadata as submitted to the source repositories.
This includes experiment titles and descriptions or abstracts.
Protocol information, which generally contains the type of sample preparation, is handled differently by different databases.
We push this information down to the sample level and provide it that way.
In the case of data imported from ArrayExpress, the protocol _may_ be the information from the full experiment and not just the sample in question.
Sample metadata in their original, unharmonized form are available as part of [refine.bio downloads](#downloadable-files).

# Processing Information

## refine.bio processed  ![refinebio-processedibadge](https://user-images.githubusercontent.com/15315514/44549308-b2621600-a6ee-11e8-897b-5cdcb5d1ed9d.png)

Because refine.bio is designed to be consistently updated, we use processing and normalization methods that operate on a single sample wherever possible.
Processing and normalization methods that require multiple samples (e.g., Robust Multi-array Average or RMA) generally have outputs that are influenced by whatever samples are included and rerunning these methods whenever a new sample is added to the system would be impractical.

### Microarray pipelines

![microarray-pipeline](https://user-images.githubusercontent.com/15315514/44549355-d3c30200-a6ee-11e8-9e1a-d5140df1f7d6.png)

#### Affymetrix

SCAN (Single Channel Array Normalization) is a normalization method for develop for single channel Affymetrix microarrays that allows us to process individual samples.
SCAN models and corrects for the effect of technical bias, such as GC content, using a mixture-modeling approach.
For more information about this approach, see the primary publication (<a href = "http://dx.doi.org/10.1016/j.ygeno.2012.08.003" target = "blank">Piccolo, et al. _Genomics._ 2012.</a>) and the <a href = "https://www.bioconductor.org/packages/release/bioc/html/SCAN.UPC.html" target = "blank">SCAN.UPC Bioconductor package</a> documentation.
We specifically use the `SCANfast` implementation of SCAN and the Brainarray packages as probe-summary packages when available.

##### Platform detection

We have encountered instances where the platform label from the source repository and the metadata included in the sample's raw data file (`.CEL` file) itself do not match.
In these cases, we take the platform information included in the raw data (`.CEL`) file header to be the true platform label.

#### Illumina BeadArrays

Dr. Stephen Piccolo, the developer of SCAN, has adapted the algorithm for use with Illumina BeadArrays for refine.bio.
Because this Illumina SCAN methodology is not yet incorporated into the SCAN.UPC package, we briefly summarize the methods below.

We require that non-normalized or raw expression values and detection p-values to be present in Illumina non-normalized data.
If we infer that background correction has not occurred in the non-normalized data (e.g., there are no negative expression values), the data are background corrected using the <a href = "http://web.mit.edu/%7Er/current/arch/i386_linux26/lib/R/library/limma/html/nec.html" target = "blank">`limma::nec`</a>  function (<a href = "https://doi.org/10.1093/nar/gkq871" target = "blank">Shi, Oshlack, and Smyth. _Nucleic Acids Research._ 2010.</a>).
Following background correction -- either upstream presumably in the Illumina BeadStudio software or in our processor, arrays are normalized with SCAN.
SCAN requires probe sequence information obtained from the <a href = "https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip" target = "blank">Illumina BeadArray Bioconductor annotation packages</a> (e.g., <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html" target = "blank">`illuminaHumanv1.db`</a>).
We only retain probes that have a "Good" or "Perfect" rating in these packages; this quality rating is in reference to how well a probe is likely to measure its target transcript.

##### Platform detection

We infer the Illumina BeadArray platform that a sample is likely to be run on by comparing the probe identifiers in the unprocessed file to probes for each of the Illumina expression arrays for a given organism.
We again use the Illumina Bioconductor annotation packages for this step.
For instance, the overlap between the probe identifiers in a human sample and the probe identifiers in each human platform (<a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html" target = "blank">`v1`</a>, <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv2.db.html" target = "blank">`v2`</a>, <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv3.db.html" target = "blank">`v3`</a>[](https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv3.db.html), and <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html" target = "blank">`v4`</a>) is calculated.
The platform with the highest overlap (provided it is >75%) is inferred to be the true platform.
Some analyses around this platform detection procedure can be found in <a href = "https://github.com/jaclyn-taroni/beadarray-platform-detection" target = "blank">this repository</a>.

### RNA-seq pipelines

![rna-seq-pipeline](https://user-images.githubusercontent.com/15315514/44549339-c86fd680-a6ee-11e8-8d62-419ae7f10a94.png)

We use <a href = "https://combine-lab.github.io/salmon/" target = "blank">Salmon</a> and <a href = "https://bioconductor.org/packages/release/bioc/html/tximport.html" target = "blank">tximport</a> to process all RNA-seq data in refine.bio.
We obtain sra files run on our <a href = "https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">supported short-read platforms</a> from NCBI Sequence Read Archive.
We use <a href = "https://ncbi.github.io/sra-tools/fastq-dump.html" target = "blank">`fastq-dump`</a> to named pipes, which allows us to support paired-end experiments, and pass these to Salmon.
Note that any unmated reads from paired experiments are discarded.

We use the library strategy and library source metadata fields to identify RNA-seq experiments.
It's possible that experiments that are inappropriate for use with Salmon will still appear in refine.bio (e.g., long-read platforms that are labeled incorrectly in the source repository).
If you find an experiment that you believe is inappropriate for use with Salmon, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

#### Salmon

Salmon is an alignment-free method for estimating transcript abundances from RNA-seq data (<a href = "http://dx.doi.org/10.1038/nmeth.4197" target = "blank">Patro, et al. _Nature Methods_. 2017.</a>).
We use it in <a href = "http://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-quasi-index-and-fmd-index-based-modes" target = "blank">quasi-mapping mode</a>, which is significantly faster than alignment-based approaches and requires us to build a Salmon transcriptome index.

##### Transcriptome index

We build a custom reference transcriptome (using <a href = "https://github.com/deweylab/RSEM" target = "blank">RSEM</a> `rsem-prepare-reference`) by filtering the Ensembl genomic DNA assembly to remove _pseudogenes_, which we expect could negatively impact the quantification of protein-coding genes.
This means we're obtaining abundance estimates for coding as well as non-coding transcripts.

Building a transcriptome index with `salmon index` requires us to specify a value for the parameter `-k` that determines the size of the k-mers used for the index.
The length of a read determines what k-mer size is appropriate.
Consistent with the recommendations of the authors of Salmon, we use an index build with _k_ = 31 when quantifying samples with reads with length > 75bp.
We use _k_ = 23 for shorter read lengths.

The refine.bio processed Salmon indices are available for download.
You can make use of our API like so:

```
https://api.refine.bio/transcriptome_indices/?organism=<ORGANISM>&length=<LENGTH>
```
Where `<ORGANISM>` is the scientific name of the species in all caps separated by underscores and `<LENGTH>` is either `SHORT` or `LONG`.

To obtain the zebrafish (_Danio rerio_) index used for >75bp reads, use:

```
https://api.refine.bio/transcriptome_indices/?organism=DANIO_RERIO&length=LONG
```

The `s3_url` field will allow you to download the index.

##### Quantification with Salmon

When quantifying transcripts with `salmon quant`, we take advantage of options that allow Salmon to learn and attempt to correct for certain biases in sequencing data.
We include the flags `--seqBias` to correct for random hexamer priming and, if this is a **paired-end** experiment, `--gcBias` to correct for GC content when running salmon quant.
We set the library type parameter such that Salmon will <a href ="http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype" target = "blank">infer the sequencing library type automatically </a> for the reads it is quantifying (`-l A`).

#### tximport

Salmon quantification is at the _transcript-level_.
To better integrate with the microarray data contained in refine.bio, we summarize the transcript-level information to the _gene-level_ with `tximport` (<a href = "http://dx.doi.org/10.12688/f1000research.7563.1" target = "blank">Soneson, Love, and Robinson. _F1000 Research._ 2015.</a>).

Our tximport implementation generates <a href ="https://www.rdocumentation.org/packages/tximport/versions/1.0.3/topics/tximport" target = "blank"> "lengthScaledTPM"</a>, which are gene-level count-scale values that are generated by scaling TPM using the average transcript length across samples and to the library size.
Note that tximport is applied at the _experiment-level_ rather than to single samples.
For additional information, see the <a href= "http://bioconductor.org/packages/release/bioc/html/tximport.html" target = "blank"> tximport Bioconductor page </a>, the <a href="http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html" target ="blank">tximport tutorial _Importing transcript abundance datasets with tximport_</a>, and <a href ="http://dx.doi.org/10.12688/f1000research.7563.1" target = "blank"> Soneson, Love, and Robinson. _F1000Research._ 2015.</a>.

## Submitter processed  ![submitter-processed-badge](https://user-images.githubusercontent.com/15315514/44549307-b2621600-a6ee-11e8-9ef4-17b81d7728fd.png)

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
If you find an issue with these mappings, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub </a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

### Illumina

We support conversion from Illumina BeadArray probe IDs to Ensembl gene IDs using
<a href = "https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip" target ="blank">Bioconductor Illumina BeadArray expression packages</a>,
allowing for one-to-many mappings.

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

####  Limitations of gene identifiers when combining across platforms

We use Ensembl gene identifiers across refine.bio and, within platform, we use the same annotation to maintain consistency (e.g., all samples from the same Affymetrix platform use the same Brainarray package or identifier mapping derived from said package).
However, Brainarray packages or Bioconductor annotation packages may be assembled from different genome builds compared each other or compared to the genome build used to generate transcriptome indices.
If there tend to be considerable differences between (relatively) recent genome builds for your organism of interest or you are performing downstream analysis that would be sensitive to these differences, we do not recommend aggregating by species.

## Transformations

### Quantile normalization

refine.bio is designed to allow for the aggregation of multiple platforms and even multiple technologies.
With that in mind, we would like the distributions of samples from different platforms/technologies to be as similar as possible.
We use <a href ="https://en.wikipedia.org/wiki/Quantile_normalization" target = "blank"> quantile normalization</a> to accomplish this.
Specifically, we generate a reference distribution for each organism from a large body of data with the `normalize.quantiles.determine.target` function from the <a href = "http://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html" target ="blank">preprocessCore</a> R package and quantile normalize samples that a user selects for download with this target (using the `normalize.quantiles.use.target` function of `preprocessCore`).
We go into more detail below.

#### Reference distribution

By performing quantile normalization, we assume that the differences in expression values between samples arise solely from technical differences.
This is not always the case; for instance, samples included in refine.bio are from multiple tissues.
We'll use as many samples as possible to generate the reference or target distribution.
By including as diverse biological conditions as we have available to us to inform the reference distribution, we attempt to generate a tissue-agnostic consensus.
To that end, we use the Affymetrix microarray platform with the largest number of samples for a given organism (e.g., `hgu133plus2` in humans) and only samples we have processed from raw as shown below.

![docs-ref-dist](https://user-images.githubusercontent.com/15315514/45969124-cb692a00-c000-11e8-9cfc-6317c92202c8.png)

##### Quantile normalizing your own data with refine.bio reference distribution

refine.bio quantile normalization reference distribution or "targets" are available for download.
You may wish to use these to normalize your own data to make it more comparable to data you obtain from refine.bio.

Quantile normalization targets can be obtained by first querying the API like so:

```
https://api.refine.bio/qn_targets/?organism=<ORGANISM>
```

Where `<ORGANISM>` is the scientific name of the species in all caps separated by underscores.

To obtain the zebrafish (_Danio rerio_) reference distribution, use:

```
https://api.refine.bio/qn_targets/?organism=DANIO_RERIO
```

The `s3_url` field will allow you to download the index.
We provide a full example for obtaining and using a refine.bio reference distribution in R <a href="https://github.com/AlexsLemonade/refinebio-examples/tree/master/normalize-own-data" target = "blank">here</a>.

#### Quantile normalizing samples for delivery

Once we have a reference/target distribution for a given organism, we use it to quantile normalize any samples that a user has selected for download.
This quantile normalization step takes place _after_ the summarization and inner join steps described above and illustrated below.

![docs-normalization](https://user-images.githubusercontent.com/15315514/46034575-f9b53b00-c0ce-11e8-893d-868a2aa98520.png)

As a result of the quantile normalization shown above, Sample 1 now has the same underlying distribution as the reference for that organism.

Note that only gene expression matrices that we are able to successfully quantile normalize will be available for download.

#### Limitations of quantile normalization across platforms with many zeroes

Quantile normalization is a strategy that can address many technical effects, generally at the cost of retaining certain sources of biological variability.
However, in cases where there are many ties that are different between samples, the transformation can produce outputs with different distributions.
This arises in data processed by refine.bio when RNA-seq and microarray data are combined.
To confirm that we have quantile normalized data correctly before returning results to the user, we evaluate the top half of expression values and confirm that a KS test produces a non-significant p-value.
Users who seek to analyze RNA-seq and microarray data together should be aware that the low-expressing genes may not be comparable across the sets.

### Gene transformations

In some cases, it may be useful to row-normalize or transform the gene expression values in a matrix (e.g., following aggregation and quantile normalization).
We offer the following options for transformations:

* **None:** No row-wise transformation is performed.

* **Z-score:** Row values are <a href = "https://en.wikipedia.org/wiki/Standard_score" target = "blank">z-scored</a> using the <a href = "http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html" target = "blank">`StandardScaler`</a> from <a href = "http://scikit-learn.org/stable/index.html" target = "blank">`scikit-learn`</a>.
This transformation is useful for examining samples' gene expression values relative to the rest of the samples in the  expression matrix (either all selected samples from that _species_ when aggregating by species or all selected samples in an _experiment_ when aggregating by experiment).
If a sample has a positive value for a gene, that gene is more highly expressed in that sample compared to the mean of all samples; if that value is negative, that gene is less expressed compared to the population.
It assumes that the data are normally distributed.

* **Zero to one:** Rows are scaled to values `[0,1]` using the <a href = "http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html" target = "blank">`MinMaxScaler`</a> from <a href = "http://scikit-learn.org/stable/index.html" target = "blank">`scikit-learn`</a>.
We expect this transformation to be most useful for certain machine learning applications (e.g., those using cross-entropy as a loss function).

In the plot below, we demonstrate the effect of different scaling options on gene expression values (using a randomly selected human dataset, microarray platform, and gene):

<img src="https://user-images.githubusercontent.com/19534205/46036353-50247880-c0d3-11e8-9b7f-07818e545d68.png" width="480">

Note that the distributions retain the same general _shape_, but the range of values and the density are altered by the transformations.

# Downloadable Files

Users can download gene expression data and associated sample and experiment metadata from refine.bio.
These files are delivered as a zip file.
The folder structure within the zip file is determined by whether a user selected to aggregate by **experiment** or by **species**.

### The download folder structure for data aggregated by experiment:

![docs-downloads-experiment-agg](https://user-images.githubusercontent.com/15315514/45906716-2f9eaa80-bdc3-11e8-9855-2aaeb74e588d.png)

In this example, two experiments were selected.
There will be as many folders as there are selected experiments.

### The download folder structure for data aggregated by species:

![docs-downloads-species-agg](https://user-images.githubusercontent.com/15315514/45906715-2f9eaa80-bdc3-11e8-8ab3-90ccc40cfa11.png)

In this example, samples from two species were selected.
There will be as many folders as there are selected experiments and this will be the case regardless of how many individual experiments were included.

In both cases, `aggregated_metadata.json` contains metadata, including both _experiment_ metadata (e.g., experiment description and title) and _sample_ metadata for everything included in the download.
Below we describe the files included in the delivered zip file.

## Gene Expression Matrix

Gene expression matrices are delived in <a href = "https://en.wikipedia.org/wiki/Tab-separated_values" target = "blank">tab-separated value</a> (TSV) format.
In these matrices, rows are _genes_ or _features_ and columns are _samples_.
Note that this format is consistent with the input expected by many programs specifically designed for working with gene expression matrices, but some machine learning libraries will expect this to be transposed.
The column names or header will contain values corresponding to sample accessions (denoted `refinebio_accession_code` in metadata files).
You can use these values in the header to map between a sample's gene expression data and its metadata (e.g., disease label or age). See also [Use Cases for Downstream Analysis](#use-cases-for-downstream-analysis).

## Sample Metadata

Sample metadata is delivered in the `metadata_<experiment-accession-id>.tsv`, `metadata_<species>.json`, `metadata_<species>.tsv`, and `aggregated_metadata.json` files.

The primary way we identify samples is by using the sample accession, denoted by `refinebio_accession_code`.
Harmonized metadata fields (see the [section on harmonized metadata](#refine.bio-harmonized-metadata)) are noted with a `refinebio_` prefix.
The `refinebio_source_archive_url` and `refinebio_source_database` fields indicate where the sample was obtained from.
If there are no keys from the source data associated with a harmonized key, the harmonized metadata field will be empty.
We also deliver submitter-supplied data; see below for more details.
**We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata.**
If you find that refine.bio metadata does not accurately reflect the metadata supplied by the submitter, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [ccdl@alexslemonade.org](mailto:ccdl@alexslemonade.org).

### TSV files

In metadata TSV files, samples are represented as rows.
The first column contains the `refinebio_accession_code` field, which match the header/column names in the gene expression matrix, followed by refine.bio-harmonized fields (e.g., `refinebio_`), and finally submitter-supplied values.
Some information from source repositories comes in the form of nested values, which we attempt to "flatten" where possible.
Note that some information from source repositories is redundant--ArrayExpress samples often have the same information in `characteristic` and `variable` fields--and we _assume_ that if a field appears in both, the values are identical.
For samples run on Illumina BeadArray platforms, information about what platform we detected and the metrics used to make that determination will also be included.

Columns in these files will often have missing values, as not all fields will be available for every sample included in a file.
This will be particularly evident when aggregating by experiments that have different submitter-supplied information associated with them (e.g., one experiment contains a `imatinib` key and all others do not).

### JSON files

Submitter supplied metadata and source urls are delivered in `refinebio_annotations`.
As described above, harmonized values are noted with a `refinebio_` prefix.

## Experiment Metadata

Experiment metadata (e.g., experiment description and title) is delivered in the `metadata_<species>.json` and `aggregated_metadata.json` files.


# Species compendia

We periodically release compendia comprised of all the samples from a species that we were able to process.
We refer to these as **species compendia**.
We process these compendia in a manner that is different from the options that are available via the web user interface.
These species compendia provide a snapshot of the most complete collection of gene expression that refine.bio can produce for each supported organism.

The refine.bio web interface does an inner join when datasets are combined, so only genes present in all datasets are included in the final matrix.
For compendia, we take the union of all genes, filling in any missing values with `NA`.
This is a "full outer join" as illustrated below.
We use a full outer join because it allows us to retain more genes in a compendium and we impute missing values during compendia creation.

![outer join](https://user-images.githubusercontent.com/15315514/44534241-4dde9100-a6c5-11e8-8a9c-aa147e294e81.png)

We perform an outer join each time samples are combined in the process of building species compendia.

![docs-species-compendia](https://user-images.githubusercontent.com/15315514/48498088-72995f00-e803-11e8-8832-3a9024748431.png)

Samples from each technology—microarray and RNA-seq—are combined separately.
In RNA-seq samples, we filter out genes with low total counts and then `log2(x + 1)` the data.
We join samples from both technologies.
We then drop genes that have missing values in greater than 30% of samples.
Finally, we drop samples that have missing values in greater than 50% of genes.
We impute the remaining missing values with IterativeSVD from <a href = "https://pypi.org/project/fancyimpute/" target = "blank">fancyimpute</a>.
We then quantile normalize all samples as described above.

We've made our analyses underlying processing choices and exploring test compendia available at our <a href = "https://github.com/AlexsLemonade/compendium-processing" target = "blank">`compendium-processing`</a> repository.


# Use Cases for Downstream Analysis

Our <a href = "https://github.com/AlexsLemonade/refinebio-examples" target = "blank">`refinebio-examples`</a> repo includes a number of different analyses you can perform with data from refine.bio. We include examples in the R programming language, and where applicable, <a href = "http://genepattern-notebook.org/example-notebooks/" target = "blank">GenePattern Notebooks</a> and scripts to prepare refine.bio data for use with GenePattern. The following examples are included:

* Differential expression analysis [<a href = "https://github.com/AlexsLemonade/refinebio-examples/tree/master/differential-expression" target = "blank">README</a>, <a href = "https://alexslemonade.github.io/refinebio-examples/differential-expression/gene_DE.nb.html" target = "blank">notebook</a>]
* Converting between different gene identifiers [<a href = "https://github.com/AlexsLemonade/refinebio-examples/tree/master/ensembl-id-convert" target = "blank">README</a>, <a href = "https://alexslemonade.github.io/refinebio-examples/ensembl-id-convert/ensembl_id_convert.nb.html" target = "blank">notebook</a>]
* Ortholog mapping [<a href = "https://github.com/AlexsLemonade/refinebio-examples/tree/master/ortholog-mapping" target = "blank">README</a>, <a href = "https://alexslemonade.github.io/refinebio-examples/ortholog-mapping/ortholog_mapping_example.nb.html" target = "blank">notebook</a>]
* Clustering/heatmap generation [<a href = "https://github.com/AlexsLemonade/refinebio-examples/tree/master/clustering" target = "blank">README</a>, <a href = "https://alexslemonade.github.io/refinebio-examples/clustering/clustering_example.nb.html" target = "blank">notebook</a>] 
