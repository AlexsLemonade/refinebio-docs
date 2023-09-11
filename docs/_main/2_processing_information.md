
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
When available, we use <a href = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp" target = "blank">BrainArray Custom CDFs</a> during processing with SCAN.

##### Affymetrix platform detection

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

##### Illumina platform detection

We infer the Illumina BeadArray platform that a sample is likely to be run on by comparing the probe identifiers in the unprocessed file to probes for each of the Illumina expression arrays for a given organism.
We again use the Illumina Bioconductor annotation packages for this step.
For instance, the overlap between the probe identifiers in a human sample and the probe identifiers in each human platform (<a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html" target = "blank">`v1`</a>, <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv2.db.html" target = "blank">`v2`</a>, <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv3.db.html" target = "blank">`v3`</a>, and <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html" target = "blank">`v4`</a>) is calculated.
The platform with the highest overlap (provided it is >75%) is inferred to be the true platform.
Some analyses around this platform detection procedure can be found in <a href = "https://github.com/jaclyn-taroni/beadarray-platform-detection" target = "blank">this repository</a>.

##### Handling Illumina probes that map to multiple Ensembl gene identifiers

Illumina probes sometimes map to multiple Ensembl gene identifiers when using the annotation in <a href = "https://www.bioconductor.org/packages/release/BiocViews.html#___IlluminaChip" target = "blank">Illumina BeadArray Bioconductor annotation packages</a> (e.g., <a href = "https://www.bioconductor.org/packages/release/data/annotation/html/illuminaHumanv1.db.html" target = "blank">`illuminaHumanv1.db`</a>).
For human platforms in particular, these genes tend to be from highly polymorphic loci, e.g., Killer-cell immunoglobulin-like receptors.
Because refine.bio allows users to combine samples from multiple platforms, we prioritize Ensembl gene identifiers that maximize compatibility with other platforms.
Specifically, we select genes in order of priority as follows:

- Pick the gene ID with the most appearances in BrainArray packages for Affymetrix platforms of the same species as the input Illumina platform
- If two or more of the associated gene IDs appear an equal number of times in BrainArray packages, or if none of the associated gene IDs appear in any BrainArray package, we break ties as follows:
  - First, we check Ensembl and filter out any gene IDs that are no longer valid
  - Next, if there are any Ensembl genes on the primary assembly, we take only the genes on the primary assembly and discard genes in <a href"https://www.ensembl.org/info/genome/genebuild/haplotypes_patches.html" target = "blank">haplotypes (alternative versions of the genome) or patches</a>.
  - If there are still two or more genes left, we pick the gene with the lowest Ensembl gene identifier to break the tie. This is an arbitrary but consistent way to break ties.

This is implemented in <a href="https://github.com/AlexsLemonade/illumina-refinery" target = "blank">`AlexsLemonade/illumina-refinery`</a>.

### RNA-seq pipelines

![rna-seq-pipeline](https://user-images.githubusercontent.com/15315514/44549339-c86fd680-a6ee-11e8-8d62-419ae7f10a94.png)

We use <a href = "https://combine-lab.github.io/salmon/" target = "blank">Salmon</a> and <a href = "https://bioconductor.org/packages/release/bioc/html/tximport.html" target = "blank">tximport</a> to process all RNA-seq data in refine.bio.
We obtain sra files run on our <a href = "https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">supported short-read platforms</a> from NCBI Sequence Read Archive.
We use <a href = "https://ncbi.github.io/sra-tools/fastq-dump.html" target = "blank">`fastq-dump`</a> to named pipes, which allows us to support paired-end experiments, and pass these to Salmon.
Note that any unmated reads from paired experiments are discarded.

We use the **library strategy** and **library source** metadata fields to identify RNA-seq experiments.
It's possible that experiments that are inappropriate for use with Salmon will still appear in refine.bio (e.g., long-read platforms that are labeled incorrectly in the source repository).
We also encounter trouble distinguishing single-cell and bulk RNA-seq data from these fields.
We strongly recommend exercising caution when using single-cell data from refine.bio as the pipeline we use may be inappropriate (e.g., correcting for gene length in 3' tagged RNA-seq data induces bias [<a href = "https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#tagged-rna-seq" target = "blank">ref</a>], Salmon TPM may overcorrect expression of long genes [<a href = "http://hemberg-lab.github.io/scRNA.seq.course/construction-of-expression-matrix.html#reads-alignment" target = "blank">ref</a>]).
If you find an experiment that you believe is inappropriate for use with our pipeline, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

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
https://api.refine.bio/v1/transcriptome_indices/?organism__name=<ORGANISM>&length=<LENGTH>
```
Where `<ORGANISM>` is the scientific name of the species in all caps separated by underscores and `<LENGTH>` is either `SHORT` or `LONG`.

To obtain the zebrafish (_Danio rerio_) index used for >75bp reads, use:

```
https://api.refine.bio/v1/transcriptome_indices/?organism__name=DANIO_RERIO&length=LONG
```

The `download_url` field will allow you to download the index.

###### Microbial transcriptome indices

For an individual microbial species, there are often multiple genome assemblies available.
Multiple assemblies from a species often reflect that multiple strains from that species have been characterized or established as laboratory strains.
We use a single genome assembly, selected by domain experts, when we generate the transcriptome index used for each microbial species.
The single genome assembly is typically from a common laboratory strain.
This strain's transcriptome index is used to quantify _all_ RNA-seq data from that species, regardless of a sample's reported strain of origin.
Quantifying using a single strain background allows us to generate compendia using a single transcriptome index, but most likely results in some loss of information for samples that do not originate from that strain.

For a list of supported microbial species and the assemblies used, please see the [`config/organism_strain_mapping.csv`](https://github.com/AlexsLemonade/refinebio/blob/master/config/organism_strain_mapping.csv) file.
For more context, please also visit the original GitHub issue at [`AlexsLemonade/refinebio#1722`](https://github.com/AlexsLemonade/refinebio/issues/1722).

##### Quantification with Salmon

When quantifying transcripts with `salmon quant`, we take advantage of options that allow Salmon to learn and attempt to correct for certain biases in sequencing data.
We include the flags `--seqBias` to correct for random hexamer priming and, if this is a **paired-end** experiment, `--gcBias` to correct for GC content when running salmon quant.
We set the library type parameter such that Salmon will <a href ="http://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype" target = "blank">infer the sequencing library type automatically </a> for the reads it is quantifying (`-l A`).

#### tximport

Salmon quantification is at the _transcript-level_.
To better integrate with the microarray data contained in refine.bio, we summarize the transcript-level information to the _gene-level_ with `tximport` (<a href = "http://dx.doi.org/10.12688/f1000research.7563.1" target = "blank">Soneson, Love, and Robinson. _F1000 Research._ 2015.</a>).

Our tximport implementation generates <a href ="https://www.rdocumentation.org/packages/tximport/versions/1.0.3/topics/tximport" target = "blank"> "lengthScaledTPM"</a>, which are gene-level count-scale values that are generated by scaling TPM using the average transcript length across samples and the library size.
Note that tximport is applied at the _experiment-level_ rather than to single samples.
For additional information, see the <a href= "http://bioconductor.org/packages/release/bioc/html/tximport.html" target = "blank"> tximport Bioconductor page </a>, the <a href="http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html" target ="blank">tximport tutorial _Importing transcript abundance datasets with tximport_</a>, and <a href ="http://dx.doi.org/10.12688/f1000research.7563.1" target = "blank"> Soneson, Love, and Robinson. _F1000Research._ 2015.</a>.

In some cases, all samples in an experiment can not be processing using Salmon (e.g., the file available for one sample is malformed).
When experiments are >80% complete and contain more than 20 samples, we may run tximport on all available samples at the time.
We've found the effect of running tximport "early" on the resulting values to be small under these conditions.
As a result, the tximport may be run on the same example multiple times; the most recent values will be available from refine.bio.

## Submitter processed  ![submitter-processed-badge](https://user-images.githubusercontent.com/15315514/44549307-b2621600-a6ee-11e8-9ef4-17b81d7728fd.png)

Sometimes raw data for a sample is either unavailable at the source repository or exists in a form that we can not process.
For microarray platforms that we support, we obtain the submitter processed expression data and use these values in refine.bio with some modification (e.g., log2-transformation where we detect it has not been performed).

As noted above, we use Ensembl gene identifiers throughout refine.bio.
Submitter processed data may use other gene (or probe) identifiers that we must convert to Ensembl gene identifiers.
We describe the processes for Affymetrix and Illumina data below.
Note in the case of one-to-many mappings when going from the ID used by the submitter to the Ensembl gene ID, expression values are duplicated:
if a probe maps to two Ensembl gene ids, those two Ensembl gene ids will both have the probe's expression value following conversion.

### Affymetrix identifier conversion

We have created custom gene mapping files for most of the Affymetrix platforms we support.
Briefly, for Brainarray supported platforms, we use the Brainarray (e.g., `hgu133plus2hsensgprobe`) and the platform-specific annotation package from Bioconductor (e.g., `hgu133plus2.db`) to generate a platform-specific mapping file that includes probe IDs, Ensembl gene IDs, gene symbols, Entrez IDs, RefSeq and Unigene identifiers.
The rationale for only using probes or IDs that are accounted for in the Brainarray package is two-fold: 1) Brainarray packages are updated as we learn more about the genome and 2) it allows for these submitter processed data to be more consistent with refine.bio processed data.
We support identifier conversion for a limited number of platforms that either do not have a Brainarray or Bioconductory annotation packages.

The code for deriving these mappings and more details are available at https://github.com/AlexsLemonade/identifier-refinery.
If you find an issue with these mappings, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub </a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

### Illumina identifier conversion

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
There is a single reference distribution per species, used to normalize all samples from that species regardless of platform or technology.
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
https://api.refine.bio/v1/qn_targets/<ORGANISM>
```

Where `<ORGANISM>` is the scientific name of the species separated by underscores.

To obtain the zebrafish (_Danio rerio_) reference distribution, use:

```
https://api.refine.bio/v1/qn_targets/danio_rerio
```

The `s3_url` field will allow you to download the index.

#### Quantile normalizing samples for delivery

Once we have a reference/target distribution for a given organism, we use it to quantile normalize any samples that a user has selected for download.
This quantile normalization step takes place _after_ the summarization and inner join steps described above and illustrated below.

![docs-normalization](https://user-images.githubusercontent.com/15315514/46034575-f9b53b00-c0ce-11e8-893d-868a2aa98520.png)

As a result of the quantile normalization shown above, Sample 1 now has the same underlying distribution as the reference for that organism.

Note that only gene expression matrices that we are able to successfully quantile normalize will be available for download.

#### Limitations of quantile normalization across platforms with many zeroes

Quantile normalization is a strategy that can address many technical effects, generally at the cost of retaining certain sources of biological variability.
We use a single reference distribution per organism, generated from the Affymetrix microarray platform with the largest number of samples we were able to process from raw data (see [_Reference distribution_](#reference=distribution)).
In cases where the unnormalized data contains many ties within samples (and the ties are different between samples) the transformation can produce outputs with somewhat different distributions.
This situation arises most often when RNA-seq data and microarray data are combined into a single dataset or matrix.
To confirm that we have quantile normalized data correctly before returning results to the user, we evaluate the top half of expression values and confirm that a KS test produces a non-significant p-value.
Users who seek to analyze RNA-seq and microarray data together should be aware that the low-expressing genes may not be comparable across the sets.

#### Skipping quantile normalization for RNA-seq experiments

When selecting RNA-seq samples for download and to aggregate by experiment, users have the option to skip quantile normalization by first selecting Advanced Options and checking the "Skip quantile normalization for RNA-seq samples" box.
In this case, the output of tximport will be delivered in TSV files (see [our section on RNA-seq data processing with tximport](#tximport)).
These data can be used for differential expression analysis as "bias corrected counts without an offset" as described in the <a href = "https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#use-with-downstream-bioconductor-dge-packages" target = "blank">_Use with downstream Bioconductor DGE packages_ section of tximport vignette</a>.
Note that these data will be less comparable to other datasets from refine.bio because this step has been skipped.

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
