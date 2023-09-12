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
You can use these values in the header to map between a sample's gene expression data and its metadata (e.g., disease label or age). See also [Downstream Analysis with refine.bio Examples](#downstream-analysis-with-refinebio-examples).

## Sample Metadata

Sample metadata is delivered in the `metadata_<experiment-accession-id>.tsv`, `metadata_<species>.json`, `metadata_<species>.tsv`, and `aggregated_metadata.json` files.

The primary way we identify samples is by using the sample accession, denoted by `refinebio_accession_code`.
Harmonized metadata fields (see the [section on harmonized metadata](#refine.bio-harmonized-metadata)) are noted with a `refinebio_` prefix.
The `refinebio_source_archive_url` and `refinebio_source_database` fields indicate where the sample was obtained from.
If there are no keys from the source data associated with a harmonized key, the harmonized metadata field will be empty.
We also deliver submitter-supplied data; see below for more details.
**We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata.**
If you find that refine.bio metadata does not accurately reflect the metadata supplied by the submitter, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

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


The `aggregated_metadata.json` file contains additional information regarding the processing of your dataset.
Specifically, the `aggregate_by` and `scale_by` fields note how the samples are grouped into gene expression matrices and how the gene expression data values were transformed, respectively.
The `quantile_normalized` fields notes whether or not quantile normalization was performed.
Currently, we only support skipping quantile normalization for RNA-seq experiments when aggregating by experiment on the web interface.
