# Source Data

## Types of Data

The current version of <a href ="https://www.refine.bio" target = "blank">refine.bio </a>is designed to process gene expression data.
This includes both microarray data and RNA-seq data.
We normalize data to Ensembl gene identifiers and provide abundance estimates.

More precisely, we support microarray platforms based on their GEO or ArrayExpress accessions.
We currently support Affymetrix microarrays and Illumina BeadArrays, and we are continuing to evaluate and add support for more platforms.
This <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_microarray_platforms.csv" target = "blank">table </a> contains the microarray platforms that we support.
We process a subset of Affymetrix platforms using the <a href = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp" target = "blank">BrainArray Custom CDFs</a>, which are denoted by a `y` in the `is_brainarray` column.
We also support RNA-seq experiments performed on <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">these</a>  short-read platforms.
For more information on how data are processed, see the [Processing Information](#processing-information) section of this document.
If there is a platform that you would like to see processed, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a>.
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

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
If the values can not be type-cast to doubles (e.g., "9yrs 2mos"), these are not added to the harmonized field.
We do not attempt to normalize differences in units (e.g., months, years, days) for the harmonized age key.
Users should consult the submitter-supplied information to determine what unit is used.

Sex is a special case; we map to `female` and `male` values if the values are one of the following:

| Harmonized `sex` value | Values |
|:-----------------------:|-------|
| `female` | `f`, `female`, `woman`|
| `male` | `m`, `male`, `man` |

Only harmonized values are displayed in the sample table on the web interface.
When downloading refine.bio data, these harmonized metadata are denoted with the `refinebio_` prefix.

We recommend that users confirm metadata fields that are particularly important via the submitter-supplied metadata. If you find that the harmonized metadata does not accurately reflect the metadata supplied by the submitter, please <a href ="https://github.com/AlexsLemonade/refinebio/issues" target = "blank">file an issue on GitHub</a> so that we can resolve it.
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

### Submitter Supplied Metadata

We also capture the metadata as submitted to the source repositories.
This includes experiment titles and descriptions or abstracts.
Protocol information, which generally contains the type of sample preparation, is handled differently by different databases.
We push this information down to the sample level and provide it that way.
In the case of data imported from ArrayExpress, the protocol _may_ be the information from the full experiment and not just the sample in question.
Sample metadata in their original, unharmonized form are available as part of [refine.bio downloads](#downloadable-files).
