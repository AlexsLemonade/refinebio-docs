# Frequently Asked Questions

#### What is the difference between refine.bio-processed and submitter-processed datasets?

Samples and the datasets they comprise are designated as refine.bio-processed if we were able to obtain raw data in a suitable format for one of our processing pipelines.
If no suitable raw data is available for a sample on a supported platform, we obtain the _processed_ data available from the source repository and modify it to be more consistent with refine.bio-processed data; these samples are termed submitter-processed.
See the **Source Data** section for more information.

#### How do you process the data?

We process samples run on <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_microarray_platforms.csv" target = "blank">supported microarray platforms</a> with Single-Channel Array Normalization (SCAN) and transcriptomic samples run on <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">supported sequencing platforms</a> with Salmon and tximport.
Please see the **Processing Information** section for more details.

#### What type of data does refine.bio support?

refine.bio currently supports gene expression data, specifically genome-scale microarray and RNA-seq data.
See our supported <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_microarray_platforms.csv" target = "blank">microarray</a> and <a href ="https://github.com/AlexsLemonade/refinebio/blob/dev/config/supported_rnaseq_platforms.txt" target = "blank">RNA-seq</a> platforms.

#### What does "corrected" metadata mean?

Scientists will often use different terminology to refer to a similar sample metadata field or _key_.
For example, `treatment` and `treatment protocol` may make reference to the same kinds of information.
We attempt to perform some mapping between keys to aid in searches.
See the **refine.bio-harmonized Metadata** section for more information, including a full list of the mappings we perform.

#### Why do the values differ a little bit if I download different datasets?

We've prioritized keeping expression values consistent _within_ a dataset based on the samples it contains.
Specifically, we remove any genes that are not measured in every sample in a dataset and we do that prior to performing quantile normalization.
When using quantile normalization, the expression value a gene is assigned in a particular sample depends on the _rank_ of that gene.
If a user download different datasets, which may have different numbers of genes, it's possible then that the same gene in the same sample would have a different expression value between them.

#### Why do I get a limited number of genes back when I aggregate samples from different platforms?

Different platforms will often measure different sets of genes.
These differences can be particularly pronounced when comparing older microarray platforms to more recent platforms.
When aggregating samples, we retain _only_ the genes present in _every sample_.
If the dataset delivered to you has fewer genes than you were expecting for that organism, it could be the result of combining multiple platforms or it may be from an older microarray platform.

#### Why can't I add certain samples to my dataset?

refine.bio will sometimes obtain the metadata (e.g., sample title or experimental protocol) associated with a sample but the raw or submitter processed expression data files are in a format that we can not process.
We do not allow you to add these samples to your dataset because we can not deliver expression values.

#### How can I find out what versions of software/packages were used to process the data?
Version information for the packages we think are _most important_ for data processing is available on the pop-up displayed when you click a sample's processing information link.

<img src="https://user-images.githubusercontent.com/15315514/46314376-5c5b7a80-c598-11e8-9bf3-c63a6a1a9696.gif" width="560">

The same package information is in the processor list available via our API:

```
https://api.refine.bio/processors/
```

In addition, you may wish to obtain <a href ="https://hub.docker.com/u/ccdl/" target = "blank">our Docker images</a> (prefixed with `dr_`) which will allow you to access version information for every dependency.

#### Are refine.bio datasets I download batch corrected?

We apply quantile normalization to mitigate issues caused by differences in the underlying distributions of gene expression values in samples.
This makes the gene expression values broadly comparable, but doesn't explicitly correct for batch, dataset, or platform.
If the scientific question and analysis methods require datasets to be batch corrected, users should first investigate the existence of batch effects using methods such as Principal Components Analysis.
If the source dataset is associated with major sources of variability in the data, users may wish to use a meta-analysis framework considering each dataset independently or to apply a batch correction tool.
For certain analyses it may be sufficient to include batch, dataset, or platform as covariates.
We use quantile normalization to make samples more comparable to one another, but this is unlikely to account for all batch effects, dataset-specific, or platform-specific biases in all cases.
