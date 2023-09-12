# refine.bio Compendia

We periodically release compendia comprised of all the samples from a species that we were able to process.
We refer to these as **refine.bio compendia**.
We offer two kinds of refine.bio compendia: [normalized compendia](#normalized-compendia) and [RNA-seq sample compendia](#rna-seq-sample-compendia).

## Normalized compendia

refine.bio normalized compendia are comprised of all the samples from a species that we were able to process, aggregate, and normalize.
Normalized compendia provide a snapshot of the most complete collection of gene expression that refine.bio can produce for each supported organism.
We process these compendia in a manner that is different from the options that are available via the web user interface.
Note that [submitter processed](#submitter-processed--) samples that are available through the web user interface are omitted from normalized compendia because [these samples can introduce unwanted technical variation](https://github.com/AlexsLemonade/refinebio/issues/2114).

The refine.bio web interface does an inner join when datasets are combined, so only genes present in all datasets are included in the final matrix.
For compendia, we take the union of all genes, filling in any missing values with `NA`.
This is a "full outer join" as illustrated below.
We use a full outer join because it allows us to retain more genes in a compendium and we impute missing values during compendia creation.

![outer join](https://user-images.githubusercontent.com/15315514/44534241-4dde9100-a6c5-11e8-8a9c-aa147e294e81.png)

We perform an outer join each time samples are combined in the process of building normalized compendia.

![docs-normalized-compendia](https://user-images.githubusercontent.com/15315514/65698014-ddcdd780-e049-11e9-8ed6-f1d2f8ac2ee7.png)

Samples from each technology—microarray and RNA-seq—are combined separately.
In RNA-seq samples, we filter out genes with low total counts and then `log2(x + 1)` the data.
We join samples from both technologies.
We then drop genes that have missing values in greater than 30% of samples.
Finally, we drop samples that have missing values in greater than 50% of genes.
We impute the remaining missing values with IterativeSVD from <a href = "https://pypi.org/project/fancyimpute/" target = "blank">fancyimpute</a>.
We then quantile normalize all samples as described above.

We've made our analyses underlying processing choices and exploring test compendia available at our <a href = "https://github.com/AlexsLemonade/compendium-processing" target = "blank">`compendium-processing`</a> repository.

### Collapsing by genus

Microarray platforms are generally designed to assay samples from a specific species.
In some cases, publicly available data surveyed by refine.bio may include samples where the microarray platform used was not specifically designed for the species as described (e.g., samples labeled _Bos indicus_ were run on _Bos taurus_ microarrays or mouse crosses that are not labeled _Mus musculus_ were run on _Mus musculus_ microarrays).
When we encounter this in refine.bio, we will include samples in a compendium from species that differ from the primary platform species when the two species share a genus (e.g., _Bos indicus_ samples run on _Bos taurus_ microarrays are included in the _Bos taurus_ normalized compendium, and _Mus_ crosses are included in the _Mus musculus_ normalized compendium).
Such non-primary species samples generally account for a small fraction of the total samples included in a normalized compendium.
If you would like to filter a normalized compendium based on a sample's species label, you can use the `refinebio_organism` column in the metadata TSV file or the `.samples[].refinebio_organism` field in the metadata JSON file included as part of the download.

Note that non-primary species samples from species that are outside the genus of the primary platform species are not currently available in any normalized compendium (e.g., _Pan troglodytes_ samples assayed on _Homo sapiens_ microarrays are not included in the _Pan troglodytes_ or _Homo sapiens_ compendia), but can be included in datasets from refine.bio.

Below is the list of organisms and their primary organisms:

| Primary Organism | Organisms included in compendium |
|:---------:|-------------------|
|`Anopheles gambiae`|`Anopheles gambiae`|
|`Arabidopsis thaliana`|`Arabidopsis thaliana`, `Arabidopsis halleri`, `Arabidopsis thaliana x arabidopsis halleri subsp. gemmifera`, `Arabidopsis lyrata subsp. petraea`, `Arabidopsis lyrata subsp. lyrata`, `Arabidopsis thaliana x arabidopsis lyrata`, `Arabidopsis halleri subsp. gemmifera`, `Arabidopsis lyrata`|
|`Bos indicus`|`Bos indicus`|
|`Bos taurus`|`Bos taurus`, `Bos indicus`, `Bos grunniens`|
|`Caenorhabditis elegans`|`Caenorhabditis elegans`|
|`Citrus sinensis`|`Citrus x paradisi`, `Citrus reticulata`, `Citrus sinensis`, `Citrus limon`, `Citrus reticulata x citrus trifoliata`, `Citrus clementina`, `Citrus unshiu`, `Citrus x tangelo`, `Citrus maxima`|
|`Danio rerio`|`Danio rerio`|
|`Drosophila melanogaster`|`Drosophila melanogaster`, `Drosophila simulans`, `Drosophila mauritiana`, `Drosophila sechellia`, `Drosophila teissieri`, `Drosophila santomea`, `Drosophila yakuba`|
|`Escherichia coli`|`Escherichia coli`, `Escherichia coli str. k-12 substr. mg1655`, `Escherichia coli k-12`, `Escherichia coli cft073`, `Escherichia coli str. k-12 substr. w3110`, `Escherichia coli uti89`, `Escherichia coli b str. rel606`, `Escherichia coli o157`, `Escherichia coli 8624`, `Escherichia coli sci-07`, `Escherichia coli bw25113`, `Escherichia coli apec o2`, `Escherichia coli str. k-12 substr. mc4100`, `Escherichia coli o08`, `Escherichia coli str. k-12 substr. dh10b`|
|`Escherichia coli k-12`|`Escherichia coli k-12`|
|`Escherichia coli str. k-12 substr. mg1655`|`Escherichia coli str. k-12 substr. mg1655`|
|`Gallus gallus`|`Gallus gallus`|
|`Glycine max`|`Glycine max`, `Glycine soja`|
|`Gossypium hirsutum`|`Gossypium herbaceum`, `Gossypium hirsutum`, `Gossypium barbadense`, `Gossypium arboreum`|
|`Homo sapiens`|`Homo sapiens`|
|`Hordeum vulgare`|`Hordeum vulgare`, `Hordeum vulgare subsp. spontaneum`|
|`Lepidium sativum`|`Lepidium sativum`|
|`Macaca fascicularis`|`Macaca fascicularis`|
|`Macaca mulatta`|`Macaca mulatta`|
|`Musculus`|`Musculus`|
|`Mus musculus`|`Mus musculus`, `Mus spretus`, `Mus caroli`, `Mus musculus musculus x m. m. domesticus`, `Mus musculus domesticus`, `Mus musculus x mus spretus`, `Mus musculus musculus x m. m. castaneus`, `Mus musculus musculus`, `Mus musculus castaneus`, `Mus sp.`|
|`Mustela putorius furo`|`Mustela putorius furo`|
|`Oryza sativa`|`Oryza sativa japonica`, `Oryza sativa`, `Oryza sativa indica group`, `Oryza longistaminata`|
|`Oryza sativa indica group`|`Oryza sativa indica group`|
|`Plasmodium falciparum`|`Plasmodium falciparum 3d7`, `Plasmodium falciparum`|
|`Populus tremula x populus alba`|`Populus tremula x populus alba`|
|`Populus trichocarpa`|`Populus trichocarpa`|
|`Populus x canadensis`|`Populus x canadensis`|
|`Pseudomonas aeruginosa`|`Pseudomonas aeruginosa`, `Pseudomonas aeruginosa pao1`, `Pseudomonas putida`, `Pseudomonas aeruginosa ucbpp-pa14`, `Pseudomonas aeruginosa pa14`, `Pseudomonas aeruginosa tbcf10839`, `Pseudomonas aeruginosa pahm4`|
|`Pseudomonas aeruginosa pao1`|`Pseudomonas aeruginosa pao1`|
|`Rattus norvegicus`|`Rattus norvegicus`, `Rattus rattus`, `Rattus norvegicus albus`|
|`Saccharomyces cerevisiae`|`Saccharomyces cerevisiae`, `Saccharomyces cerevisiae s288c`, `Saccharomyces pastorianus`, `Saccharomyces pastorianus weihenstephan 34/70`, `Saccharomyces cerevisiae vin13`, `Saccharomyces uvarum`, `Saccharomyces cerevisiae ec1118`, `Saccharomyces cerevisiae cen.pk113-7d`, `Saccharomyces cerevisiae by4741`, `Saccharomyces cerevisiae sk1`, `Saccharomyces bayanus`, `Saccharomyces cerevisiae x saccharomyces kudriavzevii`, `Saccharomyces boulardii`|
|`Schizosaccharomyces pombe`|`Schizosaccharomyces pombe`, `Schizosaccharomyces pombe 972h-`|
|`Staphylococcus aureus`|`Staphylococcus aureus`, `Staphylococcus aureus subsp. aureus rn4220`, `Staphylococcus aureus subsp. aureus n315`, `Staphylococcus aureus subsp. aureus usa300`, `Staphylococcus aureus subsp. aureus mu50`, `Staphylococcus aureus subsp. aureus str. newman`|
|`Sus scrofa`|`Sus scrofa domesticus`, `Sus scrofa`|
|`Triticum aestivum`|`Triticum aestivum`, `Triticum turgidum subsp. dicoccoides`, `Triticum turgidum subsp. durum`, `Triticum turgidum`, `Triticum carthlicum`, `Triticum monococcum`|
|`Vitis hybrid cultivar`|`Vitis hybrid cultivar`|
|`Vitis riparia`|`Vitis riparia`|
|`Vitis vinifera`|`Vitis vinifera`, `Vitis rotundifolia`, `Vitis hybrid cultivar`, `Vitis aestivalis`, `Vitis riparia`, `Vitis cinerea var. helleri x vitis vinifera`, `Vitis cinerea var. helleri x vitis rupestris`, `Vitis cinerea var. helleri x vitis riparia`|
|`Xenopus laevis`|`Xenopus muelleri`, `Xenopus laevis x xenopus muelleri`, `Xenopus laevis`, `Xenopus laevis x xenopus borealis`, `Xenopus borealis`|
|`Zea mays`|`Zea mays`|

### Normalized Compendium Download Folder

Users will receive a zipped folder with a gene expression matrix aggregated by species, along with associated metadata.
Below is the detailed folder structure:

![docs-downloads-species-compendia](https://user-images.githubusercontent.com/15315514/65180873-ddbb4f80-da2b-11e9-97e9-127c68106182.png)

## RNA-Seq Sample Compendia

refine.bio RNA-seq sample compendia are comprised of the Salmon output for the collection of RNA-seq samples from an organism that we have processed with refine.bio.
Each individual sample has its own `quant.sf` file; the samples have not been aggregated and normalized.
RNA-seq sample compendia are designed to allow users that are comfortable handling these files to generate output that is most useful for their downstream applications.
Please see the [Salmon documentation on the `quant.sf` output format](https://salmon.readthedocs.io/en/latest/file_formats.html#quantification-file) for more information.

### RNA-Seq Sample Compendium Download Folder

Users will receive a zipped folder with individual `quant.sf` files for each sample that we were able to process with Salmon, grouped into folders based on the experiment those samples come from, along with any associated metadata in refine.bio.
Please note that our RNA-seq sample metadata is limited at this time and in some cases, we could not successfully run Salmon on every sample within an experiment (e.g., our processing infrastructure encountered an error with the sample, the sequencing files were malformed).
In addition, we use the terms "sample" and "experiment" to be consistent with the rest of refine.bio, but files will use run identifiers (e.g., SRR, ERR, DRR) and project identifiers (e.g., SRP, ERP, DRP), respectively.
Below is the detailed folder structure:

![docs-downloads-quantpendia](https://user-images.githubusercontent.com/15315514/65271488-4289af00-daeb-11e9-9006-2d7536c4a103.png)
