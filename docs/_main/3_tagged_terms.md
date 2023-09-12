# Tagged Terms

A lack of uniformity in the values used in submitter-supplied metadata can make it challenging for users to find and use experiments that are relevant to their research.
For example, researchers submitting different experiments might use the labels "cytotoxic T cell" and "CD8 T cell" to refer to the same population of CD8+ T cells.
Although refine.bio metadata processing attempts to harmonize metadata keys, methods that tag or label samples with biomedical ontology terms offer a more rigorous approach to standardizing metadata values.
Biomedical ontologies provide principled ways to uniformly define and name entities such as cell types and represent relationships between entities (e.g., cytotoxic T cells are a subset of effector T cells) ([Bodenreider, Mitchell, and McCray. _Pac Symp Biocomput._ 2005.](https://doi.org/10.1142/9789812704856_0016)).
In the CD8+ T cell example above, an approach that tags samples with ontology terms may standardize samples from both experiments to `cytotoxic T cell` using the [Cell Ontology](https://doi.org/10.1186/s13326-016-0088-7) term [`CL:0000910`](http://purl.obolibrary.org/obo/CL_0000910).
We use the language "tagged term" to refer to any ontology term assigned to a sample via a prediction or labeling method, as detailed below.

<!--
## Assigning terms to samples

## Ontologies in use
-->
