# URD

**URD** is an R package designed for reconstructing transcriptional trajectories underlying specification or differentiation processes in the form of a branching tree, using single cell RNA-sequencing data. It is named after the [Norse mythological figure](https://en.wikipedia.org/wiki/Urdr) who nurtures the world tree and decides all fates.

**URD** is described in:
[Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis.](https://www.ncbi.nlm.nih.gov/pubmed/29700225)
Farrell JA & Wang Y (equal contribution), Riesenfeld SJ, Shekhar K, Regev A & Schier AF (equal contribution).
Science 26 Apr 2018. doi: 10.1126/science.aar3131

## Installation

- To install, follow the instructions in [INSTALL.md](INSTALL.md)

## Example code

- To get started, check out the Quick Start Tutorial: ([View](Analyses/QuickStart/URD-QuickStart-AxialMesoderm.md)) ([Download](Analyses/QuickStart/URD-QuickStart-AxialMesoderm.Rmd))

- The full supplementary analysis from Farrell & Wang, *et al.* that describes the reconstruction of developmental trajectories during zebrafish embryogenesis (3.3-12 hours post-fertilization) is also [available](Analyses/SupplementaryAnalysis)

- The full supplementary analysis from Siebert, Farrell, *et al.* that describes the reconstruction of developmental trajectories in adult Hydra is also [available](https://github.com/cejuliano/hydra_single_cell)

## Data

- The processed URD object from zebrafish embryogenesis (3.3 - 12 hours post-fertilization) from Farrell & Wang, *et al.* can be downloaded from the [Broad Single-cell Portal](https://portals.broadinstitute.org/single_cell/data/public/single-cell-reconstruction-of-developmental-trajectories-during-zebrafish-embryogenesis?filename=URD_Zebrafish_Object.rds) (requires log-in).

- The processed URD objects from adult Hydra from Siebert, Farrell, et al. can be downloaded from the [Broad Single-cell Portal](https://portals.broadinstitute.org/single_cell/study/SCP260/stem-cell-differentiation-trajectories-in-hydra-resolved-at-single-cell-resolution) (requires log-in) or from [Data Dryad](https://datadryad.org/resource/doi:10.5061/dryad.v5r6077).

## Version History

Detailed changes and minor versions are listed in [NEWS.md](NEWS.md)

#### Most recent bugfix: 1.1.0 (July 25, 2019)

#### 2019/07/25:
Version 1.1 released! This includes new functions for using NMF module analysis to detect and remove doublets, new functions for looking at gene expression relative to pseudotime, and new plotting functions for inspecting data. It accompanies the publication [Stem cell differentiation trajectories in Hydra resolved at single-cell resolution.](https://science.sciencemag.org/content/365/6451/eaav9314)

#### 2018/04/30:
Version 1.0 released! It accompanies our manuscript [Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis!](https://www.ncbi.nlm.nih.gov/pubmed/29700225)


	