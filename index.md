# Welcome to eCOMET: A Tool for Mass Spectrometry–Based Ecological Metabolomics

<!-- badges: start -->
[![R-CMD-check](https://github.com/Phytoecia/eCOMET/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Phytoecia/eCOMET/actions/workflows/R-CMD-check.yaml)
[![License: AGPL-3](https://img.shields.io/badge/License-AGPL--3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![GitHub issues](https://img.shields.io/github/issues/phytoecia/eCOMET)](https://github.com/phytoecia/eCOMET/issues)
[![GitHub stars](https://img.shields.io/github/stars/phytoecia/eCOMET)](https://github.com/phytoecia/eCOMET/stargazers)
<!-- badges: end -->

**eCOMET** is an R package designed specifically for ecological metabolomics (ecometabolomics), in which researchers apply untargeted mass spectrometry-based metabolomics to characterize the chemical phenotypes of organisms in their ecological and evolutionary context. Unlike metabolomics in biomedical or model-organism research, ecometabolomics frequently involve non-model, non-animal systems where genomic resources are scarce, reference spectral libraries provide limited coverage, and metabolic pathway annotations remain incomplete. Never the less, the goal of the package is to help chemical ecologists gain ecological insights from metabolomics data despite these contraints.

It has two core goals:

1. **Standardize ecometabolomics workflows** by introducing the `mmo` object, a unified data container that keeps feature abundances, sample metadata, compound annotations, and chemical similarity matrices aligned throughout an analysis. eCOMET integrates directly with the outputs of leading metabolomics tools — [MZmine](https://mzmine.github.io/), [SIRIUS/CANOPUS](https://v6.docs.sirius-ms.io/), [DreaMS](https://dreams-docs.readthedocs.io/), and [GNPS](https://gnps.ucsd.edu/) — making it easier to integrate data steams and improves reproducibility. 

2. **Provide educational resources for chemical ecologists** through annotated tutorials covering common study designs, from treatment-based experiments to multi-species field comparisons. Unlike biomedical metabolomics toolkits, eCOMET is built around the questions chemical ecologists actually ask: How does chemical diversity vary across species or habitats? Which compounds are differentially accumulated under stress? How does metabolomic composition relate to trait or performance data?

---

## Key Features

**Data ingestion and the `mmo` object**
- Import MZmine feature tables and link them to sample metadata
- Add SIRIUS/CANOPUS compound class annotations
- Add pairwise chemical similarity from DreaMS, MS2DeepScore, or cosine similarity
- Add custom annotations from in-house compound databases
- Filter the `mmo` object by sample, group, or feature list — all linked tables update together
- Match and filter associated MGF spectral files to keep MS2 data in sync

**Statistics and differential analysis**
- Differential accumulation analysis (DAMs) with pairwise group comparisons
- Volcano plots, ANOVA with post-hoc tests, and log2 fold-change summaries
- Metabolite set enrichment analysis (MSEA) using CANOPUS class annotations
- Chemical class enrichment analysis across multiple ClassyFire and NPC levels
- Phenotype association screening — correlate individual features with continuous ecological variables (height, herbivory, temperature, etc.)

**Chemical diversity**
- Alpha diversity: richness, Faith's PD, Hill numbers (q = 0, 1, 2), functional Hill numbers weighted by chemical distance, rarefaction curves
- Beta diversity: Bray-Curtis, Jaccard, CSCS (Chemical Structural and Compositional Similarity), and Generalized UniFrac using spectral similarity distances
- Ordination: PCA, PLS-DA, NMDS, PCoA, hierarchical clustering
- Specialization index for comparing chemical breadth across species or treatments

**Compound networks and dendrograms**
- Build feature dendrograms from pairwise chemical distances (DreaMS, MS2DeepScore, cosine, or custom)
- Incorporate Ion Identity Networking (IIN) and feature correlation constraints into dendrogram topology
- Visualize dendrograms with `ggtree` — circular or rectangular layouts, class coloring, fully layerable with ggplot2
- Export to iTOL for interactive annotation and publication-quality circular trees
- Export node and edge tables for Cytoscape molecular network visualization

**Visualization**
- All plot functions return `ggplot`/`ggtree` objects — add layers, themes, and annotations after the call
- Stacked bar plots of NPC compound class composition by group
- Heatmaps with flexible normalization and annotation tracks

---

<img width="2450" height="3500" alt="eCOMET workflow overview" src="./images/workflows/ecomet_workflow_total.png" />

---

## Installation

Install eCOMET from GitHub using `pak` (recommended — handles all dependencies automatically):

```r
# install.packages("pak")
pak::pak("phytoecia/eCOMET")
library(ecomet)
```

Or using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("phytoecia/eCOMET")
```


---

## Getting started

**New to eCOMET?** Start with the [Introduction](./articles/eCOMET-intro.html), which explains the `mmo` object, required input file formats, and the core preprocessing steps.

We currently support feature tables exported from **MZmine version 4 or later**. See the [MZmine documentation](https://mzmine.github.io/mzmine_documentation/module_docs/io/feat-list-export.html) for export instructions.

## Tutorials

eCOMET includes worked case studies covering two common study designs in chemical ecology:

1. [Treatment-based study](./articles/Tutorial_1_Treatment_based_study.html) — comparing chemical profiles between experimental treatments; covers normalization, DAM analysis, volcano plots, heatmaps, and enrichment analysis
2. [Interspecific comparisons](./articles/Tutorial_2_Interspecific_Comparisons.html) — comparing chemical diversity across co-occurring species; covers feature-based and structure-aware beta diversity, NMDS, PCoA, and hierarchical clustering
3. [Compound dendrograms and molecular networks](./articles/Tutorial_3_Dendrogram_and_Molecular_Networks.html) — building and exporting chemical similarity trees and networks for iTOL and Cytoscape

---

If you have questions, find a bug, or want to suggest a new feature, visit the [GitHub repository](https://github.com/phytoecia/eCOMET). Contributions are welcome — see the [Contributing](./articles/contributing.html) page for how to get involved.