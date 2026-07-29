# Package index

## Core dataset processing

Build and populate the mmo object from MZmine outputs and annotation
tools.

<!-- end list -->

  - `GetMZmineFeature()` : Import MZmine feature table and metadata to
    create a mmo (mass-spectrometry metabolomics object)

  - `AddSiriusAnnot()` :
    
    Adding annotation from SIRIUS to the `mmo`

  - `AddChemSim()` :
    
    Add chemical similarity matrices to the `mmo`

  - `AddCustomSim()` :
    
    Add a custom feature similarity matrix to the `mmo`

  - `AddChemDist()` :
    
    Add chemical distance matrices to the `mmo`

  - `AddCustomDist()` :
    
    Add a custom feature distance matrix to the `mmo`

  - `AddCustomAnnot()` :
    
    Add custom annotations to a `mmo`

  - `ReplaceZero()` :
    
    \#' Replace zero and NA values in the `mmo`

  - `MassNormalization()` : Use sample mass in the metadata file to
    normalize the peak area

  - `LogNormalization()` :
    
    Log-normalize the peak area in the `mmo`

  - `MeancenterNormalization()` :
    
    Mean-center the peak area in the `mmo`

  - `ZNormalization()` :
    
    Z-normalize the peak area in the `mmo`

  - `ReorderGroups()` :
    
    Reorder samples in the `mmo` based on group order

  - `SwitchGroup()` :
    
    Switch the group column in the `mmo`

  - `FeaturePresence()` : Convert feature abundances to presence /
    absence

## Filtering

Subset the mmo object by samples, groups, or features, and keep
associated MGF spectral files in sync.

<!-- end list -->

  - `filter_mmo()` :
    
    Filter a `mmo` by samples, groups, and/or features

  - `filter_mgf_to_mmo()` :
    
    Filter an MGF file to keep only spectra for features present in
    `mmo$feature_data$id`

  - `annotate_feature_info_ms2_from_mgf()` :
    
    Annotate `mmo$feature_info` with MS2 presence and MS2 block counts
    from an MGF

  - `filter_canopus_annotations()` :
    
    Filter CANOPUS / SIRIUS annotations in an eCOMET `mmo` by
    probability threshold

  - `filter_cosmic_structure()` : Filter SIRIUS structure (CSI:FingerID)
    annotations by COSMIC confidence score

## Differential analysis & visualization

Identify and visualize differentially accumulated metabolites between
groups.

<!-- end list -->

  - `PairwiseComp()` :
    
    Perform pairwise comparison between two groups in the `mmo`

  - `GetDAMs()` :
    
    Generates lists of DAMs (Differentially Accumulated Metabolites) for
    each comparison in the `mmo`

  - `GetLog2FoldChange()` : Calculate log2 fold change for a given
    control group

  - `GetGroupMeans()` :
    
    Calculate group means from the `mmo`

  - `VolcanoPlot()` : Volcano plot for visualizing differential
    metabolite analysis results

  - `AnovaBarPlot()` : Generate barplots for each feature and perform
    ANOVA

  - `GenerateHeatmapInputs()` : GenerateHeatmapInputs (similarity-based)

## Ordination & clustering

Multivariate methods for comparing samples and groups.

<!-- end list -->

  - `PCAplot()` : Plots PCA and performs PERMANOVA
  - `PLSDAplot()` : PLS-DA plot with feature loadings
  - `NMDSplot()` : NMDSplot
  - `PCoAplot()` : PCoAplot
  - `HCplot()` : HCplot

## Chemical diversity

Alpha and beta diversity metrics incorporating feature abundance and
chemical distances.

<!-- end list -->

  - `GetAlphaDiversity()` : GetAlphaDiversity (similarity-based)
  - `GetBetaDiversity()` : GetBetaDiversity (similarity-based)
  - `GetRichness()` : GetRichness
  - `GetHillNumbers()` : GetHillNumbers
  - `GetFunctionalHillNumber()` : GetFunctionalHillNumber
    (similarity-based)
  - `GetFaithPD()` : GetFaithPD (similarity-based)
  - `RarefactionAUC()` : RarefactionAUC
  - `CalculateGroupBetaDistance()` : CalculateGroupBetaDistance
  - `GetSpecializationIndex()` : GetSpecializationIndex

## Chemical class analysis

Enrichment and visualization of compound classes using SIRIUS/CANOPUS
annotations.

<!-- end list -->

  - `CanopusLevelEnrichmentAnal()` : Enrichment analysis for
    Canopus-predicted terms
  - `CanopusListEnrichmentPlot()` : Generate a plot for enrichment
    analysis of Canopus-predicted terms
  - `CanopusListEnrichmentPlot_2()` : Generate a plot for enrichment
    analysis of Canopus-predicted terms across multiple levels
  - `CanopusLevelEnrichmentPlot()` : Generate a plot for enrichment
    analysis of Canopus-predicted terms at a specific level using a list
    of vectors of features
  - `CanopusAllLevelEnrichmentPlot()` : Generate a plot for enrichment
    analysis of Canopus-predicted terms across all levels
  - `PlotNPCStackedBar()` : PlotNPCStackedBar
  - `MSEA()` : Metabolite Set Enrichment Analysis (MSEA)

## Compound networks & dendrograms

Build chemical similarity trees and export to iTOL or Cytoscape.

<!-- end list -->

  - `FeatureDendrogram()` : FeatureDendrogram (similarity-based)
  - `PlotFeatureDendrogram()` : PlotFeatureDendrogram
  - `ExportITOL()` : ExportITOL
  - `ExportCytoscape()` : ExportCytoscape

## Phenotype association

Correlate individual features or metabolite sets with continuous
ecological variables.

<!-- end list -->

  - `FeaturePhenotypeCorrelation()` : FeaturePhenotypeCorrelation
  - `ScreenFeaturePhenotypeCorrelation()` : Screen feature-phenotype
    correlation
  - `GetPerformanceFeatureCorrelation()` :
    GetPerformanceFeatureCorrelation
  - `GetPerformanceFeatureLMM()` : GetPerformanceFeatureLMM
  - `GetPerformanceFeatureRegression()` :
    GetPerformanceFeatureRegression
  - `PlotFoldchangeResistanceQuad()` : PlotFoldchangeResistanceQuad
  - `PlotFoldchangeResistanceRegression()` :
    PlotFoldchangeResistanceRegression
  - `PlotFoldchangeResistanceRegression_t()` :
    PlotFoldchangeResistanceRegression\_t

## Save & export

Save mmo objects and export results.

<!-- end list -->

  - `SaveMMO()` : Save entire mmo object to a file (RDS)

  - `LoadMMO()` :
    
    Load a `mmo` previously saved with SaveMMO

  - `ExportFeaturesToCSV()` : ExportFeaturesToCSV

  - `pool_mmo_by_group()` : pool\_mmo\_by\_group

## Utilities

Helper functions for data access, formatting, and internal calculations.

<!-- end list -->

  - `GetSimMat()` :
    
    Retrieve a feature similarity matrix from the `mmo`

  - `GetDistanceMat()` :
    
    Get the distance matrix from the `mmo` based on the specified
    distance metric

  - `GetNormFeature()` :
    
    Retrieve feature data from the `mmo`, with normalization options

  - `FeatureToID()` :
    
    Convert feature names to IDs in the `mmo`

  - `IDToFeature()` :
    
    Convert feature IDs to names in the `mmo`

  - `print(<mmo>)` :
    
    Print method for `mmo`s Provides a clean, human-readable overview of
    an `mmo` list object instead of dumping the entire list when the
    object is printed in the console.

  - `anova_tukey_dunnett()` :
    
    Perform ANOVA and Tukey's HSD test on the `mmo`

  - `write_anova()` : Write results of anova\_tukey\_dunnett to a CSV
    file

  - `permanova_stat()` : Perform PERMANOVA and pairwise comparisons

## Dissimilarity-based versions (superseded)

Retained versions of the diversity and dendrogram functions that read
the older `.dissim` slots written by `AddChemDist()`. The functions
without the `_derep` suffix read the `.sim` slots written by
`AddChemSim()` and are the recommended entry points.

<!-- end list -->

  - `GetAlphaDiversity_derep()` : GetAlphaDiversity

  - `GetBetaDiversity_derep()` : GetBetaDiversity

  - `GetFaithPD_derep()` : GetFaithPD

  - `GetFunctionalHillNumber_derep()` : GetFunctionalHillNumber

  - `FeatureDendrogram_derep()` : FeatureDendrogram

  - `GenerateHeatmapInputs_derep()` :
    
    Generate input files to be used for pheatmap from the `mmo`
