library('readr')
library('iheatmapr') 
library('dplyr')
library('UpSetR')
library('tidyr')
library('stringr')
library('reshape')
library('readxl')
library('ggrepel')
library('ggbeeswarm')
library('caret') # for PLS-DA
library('pls') # for PLS-DA
library('entropy') # for shannon
library('DescTools') #Dunnet test
library('rlang')
library('pheatmap') # for heatmap
#library('lme4') # for mixed effect model
library('lmerTest') # for mixed effect model
library('ggbeeswarm') # for beeswarm plot
# library('ggvenn') # for venn diagram
library('RColorBrewer') # for color palette
library(viridis)
library('vegan') # for diversity index
library('ape') # for betadiv
library('GUniFrac') # for betadiv

################################################
setwd('/home/minsoo/software/mmo/250421_demo')



########################################################################################
# Define functions for data import and normalization
########################################################################################
# Load mzmine_feature.csv and metadata to create a mmo object
# In metadata file, filename should be datafile.___.mzML, then group, then mass
GetMZmineFeature <- function(mzmine_dir, metadata_dir){
  mmo <- list()
  data <- read.csv(mzmine_dir)
  #create feature column
  data <- data %>% mutate(feature = paste(mz, rt, sep = '_'))
  area_columns <- grep("datafile.+\\.mzML.area", names(data), value = TRUE)
  feature_df <- data %>% select(id, feature, all_of(area_columns))
  feature_df$id <- gsub(" ", "", feature_df$id)
  mmo$feature_data <- feature_df
  metadata <- read.csv(metadata_dir)
  metadata$sample <- paste('datafile.', metadata$sample, '.area', sep = '')
  mmo$metadata <- metadata
  mmo$pairwise <- data.frame(feature = mmo$feature_data$feature, id = mmo$feature_data$id)
  return(mmo)
}

# Add annotation from SIRIUS to the mmo object
AddSiriusAnnot <- function(mmo, canopus_structuredir, canopus_formuladir){
  structure_identifications <- read_tsv(canopus_structuredir, show_col_types = FALSE)
  structure_identifications$mappingFeatureId <- gsub(" ", "", structure_identifications$mappingFeatureId)
  canopus_formula_summary <- read_tsv(canopus_formuladir, show_col_types = FALSE)
  canopus_formula_summary$mappingFeatureId <- gsub(" ", "", canopus_formula_summary$mappingFeatureId)
  sirius_df <- mmo$feature_data %>% select(id, feature)
  sirius_df <- sirius_df %>%
  left_join(structure_identifications, by = c("id" = "mappingFeatureId")) %>%
  left_join(canopus_formula_summary, by = c("id" = "mappingFeatureId"))
  mmo$sirius_annot <- sirius_df
  return(mmo)
}
AddCustomAnnot <- function(mmo, DB, mztol = 5, rttol = 0.5) {
  DB <- read.csv(DB)
  DB <- DB %>%
    mutate(mz = as.numeric(mz), rt = as.numeric(rt))
  
  feature_annot <- mmo$feature_data %>%
    mutate(mz = as.numeric(sapply(str_split(feature, "_"), function(x) x[1])),
           rt = as.numeric(sapply(str_split(feature, "_"), function(x) x[2])))
  
  annotated_features <- feature_annot %>%
    mutate(custom_annot = purrr::map2(mz, rt, ~ DB %>%
                                   filter(abs(.x - mz) / .x * 1e6 <= mztol,
                                          abs(.y - rt) <= rttol) %>%
                                   pull(compound)))
  
  mmo$custom_annot <- annotated_features %>%
    select(id, feature, custom_annot)
  
  return(mmo)
}

#Replace 0 and NA to half of the smallest value or 1
ReplaceZero <- function(mmo, method = 'one') {
  df <- mmo$feature_data
  df[] <- apply(df, 1, function(row) {
    # Convert the row to numeric, ignoring non-numeric columns
    numeric_row <- as.numeric(row[-c(1, 2)])  # Skip 'id' and 'feature' columns
    # Get the smallest non-zero, non-NA value in the row
    smallest_value <- min(numeric_row[numeric_row > 0], na.rm = TRUE)
    # Replace 0 and NA with half of the smallest_value
    row[-c(1, 2)] <- sapply(numeric_row, function(x) {
      if (is.na(x) || x == 0) {
        if (method == 'one') {
          return(1)
        } else if (method == 'half_mean') {
          return(smallest_value / 2)
        }
      } else {
        return(x)
      }
    })
    
    return(row)
  }) %>%
    t() %>%
    as.data.frame()  # Convert back to dataframe
  mmo$feature_data <- df
  return(mmo)
  print('Missing values were filled')
}
# Use mass in the metadata file to normalize the peak area
MassNormalization <- function(mmo){
  normalized_df <- mmo$feature_data
  metadata <- mmo$metadata
  for (sample_col in colnames(mmo$feature_data)[-c(1,2)]) {
    sample_metadata <- metadata[metadata$sample == sample_col, ]
    mass <- sample_metadata$mass
    normalized_df[[sample_col]] <- as.numeric(mmo$feature_data[[sample_col]])/mass
  }
  mmo$feature_data <- normalized_df
  print("Peak area normalized by sample mass!")
  return(mmo)
}
# Make a df of log-transformed peak area in the mmo object (mmo$log)
LogNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  log_data <- log2(feature_data_only)
  log_df <- cbind(mmo$feature_data[, 1:2], log_data)
  mmo$log <- log_df
  return(mmo)
}
# Make a df of mean-centered peak area in the mmo object (mmo$meancentered)
MeancenterNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  mean_centered_data <- t(apply(feature_data_only, 1, function(x) x - mean(x, na.rm = TRUE)))
  mean_centered_df <- cbind(mmo$feature_data[, 1:2], mean_centered_data)
  mmo$meancentered <- mean_centered_df
  return(mmo)
}
# Make a df of Z-score of peak area in the mmo object (mmo$zscore)
ZNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  zscore_df <- t(apply(feature_data_only, 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  zscore_df <- cbind(mmo$feature_data[, 1:2], zscore_df)
  mmo$zscore <- zscore_df
  return(mmo)
}
GetCosSimToDissim <- function(mmo, cos.sim.dir){
  cos.sim <- read.csv(cos.sim.dir, col.names = c('cluster1', 'cluster2', 'metric', 'cosine_similarity', 'etc'))
  cos.sim$cosine_dissimilarity <- 1 - cos.sim$cosine_similarity
  # Create a pairwise matrix of cosine dissimilarity
  clusters <- unique(c(cos.sim$cluster1, cos.sim$cluster2))
  cos_dissim_mat <- matrix(1, nrow = length(clusters), ncol = length(clusters))
  rownames(cos_dissim_mat) <- clusters
  colnames(cos_dissim_mat) <- clusters

  # Fill the matrix with cosine dissimilarity values
  for (i in 1:nrow(cos.sim)) {
    cluster1 <- as.character(cos.sim$cluster1[i])
    cluster2 <- as.character(cos.sim$cluster2[i])
    cosine_dissimilarity <- cos.sim$cosine_dissimilarity[i]
    cos_dissim_mat[cluster1, cluster2] <- cosine_dissimilarity
    cos_dissim_mat[cluster2, cluster1] <- cosine_dissimilarity
  }
  mmo$cos.dissim <- cos_dissim_mat
  return(mmo)
}
########################################################################################
# Define functions for supporting analysis
########################################################################################
# Get feature data from the mmo object, with normalization options
GetNormFeature <- function(mmo, normalization = 'None'){
  if (normalization == 'None'){
    feature <- mmo$feature_data
  } else if (normalization == 'Log'){
    feature <- mmo$log
  } else if (normalization == 'Meancentered'){
    feature <- mmo$meancentered
  } else if (normalization == 'Z'){
    feature <- mmo$zscore
  }
  return(feature)
}

# Function to get list of IDs corresponding to a list of feature names
FeatureToID <- function(mmo, feature_names) {
  feature_data <- mmo$feature_data
  feature_ids <- feature_data %>%
    filter(feature %in% feature_names) %>%
    pull(id)
  return(feature_ids)
}

IDToFeature <- function(mmo, feature_ids) {
  feature_data <- mmo$feature_data
  feature_names <- feature_data %>%
    filter(id %in% feature_ids) %>%
    pull(feature)
  return(feature_names)
}

# Calculate group means for each feature
GetGroupMeans <- function(mmo, normalization = 'None', filter = FALSE, filter_groups) {
  feature_data <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  
  # Melt the feature data to long format
  long_feature_data <- melt(feature_data, id.vars = c('id', 'feature'), variable.name = 'sample', value.name = 'feature_value')
  colnames(long_feature_data) <- c('id', 'feature', 'sample', 'feature_value')
  
  # Merge with metadata to get group information
  merged_data <- merge(long_feature_data, metadata[, c('sample', 'group')], by = 'sample')
  
  # Calculate group means
  group_means <- merged_data %>%
    group_by(group, id) %>%
    summarise(mean_value = mean(feature_value, na.rm = TRUE)) %>%
    spread(key = group, value = mean_value)
  if (filter == TRUE){
    group_means <- group_means %>% select(id, all_of(filter_groups))
  }
  
  return(group_means)
}

# Function to calculate fold change from gruop mean data
GetLog2FoldChange <- function(group_means, control_group) {
  control_means <- group_means[[control_group]]
  fold_change <- group_means %>%
    mutate(across(-id, ~ log2(. / control_means)))
  
  return(fold_change)
}


########################################################################################
# Define functions for pairwise comparison and visualization
########################################################################################

# Calculate log2FC and adjusted p-value for the designated groups
# FC is calculated by group2/group1
# correction : multiple comparison correction method. See p.adjust() 
# Adds columns to the mmo$pairwise
PairwiseComp <- function(mmo, group1, group2, correction = 'BH'){
  feature <- mmo$feature_data
  metadata <- mmo$metadata
  #Get sample names
  group1_samples <- metadata %>% filter(group == group1) %>% pull(sample)
  group2_samples <- metadata %>% filter(group == group2) %>% pull(sample)
  #Get data from the samples
  group1_data <- feature %>% select(id, feature, all_of(group1_samples))
  group2_data <- feature %>% select(id, feature, all_of(group2_samples))
  #Make empty column
  log2FC <- numeric(nrow(feature))
  pval <- numeric(nrow(feature))
  #Pairwise comparison
  for (i in 1:nrow(feature)){
    group1_value <- as.numeric(group1_data[i, -c(1,2)])
    group2_value <- as.numeric(group2_data[i, -c(1,2)])

    group1_mean <- mean(group1_value, na.rm = TRUE)
    group2_mean <- mean(group2_value, na.rm = TRUE)
    log2FC[i] <- log2(group2_mean/group1_mean)

    pval[i] <- tryCatch(
      expr = {
        p <- t.test(group1_value, group2_value, na.rm = TRUE)$p.value
        p
      }, 
      error = function(e) 
      {
        return(1)
      }
    )
    
    # ttest <- t.test(group1_value, group2_value, na.rm = TRUE)
    
    # pval[i] <- ttest$p.value
  }
  padj <- p.adjust(pval, method = correction)
  #Store in results
  results <- data.frame(
    log2FC= log2FC,
    padj = padj
  )
  names(results) <- c(paste(group1, "vs", group2, "log2FC", sep = "_"), paste(group1, "vs", group2, "padj", sep = "_"))
  #Add pairwise results to the mmo object
  mmo$pairwise <- cbind(mmo$pairwise, results)
  print(paste(group2, '/', group1, 'comparison was completed'))
  return(mmo) 
}

# Plot log2FC and -log(pval) from pairwise comparison.  PairwiseComp(mmo, 'group1', 'group2') should be precended
# mmo : mmo object with pairwise comparison matrix
# topk : the number of features to sshow labels
# pthr : the threshold of adjusted p-value to be considered significant
# outdir : plot name
# height : plot height (<20)
# width  : plot width (<20)
VolcanoPlot <- function(mmo, comp, topk = 10, pthr = 0.1, outdir = 'volcano.png', height = 5, width = 5){
  VolData <- mmo$pairwise %>% select(feature,all_of(c(paste(comp, 'log2FC', sep = '_'), paste(comp, 'padj', sep = '_'))))
  colnames(VolData) <- c('feature', 'log2FC', 'padj')
  VolData <- VolData %>% 
    mutate(
      Expression = case_when(log2FC >= 1 & padj <= pthr ~ "Up-regulated",
                            log2FC <= -1 & padj <= pthr ~ "Down-regulated",
                            TRUE ~ "Not significant")
      )

  top_features <- bind_rows(
    VolData %>%
      filter(Expression =='Up-regulated')  %>%
      arrange(desc(abs(log2FC)), padj) %>%
      head(topk),
    VolData %>%
      filter(Expression == 'Down-regulated') %>%
      arrange(desc(abs(log2FC)), padj) %>%
      head(topk)
  )

  volcano <- ggplot(VolData, aes(x = log2FC, y = -log(padj, 10))) +
    geom_point(aes(color = Expression), size = 0.4)+
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"FDR"))+
    scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    theme_classic()+
    geom_label_repel(data = top_features,
                    mapping = aes(log2FC, -log(padj,10), label = feature),
                    size = 2)
  
  volcano
  ggsave(outdir, height = height, width = width)
}

########################################################################################
# Define functions for multivariate analysis
########################################################################################

# PLS-DA plot
# colors : for scale_color_manual
# topk : number of features to show loading
# Choose which feature area to use for calculation among 'None' (mass-normalized), 'Log','Meancentered', and 'Z'
PLSDAplot <- function(mmo, color, topk = 10, outdir, normalization = 'Z', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL) {
  metadata <- mmo$metadata
  #Get appropriate feature by normalization parameter
  feature <- GetNormFeature(mmo, normalization)

  # All feature or filtered feature
  if (filter_feature == TRUE){
    feature <- feature %>% filter(feature %in% feature_list)
  }
  
  X <- t(as.matrix(feature[, -(1:2)]))
  Y <- c()
  for (col in colnames(feature)[-c(1,2)]){
    Y <- append(Y, metadata[metadata$sample == col, ]$group)
  }
  Y <- as.factor(Y)

  plsda_model <- plsda(X, Y, ncomp = 2)
  scores <- plsda_model$scores[, 1:2]
  plsda_df <- data.frame(Comp1 = scores[, 1], Comp2 = scores[, 2], Group = Y)
  loadings <- plsda_model$loadings
  loadings_comp1 <- loadings[, 1]
  loadings_comp2 <- loadings[, 2]  
  if (filter_feature == FALSE){
    loadings_df <- data.frame(Feature = mmo$feature_data$feature, 
                            Comp1_Loading = loadings_comp1, 
                            Comp2_Loading = loadings_comp2)
  } else {
    loadings_df <- data.frame(Feature = feature_list, 
                            Comp1_Loading = loadings_comp1, 
                            Comp2_Loading = loadings_comp2)
  }


  top_features <- loadings_df %>%
  mutate(abs_loading_comp1 = abs(Comp1_Loading),
         abs_loading_comp2 = abs(Comp2_Loading)) %>%
  arrange(desc(abs_loading_comp1 + abs_loading_comp2)) %>%
  head(topk)
                          

  if (filter_group == TRUE){
    plsda_df <- plsda_df %>% filter(Group %in% group_list)
  }

  ggplot(plsda_df, aes(x = Comp1, y = Comp2, color = Group)) +
  geom_point(size = 3) +
  theme_classic() +
  stat_ellipse(level = 0.90) +
  ggtitle("PLS-DA Plot") +
  labs(x = "Component 1", y = "Component 2") +
  scale_color_manual(values = custom_colors) +
  theme(legend.position = "right")+
  geom_segment(data = top_features,
               aes(x = 0, y = 0, xend = Comp1_Loading * 500, yend = Comp2_Loading * 500),  # Scale the arrows
               arrow = arrow(length = unit(0.3, "cm")), color = "grey", linewidth = 1) +
  # Add labels for the top 10 features
  geom_text_repel(data = top_features,
            aes(x = Comp1_Loading * 500, y = Comp2_Loading * 500, label = Feature),
            color = "black", vjust = 1.5, size = 3)

  ggsave(outdir, height = 6, width = 6)
  write.csv(loadings_df, 'PLSDA_loadings.csv')
  print(paste(normalization, '-normalized feature was used'))
}


################################################
#Define enrichment analysis using Canopus-predicted terms
# list_test : a vector containing names of features to analyze
# mmo : the mmo object with sirius annotation and normalized
# pthr : the threshold for adjusted p-value to be considered significant
# sig : a logical vaue to show only significant terms or not
# representation : 'greater' for overrepresentation
CanopusEnrichmentAnal <- function(mmo,list_test, pthr = 0.1, sig=TRUE, representation = 'greater'){
  all_feature <- mmo$sirius_annot
  subset_feature <- mmo$sirius_annot %>% filter(feature %in% list_test)
  # Split the classifications into separate terms
  all_feature$classifications_split <- str_split(all_feature[,46], ";")
  subset_feature$classifications_split <- str_split(subset_feature[,46], ";")
  # make a single list of all the terms are stored
  expanded_all <- all_feature %>%
    unnest(classifications_split) %>% 
    mutate(classifications_split = str_trim(classifications_split))
  expanded_subset <- subset_feature %>%
    unnest(classifications_split) %>%
    mutate(classifications_split = str_trim(classifications_split))
  # Create contingency tables
  total_term_counts <- table(expanded_all$classifications_split)
  subset_term_counts <- table(expanded_subset$classifications_split)  
  # Perform enrichment analysis using Fisher's exact test
  enrichment_results <- sapply(names(subset_term_counts), function(term) {
    contingency_matrix <- matrix(c(
      subset_term_counts[[term]],
      sum(subset_term_counts) - subset_term_counts[[term]],
      total_term_counts[[term]],
      sum(total_term_counts) - total_term_counts[[term]]
    ), nrow = 2, byrow = TRUE)
    fisher.test(contingency_matrix, alternative = representation)$p.value
  })
  # Adjust p-values for multiple testing
  adjusted_pvalues <- p.adjust(enrichment_results, method = "fdr")
  # Create a results dataframe
  results <- data.frame(
    term = names(enrichment_results),
    subsetcount = as.numeric(subset_term_counts[names(enrichment_results)]),
    totalcount = as.numeric(total_term_counts[names(enrichment_results)]),
    pval = enrichment_results,
    fdr = adjusted_pvalues
  )
  # Filter for significantly enriched terms
  significant_terms <- results %>%
    filter(fdr < pthr) %>%
    arrange(fdr)
  if(sig==TRUE){
    return(significant_terms)
  }else{
    return(results)
  }
}
################################################
#Define enrichment analysis using Canopus-predicted terms
# mmo : the mmo object with sirius annotation and normalized
# list_test : a vector containing names of features to analyze
# pthr : the threshold for adjusted p-value to be considered significant
# sig : a logical vaue to show only significant terms or not
# term_level : the level of term to use for enrichment analysis
# representation : 'greater' for overrepresentation
CanopusLevelEnrichmentAnal <- function(mmo,list_test, pthr = 0.1, sig=TRUE, term_level = 'NPC_pathway', representation = 'greater'){
  all_feature <- mmo$sirius_annot 
  subset_feature <- mmo$sirius_annot %>% filter(feature %in% list_test)
  print(paste('total features:', nrow(all_feature), 'list_test features:', nrow(subset_feature)))
  # Select the appropriate term level for enrichment analysis
  if (term_level == "NPC_pathway") {
    all_feature$classifications_split <- all_feature[[30]]
    subset_feature$classifications_split <- subset_feature[[30]]    
  } else if (term_level == "NPC_superclass") {
    all_feature$classifications_split <- all_feature[[32]]
    subset_feature$classifications_split <- subset_feature[[32]]
  } else if (term_level == "NPC_class") {
    all_feature$classifications_split <- all_feature[[34]]
    subset_feature$classifications_split <- subset_feature[[34]]
  } else if (term_level == "ClassyFire_superclass") {
    all_feature$classifications_split <- all_feature[[36]]
    subset_feature$classifications_split <- subset_feature[[36]]
  } else if (term_level == "ClassyFire_class") {
    all_feature$classifications_split <- all_feature[[38]]
    subset_feature$classifications_split <- subset_feature[[38]]
  } else if (term_level == "ClassyFire_subclass") {
    all_feature$classifications_split <- all_feature[[40]]
    subset_feature$classifications_split <- subset_feature[[40]]
  } else if (term_level == "ClassyFire_level5") {
    all_feature$classifications_split <- all_feature[[42]]
    subset_feature$classifications_split <- subset_feature[[42]]
  } else if (term_level == "ClassyFire_most_specific") {
    all_feature$classifications_split <- all_feature[[44]]
    subset_feature$classifications_split <- subset_feature[[44]]
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  } 

  total_term_counts <- table(unlist(all_feature$classifications_split))
  subset_term_counts <- table(unlist(subset_feature$classifications_split))

  total_term_counts['None'] <- sum(is.na(all_feature$classifications_split))
  subset_term_counts['None'] <- sum(is.na(subset_feature$classifications_split))

  # Perform enrichment analysis using Fisher's exact test
  enrichment_results <- sapply(names(subset_term_counts), function(term) {
    contingency_matrix <- matrix(c(
      subset_term_counts[[term]],
      sum(subset_term_counts) - subset_term_counts[[term]],
      total_term_counts[[term]],
      sum(total_term_counts) - total_term_counts[[term]]
    ), nrow = 2, byrow = TRUE)
    fisher.test(contingency_matrix, alternative = representation)$p.value
  })
  # Adjust p-values for multiple testing
  adjusted_pvalues <- p.adjust(enrichment_results, method = "fdr")
  # Create a results dataframe
  results <- data.frame(
    term_level = term_level,
    term = names(enrichment_results),
    subsetcount = as.numeric(subset_term_counts[names(enrichment_results)]),
    totalcount = as.numeric(total_term_counts[names(enrichment_results)]),
    foldenrichment = (as.numeric(subset_term_counts[names(enrichment_results)]) / length(list_test))/(as.numeric(total_term_counts[names(enrichment_results)]) / nrow(all_feature)),
    pval = enrichment_results,
    fdr = adjusted_pvalues
  )
  # Filter for significantly enriched terms
  significant_terms <- results %>%
    filter(pval < pthr) %>%
    arrange(pval)
  if(sig==TRUE){
    return(significant_terms)
  }else{
    return(results)
  }
}

# Generate plot and result files for the list of pairwise comparisons
# Here is an example of comp.list:
# comp.list <- list(sl1.up, sl2.up, sl3.up, px1.up, px2.up, px3.up, le1.up, le2.up, le3.up)
# names(comp.list) <- c('sl1.up', 'sl2.up', 'sl3.up', 'px1.up', 'px2.up', 'px3.up', 'le1.up', 'le2.up', 'le3.up')
CanopusEnrichmentPlot <- function(mmo = mmo, comp.list, pthr = 0.1, representation = 'greater', prefix = 'enrichment'){
  df.EA <- data.frame()
  sig.terms <- c()
  for(list in names(comp.list)){
    # Calculate enrichment score for all terms
    res <- CanopusEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation)
    res <- res %>% mutate(comp = list)
    df.EA <- bind_rows(df.EA, res)
    # get terms that are at least once enriched in one comparison
    res.sig <- CanopusEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]], sig=TRUE, pthr = pthr, representation = representation)
    sig.terms <- append(sig.terms, res.sig$term)
  }
  sig.terms <- unique(sig.terms)
  df.EA.sig <- df.EA %>% filter(term %in% sig.terms)
  df.EA.sig <- df.EA.sig %>% 
    mutate(label = cut(
        fdr,
        breaks = c(0,0.001, 0.01, 0.05, 1),
        labels = c("***", "**", "*", "")
    ))

  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = comp, y = term, label = label))+
    geom_point(aes(size = subsetcount, color = fdr))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')
  enrichment_plot
  write.csv(df.EA, paste(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste(prefix, '.pdf'), width = 6, height = 8)
}

CanopusLevelEnrichmentPlot <- function(mmo = mmo, comp.list, term_level = 'NPC_pathway',pthr = 0.1, representation = 'greater', prefix = 'enrichment'){
  df.EA <- data.frame()
  sig.terms <- c()
  for(list in names(comp.list)){
    # Calculate enrichment score for all terms
    res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level)
    res <- res %>% mutate(comp = list)
    df.EA <- bind_rows(df.EA, res)
    # get terms that are at least once enriched in one comparison
    res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level)
    sig.terms <- append(sig.terms, res.sig$term)
  }
  sig.terms <- unique(sig.terms)
  df.EA.sig <- df.EA %>% filter(term %in% sig.terms)
  df.EA.sig <- df.EA.sig %>% 
    mutate(label = cut(
        pval,
        breaks = c(0,0.001, 0.01, 0.05, 1),
        labels = c("***", "**", "*", "")
    ))

  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = comp, y = term, label = label))+
    geom_point(aes(size = subsetcount, color = pval))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')
  enrichment_plot
  write.csv(df.EA, paste(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste(prefix, '.pdf'), width = 6, height = 8)
}
##
# Generate plot and result files for the list of pairwise comparisons
CanopusAllLevelEnrichmentPlot <- function(mmo = mmo, comp.list, terms = 'all_terms',pthr = 0.1, representation = 'greater', prefix = 'enrichment'){
  df.EA <- data.frame()
  sig.terms <- c()
  if(terms == 'all_terms'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'NPC'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway')
  } else if (terms == 'ClassyFire'){
    term_levels = c('ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  }
  for(term_level in term_levels){
    for(list in names(comp.list)){
      # Calculate enrichment score for all terms
      res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level)
      res <- res %>% mutate(comp = list)
      df.EA <- bind_rows(df.EA, res)
      # get terms that are at least once enriched in one comparison
      res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level)
      sig.terms <- append(sig.terms, res.sig$term)
    }
    sig.terms <- unique(sig.terms)
    df.EA.sig <- df.EA %>% filter(term %in% sig.terms)
    df.EA.sig <- df.EA.sig %>% 
      mutate(label = cut(
          pval,
          breaks = c(0,0.001, 0.01, 0.05, 1),
          labels = c("***", "**", "*", "")
      ))
  }
  enrichment_plot <- ggplot(data = df.EA.sig, aes(x = comp, y = term, label = label))+
    geom_point(aes(size = subsetcount, color = pval))+
    geom_text()+
    scale_size_area(name = 'Count', max_size = 10)+
    scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    xlab('Comparisons')+
    ylab('Chemical classes')
  enrichment_plot
  write.csv(df.EA, paste(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste(prefix, '.pdf'), width = 6, height = 8)
}

# Function to get performance individual feature regression
# mmo : mmo object
# target : the name of feature 
# herbivore : the name of herbivore performance in metadata
# groups : the group from metadata containing performance data
# DAM.list : a list of DAMs
# normalization : 'None', 'Log', 'Meancentered', 'Z'
# model is lmm or lm
# output : output file name for the plot
FeaturePerformanceRegression <- function(mmo, target, herbivore, groups, DAM.list, model = 'lmm', normalization = 'Z', output){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # Get herbivore performance from the metadata, get the feature value from the feature matrix, then combine
  herbivore.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,herbivore]) %>% filter(group %in% groups)
  feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[feature$feature == target, -(1:2)]))
  combined_df <- merge(herbivore.df, feature_df, by='sample')
  
  # Perform linear mixed model or simple linear regression
  if (model == 'lmm'){
    fit <- lmer(combined_df$performance ~ combined_df$feature_value + (1|combined_df$group))
    p_value <- summary(fit)$coefficients[2, 5]
  } else if (model == 'lm'){
    fit <- lm(combined_df$performance ~ combined_df$feature_value)
    p_value <- summary(fit)$coefficients[2, 4]
  } else {
    stop("Invalid model type. Please use 'lmm' or 'lm'")
  }

  # Plot the fit using ggplot
  ggplot(combined_df, aes(x = feature_value, y = performance, color = group)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    theme_classic() +
    labs(title = paste("Regression of", target, "against", herbivore, "performance"),
         x = "Feature Value",
         y = "Performance") +
    theme(legend.position = "right") +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 3)), 
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
  ggsave(output, height = 6, width = 6)
}

# Function to get performance feature regression
# From the metadata, the herbivore performance is fetched
# From the feature data, the performance is regressed against each feature
# mmo : mmo object
# herbivore : the name of herbivore performance in metadata
# groups : the group from metadata containing performance data
# DAM.list : a list of DAMs
# comparisons : a list of pairwise comparisons
# normalization : 'None', 'Log', 'Meancentered', 'Z'
GetPerformanceFeatureRegression <- function(mmo, herbivore, groups, DAM.list, comparisons, normalization = 'Z'){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata
  
  # herbivore.sample <- metadata %>% filter(group %in% groups) %>% pull(sample)
  # herbivore.area <- feature %>% select(id, feature, all_of(herbivore.sample))

  performance.linreg <- data.frame(pval = double(), effect.size = double())
  herbivore.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,herbivore]) %>% filter(group %in% groups)

  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)  
  for (i in 1:nrow(feature)) {
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(herbivore.df, feature_df, by='sample')

    fit <- lm(combined_df$performance ~ combined_df$feature_value)
    effect.size <- coef(fit)[2]
    p_value <- summary(fit)$coefficients[2, 4]
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }
    #is.Spec <- feature_name %in% target

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = effect.size, p_value = p_value, tag = tag
    ))
  }
  # Add FC columns to regression_results
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

GetPerformanceFeatureLMM <- function(mmo, herbivore, groups, DAM.list, comparisons, normalization = 'Z'){
  feature <- GetNormFeature(mmo, normalization)
  # if (normalization == 'None'){
  #   feature <- mmo$feature_data
  # } else if (normalization == 'Log'){
  #   feature <- mmo$log
  # } else if (normalization == 'Meancentered'){
  #   feature <- mmo$meancentered
  # } else if (normalization == 'Z'){
  #   feature <- mmo$zscore
  # }
  metadata <- mmo$metadata

  # get herbivore performance data
  herbivore.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,herbivore]) %>% filter(group %in% groups)
  #create an empty dataframe to store regression results
  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)  
  # iterate regression analysis
  for (i in 1:nrow(feature)) {
    # for each feature, generate herbivore performance X feature value data
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(herbivore.df, feature_df, by='sample')

    # linear mixed model
    lmm_fit <- lmer(performance ~ feature_value + (1|group), data = combined_df)
    fixed_effects <- fixef(lmm_fit)
    effect.size <- fixed_effects[2]
    p_value <- summary(lmm_fit)$coefficients[2, 5]
    
    #tag using DAM.list
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = effect.size, p_value = p_value, tag = tag
    ))
  }
  # Add FC columns to regression_results
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

# Plot the regression results from GetPerformanceFeatureRegression or GetPerformanceFeatureLMM
# performance_regression : the regression results from GetPerformanceFeatureRegression or GetPerformanceFeatureLMM
# fold_change : the name of the fold change column in mmo$pairwise
# color : the color palette for the plot
# output_dir : the output directory for the plot
PlotFoldchangeResistanceRegression <- function(performance_regression, fold_change, color, output_dir){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste("-effect.size ~", fold_change)))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[2, 4]
  r_squared <- summary_fit$r.squared

  ggplot(performance_regression, aes(x = !!sym(fold_change), y = -effect.size)) +
    geom_point(size = 0.5, aes(color = tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab(fold_change) +
    ylab('-effect.size') +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}

PlotFoldchangeResistanceRegression_t <- function(performance_regression, fold_change, color, output_dir){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste(fold_change, "~ -effect.size")))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[4]
  r_squared <- summary_fit$r.squared

  ggplot(performance_regression, aes(x = -effect.size, y = !!sym(fold_change))) +
    geom_point(size = 0.5, aes(color = tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}

# Function to generate CSV files containing peak area 
# input_file : a CSV file containing compound name, mz, and rt
# mztol : mz tolerance for feature matching
# rtwin : rt window for feature matching
# normalization : 'None', 'Log', 'Meancentered', 'Z'
# output_dir : the directory to save the CSV file
QuantToCSV <- function(mmo, input_file, normalization = 'None', output_dir) {
  # Read the input file
  input_data <- read.csv(input_file)
  
  # Select the appropriate feature data based on normalization parameter
  feature <- GetNormFeature(mmo, normalization)
  
  # Initialize an empty list to store results
  results_list <- list()
  
  # Iterate through each row in the input data
  for (i in 1:nrow(input_data)) {
    compound <- input_data$compound[i]
    mz <- input_data$mz[i]
    rt <- input_data$rt[i]
    mztol <- input_data$mztol[i]
    rtwin <- input_data$rtwin[i]
    
    # Filter the feature data based on mz and rt tolerance
    filtered_data <- feature %>%
      filter(abs(as.numeric(sapply(str_split(feature, "_"), function(x) x[1])) - mz) / mz * 1e6 <= mztol,
             abs(as.numeric(sapply(str_split(feature, "_"), function(x) x[2])) - rt) <= rtwin)
    
    # If no matching features are found, skip to the next iteration
    if (nrow(filtered_data) == 0) {
      next
    }
    
    # Extract the peak area for each group
    peak_area <- filtered_data %>%
      select(-id, -feature)
    
    # Add the compound name to the peak area data
    peak_area <- cbind(compound = compound, peak_area)
    
    # Append the peak area data to the results list
    results_list[[i]] <- peak_area
  }
  
  # Combine all results into a single data frame
  results_df <- do.call(rbind, results_list)
  
  # Write the results to a CSV file
  write.csv(results_df, file = output_dir, row.names = FALSE)
  
  print(paste("Peak area CSV file generated at:", output_dir))

  return(results_df)
}

# Draw barplots for each row in quan
DrawBarplots <- function(mmo, quan) {
  metadata <- mmo$metadata
  # Merge quan with metadata to get group information
  quan_long <- melt(quan, id.vars = "compound", variable.name = "sample", value.name = "peak_area")
  colnames(quan_long) <- c("compound", "sample", "peak_area")
  quan_long <- merge(quan_long, metadata[, c("sample", "group")], by = "sample")
  
  # Iterate through each compound and draw barplot
  for (compound in unique(quan_long$compound)) {
    compound_data <- quan_long %>% filter(compound == !!compound)
    
    p <- ggplot(compound_data, aes(x = group, y = peak_area, fill = group)) +
      geom_bar(stat = "summary", fun = "mean", position = "dodge") +
      geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.9), width = 0.2) +
      geom_beeswarm() +
      theme_classic() +
      labs(title = paste("Peak Area for", compound), x = "Group", y = "Peak Area") +
      theme(legend.position = "none")
    
    # Save the plot
    ggsave(paste0(compound, "_barplot.png"), plot = p, width = 6, height = 4)
  }
}

# Function to export features to CSV file
ExportFeaturesToCSV <- function(mmo, feature_list, normalization = 'None', output_dir){
  feature <- GetNormFeature(mmo, normalization = normalization) # Get normalized feature data
  # Filter the feature data, annotation, and DA analysis for the list provided
  selected_annotations <- mmo$sirius_annot %>% filter(feature %in% feature_list) %>% 
    select(id = 1, feature = 2, formula = 12, compound = 15, pubchem = 18, ionmass = 21, NPC_pathway = 30, NPC_superclass = 32, NPC_class = 34, ClassyFire_superclass = 36, ClassyFire_class = 38, ClassyFire_subclass = 40, ClassyFire_level5 = 42, ClassyFire_most_specific = 44)
  selected_feature <- feature %>% filter(feature %in% feature_list)
  selected_pairwise <- mmo$pairwise %>% filter(feature %in% feature_list)
  # Merge all
  merged_df <- merge(selected_annotations, selected_feature, by = 'feature')
  merged_df <- merge(merged_df, selected_pairwise, by = 'feature')

  write.csv(merged_df, output_dir)
}

##### Chemical Diversity indices #####
# Function to calculate unweighted Hill numbers
GetHillNumbers <- function(mmo, normalization = 'None', q = 0) {
  feature <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  
  hill_numbers <- apply(feature[, -(1:2)], 2, function(x) {
    p <- x / sum(x)
    if (q == 0) {
      return(length(p))
    } else if (q == 1) {
      return(exp(-sum(p * log(p))))
    } else {
      return((sum(p^q))^(1 / (1 - q)))
    }
  })
  
  groups <- c()
  for (col in colnames(feature)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }
  
  hill_df <- data.frame(group = groups, hill_number = hill_numbers)
  
  
  return(hill_df)
}
# Function to calculate functional Hill numbers weighted by qemistree distance
GetFunctionalHillNumber <- function(mmo, normalization = 'None',q = 1){
  feature <- GetNormFeature(mmo, normalization = 'None')
  metadata <- mmo$metadata
  # Scale the qemistree distance matrix to be between 0 and 1
  scaled_dissimilarity <- mmo$qemistree / max(mmo$qemistree)
  # Calculate the relative proportions of each feature and reorder them to match the order of the distance matrix
  q.feature <- feature %>% filter(id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ]

  # Calculate Hill
  functional_hill_number <- c()
  if (q == 1){
    for (column in 1:ncol(relative_proportions)){
      val = 0
      for (i in 1:nrow(scaled_dissimilarity)){
        for (j in 1:nrow(scaled_dissimilarity)){
          val = val + scaled_dissimilarity[i,j]*relative_proportions[i,column] * relative_proportions[j,column] * log(relative_proportions[i,column] * relative_proportions[j,column])
        }
      }
      functional_hill_number[column] <- exp(-val)
      print(paste(column,'per',ncol(relative_proportions), 'columns processed'))
    }    
  } else {
    for (column in 1:ncol(relative_proportions)){
      val = 0
      for (i in 1:nrow(scaled_dissimilarity)){
        for (j in 1:nrow(scaled_dissimilarity)){
          val = val + scaled_dissimilarity[i,j] * (relative_proportions[i,column] * relative_proportions[j,column])^q
        }
      }
      functional_hill_number[column] <- val^(1/(1-q))
      print(paste(column,'per',ncol(relative_proportions), 'columns processed'))
    }
  }
  names(functional_hill_number) <- colnames(relative_proportions)
  # Get the group information
  groups <- c()
  for (col in colnames(feature)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }
  
  hill_df <- data.frame(group = groups, hill_number = functional_hill_number)
  return(hill_df)
}
GetFunctionalHillNumberwithRaoQ <- function(mmo, normalization = 'None',q = 1){
  feature <- GetNormFeature(mmo, normalization = 'None')
  metadata <- mmo$metadata
  # Scale the qemistree distance matrix to be between 0 and 1
  scaled_dissimilarity <- mmo$qemistree / max(mmo$qemistree)
  # Calculate the relative proportions of each feature and reorder them to match the order of the distance matrix
  q.feature <- feature %>% filter(id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ]

  functional_hill_number <- c()
  # Calculate Hill
  if (q == 1){
    for (column in 1:ncol(relative_proportions)){
      raoQ <- sum(as.matrix(scaled_dissimilarity) %*% relative_proportions[,column] %*% t(relative_proportions[,column]))
      val = 0
      for (i in 1:nrow(scaled_dissimilarity)){
        for (j in 1:nrow(scaled_dissimilarity)){
          val = val + scaled_dissimilarity[i,j]*relative_proportions[i,column] * relative_proportions[j,column]/raoQ * log(relative_proportions[i,column] * relative_proportions[j,column]/raoQ)
        }
      }
      functional_hill_number[column] <- exp(-val)
      print(paste(column,'per',ncol(relative_proportions), 'columns processed'))
    }    
  } else {
    for (column in 1:ncol(relative_proportions)){
      val = 0
      for (i in 1:nrow(scaled_dissimilarity)){
        for (j in 1:nrow(scaled_dissimilarity)){
          val = val + scaled_dissimilarity[i,j] * (relative_proportions[i,column] * relative_proportions[j,column]/raoQ)^q
        }
      }
      functional_hill_number[column] <- val^(1/(1-q))
      print(paste(column,'per',ncol(relative_proportions), 'columns processed'))
    }
  }
  names(functional_hill_number) <- colnames(relative_proportions)
  # Get the group information
  groups <- c()
  for (col in colnames(feature)[-c(1, 2)]) {
    groups <- append(groups, metadata[metadata$sample == col, ]$group)
  }
  
  hill_df <- data.frame(group = groups, hill_number = functional_hill_number)
  return(hill_df)
}
GetBetaDiversity <- function(mmo, method = 'Gen.Uni', normalization = 'None') {    
  # Get compound distance and build tree for UniFrac
  scaled_dissimilarity <- mmo$qemistree / max(mmo$qemistree)
  compound_tree <- as.phylo(hclust(as.dist(scaled_dissimilarity), method = "average"))

  # Get feature matrix of relative proportion
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata
  feature <- feature %>% filter(id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x)) 
  rownames(relative_proportions) <- feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ] #reorder
  relative_proportions <- t(relative_proportions)
  # Calculate Generalized UniFrac
  if (method == 'Gen.Uni') {
    guni <- GUniFrac(relative_proportions, compound_tree, alpha = c(0, 0.5, 1), verbose = TRUE)
    beta_div <- guni$unifracs
  } else if (method == 'bray') {
    beta_div <- as.matrix(vegdist(relative_proportions, method = 'bray'))
  } else if (method == 'jaccard') {
    beta_div <- as.matrix(vegdist(relative_proportions, method = 'jaccard'))
  } else {
    stop("Invalid method. Please use 'Gen.Uni', 'bray' or 'jaccard'")
  }


  return(beta_div)
}


########################################################################################
# 1. Load Data
########################################################################################
# 1.1. Give directories
mzmine_featuredir <- 'raw_data/mzmine_feature.csv'
metadatadir <- "raw_data/metadata.csv"
canopus_formuladir <- "raw_data/canopus_formula_summary.tsv"
canopus_structuredir <- "raw_data/structure_identifications.tsv"
qemistry_dir <- 'raw_data/fingerprint_distance.csv'
cos.sim.dir <- 'raw_data/modified_cos_similarity.csv'

# 1.2. Initiate object
mmo <- GetMZmineFeature(mzmine_featuredir, metadatadir) # load mzmine featurelist and metadata
mmo <- AddSiriusAnnot(mmo, canopus_structuredir, canopus_formuladir) # Add annotation from sirius
mmo <- ReplaceZero(mmo) # Replace 0 and NA values by half of the minimum value
mmo <- MassNormalization(mmo) # Normalize peak area by sample mass in metadata
mmo <- MeancenterNormalization(mmo) # Add mean-centered area
mmo <- LogNormalization(mmo) # Add log-transformed area
mmo <- ZNormalization(mmo) # Add Zscore
qemistree <- read.csv(qemistry_dir, row.names = 1, header = TRUE, check.names = FALSE) # Load qemistree distance matrix
mmo$qemistree <- qemistree # Add qemistree distance matrix to mmo object
mmo <- GetCosSimToDissim(mmo, cos.sim.dir) # Add cosine dissimilarity matrix to mmo object
#some features are in the qemistree (i.e., sirius), but not in the mzmine feature list
siriused_ids <- rownames(mmo$qemistree)
siriused_features <- IDToFeature(mmo, siriused_ids)



# 1.3. Get annotations from sirius
features <- mmo$sirius_annot
# Get list of glucosinolated using sirius annotation
GLSs <- features %>% filter(str_detect(features[[46]], "Glucosinolate")) %>% pull(feature)
# Get list of flavonoids using sirius annotation
FLVs <- features %>% filter(str_detect(features[[46]], "Flavonoid")) %>% pull(feature)



########################################################################################
# 2. PLS-DA plot
#######################################################################################
# Define your custom colors for each group
custom_colors <- c("con" = "#999999", "sl1" = "#c89c7d", "px1" = "#55a5d3", "px2" = "#2575a4", 
                   "px3" = "#00446c", "sl2" = "#eeb696", "sl3" = "#bc6129", "le1" = "#a4d069", 
                   "le2" = "#618a29", "le3" = "#386200")
#Or automatically give colors
custom_colors <- setNames(brewer.pal(length(unique(mmo$metadata$group)), "Set3"), unique(mmo$metadata$group))

# Plot pls-da for all possible normalizations 
normalizations <- c('None', 'Log', 'Meancentered', 'Z')
topks <- c(0, 10)
for (normalization in normalizations){
  for (topk in topks){
    PLSDAplot(mmo, color = custom_colors, topk = topk, outdir = paste('plots/plsda/PLSDA_', normalization, '_top', topk, '.pdf'), normalization = normalization)
  }
}

# All features for all group
PLSDAplot(mmo, custom_colors, topk = 0, outdir = 'plots/plsda.pdf', normalization = 'Z', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL)
# Glucosinolates for onlt sl groups
PLSDAplot(mmo, custom_colors, topk = 0, outdir = 'plots/plsda_GLS_sl.pdf', normalization = 'Z', filter_feature = TRUE, feature_list = GLSs, filter_group = TRUE, group_list = c('con', 'sl1', 'sl2', 'sl3'))


########################################################################################
# 3. Pairwise comparison
########################################################################################

# 3.1. Add pairwise comparisons
mmo <- PairwiseComp(mmo, 'con', 'sl1')
mmo <- PairwiseComp(mmo, 'con', 'sl2')
mmo <- PairwiseComp(mmo, 'con', 'sl3')
mmo <- PairwiseComp(mmo, 'con', 'px1')
mmo <- PairwiseComp(mmo, 'con', 'px2')
mmo <- PairwiseComp(mmo, 'con', 'px3')
mmo <- PairwiseComp(mmo, 'con', 'le1')
mmo <- PairwiseComp(mmo, 'con', 'le2')
mmo <- PairwiseComp(mmo, 'con', 'le3')


# 3.2. Get lists of DAMs for each comparison
# Generate the list of comparisons automatically by looking up mmo$pairwise
comparison_columns <- colnames(mmo$pairwise)
log2FC_columns <- grep("log2FC", comparison_columns, value = TRUE)
comparisons <- unique(sub("log2FC", "", log2FC_columns))
comparisons <- sub("_$", "", comparisons)  # Remove trailing underscore from comparisons

# 3.2.1. FC > 2
# Make list of DAMs for up and downregulation for each comparison
DAMs_up <- list()
DAMs_down <- list()
for (comp in comparisons) {
  group1 <- strsplit(comp, "_vs_")[[1]][1]
  group2 <- strsplit(comp, "_vs_")[[1]][2]
  DAMs_up[[paste(comp, "up", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) > 1 & get(paste(comp, "padj", sep = "_")) < 0.1)$feature
  DAMs_down[[paste(comp, "down", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) < -1 & get(paste(comp, "padj", sep = "_")) < 0.1)$feature
}          
names(DAMs_up) <- paste(comparisons, "up", sep = ".")
names(DAMs_down) <- paste(comparisons, "down", sep = ".")


########################################################################################
########################################################################################
# 4. Upset
########################################################################################
# 4.1. Define groups

UpsetInput1 <- list(
  sl1.up = DAMs_up$con_vs_sl1.up,
  sl2.up = DAMs_up$con_vs_sl2.up,
  sl3.up = DAMs_up$con_vs_sl3.up,
  px1.up = DAMs_up$con_vs_px1.up,
  px2.up = DAMs_up$con_vs_px2.up,
  px3.up = DAMs_up$con_vs_px3.up,
  le1.up = DAMs_up$con_vs_le1.up,
  le2.up = DAMs_up$con_vs_le2.up,
  le3.up = DAMs_up$con_vs_le3.up
)
UpsetInput2 <- list(
  sl3.up = DAMs_up$con_vs_sl3.up,
  px3.up = DAMs_up$con_vs_px3.up,
  le3.up = DAMs_up$con_vs_le3.up
)

#4.2. Plot
# Upset
pdf("plots/upset/Upset_gall.pdf",7, 6)
upset(fromList(UpsetInput), nsets=10, nintersects=20,order.by='freq', mainbar.y.label='Features in Set', line.size=1, point.size=4, shade.color='white', text.scale=1, show.numbers=FALSE)
dev.off()
# Venn
ggvenn(UpsetInput2)

########################################################################################
# 5. Volcano Plot
#######################################################################################
for (comp in comparisons){
  VolcanoPlot(mmo, comp = comp, topk = 10, outdir = paste('plots/volcano_', comp, '.pdf', sep = ''))
}


########################################################################################
# 6. Canopus Enrichment analysis 
#######################################################################################
#6.1. Enrichment analysis for all terms
for (term in c('NPC', 'ClassyFire', 'all_terms')){
  CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, pthr = 0.1, terms = term, representation = 'greater', prefix = paste('plots/enrichment/All_levels_', term, '_degree_merged'))
}
#6.2. Enrichment analysis for each term level
term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
for (term_level in term_levels){
  CanopusLevelEnrichmentPlot(mmo, DAMs_up, term_level = term_level, pthr = 0.1, prefix = paste('plots/enrichment/DAM_merged_upreg_', term_level, sep = ''))
}



########################################################################################
# 7. QEMISTREE and heatmap
#######################################################################################
#Alter this part for heatmap
filter_feature <- FALSE
feature_list <- GLSs # list of features to show in heatmap
filter_group <- FALSE
group_list <- c('con', 'sl1', 'sl2', 'sl3') # list of groups to show in heatmap
normalization <- 'Z' #None, Log, Meancentered, Z
distance <- 'cosine' #cosine or qemistree or combined

# 7.1. Process data for heatmap
# 7.1.1. Get summarized data (group mean or FC)
if (filter_group){
  group_means <- GetGroupMeans(mmo, normalization = normalization, filter = TRUE, filter_groups = group_list)
} else {
  group_means <- GetGroupMeans(mmo, normalization = normalization)
}
heatmap_data <- group_means

head(heatmap_data)
# 7.1.2. Filter features
# Determine distance metric
if (distance == 'qemistree'){
  distance_matrix <- mmo$qemistree
} else if (distance == 'cosine'){
  distance_matrix <- mmo$cos.dissim
} else if (distance == 'combined'){
  distance_matrix <- mmo$qemistree * mmo$cos.dissim
}
heatmap_data <- heatmap_data %>% filter(id %in% rownames(distance_matrix)) # get features with fingerprints
# make matrix for heatmap
FC_matrix <- as.matrix(heatmap_data[,-1])
rownames(FC_matrix) <- heatmap_data$id
# Reorder the rows of distance_matrix to match the order of FC_matrix_
distance_matrix <- distance_matrix[rownames(FC_matrix), rownames(FC_matrix)]
dist_matrix <- as.dist(distance_matrix)



if (filter_feature){
  filter_list <- feature_list #GLSs
  filter_id <- FeatureToID(mmo, filter_list)
  filter_id <- filter_id[filter_id %in% rownames(distance_matrix)] # remove custom-annotated but not siriused features
  filter_distance <- mmo$qemistree[filter_id, filter_id]
  heatmap_data <- heatmap_data %>% filter(id %in% filter_id)


  heatmap_data$con <- NULL # as we are looking at FC, all con are 0
  # make matrix for heatmap
  FC_matrix <- as.matrix(heatmap_data[,-1])
  rownames(FC_matrix) <- heatmap_data$id


  # Reorder the rows of distance_matrix to match the order of FC_matrix_
  filter_distance <- filter_distance[rownames(FC_matrix), rownames(FC_matrix)]
  dist_matrix <- as.dist(filter_distance)
  #Label custm-annotated features
  row_label <- rownames(FC_matrix)
  for (i in 1:length(rownames(FC_matrix))){
    id <- rownames(FC_matrix)[i]
    custom_annot <- mmo$custom_annot$custom_annot[mmo$custom_annot$id == id]
    if (length(custom_annot[[1]]) > 0) {
      row_label[i] <- custom_annot[[1]]
    }
  }
}

# 7.2.1.2. Get annotation for Qemistree
#SIRIUS superclass and specific classes
sirius_annot <- mmo$sirius_annot
# Get NPC Annotations
NPC_pathway <- unique(sirius_annot[[30]])[-1] #remove NA
NPC_superclass <- unique(sirius_annot[[32]])[-1] #remove NA
NPC_class <- unique(sirius_annot[[34]])[-1] #remove NA

sirius_annot_qemistry <- sirius_annot %>%
  select(id = 1, NPC_class = 34, NPC_superclass = 32, NPC_pathway = 30) %>%
  filter(id %in% rownames(distance_matrix)) # get features with fingerprints

rownames(sirius_annot_qemistry) <- sirius_annot_qemistry$id
sirius_annot_qemistry$id <- NULL

# 7.2.1.3. Set colors to show in heatmap
ann_colors = list(
    NPC_pathway = c('Alkaloids' = 'purple', 'Carbohydrates' = 'green', 'Polyketides' = 'yellow', 'Terpenoids' = 'orange', 
        'Amino acids and Peptides' = 'blue', 'Shikimates and Phenylpropanoids' = 'red', 'Fatty acids' = 'brown'),
    NPC_class = setNames(viridis(length(NPC_class)), NPC_class),
    NPC_superclass = setNames(viridis(length(NPC_superclass)), NPC_superclass)
    )

##### Specific colors for specific classes of interest
ann_colors$NPC_class['Glucosinolates'] <- 'red'

# 7.2.1.4. Plot
pdf("plots/test.pdf", width = 20, height = 40)
pheatmap(mat = FC_matrix, 
     cluster_rows = TRUE, 
     clustering_distance_rows = dist_matrix, 
     cluster_cols = TRUE, 
     clustering_method = "average", #UPGMA
     show_rownames = TRUE, 
     show_colnames = TRUE,
     annotation_row = sirius_annot_qemistry, # Dataframe with the rownames are identical with 'mat' and gives annotation
     annotation_colors = ann_colors,
     cellwidth = 25,
     cellheight = 0.5,
     treeheight_row = 100,
     fontsize_row = 5,
     fontsize_col = 15,
     scale = 'none',
     annotation_names_row = TRUE,
     color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

########################################################################################
# 8. Chemical Diversity
#######################################################################################
# 8.1. Alpha diversity
# 8.1.1. unweighted Hill
unweighted_q1 <- GetHillNumbers(mmo, q = 1)
unweighted_q2 <- GetHillNumbers(mmo, q = 2)
# 8.1.2. Functional Hill number without Rao's Q
functional.hill.q1 <- GetFunctionalHillNumber(mmo, normalization = 'None', q = 1)
functional.hill.q2 <- GetFunctionalHillNumber(mmo, normalization = 'None', q = 2)
# 8.1.3. Functional Hill number with Rao's Q standardization
functional.hill.q1.raoQ <- GetFunctionalHillNumberwithRaoQ(mmo, normalization = 'None', q = 1)
functional.hill.q2.raoQ <- GetFunctionalHillNumberwithRaoQ(mmo, normalization = 'None', q = 2)
#Plot
ggplot(functional.hill.q1, aes(x = group, y = hill_number)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 0.5) +
  theme_classic() 
ggsave('plots/diversity/functional_q1.png', width = 8,height = 8)

# 8.2. Beta diversity
# Generalized UniFrac distance using fingerprint distance. The alpha value can be adjusted
feature <- mmo$zscore
metadata <- mmo$metadata
guni <- GetBetaDiversity(mmo, method = 'Gen.Uni', normalization = 'Log')
guni.0 <- guni[,,'d_0'] # GUniFrac with alpha 0
guni.05 <- guni[,,'d_0.5'] # GUniFrac with alpha 0.5
guni.1 <- guni[,,'d_1'] # Weighted UniFrac

# Perform NMDS
nmds <- metaMDS(guni.1, k = 2, try = 50, trymax = 100)
nmds_coords <- as.data.frame(scores(nmds))
groups <- c()
for (col in colnames(feature)[-c(1, 2)]) {
  groups <- append(groups, metadata[metadata$sample == col, ]$group)
}
nmds_coords$group <- groups
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = group)) +
      geom_point(size = 3) +
      #geom_text_repel(aes(label = group), size = 3) +
      theme_classic() +
      stat_ellipse(level = 0.90) +
      labs(x = "NMDS1", y = "NMDS2") +
      theme(legend.position = "right")
    
ggsave(paste('plots/nmds.pdf'), height = 6, width = 6)

########################################################################################
# 9. Performance regression
#######################################################################################
# fetch pairwise comparisons from mmo$pairwise
comparisons <- c('con_vs_sl1', 'con_vs_sl2', 'con_vs_sl3', 'con_vs_px1', 'con_vs_px2', 'con_vs_px3', 'con_vs_le1', 'con_vs_le2', 'con_vs_le3')
# Define the DAMs
sl.up <- Reduce(union, c(
  DAMs_up$con_vs_sl1.up, 
  DAMs_up$con_vs_sl2.up, 
  DAMs_up$con_vs_sl3.up
))
sl.down <- Reduce(union, c(
  DAMs_down$con_vs_sl1.down, 
  DAMs_down$con_vs_sl2.down, 
  DAMs_down$con_vs_sl3.down
))

# Fit LMM and find significant features
sl.lmm <- GetPerformanceFeatureLMM(mmo, herbivore = 'Sl', groups = c('sl1', 'sl2', 'sl3'), DAM.list = list(sl.up = sl.up, sl.down = sl.down), comparisons)
sl.lmm.sig <- sl.lmm %>% filter(p_value < 0.05)
sl.lmm.sig.neg <- sl.lmm.sig %>% filter(effect.size < 0)
# Plot the association between fold change and effect size
PlotFoldchangeResistanceRegression(performance_regression = sl.lmm.sig, 
  fold_change = 'con_vs_sl3_log2FC', 
  color = c('sl.up' = '#ba3b3c', 'sl.down' = '#1a3f9e', 'else' = 'grey'), 
  output_dir = 'plots/sl_lmm.png')
