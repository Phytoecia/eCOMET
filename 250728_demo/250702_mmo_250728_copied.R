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
library('ggvenn')
library('multcomp')
library('ggsignif')
library('readr')
library('multcompView')
select <- dplyr::select 
################################################
setwd('/home/minsoo/software/mmo/250728_demo')



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
  print('MMO object created.')
  print(paste0('Feature number: ', nrow(mmo$feature_data)))
  print(paste0(nrow(mmo$metadata), ' samples in ', length(unique(mmo$metadata$group)), ' groups'))
  return(mmo)
}

# Add annotation from SIRIUS to the mmo object
AddSiriusAnnot <- function(mmo, canopus_structuredir, canopus_formuladir){
  structure_identifications <- read_tsv(canopus_structuredir, show_col_types = FALSE)
  structure_identifications$mappingFeatureId <- gsub(" ", "", structure_identifications$mappingFeatureId)
  canopus_formula_summary <- read_tsv(canopus_formuladir, show_col_types = FALSE)
  canopus_formula_summary$mappingFeatureId <- gsub(" ", "", canopus_formula_summary$mappingFeatureId)
  siriused_ids <- unique(union(structure_identifications$mappingFeatureId, canopus_formula_summary$mappingFeatureId))
  sirius_df <- mmo$feature_data %>% select(id, feature)
  sirius_df <- sirius_df %>%
  left_join(structure_identifications, by = c("id" = "mappingFeatureId")) %>%
  left_join(canopus_formula_summary, by = c("id" = "mappingFeatureId"))
  mmo$sirius_annot <- sirius_df
  print('SIRIUS annotation added to mmo$sirius_annot')
  return(mmo)
}
AddCustomAnnot <- function(mmo, DB_file, mztol = 5, rttol = 0.5) {
  DB <- read.csv(DB_file)
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
  print(paste('Custom annotation added to mmo$custom_annot using ', DB_file))
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
  print(paste('Missing values were filled with', method))
  return(mmo)
}
# Use mass in the metadata file to normalize the peak area
MassNormalization <- function(mmo){
  normalized_df <- mmo$feature_data
  metadata <- mmo$metadata
  mean_mass <- mean(mmo$metadata$mass)
  for (sample_col in colnames(mmo$feature_data)[-c(1,2)]) {
    sample_metadata <- metadata[metadata$sample == sample_col, ]
    mass <- sample_metadata$mass
    normalized_df[[sample_col]] <- as.numeric(mmo$feature_data[[sample_col]])*mean_mass/mass
  }
  mmo$feature_data <- normalized_df
  print("Peak area are normalized by sample mass")
  return(mmo)
}
# Make a df of log-transformed peak area in the mmo object (mmo$log)
LogNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  log_data <- log2(feature_data_only)
  log_df <- cbind(mmo$feature_data[, 1:2], log_data)
  mmo$log <- log_df
  print('Log-normalized values were added to mmo$log')
  return(mmo)
}
# Make a df of mean-centered peak area in the mmo object (mmo$meancentered)
MeancenterNormalization <- function(mmo){
  feature_data_only <- mmo$feature_data[,-(1:2)]
  mean_centered_data <- t(apply(feature_data_only, 1, function(x) x - mean(x, na.rm = TRUE)))
  mean_centered_df <- cbind(mmo$feature_data[, 1:2], mean_centered_data)
  mmo$meancentered <- mean_centered_df
  print('Meancentered values were added to mmo$meancentered')
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
  print('Z-score values were added to mmo$zscore')
  return(mmo)
}

# Add cosine, DREAMS and MS2DeepScore dissimilarity matrices to the mmo object
AddChemDist <- function(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL) {
  add_dissim_matrix <- function(mmo, sim_dir, slot_name) {
    sim <- read.csv(sim_dir, col.names = c('cluster1', 'cluster2', 'metric', 'similarity', 'etc'))
    sim$dissimilarity <- 1 - sim$similarity
    clusters <- unique(c(sim$cluster1, sim$cluster2))
    dissim_mat <- matrix(1, nrow = length(clusters), ncol = length(clusters))
    rownames(dissim_mat) <- clusters
    colnames(dissim_mat) <- clusters
    diag(dissim_mat) <- 0  # Set diagonal to 0
    for (i in 1:nrow(sim)) {
      c1 <- as.character(sim$cluster1[i])
      c2 <- as.character(sim$cluster2[i])
      dissim <- sim$dissimilarity[i]
      dissim_mat[c1, c2] <- dissim
      dissim_mat[c2, c1] <- dissim
    }
    mmo[[slot_name]] <- dissim_mat
    print(paste(slot_name, "added to mmo"))
    return(mmo)
  }
  if (!is.null(cos_dir)) mmo <- add_dissim_matrix(mmo, cos_dir, "cos.dissim")
  if (!is.null(dreams_dir)) mmo <- add_dissim_matrix(mmo, dreams_dir, "dreams.dissim")
  if (!is.null(m2ds_dir)) mmo <- add_dissim_matrix(mmo, m2ds_dir, "m2ds.dissim")
  if (is.null(cos_dir) && is.null(dreams_dir) && is.null(m2ds_dir)) {
    stop("Please provide at least one valid directory.")
  }
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
  print('Cosine similarity were transformed to dissimilarity and added to mmo$cos.dissim')
  return(mmo)
}
# Add DREAMS dissimilarity matrix to the mmo object
AddDreamsDissim <- function(mmo, dreams_dir) {
  # Load DREAMS similarity matrix
  dreams_sim <- read.csv(dreams_dir, row.names = 1, header = TRUE, check.names = FALSE)
  # Convert similarity to dissimilarity
  dreams_dissim <- 1 - dreams_sim
  # Add DREAMS dissimilarity matrix to mmo object
  mmo$dreams_dissim <- dreams_dissim
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
  } else {
    print('The normalization should be None, Log, Meancentered, or Z')
  }
  return(feature)
}
GetDistanceMat <- function(mmo, distance = 'dreams'){
  if (distance == 'dreams'){
    distance_matrix <- mmo$dreams.dissim
  } else if (distance == 'cosine'){
    distance_matrix <- mmo$cos.dissim
  } else if (distance == 'm2ds'){
    distance_matrix <- mmo$m2ds.dissim
  }
  return(distance_matrix)
}
# Function to get list of IDs corresponding to a list of feature names
FeatureToID <- function(mmo, feature_names) {
  feature_data <- mmo$feature_data
  feature_ids <- feature_data %>%
    filter(feature %in% feature_names) %>%
    select(feature, id)
  # Match the order of feature_names
  feature_ids <- feature_ids$id[match(feature_names, feature_ids$feature)]
  return(feature_ids)
}
# Function to convert ID to feature
IDToFeature <- function(mmo, feature_ids) {
  feature_data <- mmo$feature_data
  feature_names <- feature_data %>%
    filter(id %in% feature_ids) %>%
    select(id, feature)
  # Match the order of feature_ids
  feature_names <- feature_names$feature[match(feature_ids, feature_names$id)]
  return(feature_names)
}

# returns group-averaged feature dataframe 
GetGroupMeans <- function(mmo, normalization = 'None', filter_feature = FALSE, feature_list = NULL, filter_group = FALSE, group_list = NULL) {
  feature_data <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  
  # Melt the feature data to long format
  long_feature_data <- melt(feature_data, id.vars = c('id', 'feature'), variable.name = 'sample', value.name = 'feature_value')
  colnames(long_feature_data) <- c('id', 'feature', 'sample', 'feature_value')
  
  # Merge with metadata to get group information
  merged_data <- merge(long_feature_data, metadata[, c('sample', 'group')], by = 'sample')
  if (filter_group == TRUE){
    merged_data <- merged_data %>% filter(group %in% group_list)
  }
  # Calculate group means
  group_means <- merged_data %>%
    group_by(group, id) %>%
    summarise(mean_value = mean(feature_value, na.rm = TRUE)) %>%
    spread(key = group, value = mean_value)
  if (filter_feature == TRUE){
    group_means <- group_means %>% filter(id %in% FeatureToID(mmo, feature_list))
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

res <- DunnettTest(hill_number ~ group, data = alphadiv)

anova_tukey_dunnett <- function(df, formula) {
  aov_res <- aov(as.formula(formula), data = df)
  tukey_res <- TukeyHSD(aov_res)
  tukey_sig <- multcompLetters4(aov_res, tukey_res)
  dunnett_res <- DunnettTest(as.formula(formula), data = df)
  return(list(aov_res = aov_res, tukey_res = tukey_res, tukey_sig = tukey_sig, dunnett_res = dunnett_res))
}
write_anova <- function(anova_data, outdir, way='oneway'){
  way_num <- switch(way, oneway = 1, twoway = 3)
  # Perform ANOVA and Tukey HSD
  aov_res <- anova_data$aov_res
  tukey_res <- anova_data$tukey_res
  tukey_sig <- anova_data$tukey_sig[way_num]
  dunnett_res <- anova_data$dunnett_res

  # Save ANOVA and Tukey HSD results
  anova_df <- as.data.frame(summary(aov_res)[[1]])
  anova_df$Comparison <- rownames(anova_df)
  tukey_df <- as.data.frame(tukey_res[way_num])
  tukey_df$Comparison <- rownames(tukey_df)
  sig_letter <- as.data.frame(unlist(tukey_sig))
  sig_letter$Comparison <- rownames(sig_letter)
  dunnett_df <- as.data.frame(dunnett_res[[1]])
  dunnett_df$comp <- rownames(dunnett_df)
  
  # Create a combined results data frame
  combined_df <- bind_rows(
    tibble(Test = "ANOVA", anova_df),
    tibble(Test = "Tukey", tukey_df),
    tibble(Test = 'sig', sig_letter),
    tibble(Test = 'Dunnett', dunnett_df)
  )
  write_csv(combined_df, file = outdir)
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

# Get lists of DAMs (Differentially Accumulated Metabolites) for each comparison
# mmo : mmo object with pairwise comparison matrix
# fc_cutoff : the threshold of log2FC to be considered significant
# pval_cutoff : the threshold of adjusted p-value to be considered significant
GetDAMs <- function(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.1) {
  # Generate the list of comparisons automatically by looking up mmo$pairwise
  comparison_columns <- colnames(mmo$pairwise)
  log2FC_columns <- grep("log2FC", comparison_columns, value = TRUE)
  comparisons <- unique(sub("log2FC", "", log2FC_columns))
  comparisons <- sub("_$", "", comparisons)  # Remove trailing underscore from comparisons
  # Make list of DAMs for up and downregulation for each comparison
  DAMs_up <- list()
  DAMs_down <- list()
  for (comp in comparisons) {
    group1 <- strsplit(comp, "_vs_")[[1]][1]
    group2 <- strsplit(comp, "_vs_")[[1]][2]
    DAMs_up[[paste(comp, "up", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) > 0.5849625 & get(paste(comp, "padj", sep = "_")) < 0.1)$feature
    DAMs_down[[paste(comp, "down", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) < -0.5849625 & get(paste(comp, "padj", sep = "_")) < 0.1)$feature
  }          
  names(DAMs_up) <- paste(comparisons, "up", sep = ".")
  names(DAMs_down) <- paste(comparisons, "down", sep = ".")
  return(list(DAMs_up = DAMs_up, DAMs_down = DAMs_down))
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
  if (topk > 0){
  loading_scale <- max(abs(scores))/(4*max(abs(top_features$Comp1_Loading)))
  }
  
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
               aes(x = 0, y = 0, xend = Comp1_Loading * loading_scale, yend = Comp2_Loading * loading_scale),  # Scale the arrows
               arrow = arrow(length = unit(0.3, "cm")), color = "grey", linewidth = 1) +
  # Add labels for the top 10 features
  geom_text_repel(data = top_features,
            aes(x = Comp1_Loading * loading_scale, y = Comp2_Loading * loading_scale, label = Feature),
            color = "black", vjust = 1.5, size = 3)

  ggsave(outdir, height = 6, width = 6)
  write.csv(loadings_df, 'PLSDA_loadings.csv')
  print(paste(normalization, '-normalized feature was used'))
}

# Generate heatmap inputs
# mmo : mmo object with sirius annotation and normalized
# filter_feature : whether to filter features by feature_list
# feature_list : a vector of feature names to filter
# filter_group : whether to filter groups by group_list
# group_list : a vector of group names to filter
# summarize : 'fold_change' or 'mean'. If 'fold_change', control_group is used to calculate fold change
# control_group : the group to use as control for fold change calculation
# normalization : 'None', 'Log', 'Meancentered', or 'Z'
# distance : 'dreams', 'cosine', or 'm2ds'
GenerateHeatmapInputs <- function(mmo, filter_feature = FALSE, feature_list = NULL, 
                                filter_group = FALSE, group_list = NULL, 
                                summarize = 'fold_change', control_group = 'ctrl', 
                                normalization = 'None', distance = 'dreams') {
  # 12.1.1. Get summarized data (group mean or FC)
  if (filter_group){
    group_means <- GetGroupMeans(mmo, normalization = normalization, filter_group = TRUE, group_list = group_list)
  } else {
    group_means <- GetGroupMeans(mmo, normalization = normalization)
  }
  if (summarize == 'fold_change'){
    fold_change <- GetLog2FoldChange(group_means, control_group = control_group)
    heatmap_data <- fold_change
    heatmap_data[[control_group]] <- NULL
  } else if(summarize == 'mean'){ 
    heatmap_data <- group_means
  }
  # 12.1.2. Filter features
  # Determine distance metric
  distance_matrix <- GetDistanceMat(mmo, distance = distance)
  heatmap_data <- heatmap_data %>% filter(id %in% rownames(distance_matrix)) # remove features not in distance matrix

  # make matrix for heatmap
  FC_matrix <- as.matrix(heatmap_data[,-1])
  rownames(FC_matrix) <- heatmap_data$id
  # Reorder the rows of distance_matrix to match the order of FC_matrix_
  distance_matrix <- distance_matrix[rownames(FC_matrix), rownames(FC_matrix)]
  dist_matrix <- as.dist(distance_matrix)

  row_label <- rownames(FC_matrix)
  if (filter_feature){
    filter_list <- feature_list 
    filter_id <- FeatureToID(mmo, filter_list)
    filter_id <- filter_id[filter_id %in% rownames(distance_matrix)] # remove custom-annotated but not in the distance matrix
    filter_distance <- distance_matrix[filter_id, filter_id]
    heatmap_data <- heatmap_data %>% filter(id %in% filter_id)

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
  return(list(FC_matrix = FC_matrix, dist_matrix = dist_matrix, row_label = row_label, heatmap_data = heatmap_data))
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
  # print(paste('total features:', nrow(all_feature), 'list_test features:', nrow(subset_feature)))
  # Select the appropriate term level for enrichment analysis
  if (term_level == "NPC_pathway") {
    all_feature$classifications_split <- all_feature[[32]]
    subset_feature$classifications_split <- subset_feature[[32]]    
  } else if (term_level == "NPC_superclass") {
    all_feature$classifications_split <- all_feature[[34]]
    subset_feature$classifications_split <- subset_feature[[34]]
  } else if (term_level == "NPC_class") {
    all_feature$classifications_split <- all_feature[[36]]
    subset_feature$classifications_split <- subset_feature[[36]]
  } else if (term_level == "ClassyFire_superclass") {
    all_feature$classifications_split <- all_feature[[38]]
    subset_feature$classifications_split <- subset_feature[[38]]
  } else if (term_level == "ClassyFire_class") {
    all_feature$classifications_split <- all_feature[[40]]
    subset_feature$classifications_split <- subset_feature[[40]]
  } else if (term_level == "ClassyFire_subclass") {
    all_feature$classifications_split <- all_feature[[42]]
    subset_feature$classifications_split <- subset_feature[[42]]
  } else if (term_level == "ClassyFire_level5") {
    all_feature$classifications_split <- all_feature[[44]]
    subset_feature$classifications_split <- subset_feature[[44]]
  } else if (term_level == "ClassyFire_most_specific") {
    all_feature$classifications_split <- all_feature[[46]]
    subset_feature$classifications_split <- subset_feature[[46]]
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
    foldenrichment = (as.numeric(subset_term_counts[names(enrichment_results)]) / length(subset_feature))/(as.numeric(total_term_counts[names(enrichment_results)]) / nrow(all_feature)),
    pval = enrichment_results,
    fdr = adjusted_pvalues
  )
  results <- results %>% filter(term != 'None')
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


CanopusListEnrichmentPlot <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus <- sig.canopus %>% arrange(desc(foldenrichment))
  ggplot(sig.canopus, aes(x = foldenrichment, y = reorder(term, foldenrichment), color = -log(pval), size = subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')+
    xlim(0,max(sig.canopus$foldenrichment+1))
    #facet_wrap(~term_level, ncol = 1, scales = 'free_y', strip.position = 'right', shrink = TRUE)
    
  ggsave(outdir, height = height, width = width)
}
CanopusListEnrichmentPlot_2 <- function(mmo, feature_list, pthr = 0.05, outdir, height = 5, width = 5, topn = 5){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, feature_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater'))
  }
  sig.canopus$term <- paste(sig.canopus$term, ';', sig.canopus$term_level)
  sig.canopus <- sig.canopus %>% slice_max(order_by = -pval, n = topn)
  sig.canopus <- sig.canopus %>% arrange(desc(foldenrichment))
  ggplot(sig.canopus, aes(x = foldenrichment, y = reorder(term, foldenrichment), color = -log(pval), size = subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    xlim(0,max(sig.canopus$foldenrichment+1))+
    ylab('Chemical Class')
    
  ggsave(outdir, height = height, width = width)
}

CanopusLevelEnrichmentPlot <- function(mmo = mmo, comp.list, term_level = 'NPC_pathway',pthr = 0.1, representation = 'greater', prefix = 'enrichment', height = 5, width = 5){
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
  ggsave(paste(prefix, '.pdf'), width = width, height = height)
}
##
# Generate plot and result files for the list of pairwise comparisons
CanopusAllLevelEnrichmentPlot <- function(mmo = mmo, comp.list, terms = 'all_terms', term_levels = NULL, pthr = 0.1, representation = 'greater', prefix = 'enrichment', height = 10, width = 8){
  df.EA <- data.frame()
  sig.terms <- c()
  if(terms == 'all_terms'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'NPC'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway')
  } else if (terms == 'ClassyFire'){
    term_levels = c('ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'custom'){
    term_levels = term_levels
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
    ylab('Chemical classes')+
    facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')
  enrichment_plot
  write.csv(df.EA, paste0(prefix, '.csv'), row.names = FALSE)
  write.csv(df.EA.sig, paste0(prefix, '_sig.csv'), row.names = FALSE)
  ggsave(paste0(prefix, '.pdf'), width = width, height = height)
}

# Function to get performance individual feature regression
# mmo : mmo object
# target : the name of feature 
# phenotype : the name of phenotype performance in metadata
# groups : the group from metadata containing performance data
# normalization : 'None', 'Log', 'Meancentered', 'Z'
# model is lmm or lm
# output : output file name for the plot
FeaturePerformanceRegression <- function(mmo, target, phenotype, groups, model = 'lmm', normalization = 'Z', output){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # Get phenotype performance from the metadata, get the feature value from the feature matrix, then combine
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) %>% filter(group %in% groups)
  feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[feature$feature == target, -(1:2)]))
  combined_df <- merge(phenotype.df, feature_df, by='sample')
  
  # Perform linear mixed model or simple linear regression
  if (model == 'lmm'){
    fit <- lmer(combined_df$performance ~ combined_df$feature_value + (1|combined_df$group))
    p_value <- summary(fit)$coefficients[2, 5]
  } else if (model == 'lm'){
    fit <- lm(combined_df$performance ~ combined_df$feature_value)
    p_value <- summary(fit)$coefficients[2, 4]
  } else if (model == 'pearson'){
    pearson <- cor.test(combined_df$performance, combined_df$feature_value)
    p_value <- pearson[[3]]
  } else {
    stop("Invalid model type. Please use 'lmm' or 'lm' or 'pearson")
  }

  # Plot the fit using ggplot
  ggplot(combined_df, aes(x = feature_value, y = performance, color = group)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE) +
    theme_classic() +
    labs(title = paste("Regression of", target, "against", phenotype, "performance"),
         x = "Feature Value",
         y = "Performance") +
    theme(legend.position = "right") +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 4)), 
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
  ggsave(output, height = 6, width = 6)
}


# Function to get performance feature regression
# From the metadata, the phenotype performance is fetched
# From the feature data, the performance is regressed against each feature
# mmo : mmo object
# phenotype : the name of phenotype performance in metadata
# groups : the group from metadata containing performance data
# DAM.list : a list of DAMs
# comparisons : a list of pairwise comparisons
# normalization : 'None', 'Log', 'Meancentered', 'Z'
GetPerformanceFeatureRegression <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'None'){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata
  
  # phenotype.sample <- metadata %>% filter(group %in% groups) %>% pull(sample)
  # phenotype.area <- feature %>% select(id, feature, all_of(phenotype.sample))

  performance.linreg <- data.frame(pval = double(), effect.size = double())
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) %>% filter(group %in% groups)

  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)  
  for (i in 1:nrow(feature)) {
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')

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

GetPerformanceFeatureLMM <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'Z'){
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

  # get phenotype performance data
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) %>% filter(group %in% groups)
  #create an empty dataframe to store regression results
  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)  
  # iterate regression analysis
  for (i in 1:nrow(feature)) {
    # for each feature, generate phenotype performance X feature value data
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')

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
GetPerformanceFeatureCorrelation <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'None'){
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata
  
  # phenotype.sample <- metadata %>% filter(group %in% groups) %>% pull(sample)
  # phenotype.area <- feature %>% select(id, feature, all_of(phenotype.sample))

  performance.linreg <- data.frame(pval = double(), effect.size = double())
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, performance = metadata[,phenotype]) %>% filter(group %in% groups)

  regression_results <- data.frame(feature = character(), effect.size = numeric(), p_value = numeric(), is.Spec = logical(), stringsAsFactors = FALSE)  
  for (i in 1:nrow(feature)) {
    feature_name <- feature$feature[i]
    feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[i, -(1:2)]))
    combined_df <- merge(phenotype.df, feature_df, by='sample')
    pearson <- cor.test(combined_df$performance, combined_df$feature_value)
    pval <- pearson[[3]]
    cor <- pearson[[4]]
    tag <- "else"
    for (list_name in names(DAM.list)) {
      if (feature_name %in% DAM.list[[list_name]]) {
        tag <- list_name
      }
    }
    #is.Spec <- feature_name %in% target

    regression_results <- rbind(regression_results, data.frame(
      feature = feature_name, effect.size = cor, p_value = pval, tag = tag
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
get_quadrant <- function(x, y) {
  if (x > 0 & y > 0) return("Q1")
  if (x < 0 & y > 0) return("Q2")
  if (x < 0 & y < 0) return("Q3")
  if (x > 0 & y < 0) return("Q4")
  return("Edge")  # For points on axes
}

PlotFoldchangeResistanceQuad <- function(performance_regression, fold_change, color, output_dir){
  performance_regression <- performance_regression %>%
  mutate(
    quadrant = case_when(
      -effect.size > 0 & !!sym(fold_change) > 0 ~ "Q1",  
      -effect.size < 0 & !!sym(fold_change) < 0 ~ "Q3",  
      -effect.size < 0 & !!sym(fold_change) > 0 ~ "Q2",  
      -effect.size > 0 & !!sym(fold_change) < 0 ~ "Q4",  
      TRUE ~ "Edge"             # For points on axes
    )
  )
  q_counts <- table(performance_regression$quadrant)
  q13 <- sum(q_counts[c("Q1", "Q3")], na.rm = TRUE)
  q24 <- sum(q_counts[c("Q2", "Q4")], na.rm = TRUE)
  binom_test <- binom.test(q13, q13+q24, p = 0.5, alternative = "two.sided")

  ggplot(performance_regression, aes(x = -effect.size, y = !!sym(fold_change))) +
    geom_point(size = 0.5, aes(color = tag)) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(binom_test[[3]], 500)), 
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  ggsave(output_dir, height = 6, width = 6)
}
################### Singlevariate analyses ###################
# Function to plot bar plots and perform ANOVA and Tukey's HSD
AnovaBarPlot <- function(mmo, ID_list, outdir, normalization = 'None', filter_group = FALSE, group_list = NULL) {
  # Extract metadata and feature data
  metadata <- mmo$metadata
  feature_data <- GetNormFeature(mmo, normalization)
  
  # Iterate through each feature ID
  for (target_id in ID_list) {
    # Extract feature values and merge with metadata
    feature_values <- feature_data %>%
      filter(id == target_id) %>%
      select(-id, -feature) %>%
      t() %>%
      as.data.frame()
    colnames(feature_values) <- "value"
    feature_values$sample <- rownames(feature_values)
    feature_values <- merge(feature_values, metadata, by = "sample")
    if (filter_group == TRUE){
      feature_values <- feature_values %>% filter(group %in% group_list)
    }
    # Perform ANOVA
    anova <- anova_tukey_dunnett(feature_values, 'value ~ group')
    write_anova(anova, outdir = paste0(outdir,'/', target_id, '_anova.csv'), way = 'oneway')

    
    # Generate bar plot
    p <- ggplot(feature_values, aes(x = group, y = value, fill = group)) +
      geom_bar(stat = "summary", fun = "mean", position = "dodge") +
      geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.9), width = 0.2) +
      geom_beeswarm() +
      theme_classic() +
      labs(title = paste("Feature:", target_id), x = "Group", y = "Value") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    # Save the plot
    ggsave(file.path(outdir, paste0(target_id, "_barplot.png")), plot = p, width = 6, height = 4)
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

GetFunctionalHillNumber <- function(mmo, normalization = 'None',q = 1, distance = 'dreams', filter_feature = FALSE, feature_list = NULL){
  feature <- GetNormFeature(mmo, normalization = normalization)
  metadata <- mmo$metadata
  distance_matrix <- GetDistanceMat(mmo, distance = distance)
  # Scale the  distance matrix to be between 0 and 1
  
  if (filter_feature == TRUE){
    id_list <- FeatureToID(mmo, feature_list)
    distance_matrix <- distance_matrix[id_list, id_list] 
  }
  scaled_dissimilarity <- distance_matrix / max(distance_matrix)
  # Calculate the relative proportions of each feature and reorder them to match the order of the distance matrix
  q.feature <- feature %>% filter(id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ]
  scaled_dissimilarity <- as.matrix(scaled_dissimilarity)
  # Calculate Hill
  functional_hill_number <- c()
  if (q == 1){
    mask <- relative_proportions > 0
    Plog <- ifelse(mask, relative_proportions * log(relative_proportions), 0)
    DP <- scaled_dissimilarity %*% relative_proportions
    vals <- 2 * colSums(Plog * DP)
    functional_hill_number <- exp(-vals)    
  } else {
    Pq <- relative_proportions^q
    DPq <- scaled_dissimilarity %*% Pq
    vals <- colSums(Pq*DPq)
    functional_hill_number <- vals^(1/(1-q))
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
GetHillNumbers <- function(mmo, normalization = 'None', q = 0, filter_feature = FALSE, feature_list = NULL) {
  feature <- GetNormFeature(mmo, normalization = normalization)
  if (filter_feature == TRUE) {
    feature <- feature %>% filter(feature %in% feature_list)
  }
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


GetAlphaDiversity <- function(mmo, q = 1, normalization = 'None', mode = 'weighted', distance = 'dreams', filter_feature = FALSE, feature_list = NULL){
  if (mode == 'weighted'){
    GetFunctionalHillNumber(mmo, normalization = normalization, q = q, distance = distance, filter_feature = filter_feature, feature_list = feature_list)
  } else if (mode == 'unweighted'){
    GetHillNumbers(mmo, normalization = normalization, q = q, filter_feature = filter_feature, feature_list = feature_list)
  } else{
    print('mode should be weighted or unweighted')
  }
}

GetBetaDiversity <- function(mmo, method = 'Gen.Uni', normalization = 'None', distance = 'dreams', filter_feature = FALSE, feature_list = NULL) {    
  # Get compound distance and build tree for UniFrac
  scaled_dissimilarity <- GetDistanceMat(mmo, distance = distance) / max(GetDistanceMat(mmo, distance = distance))
  if (filter_feature == TRUE) {
    id_list <- FeatureToID(mmo, feature_list)
    scaled_dissimilarity <- scaled_dissimilarity[id_list, id_list]
  }
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
  } else if (method == 'CSCS') {
    CSS <- 1-GetDistanceMat(mmo, distance = distance)
    diag(CSS)
    q.feature <- GetNormFeature(mmo, normalization = normalization) %>% filter(id %in% colnames(CSS))
    relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
    CSCS_all <- t(relative_proportions) %*% CSS %*% relative_proportions

    sample_names <- colnames(relative_proportions)
    n_samples <- length(sample_names)
    CSCS_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
    rownames(CSCS_matrix) <- sample_names
    colnames(CSCS_matrix) <- sample_names
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        CSCS_matrix[i,j] <- CSCS_all[i, j] / max(CSCS_all[i, i], CSCS_all[j, j])
      }
    }
    beta_div <- 1-CSCS_matrix
  } else {
    stop("Invalid method. Please use 'Gen.Uni', 'bray' or 'jaccard'")
  }


  return(beta_div)
}




# Function to calculate group distances against a reference group
CalculateGroupBetaDistance <- function(mmo, beta_div, reference_group, groups) {
  metadata <- mmo$metadata
  distances <- data.frame(group = character(), distance = numeric())
  
  for (group in groups) {
    if (group != reference_group) {
      group_samples <- metadata %>% filter(group == !!group) %>% pull(sample)
      reference_samples <- metadata %>% filter(group == !!reference_group) %>% pull(sample)
      
      for (sample in group_samples) {
        for (ref_sample in reference_samples) {
          distance <- beta_div[sample, ref_sample]
          distances <- rbind(distances, data.frame(group = group,sample = sample, distance = distance))
        }
      }
    }
  }
  
  return(distances)
}
########################################################################################
# 1. Load Data
########################################################################################
# 1.1. Give directories
mzmine_featuredir <- 'raw_data/250724_features_ms2.csv'
metadatadir <- "raw_data/250728_mmo_metadata.csv"
canopus_formuladir <- "raw_data/canopus_formula_summary.tsv"
canopus_structuredir <- "raw_data/structure_identifications.tsv"
gls_db <- 'raw_data/250428_GLS_Ath_simplegrad_MMO.csv'
cos_dir <- 'raw_data/250728_sim_ms2_(modified)_cosine.csv'
dreams_dir <- 'raw_data/250728_sim_dreams.csv'
m2ds_dir <- 'raw_data/250728_sim_ms2deepscore.csv'


# 1.2. Initiate object
mmo <- GetMZmineFeature(mzmine_featuredir, metadatadir) # load mzmine featurelist and metadata
head(mmo$feature_data) # check the feature data
# if any outlier samples are present, remove them
# outlier_samples <- c('datafile.sl1_2.mzML.area')
# mmo$feature_data <- mmo$feature_data %>% select(-all_of(outlier_samples))
# mmo$metadata <- mmo$metadata %>% filter(!sample %in% outlier_samples)

# Add annotations
mmo <- AddSiriusAnnot(mmo, canopus_structuredir, canopus_formuladir) # Add annotation from sirius
head(mmo$sirius_annot)
mmo <- AddCustomAnnot(mmo, DB = gls_db, mztol = 5, rttol = 0.1) # Add annotation from custom database

# Normalize data
mmo <- ReplaceZero(mmo, method = 'one') # Replace 0 and NA values by 1
mmo <- MassNormalization(mmo) # Normalize peak area by sample mass in metadata
mmo <- MeancenterNormalization(mmo) # Add mean-centered area
mmo <- LogNormalization(mmo) # Add log-transformed area
mmo <- ZNormalization(mmo) # Add Zscore

# Import chemical distance data for chemical diversity analyses
mmo <- AddChemDist(mmo, cos_dir = cos_dir, dreams_dir = dreams_dir, m2ds_dir = m2ds_dir) 
#siriused_features <- mmo$sirius_annot %>% filter(!is.na(mmo$sirius_annot[[46]]) & mmo$sirius_annot[[46]] != "") %>% pull(feature)

summary(mmo)







# write.csv(mmo$zscore, 'zscore.csv')

########################################################################################
# 2. Annotation lists
########################################################################################
features <- mmo$sirius_annot

# Get list of glucosinolated using sirius annotation
gls_hits <- mmo$custom_annot %>% filter(lengths(custom_annot) > 0)
GLSs <- gls_hits %>% pull(feature)
# Get list of flavonoids using sirius annotation
FLVs <- features %>% filter(str_detect(features[[46]], "Flavonoid")) %>% pull(feature)

GLS_id <- FeatureToID(mmo, GLSs)
GLS_feature <- IDToFeature(mmo, GLS_id)


########################################################################################
# 3. Pairwise comparison
########################################################################################

# 3.1. Add pairwise comparisons
mmo <- PairwiseComp(mmo, 'ctrl', 'sl1')
mmo <- PairwiseComp(mmo, 'ctrl', 'le1')

# 3.2. Get DAMs from the comparisons
DAMs <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.1)
DAMs_up <- DAMs$DAMs_up
DAMs_down <- DAMs$DAMs_down
head(DAMs_up)



########################################################################################
# 4. PLS-DA plot
#######################################################################################
# Define your custom colors for each group
custom_colors <- c("ctrl" = "#999999", "sl1" = "#fdcdac", "le1" = "#b3e2cd")
#Or automatically give colors
custom_colors <- setNames(brewer.pal(length(unique(mmo$metadata$group)), "Set3"), unique(mmo$metadata$group))

# Then plot the PLS-DA
PLSDAplot(mmo, color = custom_colors, outdir = 'plots/plsda/PLSDA_Z.pdf', normalization = 'Z')
PLSDAplot(mmo, color = custom_colors, topk = 0, outdir = 'plots/plsda/PLSDA_Z_noLoadings.pdf', normalization = 'Z')





########################################################################################
########################################################################################
# 5. Venn Diagram
########################################################################################
# 5.1. Define groups
VennInput <- list(
  sl1.up = DAMs_up$ctrl_vs_sl1.up,
  le1.up = DAMs_up$ctrl_vs_le1.up
)


# 5.2. Plot
# 5.2.1. Venn Diagram
ggvenn(VennInput, stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) +
  theme(legend.position = "none") 
ggsave("plots/Venn/Venn_Upreg.pdf", height = 5, width = 5)

# 5.2.2. Upset plot
pdf("plots/Venn/Upset_Upreg.pdf",7, 5)
upset(fromList(VennInput), nsets=10, nintersects=20,order.by='freq', mainbar.y.label='Features in Set', line.size=1, point.size=4, shade.color='white', text.scale=1, show.numbers=FALSE)
dev.off()

########################################################################################
# 6. Volcano Plot
#######################################################################################
comparison_columns <- colnames(mmo$pairwise)
log2FC_columns <- grep("log2FC", comparison_columns, value = TRUE)
comparisons <- sub("_$", "", unique(sub("log2FC", "", log2FC_columns)))  # Remove trailing underscore from comparisons
for (comp in comparisons){
  VolcanoPlot(mmo, comp = comp, topk = 10, outdir = paste('plots/volcano/volcano_', comp, '.pdf', sep = ''))
}

########################################################################################
# 7. Heat Map
#######################################################################################
# 7.1. Generate inputs for heatmap
HMinput_total <- GenerateHeatmapInputs(mmo, filter_feature = FALSE, feature_list = feature_list, 
                                filter_group = FALSE, group_list = group_list, 
                                summarize = 'mean', control_group = 'ctrl', 
                                normalization = 'Z', distance = 'dreams')
summary(HMinput)

HMinput_GLSs <- GenerateHeatmapInputs(mmo, filter_feature = TRUE, feature_list = GLSs, 
                                filter_group = FALSE, group_list = group_list, 
                                summarize = 'mean', control_group = 'ctrl', 
                                normalization = 'Z', distance = 'dreams')

HM_input <- HM_input_total

# 7.2. Generate NPC-based annotation table for heatmap
sirius_annot <- mmo$sirius_annot
# Get NPC Annotations
NPC_pathway <- unique(sirius_annot[[32]])
NPC_pathway <- NPC_pathway[!is.na(NPC_pathway) & NPC_pathway != ""] # remove NA and empty 
NPC_superclass <- unique(sirius_annot[[34]])
NPC_superclass <- NPC_superclass[!is.na(NPC_superclass) & NPC_superclass != ""] # remove NA and empty
NPC_class <- unique(sirius_annot[[36]])
NPC_class <- NPC_class[!is.na(NPC_class) & NPC_class != ""] # remove NA and empty

sirius_annot_filtered <- sirius_annot %>%
  # select(id = 1, NPC_pathway = 30, NPC_class = 34) %>%
  select(id = 1, NPC_class = 36, NPC_superclass = 34, NPC_pathway = 32) %>%
  # select(id = 1, NPC_pathway = 30, NPC_superclass = 32) %>%
  filter(id %in% rownames(HMinput$dist_matrix)) # get features with fingerprints

rownames(sirius_annot_filtered) <- sirius_annot_filtered$id
sirius_annot_filtered$id <- NULL

ann_colors = list(
    NPC_pathway = setNames(brewer.pal(length(NPC_pathway), "Set2"), NPC_pathway),
    NPC_class = setNames(viridis(length(NPC_class)), NPC_class),
    NPC_superclass = setNames(viridis(length(NPC_superclass)), NPC_superclass)
    )
# 7.3. Plot heatmap
pdf("plots/heatmap/dreams_total_Z.pdf", width = 15, height = 20)
pheatmap(mat = HMinput$FC_matrix, 
     cluster_rows = TRUE, #do not change
     clustering_distance_rows = HMinput$dist_matrix, # Delete this line for UPGMA clustering of rows
     cluster_cols = FALSE, 
     clustering_method = "average", #UPGMA
     show_rownames = TRUE, 
     show_colnames = TRUE,
     annotation_row = sirius_annot_filtered, # Dataframe with the rownames are identical with 'mat' and gives annotation
     annotation_colors = ann_colors,
     cellwidth = 25,
     cellheight = 0.05,
     treeheight_row = 100,
     fontsize_row = 3,
     fontsize_col = 15,
     scale = 'none',
     annotation_names_row = TRUE,
    #  labels_row = row_label,
     border_color = NA,
     color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


########################################################################################
# 8. Chemical term enrichment analysis 
#######################################################################################
# 8.1. For a single set of features, a detailed enrichment plot can be generated
# There are two plotting styles available
CanopusListEnrichmentPlot(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plots/enrichment/sl1_up.pdf', height = 6, width = 6)
CanopusListEnrichmentPlot_2(mmo, DAMs_up$ctrl_vs_sl1.up, pthr = 0.1, outdir = 'plots/enrichment/sl1_up_2.pdf', topn = 10, height = 6, width = 6)


# 8.2. For a list of sets features, a summary enrichment plot can be generated
# The summary enrichment plot can be generated for either a single level of CANOPUS classification (8.2.1) or for all levels (8.2.2)
# 8.2.1. For a single level of CANOPUS classification
term_levels <- c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_most_specific', 'ClassyFire_level5', 'ClassyFire_subclass', 'ClassyFire_class', 'ClassyFire_superclass')
CanopusLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'NPC_class', pthr = 0.1, prefix = 'plots/enrichment/DAMs_up_NPC_class.pdf')

# 8.2.2. For all levels of CANOPUS classification
# All levels, or onlt NPC or ClassyFire
terms <- c('NPC', 'ClassyFire', 'all_terms')
CanopusAllLevelEnrichmentPlot(mmo, DAMs_up, term_level = 'all_terms', pthr = 0.1, prefix = 'plots/enrichment/DAMs_up_all_terms', width = 8, height = 12)

########################################################################################
# 9. Regression with metadata
#######################################################################################
# To perform regression, the phenotype of interest should be defined in the metadata as a column
# 9.1. Regression of individual feature against a phenotype
# model can be 'lm' or 'pearson' or 'lmm'
# groups can be a vector of group names, or a single group name which has the phenotpye values in the metadata
FeaturePerformanceRegression(mmo, target = '477.0636_7.5687', phenotype = 'sl', groups = c('sl1'), 
  model = 'lm', normalization = 'Z', 
  output = paste('plots/phenotype_regression/477.0636_7.5687_sl_lm.pdf'))



# 9.2. Regression of all features against a phenotype and plotting the results 
sl.lm <- GetPerformanceFeatureRegression(mmo, phenotype = 'sl', groups = c('sl1'), 
                                        DAM.list = list(sl.up = DAMs_up$ctrl_vs_sl1.up, sl.down = DAMs_down$ctrl_vs_sl1.down), comparisons = c('ctrl_vs_sl1'))
sl.lm.sig <- sl.lm %>% filter(p_value < 0.1)
head(sl.lm.sig)

sl.cor <- GetPerformanceFeatureCorrelation(mmo, phenotype = 'sl', groups = c('sl1'), 
                                        DAM.list = list(sl.up = DAMs_up$ctrl_vs_sl1.up, sl.down = DAMs_down$ctrl_vs_sl1.down), comparisons = c('ctrl_vs_sl1'))
sl.cor.sig <- sl.cor %>% filter(p_value < 0.1)
head(sl.cor.sig)

# One may want to plot the regression results along with the fold change values
PlotFoldchangeResistanceRegression(sl.lm.sig, fold_change = 'ctrl_vs_sl1_log2FC', 
                                   color = c('sl.up' = '#d42525ff', 'sl.down' = '#281e99ff'), 
                                   output_dir = 'plots/phenotype_regression/sl_lm_sig.pdf')


########################################################################################
# 10. Chemical Diversity Measures
#######################################################################################
# 10.1. Alpha diversity by Hill numbers
# q : Hill number order, 0 for richness, 1 for Shannon, 2 for Simpson
# mode : weighted for the chemical structure (use distance argument), unweighted for no chemical weight
# filter_feature : if TRUE, only features in the feature_list are used
# feature_list : a list of features to filter the alpha diversity calculation
alphadiv <- GetAlphaDiversity(mmo, q = 3, normalization = 'Log', mode = 'weighted', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)

# 10.1.1 Plot the alpha diversity
ggplot(alphadiv, aes(x = group, y = hill_number)) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 0.5) +
  theme_classic() +
  labs(title = 'Alpha Diversity by Hill Numbers', x = 'Group', y = 'Hill Number') +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('plots/alphadiv/dreams_q3_log.pdf', width = 5, height = 5)
# Test for significant differences between groups with ANOVA
anova <- anova_tukey_dunnett(alphadiv, 'hill_number ~ group')
write_anova(anova, 'plots/alphadiv/anova_hill_number.csv')

# 10.1.2 Plot the feature distribution
# Function to plot histogram of values for each feature for given groups
groups <- c('ctrl', 'sl1', 'le1')
group.mean <- GetGroupMeans(mmo, normalization = 'Z', filter_feature = FALSE, feature_list = NULL)[,-1]
long.group.mean <- data.frame(value = double(), group = character())
for (group in groups) {
  group_data <- data.frame(value = group.mean[, group], group = group)
  colnames(group_data) <- c('value', 'group')
  long.group.mean <- rbind(long.group.mean, group_data)
}

ggplot(data = long.group.mean, aes(x = value, fill = group, color = group)) +
  geom_density(position = 'identity', alpha = 0) +
  theme_classic() +
  labs(x = "Normalized peak intensity", y = "Density") +
  # scale_fill_manual(values = custom_colors) +
  theme(legend.position = "right")
ggsave('plots/alphadiv/density_function_Z.pdf', height = 6, width = 6)  

metadata <- mmo$metadata
feature <- mmo$zscore
# The distribution should be shown!
plot_data <- data.frame(value = double(), rank = double(), sample = character(), group = character())
for (group in groups) {
  group_samples <- metadata %>% filter(group == !!group) %>% pull(sample)
  for (sample in group_samples) {
    sample_data <- unlist(as.vector(feature[, sample]))
    sorted_data <- sort(sample_data, decreasing = TRUE)
    sample_plot_data <- data.frame(value = sorted_data, rank = seq_along(sorted_data), group = group, sample = sample)
    plot_data <- rbind(plot_data, sample_plot_data)
  }
}


ggplot(plot_data, aes(x = rank, y = value, color = group)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5), alpha = 0.06) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "Sorted Feature", x = "Rank", y = "Value")
ggsave('plots/alphadiv/sorted_intensity_Z.png', height = 6, width = 6)



# 10.2. Beta diversity
# Calculate beta diversity for different normalizations and methods
# For unweighted beta diversity, Bray-Curtis or Jaccard distance can be used
bray <- GetBetaDiversity(mmo, method = 'bray', normalization = 'Log', filter_feature = FALSE, feature_list = NULL)
jaccard <- GetBetaDiversity(mmo, method = 'jaccard', normalization = 'Log', filter_feature = FALSE, feature_list = NULL)

# For weighted beta diversity, Generalized UniFrac can be used
guni <- GetBetaDiversity(mmo, method = 'Gen.Uni', normalization = 'Log', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)
guni.0 <- guni[,,'d_0'] # GUniFrac with alpha 0
guni.05 <- guni[,,'d_0.5'] # GUniFrac with alpha 0.5
guni.1 <- guni[,,'d_1'] # GUniFrac with alpha 1

CSCS <- GetBetaDiversity(mmo, method = 'CSCS', normalization = 'Log', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)

# 10.2.1. NMDS plots for beta diversity
nmds <- metaMDS(CSCS, k = 2, try = 50, trymax = 100)
nmds_coords <- as.data.frame(scores(nmds))
groups <- c()
for (col in colnames(mmo$feature_data)[-c(1, 2)]) {
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
ggsave('plots/betadiv/NDMS_guni05.pdf', height = 6, width = 6)
# 10.2.2. Quantification of beta diversity 
group_distances <- CalculateGroupBetaDistance(mmo, beta_div = guni.05, reference_group = 'ctrl', groups = c('le1', 'sl1'))
ggplot(group_distances, aes(x = group, y = distance)) +
      geom_boxplot(outlier.shape = NA) +
      geom_beeswarm(size = 0.5) +
      theme_classic() +
      labs(x = "Group", y = "Beta Diversity")
ggsave('plots/betadiv/group_dist.pdf', height = 6, width = 6)

########################################################################################
# 11. Exporting compounds of interest
#######################################################################################

ExportFeaturesToCSV(mmo, feature_list = DAMs_up$ctrl_vs_sl1.up, normalization = 'None', output_dir = 'output/sl1_up_features.csv')

distance = 'dreams'
normalization = 'None'


######
# CSCS Implementation
CSS <- 1-GetDistanceMat(mmo, distance = distance)
diag(CSS)
q.feature <- GetNormFeature(mmo, normalization = normalization) %>% filter(id %in% colnames(CSS))
relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) x / sum(x))
CSCS_all <- t(relative_proportions) %*% CSS %*% relative_proportions

sample_names <- colnames(relative_proportions)
n_samples <- length(sample_names)
CSCS_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
rownames(CSCS_matrix) <- sample_names
colnames(CSCS_matrix) <- sample_names
for (i in 1:n_samples) {
  for (j in 1:n_samples) {
    CSCS_matrix[i,j] <- CSCS_all[i, j] / max(CSCS_all[i, i], CSCS_all[j, j])
  }
}

guni <- GetBetaDiversity(mmo, method = 'Gen.Uni', normalization = 'None', distance = 'dreams', filter_feature = FALSE, feature_list = NULL)
