# author: arezou kahnamouei
# Comprehensive-ESCRT-Protein-Analysis-Pipeline
# Load required packages for analysis
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

# Define ESCRT complexes and their component genes
escrt_complexes <- list(
  "ESCRT-0" = c("HGS", "STAM1", "STAM2", "HRS"),
  "ESCRT-I" = c("TSG101", "VPS28", "VPS37A", "VPS37B", "VPS37C", "VPS37D", 
                "MVB12A", "MVB12B", "UBAP1"),
  "ESCRT-II" = c("VPS36", "VPS25", "VPS22", "SNF8", "EAP45", "EAP20", "EAP30"),
  "ESCRT-II/III" = c("CHMP7", "VPS31"),
  "ESCRT-III" = c("CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", 
                  "CHMP6", "IST1", "CHMP5", "CHMP1A", "CHMP1B", "VPS24", "VPS2", 
                  "VPS20", "VPS32", "VPS60", "DID2"),
  "VPS4-Complex" = c("VPS4A", "VPS4B", "VTA1", "LIP5", "SKD1"),
  "BRO1" = c("PDCD6IP", "ALIX", "PTPN23", "HD-PTP", "BROX")
)

# Function to analyze expression data with sample-specific patterns
analyze_expression <- function(data_file) {
  message("Reading and processing expression data...")
  data <- read.delim(data_file, stringsAsFactors = FALSE)
  
  # Identify intensity columns specifically for your samples
  # First try LFQ intensity columns
  nuclei_cols <- grep("LFQ.intensity.APEXP62_0[1-5]", colnames(data), value = TRUE)
  micronuclei_cols <- grep("LFQ.intensity.APEXONLY_0[1-5]", colnames(data), value = TRUE)
  
  # If no LFQ columns found, use regular intensity columns
  if(length(nuclei_cols) == 0 || length(micronuclei_cols) == 0) {
    message("Using regular intensity columns...")
    nuclei_cols <- grep("Intensity.APEXP62_0[1-5]", colnames(data), value = TRUE)
    micronuclei_cols <- grep("Intensity.APEXONLY_0[1-5]", colnames(data), value = TRUE)
  }
  
  # Print identified columns
  message("\nColumns being used for analysis:")
  message("Nuclei columns: ", paste(nuclei_cols, collapse = ", "))
  message("Micronuclei columns: ", paste(micronuclei_cols, collapse = ", "))
  
  if(length(nuclei_cols) == 0 || length(micronuclei_cols) == 0) {
    stop("Could not find appropriate intensity columns")
  }
  
  # Process genes
  data$Gene.List <- strsplit(as.character(data$Gene.names), ";")
  all_target_genes <- unique(unlist(escrt_complexes))
  
  # Find ESCRT genes in the data
  target_rows <- which(sapply(data$Gene.List, function(genes) {
    any(genes %in% all_target_genes)
  }))
  
  if(length(target_rows) == 0) {
    stop("No ESCRT genes found in dataset")
  }
  
  message(sprintf("\nFound %d ESCRT genes", length(target_rows)))
  
  # Extract and process target data
  target_data <- data[target_rows, ]
  
  # Calculate expression statistics
  results <- data.frame(
    Gene = sapply(target_data$Gene.List, function(genes) {
      intersect(genes, all_target_genes)[1]
    }),
    Nuclei_Mean = rowMeans(as.matrix(target_data[, nuclei_cols]), na.rm = TRUE),
    Nuclei_SD = apply(as.matrix(target_data[, nuclei_cols]), 1, sd, na.rm = TRUE),
    Micronuclei_Mean = rowMeans(as.matrix(target_data[, micronuclei_cols]), na.rm = TRUE),
    Micronuclei_SD = apply(as.matrix(target_data[, micronuclei_cols]), 1, sd, na.rm = TRUE)
  )
  
  # Add analysis metrics
  results <- results %>%
    mutate(
      log2FC = log2(Micronuclei_Mean / Nuclei_Mean),
      Complex = sapply(Gene, function(g) {
        names(escrt_complexes)[sapply(escrt_complexes, function(x) g %in% x)][1]
      })
    )
  
  # Perform statistical testing
  results$p_value <- sapply(target_rows, function(row) {
    nuclei_values <- as.numeric(data[row, nuclei_cols])
    micro_values <- as.numeric(data[row, micronuclei_cols])
    
    if(sum(!is.na(nuclei_values)) >= 2 && sum(!is.na(micro_values)) >= 2) {
      tryCatch({
        t.test(micro_values, nuclei_values)$p.value
      }, error = function(e) NA)
    } else {
      NA
    }
  })
  
  # Add multiple testing correction and significance
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results$significance <- case_when(
    is.na(results$p_adj) ~ "Not tested",
    results$p_adj < 0.001 ~ "***",
    results$p_adj < 0.01 ~ "**",
    results$p_adj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  return(results)
}

# Function to analyze complex effects
analyze_complexes <- function(results) {
  complex_stats <- results %>%
    group_by(Complex) %>%
    summarise(
      total_genes = n(),
      detected_genes = sum(!is.na(log2FC)),
      significant_genes = sum(significance %in% c("*", "**", "***"), na.rm = TRUE),
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      median_log2FC = median(log2FC, na.rm = TRUE),
      sd_log2FC = sd(log2FC, na.rm = TRUE),
      upregulated = sum(log2FC > 0 & significance != "ns", na.rm = TRUE),
      downregulated = sum(log2FC < 0 & significance != "ns", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(abs(mean_log2FC)))
  
  return(complex_stats)
}

# Create enhanced visualizations
create_visualizations <- function(results, complex_stats) {
  dir.create("plots", showWarnings = FALSE)
  
  # Enhanced volcano plot
  volcano_plot <- ggplot(results, aes(x = log2FC, y = -log10(p_adj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Complex, size = abs(log2FC)), alpha = 0.7) +
    geom_text_repel(
      data = subset(results, p_adj < 0.05),
      aes(label = paste0(Gene, " (", significance, ")")),
      size = 3,
      max.overlaps = 15
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    labs(
      title = "ESCRT Gene Expression Changes",
      subtitle = "Micronuclei vs Nuclei",
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)"
    )
  
  # Complex effect plot
  complex_plot <- ggplot(complex_stats, 
                         aes(x = reorder(Complex, abs(mean_log2FC)))) +
    geom_bar(aes(y = mean_log2FC, fill = Complex), stat = "identity") +
    geom_errorbar(
      aes(ymin = mean_log2FC - sd_log2FC,
          ymax = mean_log2FC + sd_log2FC),
      width = 0.2
    ) +
    geom_text(aes(y = ifelse(mean_log2FC >= 0,
                             mean_log2FC + sd_log2FC + 0.1,
                             mean_log2FC - sd_log2FC - 0.1),
                  label = sprintf("%d/%d sig", significant_genes, detected_genes)),
              size = 3) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "ESCRT Complex Changes",
      subtitle = "Mean log2 Fold Change with Standard Deviation",
      x = "Complex",
      y = "Mean log2 Fold Change"
    )
  
  # Save plots
  ggsave("plots/volcano_plot.pdf", volcano_plot, width = 12, height = 8)
  ggsave("plots/complex_effects.pdf", complex_plot, width = 12, height = 8)
  
  return(list(volcano = volcano_plot, complex = complex_plot))
}

# Main analysis function
main_analysis <- function() {
  message("Starting ESCRT gene expression analysis...")
  
  # Analyze expression
  results <- analyze_expression("proteinGroups_APEX_P62.txt")
  
  # Analyze complex effects
  complex_stats <- analyze_complexes(results)
  
  # Create visualizations
  plots <- create_visualizations(results, complex_stats)
  
  # Save detailed results
  write.csv(results, "gene_expression_results.csv", row.names = FALSE)
  write.csv(complex_stats, "complex_impact_summary.csv", row.names = FALSE)
  
  # Print detailed summary
  cat("\nAnalysis Summary\n")
  cat("================\n\n")
  
  cat("1. Overall Statistics\n")
  cat("-------------------\n")
  cat(sprintf("Total ESCRT genes analyzed: %d\n", nrow(results)))
  cat(sprintf("Genes with significant changes: %d\n", 
              sum(results$significance %in% c("*", "**", "***"), na.rm = TRUE)))
  
  cat("\n2. Complex-Level Analysis\n")
  cat("----------------------\n")
  print(complex_stats)
  
  cat("\n3. Most Changed Genes\n")
  cat("-----------------\n")
  top_genes <- results %>%
    arrange(desc(abs(log2FC))) %>%
    head(10)
  print(top_genes[, c("Gene", "Complex", "log2FC", "significance")])
  
  return(list(
    results = results,
    complex_stats = complex_stats,
    plots = plots
  ))
}

# Run the analysis
results <- main_analysis()


################################################################################
# Set working directory
setwd("C:/Users/ASUS/Desktop/Projects/Mrs.Kahnemouei/2/datasets/finaldata/final2")

# Load required packages for advanced visualization
library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

# Create a custom color palette for ESCRT complexes
escrt_colors <- c(
  "ESCRT-0" = "#E64B35",
  "ESCRT-I" = "#4DBBD5",
  "ESCRT-III" = "#00A087",
  "VPS4-Complex" = "#3C5488",
  "BRO1" = "#F39B7F"
)

# 1. Enhanced Volcano Plot with Custom Design
create_enhanced_volcano <- function(results) {
  ggplot(results, aes(x = log2FC, y = -log10(p_adj))) +
    # Add shaded regions for fold change thresholds
    annotate("rect", xmin = -Inf, xmax = -1, ymin = -log10(0.05), ymax = Inf,
             fill = "lightblue", alpha = 0.2) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -log10(0.05), ymax = Inf,
             fill = "lightpink", alpha = 0.2) +
    # Add threshold lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    # Add points
    geom_point(aes(color = Complex, size = abs(log2FC)), alpha = 0.8) +
    geom_text_repel(
      aes(label = Gene),
      size = 3.5,
      box.padding = 0.5,
      max.overlaps = 20,
      segment.color = "gray50",
      segment.alpha = 0.6
    ) +
    # Customize appearance
    scale_color_manual(values = escrt_colors) +
    scale_size_continuous(range = c(2, 6)) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80"),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "ESCRT Gene Expression Changes",
      subtitle = "Comparison of Micronuclei vs Nuclei",
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "ESCRT Complex",
      size = "|log2FC|"
    )
}

# 2. Complex Impact Visualization
create_complex_impact <- function(complex_stats) {
  # Calculate percentage of significant genes
  complex_stats <- complex_stats %>%
    mutate(sig_percent = significant_genes / total_genes * 100)
  
  ggplot(complex_stats, aes(x = reorder(Complex, abs(mean_log2FC)))) +
    # Add bars for mean change
    geom_bar(aes(y = mean_log2FC, fill = Complex), stat = "identity", width = 0.7) +
    # Add error bars
    geom_errorbar(
      aes(ymin = mean_log2FC - sd_log2FC,
          ymax = mean_log2FC + sd_log2FC),
      width = 0.3,
      color = "gray40"
    ) +
    # Add gene count labels
    geom_text(aes(y = ifelse(mean_log2FC >= 0,
                             mean_log2FC + sd_log2FC + 0.1,
                             mean_log2FC - sd_log2FC - 0.1),
                  label = sprintf("%d genes", total_genes)),
              size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = escrt_colors) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Impact of ESCRT Complexes",
      x = "ESCRT Complex",
      y = "Mean log2 Fold Change"
    )
}

# 3. Expression Heatmap
create_expression_heatmap <- function(results) {
  # Prepare data for heatmap
  heatmap_data <- results %>%
    select(Gene, Complex, log2FC) %>%
    arrange(Complex, log2FC)
  
  # Create matrix for heatmap
  mat <- matrix(
    heatmap_data$log2FC,
    ncol = 1,
    dimnames = list(heatmap_data$Gene, "log2FC")
  )
  
  # Create annotation
  ha <- rowAnnotation(
    Complex = heatmap_data$Complex,
    col = list(Complex = escrt_colors),
    show_legend = TRUE
  )
  
  # Create color mapping function
  col_fun <- colorRamp2(
    c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
    c("#2166AC", "white", "#B2182B")
  )
  
  # Create heatmap
  ht <- Heatmap(
    mat,
    name = "log2FC",
    col = col_fun,
    right_annotation = ha,
    row_names_gp = gpar(fontsize = 10),
    row_title = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = "Expression Changes by Gene and Complex",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    width = unit(3, "inches"),
    height = unit(8, "inches")
  )
  
  return(ht)
}

# 4. Gene Ranking Plot
create_gene_ranking <- function(results) {
  results %>%
    arrange(desc(abs(log2FC))) %>%
    ggplot(aes(x = reorder(Gene, abs(log2FC)), y = abs(log2FC))) +
    geom_bar(aes(fill = Complex), stat = "identity") +
    geom_text(aes(label = sprintf("%.2f", log2FC)),
              hjust = -0.2,
              size = 3) +
    coord_flip() +
    scale_fill_manual(values = escrt_colors) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Ranking of Gene Expression Changes",
      x = "Gene",
      y = "Absolute log2 Fold Change"
    )
}

# Create and save all plots
create_all_plots <- function(results, complex_stats) {
  # Create directory for new plots
  dir.create("advanced_plots", showWarnings = FALSE)
  
  # 1. Enhanced Volcano Plot
  volcano <- create_enhanced_volcano(results)
  ggsave("advanced_plots/enhanced_volcano.pdf", volcano, width = 12, height = 10)
  
  # 2. Complex Impact Plot
  impact <- create_complex_impact(complex_stats)
  ggsave("advanced_plots/complex_impact.pdf", impact, width = 10, height = 8)
  
  # 3. Expression Heatmap
  pdf("advanced_plots/expression_heatmap.pdf", width = 8, height = 12)
  draw(create_expression_heatmap(results))
  dev.off()
  
  # 4. Gene Ranking Plot
  ranking <- create_gene_ranking(results)
  ggsave("advanced_plots/gene_ranking.pdf", ranking, width = 10, height = 8)
  
  # Create combined overview plot
  combined_plot <- (volcano + impact) / (ranking)
  ggsave("advanced_plots/combined_overview.pdf", combined_plot, width = 16, height = 20)
  
  message("Advanced plots have been saved in the 'advanced_plots' directory")
}

# Run the visualization
create_all_plots(results$results, results$complex_stats)


##########################################################################################
# Load required packages
library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gridExtra)

# Set working directory
setwd("C:/Users/ASUS/Desktop/Projects/Mrs.Kahnemouei/2/datasets/finaldata/final2")

# Create custom color schemes
custom_colors <- c(
  "ESCRT-0" = "#FF6B6B",
  "ESCRT-I" = "#4ECDC4",
  "ESCRT-III" = "#45B7D1",
  "VPS4-Complex" = "#96CEB4",
  "BRO1" = "#FFEEAD"
)

# 1. Enhanced Molecular Functions Plot (similar to Image 1B)
create_molecular_functions_plot <- function(results) {
  # Calculate mean fold changes by function
  function_data <- results %>%
    group_by(Complex) %>%
    summarise(
      mean_change = mean(abs(log2FC), na.rm = TRUE),
      n_genes = n()
    ) %>%
    arrange(desc(mean_change))
  
  ggplot(function_data, aes(x = reorder(Complex, mean_change), y = mean_change)) +
    geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
    geom_text(aes(label = n_genes), hjust = -0.2, size = 3.5) +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "ESCRT Complex Functions",
      subtitle = "Micronuclei vs. Primary Nuclei",
      x = "",
      y = "Mean Absolute log2 Fold Change"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

# 2. Time Course-Style Plot (similar to Image 2K)
create_expression_timeline <- function(results) {
  results %>%
    arrange(log2FC) %>%
    mutate(order = row_number()) %>%
    ggplot(aes(x = order, y = log2FC)) +
    geom_line(color = "#FF6B6B", size = 1) +
    geom_point(aes(color = Complex), size = 3) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    labs(
      title = "Expression Changes Across ESCRT Genes",
      x = "Genes (ordered by expression change)",
      y = "log2 Fold Change (Micronuclei/Nuclei)"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# 3. Enhanced Heatmap with Clustering (similar to Image 3)
create_enhanced_heatmap <- function(results) {
  # Prepare matrix for heatmap
  mat <- matrix(
    results$log2FC,
    ncol = 1,
    dimnames = list(results$Gene, "log2FC")
  )
  
  # Create annotations
  ha <- rowAnnotation(
    Complex = results$Complex,
    col = list(Complex = custom_colors)
  )
  
  # Create color mapping
  col_fun <- colorRamp2(
    c(min(results$log2FC), 0, max(results$log2FC)),
    c("#3288BD", "white", "#D53E4F")
  )
  
  Heatmap(
    mat,
    name = "log2FC",
    col = col_fun,
    right_annotation = ha,
    cluster_rows = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_title = "Gene Expression Changes",
    width = unit(4, "inches"),
    height = unit(8, "inches")
  )
}

# 4. Bar Plot with Error Bars (similar to Image 4B)
create_complex_barplot <- function(results) {
  complex_stats <- results %>%
    group_by(Complex) %>%
    summarise(
      mean_fc = mean(log2FC, na.rm = TRUE),
      se = sd(log2FC, na.rm = TRUE) / sqrt(n())
    )
  
  ggplot(complex_stats, aes(x = Complex, y = mean_fc, fill = Complex)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_fc - se, ymax = mean_fc + se),
      width = 0.2
    ) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal() +
    labs(
      title = "ESCRT Complex Expression Changes",
      x = "",
      y = "Mean log2 Fold Change"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    )
}

# 5. Circular Visualization (New plot style)
create_circular_plot <- function(results) {
  # Prepare data
  results <- results %>%
    arrange(Complex, abs(log2FC))
  
  # Create circular plot
  ggplot(results, aes(x = reorder(Gene, log2FC), y = abs(log2FC))) +
    geom_segment(aes(xend = Gene, yend = 0, color = Complex), size = 1) +
    geom_point(aes(color = Complex), size = 3) +
    scale_color_manual(values = custom_colors) +
    coord_polar() +
    theme_minimal() +
    labs(
      title = "Circular View of Expression Changes",
      subtitle = "Distance from center represents magnitude of change"
    ) +
    theme(
      axis.text.x = element_text(angle = 0),
      plot.title = element_text(size = 14, face = "bold")
    )
}

# Create and save all plots
create_publication_plots <- function(results) {
  dir.create("publication_plots", showWarnings = FALSE)
  
  # Create individual plots
  p1 <- create_molecular_functions_plot(results)
  p2 <- create_expression_timeline(results)
  p4 <- create_complex_barplot(results)
  p5 <- create_circular_plot(results)
  
  # Save individual plots
  ggsave("publication_plots/molecular_functions.pdf", p1, width = 10, height = 8)
  ggsave("publication_plots/expression_timeline.pdf", p2, width = 12, height = 6)
  ggsave("publication_plots/complex_barplot.pdf", p4, width = 8, height = 6)
  ggsave("publication_plots/circular_plot.pdf", p5, width = 10, height = 10)
  
  # Save heatmap separately
  pdf("publication_plots/enhanced_heatmap.pdf", width = 8, height = 12)
  draw(create_enhanced_heatmap(results))
  dev.off()
  
  # Create combined overview
  combined_plot <- (p1 + p2) / (p4 + p5)
  ggsave("publication_plots/combined_overview.pdf", combined_plot, 
         width = 16, height = 16)
  
  message("All publication-quality plots have been saved in the 'publication_plots' directory")
}

# Run the visualization with your results
create_publication_plots(results$results)



########################################################################
# Load required packages for visualization
library(tidyverse)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(gridExtra)

# Create directory for plots
dir.create("state_comparison_plots", showWarnings = FALSE)

# Define color schemes for better visualization
state_colors <- c("MICRONUCLEIul" = "#FF6B6B", "NUCLEI" = "#4ECDC4")
escrt_colors <- c(
  "ESCRT-0" = "#E64B35",
  "ESCRT-I" = "#4DBBD5",
  "ESCRT-III" = "#00A087",
  "VPS4-Complex" = "#3C5488",
  "BRO1" = "#F39B7F"
)

# 1. Direct State Comparison Plot
create_state_comparison <- function(results) {
  # Calculate mean values and standard error for each gene
  comparison_data <- results %>%
    mutate(
      fold_change = Micronuclei_Mean / Nuclei_Mean,
      log2FC = log2(fold_change)
    ) %>%
    arrange(desc(abs(log2FC)))
  
  # Create comparison plot
  ggplot(comparison_data, aes(x = reorder(Gene, log2FC), y = log2FC)) +
    # Add bars for fold change
    geom_bar(stat = "identity", aes(fill = ifelse(log2FC > 0, "Up", "Down"))) +
    # Add significance indicators
    geom_text(aes(label = significance, y = ifelse(log2FC >= 0, log2FC + 0.1, log2FC - 0.1)), 
              size = 3) +
    # Add fold change values
    geom_text(aes(label = sprintf("%.2f", log2FC), 
                  y = ifelse(log2FC >= 0, log2FC/2, log2FC/2)),
              color = "white", size = 3) +
    scale_fill_manual(values = c("Up" = "#FF6B6B", "Down" = "#4ECDC4")) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = "Gene Expression Changes between States",
      subtitle = "MICRONUCLEIul vs NUCLEI",
      x = "Gene",
      y = "log2 Fold Change",
      fill = "Direction"
    )
}

# 2. Expression Level Comparison
create_expression_levels <- function(results) {
  results %>%
    select(Gene, Complex, Nuclei_Mean, Micronuclei_Mean) %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(State = gsub("_Mean", "", State)) %>%
    ggplot(aes(x = reorder(Gene, Expression), y = Expression, fill = State)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Complex, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = state_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 9)
    ) +
    labs(
      title = "Expression Levels by State and ESCRT Complex",
      x = "Gene",
      y = "Expression Level"
    )
}

# 3. ESCRT Complex Response Plot
create_complex_response <- function(results) {
  # Calculate complex-level statistics
  complex_stats <- results %>%
    group_by(Complex) %>%
    summarise(
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      se = sd(log2FC, na.rm = TRUE) / sqrt(n()),
      n_sig = sum(significance != "ns", na.rm = TRUE),
      n_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      response = case_when(
        mean_log2FC > 0.5 ~ "Strong Up",
        mean_log2FC > 0 ~ "Mild Up",
        mean_log2FC < -0.5 ~ "Strong Down",
        mean_log2FC < 0 ~ "Mild Down",
        TRUE ~ "No Change"
      )
    )
  
  ggplot(complex_stats, aes(x = reorder(Complex, mean_log2FC))) +
    geom_bar(aes(y = mean_log2FC, fill = response), stat = "identity") +
    geom_errorbar(aes(ymin = mean_log2FC - se, ymax = mean_log2FC + se), width = 0.2) +
    geom_text(aes(y = ifelse(mean_log2FC >= 0, 
                             mean_log2FC + se + 0.1, 
                             mean_log2FC - se - 0.1),
                  label = sprintf("%d/%d genes", n_sig, n_total)),
              size = 3) +
    scale_fill_manual(values = c(
      "Strong Up" = "#FF6B6B",
      "Mild Up" = "#FFB4B4",
      "No Change" = "#CCCCCC",
      "Mild Down" = "#B4E4DE",
      "Strong Down" = "#4ECDC4"
    )) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    ) +
    labs(
      title = "ESCRT Complex Responses to Micronuclei Formation",
      x = "Complex",
      y = "Mean log2 Fold Change",
      fill = "Response Type"
    )
}

# 4. Detailed State-Complex Interaction Plot
create_state_complex_interaction <- function(results) {
  # Calculate mean expression by state and complex
  interaction_data <- results %>%
    group_by(Complex) %>%
    mutate(
      relative_nuclei = Nuclei_Mean / mean(Nuclei_Mean, na.rm = TRUE),
      relative_micro = Micronuclei_Mean / mean(Micronuclei_Mean, na.rm = TRUE)
    ) %>%
    ungroup()
  
  ggplot(interaction_data, aes(x = relative_nuclei, y = relative_micro)) +
    geom_point(aes(color = Complex, size = abs(log2FC))) +
    geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 15) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = escrt_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "State-Complex Interaction Analysis",
      subtitle = "Relative expression levels between states",
      x = "Relative Expression in Nuclei",
      y = "Relative Expression in Micronuclei",
      color = "ESCRT Complex",
      size = "|log2FC|"
    )
}

# Create and save all plots
create_state_analysis_plots <- function(results) {
  # Create individual plots
  p1 <- create_state_comparison(results)
  p2 <- create_expression_levels(results)
  p3 <- create_complex_response(results)
  p4 <- create_state_complex_interaction(results)
  
  # Save individual plots
  ggsave("state_comparison_plots/state_comparison.pdf", p1, width = 12, height = 10)
  ggsave("state_comparison_plots/expression_levels.pdf", p2, width = 14, height = 10)
  ggsave("state_comparison_plots/complex_response.pdf", p3, width = 10, height = 8)
  ggsave("state_comparison_plots/state_complex_interaction.pdf", p4, width = 12, height = 10)
  
  # Create and save combined overview
  combined_plot <- (p1 + p3) / (p2 + p4) +
    plot_annotation(
      title = "Comprehensive Analysis of MICRONUCLEIul vs NUCLEI States",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  ggsave("state_comparison_plots/combined_analysis.pdf", 
         combined_plot, 
         width = 20, 
         height = 20)
  
  message("All state comparison plots have been saved in 'state_comparison_plots' directory")
}

# Run the analysis with error handling
tryCatch({
  create_state_analysis_plots(results$results)
}, error = function(e) {
  message("Error creating plots: ", e$message)
  message("\nPlease check your data structure and try again.")
})



#########################################################################
# Load required packages for creating beautiful visualizations
library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(scales)
library(gridExtra)

# Create a directory for our enhanced visualizations
dir.create("enhanced_plots", showWarnings = FALSE)

# Define a beautiful color palette for ESCRT complexes
escrt_colors <- c(
  "ESCRT-0" = "#FF9999",      # Soft red
  "ESCRT-I" = "#66B2FF",      # Sky blue
  "ESCRT-III" = "#99FF99",    # Soft green
  "VPS4-Complex" = "#FFCC99", # Soft orange
  "BRO1" = "#FF99FF"          # Soft purple
)

# Function to create an enhanced waterfall plot
create_waterfall_plot <- function(data) {
  # Order genes by log2FC
  data <- data %>%
    mutate(Gene = factor(Gene, levels = Gene[order(log2FC)]))
  
  ggplot(data, aes(x = Gene, y = log2FC)) +
    # Add connecting lines for waterfall effect
    geom_segment(aes(x = Gene, xend = Gene, 
                     y = 0, yend = log2FC,
                     color = Complex),
                 size = 1.2) +
    # Add points
    geom_point(aes(color = Complex), size = 4) +
    # Add gene labels
    geom_text_repel(
      aes(label = sprintf("%s\n(%.2f)", Gene, log2FC)),
      size = 3,
      direction = "y",
      hjust = -0.2
    ) +
    # Customize colors and theme
    scale_color_manual(values = escrt_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = "Expression Changes Between States",
      subtitle = "Log2 Fold Change (Micronuclei vs Nuclei)",
      y = "log2 Fold Change",
      color = "ESCRT Complex"
    ) +
    coord_flip()
}

# Function to create a bubble plot of expression levels
create_bubble_plot <- function(data) {
  # Calculate mean expression for sizing bubbles
  data <- data %>%
    mutate(mean_expression = (Nuclei_Mean + Micronuclei_Mean)/2)
  
  ggplot(data, aes(x = Nuclei_Mean/1e6, y = Micronuclei_Mean/1e6)) +
    # Add diagonal reference line
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray70") +
    # Add bubbles
    geom_point(aes(size = mean_expression/1e6, 
                   color = Complex,
                   alpha = abs(log2FC))) +
    # Add gene labels
    geom_text_repel(
      aes(label = Gene),
      size = 3,
      max.overlaps = 20
    ) +
    # Customize appearance
    scale_color_manual(values = escrt_colors) +
    scale_size_continuous(range = c(3, 15)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    ) +
    labs(
      title = "Expression Level Comparison",
      subtitle = "Micronuclei vs Nuclei States",
      x = "Nuclei Expression (millions)",
      y = "Micronuclei Expression (millions)",
      size = "Mean Expression\n(millions)",
      alpha = "|log2FC|",
      color = "ESCRT Complex"
    )
}

# Function to create a complex-specific violin plot
create_complex_violin <- function(data) {
  # Reshape data for violin plot
  long_data <- data %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(
      State = factor(State, 
                     levels = c("Nuclei_Mean", "Micronuclei_Mean"),
                     labels = c("Nuclei", "Micronuclei"))
    )
  
  ggplot(long_data, aes(x = Complex, y = Expression/1e6)) +
    # Add violin plot
    geom_violin(aes(fill = Complex), alpha = 0.6) +
    # Add individual points
    geom_point(aes(color = Complex), 
               position = position_jitter(width = 0.2),
               size = 2) +
    # Add gene labels
    geom_text_repel(
      aes(label = Gene),
      size = 3,
      max.overlaps = 15
    ) +
    # Facet by state
    facet_wrap(~State) +
    # Customize appearance
    scale_fill_manual(values = escrt_colors) +
    scale_color_manual(values = escrt_colors) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Expression Distribution by Complex and State",
      y = "Expression Level (millions)",
      x = "ESCRT Complex"
    )
}

# Function to create a heatmap-style visualization
create_expression_heatmap <- function(data) {
  # Calculate z-scores for better visualization
  data <- data %>%
    mutate(
      Nuclei_Z = scale(Nuclei_Mean)[,1],
      Micro_Z = scale(Micronuclei_Mean)[,1]
    )
  
  # Reshape for heatmap
  long_data <- data %>%
    select(Gene, Complex, Nuclei_Z, Micro_Z) %>%
    pivot_longer(
      cols = c(Nuclei_Z, Micro_Z),
      names_to = "State",
      values_to = "Z_Score"
    )
  
  ggplot(long_data, aes(x = State, y = reorder(Gene, Z_Score))) +
    # Add heatmap tiles
    geom_tile(aes(fill = Z_Score)) +
    # Add Complex annotation
    geom_tile(aes(x = 2.5, fill = Complex), width = 0.3) +
    # Customize colors
    scale_fill_viridis_c(option = "A", name = "Z-Score") +
    scale_fill_manual(values = escrt_colors, guide = "none") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      title = "Expression Patterns Across States",
      x = "State",
      y = "Gene"
    )
}

# Create and save all plots
create_all_visualizations <- function(data) {
  # Create individual plots
  p1 <- create_waterfall_plot(data)
  p2 <- create_bubble_plot(data)
  p3 <- create_complex_violin(data)
  p4 <- create_expression_heatmap(data)
  
  # Save individual plots
  ggsave("enhanced_plots/waterfall_plot.pdf", p1, width = 12, height = 8)
  ggsave("enhanced_plots/bubble_plot.pdf", p2, width = 10, height = 8)
  ggsave("enhanced_plots/violin_plot.pdf", p3, width = 14, height = 8)
  ggsave("enhanced_plots/expression_heatmap.pdf", p4, width = 10, height = 12)
  
  # Create combined visualization
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Comprehensive ESCRT Expression Analysis",
      subtitle = "Comparing Micronuclei and Nuclei States",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 14, color = "gray40")
      )
    )
  
  # Save combined plot
  ggsave("enhanced_plots/combined_visualization.pdf", 
         combined_plot, 
         width = 20, 
         height = 20)
  
  message("Enhanced visualizations have been saved in the 'enhanced_plots' directory")
}

# Read and format the data
data <- data.frame(
  Gene = c("VTA1", "HGS", "VPS37C", "MVB12A", "VPS28", "TSG101", "VPS4A", 
           "STAM2", "IST1", "BROX", "PDCD6IP", "CHMP2B", "PTPN23", "CHMP4B", 
           "VPS37B", "CHMP5"),
  Nuclei_Mean = c(68069400, 225948000, 123364600, 18563400, 0, 215460000, 
                  57438800, 656718000, 58896400, 2042600, 370674000, 22252400, 
                  14344000, 22960200, 92180800, 21866200),
  Nuclei_SD = c(14693618.1, 94581628.6, 78515199.4, 18181012.4, 0, 82396671.0,
                11783931.6, 371002968.5, 17599994.9, 4567392.5, 93323484.1,
                13278773.9, 13773545.1, 27750948.9, 54721874.9, 21668465.9),
  Micronuclei_Mean = c(88261000, 108526200, 86404400, 2980200, 0, 79990800,
                       69454600, 321476000, 46759200, 4285400, 296154000,
                       27779200, 23699800, 26504800, 13325880, 22626600),
  Micronuclei_SD = c(24599315.3, 15765194.4, 18118693.2, 6663929.8, 0,
                     43046441.8, 19303529.1, 132523360.1, 13346444.7,
                     9582445.7, 62653064.4, 4962482.8, 7670642.2, 9944288.8,
                     9824017.8, 13077937.0),
  log2FC = c(0.3748, -1.0579, -0.5138, -2.6390, NA, -1.4295, 0.2740, -1.0306,
             -0.3329, 1.0690, -0.3238, 0.3200, 0.7244, 0.2071, -2.7902, 0.0493),
  Complex = c("VPS4-Complex", "ESCRT-0", "ESCRT-I", "ESCRT-I", "ESCRT-I",
              "ESCRT-I", "VPS4-Complex", "ESCRT-0", "ESCRT-III", "BRO1",
              "BRO1", "ESCRT-III", "BRO1", "ESCRT-III", "ESCRT-I", "ESCRT-III"),
  p_value = c(0.1621, 0.0490, 0.3578, 0.1312, NA, 0.0172, 0.2757, 0.1154,
              0.2565, 0.6541, 0.1818, 0.4225, 0.2308, 0.7987, 0.0310, 0.9484),
  p_adj = c(0.3896, 0.2451, 0.4879, 0.3896, NA, 0.2324, 0.4136, 0.3896,
            0.4136, 0.7547, 0.3896, 0.5281, 0.4136, 0.8558, 0.2324, 0.9484),
  significance = c("ns", "ns", "ns", "ns", "Not tested", "ns", "ns", "ns",
                   "ns", "ns", "ns", "ns", "ns", "ns", "ns", "ns")
)

# Run the visualization
create_all_visualizations(data)


####################################################################
# Load required packages
library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(scales)

# Set custom theme for consistent, professional appearance
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Define color schemes
state_colors <- c(
  "Nuclei" = "#2C3E50",         # Deep blue
  "Micronuclei" = "#E74C3C"     # Bright red
)

complex_colors <- c(
  "ESCRT-0" = "#FF9999",      # Soft red
  "ESCRT-I" = "#66B2FF",      # Sky blue
  "ESCRT-III" = "#99FF99",    # Soft green
  "VPS4-Complex" = "#FFCC99", # Soft orange
  "BRO1" = "#FF99FF"          # Soft purple
)

# 1. Side-by-side Expression Comparison
create_state_comparison <- function(data) {
  # Reshape data for side-by-side comparison
  long_data <- data %>%
    select(Gene, Complex, Nuclei_Mean, Micronuclei_Mean) %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(
      State = str_replace(State, "_Mean", ""),
      Expression = Expression / 1e6  # Convert to millions for readability
    )
  
  # Create plot
  ggplot(long_data, aes(x = reorder(Gene, Expression), y = Expression)) +
    geom_bar(aes(fill = State), 
             stat = "identity", 
             position = "dodge",
             width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", Expression),
                  y = Expression + 1),
              position = position_dodge(width = 0.7),
              size = 3,
              angle = 90) +
    facet_wrap(~Complex, scales = "free_y") +
    scale_fill_manual(values = state_colors,
                      name = "State",
                      labels = c("Nuclei", "Micronuclei")) +
    coord_flip() +
    custom_theme +
    labs(
      title = "Expression Levels by State and ESCRT Complex",
      subtitle = "Expression shown in millions",
      x = "Gene",
      y = "Expression Level (millions)"
    )
}

# 2. Expression Ratio Plot
create_ratio_plot <- function(data) {
  # Calculate ratio and prepare data
  data <- data %>%
    mutate(ratio = Micronuclei_Mean / Nuclei_Mean)
  
  ggplot(data, aes(x = reorder(Gene, ratio), y = ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_segment(aes(x = Gene, xend = Gene,
                     y = 1, yend = ratio,
                     color = Complex),
                 size = 1.2) +
    geom_point(aes(color = Complex), size = 4) +
    geom_text_repel(
      aes(label = sprintf("%s\n(%.2fx)", Gene, ratio)),
      size = 3,
      direction = "y"
    ) +
    scale_color_manual(values = complex_colors) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    coord_flip() +
    custom_theme +
    labs(
      title = "Expression Ratio Between States",
      subtitle = "Micronuclei / Nuclei (log scale)",
      x = "Gene",
      y = "Expression Ratio",
      color = "ESCRT Complex"
    )
}

# 3. State Distribution Plot
create_distribution_plot <- function(data) {
  # Prepare data for distribution plot
  long_data <- data %>%
    select(Gene, Complex, Nuclei_Mean, Micronuclei_Mean) %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(
      State = str_replace(State, "_Mean", ""),
      Expression = Expression / 1e6
    )
  
  ggplot(long_data, aes(x = Expression, y = Complex)) +
    geom_point(aes(color = State, shape = State),
               size = 3,
               position = position_dodge(width = 0.5)) +
    geom_text_repel(
      aes(label = Gene, color = State),
      size = 3,
      position = position_dodge(width = 0.5),
      direction = "y"
    ) +
    scale_color_manual(values = state_colors) +
    custom_theme +
    labs(
      title = "Expression Distribution by Complex and State",
      subtitle = "Expression values in millions",
      x = "Expression Level (millions)",
      y = "ESCRT Complex",
      color = "State",
      shape = "State"
    )
}

# 4. State Correlation Plot
create_correlation_plot <- function(data) {
  ggplot(data, aes(x = Nuclei_Mean/1e6, y = Micronuclei_Mean/1e6)) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "gray50") +
    geom_point(aes(color = Complex), size = 4) +
    geom_text_repel(
      aes(label = Gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(values = complex_colors) +
    custom_theme +
    labs(
      title = "Expression Correlation Between States",
      subtitle = "Values shown in millions",
      x = "Nuclei Expression (millions)",
      y = "Micronuclei Expression (millions)",
      color = "ESCRT Complex"
    )
}

# Create all plots and save them
create_state_visualizations <- function(data) {
  # Create directory for plots
  dir.create("state_comparison_plots", showWarnings = FALSE)
  
  # Create individual plots
  p1 <- create_state_comparison(data)
  p2 <- create_ratio_plot(data)
  p3 <- create_distribution_plot(data)
  p4 <- create_correlation_plot(data)
  
  # Save individual plots
  ggsave("state_comparison_plots/state_comparison.pdf", p1, width = 15, height = 10)
  ggsave("state_comparison_plots/ratio_plot.pdf", p2, width = 12, height = 8)
  ggsave("state_comparison_plots/distribution_plot.pdf", p3, width = 12, height = 8)
  ggsave("state_comparison_plots/correlation_plot.pdf", p4, width = 12, height = 8)
  
  # Create and save combined visualization
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Comprehensive State Comparison Analysis",
      subtitle = "Comparing Expression in Nuclei vs Micronuclei States",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 14, color = "gray40")
      )
    )
  
  ggsave("state_comparison_plots/combined_state_analysis.pdf", 
         combined_plot, 
         width = 24, 
         height = 20)
  
  message("State comparison plots have been saved in 'state_comparison_plots' directory")
}

# Run the visualization with the data
create_state_visualizations(data)



########################################################################################
# Load required packages for advanced visualization
library(tidyverse)
library(ggrepel)
library(patchwork)
library(viridis)
library(scales)
library(RColorBrewer)

# Create directory for plots
dir.create("advanced_publication_plots", showWarnings = FALSE)

# Define professional color schemes
pub_colors <- list(
  state = c("Nuclei" = "#4169E1", "Micronuclei" = "#DC143C"),
  complex = c(
    "ESCRT-0" = "#FF69B4",    # Hot pink
    "ESCRT-I" = "#4682B4",    # Steel blue
    "ESCRT-III" = "#32CD32",  # Lime green
    "VPS4-Complex" = "#FFD700",# Gold
    "BRO1" = "#9370DB"        # Medium purple
  )
)

# Create a professional theme for consistency
pub_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"),
    strip.text = element_text(size = 11, face = "bold")
  )

# 1. Expression Profile by State
create_expression_profile <- function(data) {
  # Calculate log2 transformed expression for better visualization
  long_data <- data %>%
    select(Gene, Complex, Nuclei_Mean, Micronuclei_Mean) %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(
      State = str_replace(State, "_Mean", ""),
      log2_expr = log2(Expression + 1)
    )
  
  ggplot(long_data, aes(x = reorder(Gene, log2_expr), y = log2_expr)) +
    geom_bar(aes(fill = Complex), stat = "identity") +
    geom_text(aes(label = sprintf("%.1f", log2_expr)),
              hjust = -0.2, size = 3, angle = 0) +
    facet_wrap(~State) +
    scale_fill_manual(values = pub_colors$complex) +
    coord_flip() +
    pub_theme +
    labs(
      title = "Gene Expression Profile by State",
      subtitle = "Log2 transformed expression values",
      x = "Gene",
      y = "log2(Expression + 1)"
    )
}

# 2. Expression Changes Between States
create_state_comparison <- function(data) {
  ggplot(data, aes(x = log2FC, y = -log10(p_adj))) +
    # Add significance regions
    annotate("rect",
             xmin = -Inf, xmax = -1,
             ymin = -log10(0.05), ymax = Inf,
             fill = "#2166AC", alpha = 0.1) +
    annotate("rect",
             xmin = 1, xmax = Inf,
             ymin = -log10(0.05), ymax = Inf,
             fill = "#B2182B", alpha = 0.1) +
    # Add threshold lines
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "gray50") +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed",
               color = "gray50") +
    # Add points and labels
    geom_point(aes(color = Complex, size = abs(log2FC))) +
    geom_text_repel(
      aes(label = Gene),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5
    ) +
    scale_color_manual(values = pub_colors$complex) +
    scale_size_continuous(range = c(2, 6)) +
    pub_theme +
    labs(
      title = "Differential Expression Analysis",
      subtitle = "Micronuclei vs Nuclei comparison",
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "ESCRT Complex",
      size = "|log2FC|"
    )
}

# 3. Complex-Specific Expression Patterns
create_complex_patterns <- function(data) {
  # Calculate mean expression by complex
  complex_data <- data %>%
    group_by(Complex) %>%
    summarize(
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      sd_log2FC = sd(log2FC, na.rm = TRUE),
      genes = n()
    ) %>%
    mutate(
      direction = ifelse(mean_log2FC >= 0, "Up", "Down"),
      label = sprintf("%s\n(n=%d)", Complex, genes)
    )
  
  ggplot(complex_data, aes(x = reorder(Complex, mean_log2FC))) +
    geom_bar(aes(y = mean_log2FC, fill = Complex),
             stat = "identity") +
    geom_errorbar(aes(ymin = mean_log2FC - sd_log2FC,
                      ymax = mean_log2FC + sd_log2FC),
                  width = 0.2) +
    geom_text(aes(y = ifelse(mean_log2FC >= 0,
                             mean_log2FC + sd_log2FC,
                             mean_log2FC - sd_log2FC),
                  label = sprintf("n=%d", genes)),
              vjust = ifelse(complex_data$mean_log2FC >= 0, -0.5, 1.5)) +
    scale_fill_manual(values = pub_colors$complex) +
    coord_flip() +
    pub_theme +
    labs(
      title = "ESCRT Complex Expression Patterns",
      subtitle = "Mean log2 Fold Change with Standard Deviation",
      x = "Complex",
      y = "Mean log2 Fold Change"
    )
}

# 4. Expression Distribution Plot
create_expression_distribution <- function(data) {
  long_data <- data %>%
    select(Gene, Complex, Nuclei_Mean, Micronuclei_Mean) %>%
    pivot_longer(
      cols = c(Nuclei_Mean, Micronuclei_Mean),
      names_to = "State",
      values_to = "Expression"
    ) %>%
    mutate(
      State = str_replace(State, "_Mean", ""),
      log10_expr = log10(Expression + 1)
    )
  
  ggplot(long_data, aes(x = log10_expr, fill = State)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~Complex, scales = "free_y") +
    scale_fill_manual(values = pub_colors$state) +
    pub_theme +
    labs(
      title = "Expression Distribution by Complex",
      subtitle = "Density plot of log10-transformed expression values",
      x = "log10(Expression + 1)",
      y = "Density"
    )
}

# Create and save all plots
create_publication_plots <- function(data) {
  # Create individual plots
  p1 <- create_expression_profile(data)
  p2 <- create_state_comparison(data)
  p3 <- create_complex_patterns(data)
  p4 <- create_expression_distribution(data)
  
  # Save individual plots
  ggsave("advanced_publication_plots/expression_profile.pdf", p1,
         width = 15, height = 10, dpi = 300)
  ggsave("advanced_publication_plots/state_comparison.pdf", p2,
         width = 12, height = 10, dpi = 300)
  ggsave("advanced_publication_plots/complex_patterns.pdf", p3,
         width = 10, height = 8, dpi = 300)
  ggsave("advanced_publication_plots/expression_distribution.pdf", p4,
         width = 15, height = 10, dpi = 300)
  
  # Create and save combined visualization
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Comprehensive ESCRT Expression Analysis",
      subtitle = "Multi-dimensional analysis of gene expression changes",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 14, color = "gray40")
      )
    )
  
  ggsave("advanced_publication_plots/combined_analysis.pdf",
         combined_plot,
         width = 20, height = 24, dpi = 300)
  
  message("Advanced publication-quality plots have been saved in 'advanced_publication_plots' directory")
}

# Run the visualization
create_publication_plots(data)
