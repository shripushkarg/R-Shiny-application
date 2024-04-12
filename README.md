# Overview

Shiny app to visualize mRNA-Seq expression data profiling Huntington's Disease vs neurologically normal individuals

Link to hosted webapp: https://shripushkar.shinyapps.io/mRNA-Seq-Shiny-Visualization/

## Tab 1: Sample Information Exploration

Summary table to give an overview of each column
Displays a table of data and histograms of numeric data

## Tab 2: Counts Matrix Exploration

Provides a summary table of counts matrix table
Scatter plot to show user input filtered genes
Clustered heatmap of counts post-filter
User selected PCA plot of first 5 PC's

## Tab 3: Differential Expression

Plot of user defined x or y axis variables with threshold colored p-adj data points
Data table of resulting values after filtering by threshold

## Tab 4: Geneset Enrichment Analysis Exploration

Normalized enrichment score (NES) plot of up and down-regulated genes
Table of results with filtering options and download feature
NES scatterplot with threshold colored p-adj data points
Files to Utilize App

## File upload for each tab:
Counts.csv <br>
DE.csv <br>
Metadata.csv

Data obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810
