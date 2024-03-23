# Install necessary packages if not already installed
# install.packages(c("shiny", "DT", "ggplot2", "viridis"))

library(shiny)
library(DT)
library(ggplot2)
library(viridis)
library(tidyverse)
library(ggbeeswarm)
library(bslib)
library(dplyr)
library(colourpicker)

options(shiny.maxRequestSize = 60*1024^2)

# File Input Validation Function
validate_file <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("File not found.")
  }
  
  # Check if the file extension is either CSV or TSV
  file_extension <- tools::file_ext(file_path)
  if (!(file_extension %in% c("csv", "tsv"))) {
    stop("Invalid file format. Only CSV or TSV files are allowed.")
  }
  
  # Additional checks for well-formatted file can be added if needed
}

# Shiny App UI
ui <- fluidPage(
  theme = bslib::bs_theme(
    background = "lightgray",  # Change to your desired background color
    primary = "darkblue"  # Change to your desired primary color
  ),
  titlePanel("BF591 - Project"), ("Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Samples",
                 fileInput("sample_file", "Choose CSV File for Samples", accept = c(".csv")),
                 tags$hr(),
                 helpText("Maximum file size: 30MB")
        ),
        tabPanel("Counts",
                 fileInput("count_file", "Choose CSV File for Counts", accept = c(".csv")),
                 tags$hr(),
                 helpText("Maximum file size: 30MB")
        ),
        tabPanel("DE",
                 fileInput("de_file", "Choose CSV File for DE", accept = c(".csv")),
                 tags$hr(),
                 helpText("Maximum file size: 30MB")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Samples",
                 tabsetPanel(
                   tabPanel("Summary", DTOutput("summary_table")),
                   tabPanel("Sample Table", DTOutput("sample_table")),
                   tabPanel("Density Plot",
                            selectInput("variable", "Select Variable for Density Plot", ""),
                            plotOutput("density_plot"))
                 )
        ),
        tabPanel("Counts",
                 sliderInput("variance_threshold", "Variance Threshold",
                             min = 0, max = 100, value = 80,
                             step = 1, post = "%"),
                 sliderInput("nonzero_samples_threshold", "Non-Zero Samples Threshold",
                             min = 0, max = 100, value = 60,
                             step = 1, post = "%"),
                 tabsetPanel(
                   tabPanel("Filter Summary", DTOutput("filter_summary")),
                   tabPanel("Diagnostic Plots",
                            plotOutput("median_vs_variance"),
                            plotOutput("median_vs_zeros")),
                   tabPanel("Clustered Heatmap", plotOutput("clustered_heatmap")),
                   tabPanel("PCA Scatter Plot",
                            numericInput("num_pcs", "Number of Principal Components", value = 2, min = 2),
                            plotOutput("pca_scatter_plot"))
                 )
        ),
        tabPanel("DE",
                 tabsetPanel(
                   tabPanel("DE Table", DTOutput("de_table")),
                   tabPanel("Volcano Plot",
                            radioButtons("x_variable", "Choose the column for X-axis",
                                         choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                         selected = "log2FoldChange"),
                            radioButtons("y_variable", "Choose the column for Y-axis",
                                         choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                         selected = "padj"),
                            colourInput("color1", "Base Point Color", value = "#FFCF56"),
                            colourInput("color2", "Highlight Point Color", value = "#22577A"),
                            sliderInput("magnitude_slider", "Select the magnitude of the p adjusted coloring:",
                                        min = -35, max = 0, value = -8),
                            actionButton("plot_button", "Plot"),
                            plotOutput("volcano_plot"),
                            plotOutput("histogram_value_counts"),
                            plotOutput("histogram_log2foldchange_counts")),
                   tabPanel("Filtered Data",
                            DTOutput("filtered_data_table")
                   ),
                   tabPanel("Histogram",
                            plotOutput("histogram_log2foldchange")
                   )
                 )
        ),
        tabPanel("Gene Expression Visualization",
                 fileInput("counts_file", "Choose Gene Counts CSV File", accept = c(".csv")),
                 fileInput("info_file", "Choose Sample Information CSV File", accept = c(".csv")),
                 selectInput("categorical_var", "Select Categorical Variable", ""),
                 selectInput("gene_var", "Select Gene", ""),
                 selectInput("plot_type", "Select Plot Type",
                             choices = c("Bar Plot", "Boxplot", "Violin Plot", "Beeswarm Plot")),
                 actionButton("go_button", "Generate Plot"),
                 plotOutput("gene_plot"),  # Add this line to display the plot
                 tags$hr(),
                 helpText("Maximum file size: 30MB")
        )
      )
    )
  )
)

# Shiny App Server
server <- function(input, output, session) {
  # Read the sample information matrix
  samples_data <- reactive({
    req(input$sample_file)
    validate_file(input$sample_file$datapath)
    read.csv(input$sample_file$datapath, stringsAsFactors = FALSE)  # Ensure strings are not converted to factors
  })
  
  # Update the variable selection dropdown based on the uploaded data
  observe({
    req(samples_data())
    
    # Include only the last four columns in the choices
    variable_choices <- tail(names(samples_data()), 4)
    
    updateSelectInput(session, "variable", choices = variable_choices)
  })
  
  # Display summary of the table
  output$summary_table <- renderDT({
    summary_data <- data.frame(
      Column_Name = character(),
      Type = character(),
      Summary = character(),
      stringsAsFactors = FALSE
    )
    
    for (col_name in names(samples_data())) {
      col_type <- class(samples_data()[[col_name]])
      col_summary <- if (col_type %in% c("integer", "numeric")) {
        mean_val <- mean(samples_data()[[col_name]], na.rm = TRUE)
        sd_val <- sd(samples_data()[[col_name]], na.rm = TRUE)
        paste0(sprintf("%.2f", mean_val), " (+/- ", sprintf("%.2f", sd_val), ")")
      } else if (col_type == "factor") {
        levels_str <- paste(levels(samples_data()[[col_name]]), collapse = ", ")
        paste0(levels_str)
      } else if (col_type == "character") {
        if (col_name %in% c("tissue", "diagnosis")) {
          col_summary <- paste(unique(samples_data()[[col_name]]), collapse = ", ")
        } else {
          col_summary <- paste0("Values: ", paste(unique(samples_data()[[col_name]]), collapse = ", "))
        }
      } else {
        col_summary <- "N/A"
      }
      
      summary_data <- rbind(
        summary_data,
        data.frame(
          Column_Name = col_name,
          Type = col_type,
          Summary = col_summary,
          stringsAsFactors = FALSE
        )
      )
    }
    
    datatable(summary_data, options = list(pageLength = 10))
  })
  
  # Display the sample table
  output$sample_table <- renderDT({
    datatable(samples_data(), options = list(pageLength = 10))
  })
  
  # Display density plot based on selected variable
  output$density_plot <- renderPlot({
    req(samples_data(), input$variable)
    
    # Assuming the selected column for the density plot
    selected_variable <- input$variable
    variable_values <- as.numeric(samples_data()[[selected_variable]])
    
    # Check if the variable values are numeric
    if (!is.numeric(variable_values)) {
      return(plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "Invalid Data Type"))
    }
    
    # Add a small constant to all values before log transformation
    constant <- 1e-6
    log_variable_values <- log(variable_values + constant)
    
    # Print some information for debugging
    print(paste("Selected Variable:", selected_variable))
    print(paste("Min Log-Transformed Value:", min(log_variable_values)))
    print(paste("Max Log-Transformed Value:", max(log_variable_values)))
    
    ggplot(data.frame(Variable_Values = log_variable_values), aes(x = Variable_Values)) +
      geom_density(fill = "blue", color = "black") +
      labs(title = paste("Density Plot of Log-Transformed", selected_variable))
  })
  
  # Read the counts matrix
  counts_data <- reactive({
    req(input$count_file)
    validate_file(input$count_file$datapath)
    read.csv(input$count_file$datapath, stringsAsFactors = FALSE, row.names = 1)  # Assuming gene names are in the first column
  })
  
  # Filter genes based on variance and non-zero samples thresholds
  filtered_counts_data <- reactive({
    req(counts_data())
    
    variance_threshold <- quantile(apply(counts_data(), 1, var), input$variance_threshold / 100)
    non_zero_samples_threshold <- quantile(apply(counts_data() > 0, 1, sum), input$nonzero_samples_threshold / 100)
    
    filtered_genes <- rownames(counts_data())[apply(counts_data(), 1, function(x) var(x) >= variance_threshold & sum(x > 0) >= non_zero_samples_threshold)]
    
    counts_data()[filtered_genes, , drop = FALSE]
  })
  
  # Display filter summary table
  output$filter_summary <- renderDT({
    total_genes <- nrow(counts_data())
    total_samples <- ncol(counts_data())
    filtered_genes <- nrow(filtered_counts_data())
    remaining_genes <- total_genes - filtered_genes
    
    percentage_passed <- paste0(round(filtered_genes / total_genes * 100, 1), "%")
    percentage_not_passed <- paste0(round(remaining_genes / total_genes * 100, 1), "%")
    
    filter_summary_data <- data.frame(
      Metric = c("Total Genes", "Total Samples", "Filtered Genes", "Remaining Genes",
                 "Percentage Genes Passed", "Percentage Genes Not Passed"),
      Value = c(total_genes, total_samples, filtered_genes, remaining_genes,
                percentage_passed, percentage_not_passed)
    )
    
    datatable(filter_summary_data)
  })
  
  # Diagnostic scatter plots
  output$median_vs_variance <- renderPlot({
    # Plot median count vs variance
    plot(log(apply(counts_data(), 1, median) + 1), log(apply(counts_data(), 1, var) + 1),
         col = ifelse(rownames(counts_data()) %in% rownames(filtered_counts_data()), "blue", "red"),
         pch = 16, cex = 0.7,
         xlab = "Log(Median Count)", ylab = "Log(Variance)",
         main = "Log(Median Count) vs Log(Variance)")
    
    legend("topright", legend = c("Passed Filter", "Not Passed Filter"), col = c("blue", "red"), pch = 16)
  })
  
  output$median_vs_zeros <- renderPlot({
    # Plot median count vs number of zeros
    plot(log(apply(counts_data(), 1, median) + 1), log(apply(counts_data() == 0, 1, sum) + 1),
         col = ifelse(rownames(counts_data()) %in% rownames(filtered_counts_data()), "blue", "red"),
         pch = 16, cex = 0.7,
         xlab = "Log(Median Count)", ylab = "Log(Number of Zeros)",
         main = "Log(Median Count) vs Log(Number of Zeros)")
    
    legend("topright", legend = c("Passed Filter", "Not Passed Filter"), col = c("blue", "red"), pch = 16)
  })
  
  # Clustered heatmap
  output$clustered_heatmap <- renderPlot({
    # Plot clustered heatmap of counts matrix
    heatmap(log1p(as.matrix(filtered_counts_data())), scale = "row", col = viridis::viridis(256), main = "Clustered Heatmap",
            xlab = "Samples", ylab = "Genes", key.title = "Log(counts + 1)", key.xlab = "Expression")
  })
  
  # PCA scatter plot
  output$pca_scatter_plot <- renderPlot({
    # Perform PCA on the filtered counts matrix
    pca_result <- prcomp(t(log1p(filtered_counts_data())))
    
    # Plot PCA scatter plot
    plot(pca_result$x[, 1], pca_result$x[, 2], col = "blue", pch = 16, cex = 0.7,
         xlab = "PC1", ylab = "PC2",
         main = "PCA Scatter Plot")
  })
  
  # Read the DE data
  de_data <- reactive({
    req(input$de_file)
    validate_file(input$de_file$datapath)
    read.csv(input$de_file$datapath)
  })
  
  # Display the DE table
  output$de_table <- renderDT({
    datatable(de_data(), options = list(pageLength = 10, order = list(2, 'asc')))
  })
  
  
  # Reactive function for generating the volcano plot
  generate_volcano_plot <- reactive({
    req(de_data(), input$x_variable, input$y_variable, input$magnitude_slider, input$color1, input$color2)
    
    dataf <- de_data()
    x_name <- input$x_variable
    y_name <- input$y_variable
    slider <- input$magnitude_slider
    color1 <- input$color1
    color2 <- input$color2
    
    plot <- ggplot(dataf, aes_string(x = x_name, y = sprintf("-log10(%s) + 1e-10", y_name), color = "factor(padj < 10^slider)")) +
      geom_point(size = 1, alpha = 0.7) +
      scale_color_manual(values = c(color2, color1), name = sprintf("padj < 1 × 10^%d", slider)) +
      labs(x = x_name, y = sprintf("-log10(%s)", y_name)) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.direction = "horizontal")  
    
    return(plot)
  })
  
  # Display volcano plot
  output$volcano_plot <- renderPlot({
    generate_volcano_plot()
  })
  
  # Reactive function for generating filtered data table
  generate_filtered_data_table <- reactive({
    req(de_data(), input$magnitude_slider)
    
    # Filter data based on the padj threshold
    filtered_data <- de_data()[de_data()$padj < 10^input$magnitude_slider, ]
    
    return(filtered_data)
  })
  
  # Display filtered data table
  output$filtered_data_table <- renderDT({
    datatable(generate_filtered_data_table())
  })
  
  # Display histogram for value counts
  output$histogram_value_counts <- renderPlot({
    req(de_data())
    
    # Your code for generating the histogram for value counts goes here
    # Example: ggplot(data = de_data(), aes(x = value)) + geom_histogram()
  })
  
  # Display histogram for log2foldchange counts
  output$histogram_log2foldchange_counts <- renderPlot({
    req(de_data())
    
    # Your code for generating the histogram for log2foldchange counts goes here
    # Example: ggplot(data = de_data(), aes(x = log2FoldChange)) + geom_histogram()
  })
  
  # Display histogram for log2foldchange
  output$histogram_log2foldchange <- renderPlot({
    req(de_data())
    
    # Actual code for generating the histogram for log2foldchange
    ggplot(data = de_data(), aes(x = log2FoldChange)) + 
      geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +
      labs(title = "Log2FoldChange Histogram", x = "log2FoldChange", y = "Frequency")
  })
  
  # Read the gene counts matrix
  counts_data2 <- reactive({
    req(input$counts_file)
    validate_file(input$counts_file$datapath)
    
    # Read the counts file and transpose it
    counts_matrix <- read.csv(input$counts_file$datapath, stringsAsFactors = FALSE, row.names = 1)
    counts_matrix_transposed <- t(counts_matrix)
    
    # Convert row names to a new column named "Sample"
    counts_df <- data.frame(Sample = row.names(counts_matrix_transposed), counts_matrix_transposed)
    
    # Reset row names
    row.names(counts_df) <- NULL
    
    counts_df
  })
  
  # Read the sample information matrix
  info_data <- reactive({
    req(input$info_file)
    validate_file(input$info_file$datapath)
    read.csv(input$info_file$datapath, stringsAsFactors = FALSE)
  })
  
  # Update the categorical variable selection dropdown based on the sample information data
  observe({
    req(info_data())
    variable_choices <- names(info_data())
    updateSelectInput(session, "categorical_var", choices = variable_choices[-1])  # Exclude the first column
  })
  
  # Update the gene selection dropdown based on the gene counts data
  observe({
    req(counts_data2())
    gene_choices <- names(counts_data2())[2:length(names(counts_data2()))]  # Exclude the first column (Sample)
    updateSelectInput(session, "gene_var", choices = gene_choices)
  })
  
  # Generate and display the gene expression plot
  output$gene_plot <- renderPlot({
    req(counts_data2(), info_data(), input$gene_var, input$categorical_var, input$plot_type)
    
    # Assuming the selected column for gene counts
    selected_gene <- input$gene_var
    selected_category <- input$categorical_var
    plot_type <- input$plot_type
    
    # Prepare data for plotting
    data_to_plot <- data.frame(
      Sample = counts_data2()$Sample,
      Expression = counts_data2()[[selected_gene]],
      Category = info_data()[[selected_category]]
    )
    
    # Generate plot based on selected plot type
    if (plot_type == "Bar Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_bar(stat = "identity", fill = "blue") +
        labs(title = paste("Bar Plot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Boxplot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_boxplot(fill = "blue") +
        labs(title = paste("Boxplot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Violin Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_violin(fill = "blue") +
        labs(title = paste("Violin Plot of", selected_gene, "by", selected_category))
    } else if (plot_type == "Beeswarm Plot") {
      ggplot(data_to_plot, aes(x = Category, y = Expression)) +
        geom_beeswarm(fill = "blue") +
        labs(title = paste("Beeswarm Plot of", selected_gene, "by", selected_category))
    }
  })
}

# Run the Shiny App
shinyApp(ui, server)
