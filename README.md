# Assesssment 4: Bioinformatics final project.

## Introduction:
In this project, we aim to systematically address a series of computational biology questions by leveraging R for data analysis and visualization. Our approach involved several key steps: first, we meticulously analyzed the problem statements and identified the necessary datasets and methods required for each point. We implemented robust coding practices to ensure that our solutions not only provide accurate answers but are also reproducible and easy to understand.

To enhance collaboration and version control, we utilized a GitHub repository where all team members contributed their expertise and enhance transparency. This README file details the purpose of each script, along with their inputs and outputs, which facilitates easier navigation and understanding for future users.

Throughout the code, we prioritized clarity by incorporating sufficient comments that elucidate the functionality of each code segment. This documentation will assist others in grasping the methodology and rationale behind our solutions, promoting knowledge sharing and educational growth within the computational biology community.

### Libraries and Their Purposes:

#### 1) R.utils

```rlibrary(R.utils)
```
Purpose: Provides utility functions that enhance base R capabilities, including file handling, string manipulation, and data compression (Bengtsson, 2023). It helps streamline various operations within the R environment.

#### 2) Biostrings:

```rlibrary(Biostrings)
```
Purpose: Provides efficient methods for manipulating biological sequences, particularly DNA, RNA, and protein sequences (Pagès et al., 2023). It allows for reading, writing, and performing operations on sequences with high performance.

#### 3) seqinr:

```rlibrary(seqinr)
```
Purpose: Facilitates biological sequence analysis, allowing for functions to compute frequency distributions of k-mers (short subsequences). It is useful for various biological computations involving nucleotide and amino acid sequences (Charif and Lobry, 2007).

#### 4) dplyr:

```rlibrary(dplyr)
```
Purpose: A powerful package for data manipulation and transformation. It provides a set of functions for filtering, arranging, and summarizing data frames, making data analysis more efficient efficient (Wickham et al., 2023).

#### 5) ggplot2:

```rlibrary(ggplot2)
```
Purpose: A widely-used data visualization package that implements the Grammar of Graphics. It enables the creation of complex and customizable plots, allowing for clear data representation (Wickham and Wickham, 2016).

#### 6) readr:

```rlibrary(readr)
```
Purpose: Provides functions for reading and writing data in various formats (e.g., CSV, TSV). It is optimized for speed and ease of use, making data import and export straightforward (Wickham et al., 2018).

#### 7) tidyr:

```rlibrary(tidyr)
```
Purpose: Focuses on tidying data, ensuring datasets are in a structured format suitable for analysis. It includes functions for reshaping and organizing data frames, enhancing the clarity of data preparation steps (Wickham and Henry, 2020).

#### Inputs:

Sequence Files: FASTA or other formats containing DNA or protein sequences.
Data Parameters: Configurable settings for analysis, such as k-mer lengths, frequency thresholds, and visualization preferences.

#### Outputs:

Frequency Tables: Data frames summarizing the counts or proportions of codons or k-mers.
Visualizations: Plots generated using ggplot2, showcasing the distribution of codon usage or other metrics.
Cleaned Data Frames: Manipulated datasets ready for further analysis or export.

## Part 1. Task 1

In this task, we aim to analyze gene expression data by reading a file containing gene identifiers and their expression values. We'll display the first six genes to understand the dataset, calculate the mean expression for each gene, identify the top 10 genes with the highest mean expression, count how many genes have a mean expression below 10, and visualize the distribution of low-expression genes using a histogram plot. This analysis will help us understand gene expression patterns and identify potential candidates for further investigation.

### Download data and read it

The following code downloads a gene expression dataset from a specified GitHub repository and then reads that dataset into R as a data frame. This allows for further analysis and processing of the gene expression data in subsequent steps.

```r
URL = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
download.file(URL, destfile = "gene_expression.tsv")

gene_expression <- read.table("gene_expression.tsv", header = TRUE, sep = "\t", row.names = 1)
```

#### Inputs: 
URL: A string that specifies the location of the dataset on GitHub. In this case, it points to a TSV (tab-separated values) file containing gene expression data.

destfile: A string that specifies the name of the file to be saved locally.

#### Outputs:
Downloaded File: The code downloads the file from the specified URL and saves it in the working directory with the name "gene_expression.tsv".

Data Frame: The code reads the downloaded TSV file into R and stores it in a variable called gene_expression. The resulting data frame will:
- Have the first row of the file treated as column headers (header = TRUE).
- Use the first column of the file as row names (row.names = 1).
- Be structured with tab (\t) as the separator between values.

### Inspecting the data frame
The purpose of this code is to display the first six rows of the gene_expression data frame. This allows users to quickly inspect the structure and content of the data, ensuring that it has been loaded correctly and giving an overview of the gene expression values.

```r
head(n=6, gene_expression)
View( head (n=6, gene_expression))
```

#### Inputs:
gene_expression: This is the data frame that was previously created by reading the gene expression data from a TSV file. It contains gene expression information where rows correspond to genes and columns correspond to different conditions or samples.

n=6: This argument specifies the number of rows to be displayed, in this case, the first six rows.

#### Outputs:
Console Output: The first line head(n=6, gene_expression) will output the first six rows of the gene_expression data frame directly in the R console. This is a simple text-based representation.

Data Viewer Window: The second line View(head(n=6, gene_expression)) opens a separate data viewer window in RStudio that displays the first six rows of the gene_expression data frame allowing users to easily scroll and inspect the data.

### Calculate Mean Expression

This code calculates the mean expression level for each gene across all samples and add it as a new column to the gene_expression data frame.

```r
gene_expression <- gene_expression %>%
  mutate(mean_expression = rowMeans(select(., everything())))
  
head(n=6, gene_expression)
```

#### Inputs:
gene_expression: The data frame containing gene expression data, with rows for genes and columns for different expression values across samples.

select(., everything()): This function selects all columns in the data frame to compute the row means.

#### Outputs:
gene_expression: The modified data frame now includes a new column mean_expression that contains the mean expression values for each gene.

Console Output: The output of head(n=6, gene_expression) displays the first six rows of the updated gene_expression data frame, showing the new mean_expression column alongside the original data.

###Identify Top 10 Genes

The top 10 genes can be identified and list using the highest mean expression values from the gene_expression data frame, showed in the following code.

```r

top_genes <- gene_expression %>%
  arrange(desc(mean_expression)) %>%
  head(10)

print(top_genes)
```

#### Inputs:
gene_expression: The data frame containing gene expression data, including the mean_expression column.

#### Outputs:
top_genes: A data frame containing the top 10 genes sorted by their mean expression values in descending order.

Console Output: The print(top_genes) command displays the details of these top 10 genes in the console.

### Count and Plot the Genes with Low Expression (Mean < 10)

The purpose of this code is to identify and list the top 10 genes with the highest mean expression values from the gene_expression data frame.

```r
num_genes_below_10 <- sum(gene_expression$mean_expression < 10)

print(num_genes_below_10)

filtered_data <- gene_expression[gene_expression$mean_expression < 10, ]
ggplot(filtered_data, aes(x = mean_expression)) +
  geom_histogram( binwidth = 1, fill = "orange", color = "black") +
  labs(title = "Histogram of Mean Low Gene Expression", x = "Mean Expression", y = "Frequency") +
  theme_minimal()

ggsave("histogram_mean_low_gene_expression.png")
```

#### Inputs:
gene_expression: The data frame containing gene expression data, including the mean_expression column.

#### Outputs:
top_genes: A data frame containing the top 10 genes sorted by their mean expression values in descending order.

Console Output: The print(top_genes) command displays the details of these top 10 genes in the console.

## Part 1. Task 2

In this task, we aim to analyze tree circumference measurements collected from two sites planted 20 years ago. We will import the data from a CSV file and examine the column names to understand the dataset structure. Next, we'll calculate the mean and standard deviation of tree circumference at both sites at the start and end of the study. We will visualize the data using a box plot to compare tree circumferences at both time points. Additionally, we'll compute the mean growth over the last 10 years for each site and use a t-test to assess whether there is a statistically significant difference in growth between the two sites. This analysis will provide insights into the effects of the treatment on tree growth.

### Read Data to perfrorm a Tree Circumference Data Analysis
The code below was used to import tree circumference growth data from a CSV file in a Github repository. It also displays the column names for further analysis.

```r
URL = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
download.file(URL, destfile = "growth_data.csv")

growth_data <- read.csv("growth_data.csv")
View(growth_data)

cat("The column names are:", colnames(growth_data))
```

#### Inputs:
URL: A string containing the link to the CSV file with tree circumference measurements.

growth_data: A data frame containing the imported tree circumference data.

#### Outputs:
Console Output: The column names of the growth_data data frame are printed to the console.

View: The growth_data data frame is opened in a viewer for inspection.


