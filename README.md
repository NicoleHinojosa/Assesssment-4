# Assesssment 4: Bioinformatics final project.

## Introduction:
In this project, we aim to systematically address a series of computational biology questions by leveraging R for data analysis and visualization. Our approach involved several key steps: first, we meticulously analyzed the problem statements and identified the necessary datasets and methods required for each point. We implemented robust coding practices to ensure that our solutions not only provide accurate answers but are also reproducible and easy to understand.

To enhance collaboration and version control, we utilized a GitHub repository where all team members contributed their expertise and enhance transparency. This README file details the purpose of each script, along with their inputs and outputs, which facilitates easier navigation and understanding for future users.

Throughout the code, we prioritized clarity by incorporating sufficient comments that elucidate the functionality of each code segment. This documentation will assist others in grasping the methodology and rationale behind our solutions, promoting knowledge sharing and educational growth within the computational biology community.

### Libraries and Their Purposes:

#### 1) R.utils

```{r}library(R.utils)
```
Purpose: Provides utility functions that enhance base R capabilities, including file handling, string manipulation, and data compression (Bengtsson, 2023). It helps streamline various operations within the R environment.

#### 2) Biostrings:

```{r}library(Biostrings)
```
Purpose: Provides efficient methods for manipulating biological sequences, particularly DNA, RNA, and protein sequences (Pagès et al., 2023). It allows for reading, writing, and performing operations on sequences with high performance.

#### 3) seqinr:

```{r}library(seqinr)
```
Purpose: Facilitates biological sequence analysis, allowing for functions to compute frequency distributions of k-mers (short subsequences). It is useful for various biological computations involving nucleotide and amino acid sequences (Charif and Lobry, 2007).

#### 4) dplyr:

```{r}library(dplyr)
```
Purpose: A powerful package for data manipulation and transformation. It provides a set of functions for filtering, arranging, and summarizing data frames, making data analysis more efficient efficient (Wickham et al., 2023).

#### 5) ggplot2:

```{r}library(ggplot2)
```
Purpose: A widely-used data visualization package that implements the Grammar of Graphics. It enables the creation of complex and customizable plots, allowing for clear data representation (Wickham and Wickham, 2016).

#### 6) readr:

```{r}library(readr)
```
Purpose: Provides functions for reading and writing data in various formats (e.g., CSV, TSV). It is optimized for speed and ease of use, making data import and export straightforward (Wickham et al., 2018).

#### 7) tidyr:

```{r}library(tidyr)
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

```{r}
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

```{r}
head(n=6, gene_expression)
View( head (n=6, gene_expression))
```

#### Inputs:
gene_expression: This is the data frame that was previously created by reading the gene expression data from a TSV file. It contains gene expression information where rows correspond to genes and columns correspond to different conditions or samples.

n=6: This argument specifies the number of rows to be displayed, in this case, the first six rows.

#### Outputs:
Console Output: The first line head(n=6, gene_expression) will output the first six rows of the gene_expression data frame directly in the R console. This is a simple text-based representation.

Data Viewer Window: The second line View(head(n=6, gene_expression)) opens a separate data viewer window in RStudio that displays the first six rows of the gene_expression data frame allowing users to easily scroll and inspect the data.
