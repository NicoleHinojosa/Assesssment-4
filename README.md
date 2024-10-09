# Assesssment 4: Bioinformatics final project.
## Group 7: Nicole Hinojosa, Saqeena Thapa & Vanshika

## Introduction:

In this project, we aim to systematically address a series of computational biology questions by leveraging R for data analysis and visualization. Our approach involved several key steps: first, we meticulously analyzed the problem statements and identified the necessary datasets and methods required for each point. We implemented robust coding practices to ensure that our solutions not only provide accurate answers but are also reproducible and easy to understand.

To enhance collaboration and version control, we utilized a GitHub repository where all team members contributed their expertise and enhance transparency. This README file details the purpose of each script, along with their inputs and outputs, which facilitates easier navigation and understanding for future users.

Throughout the code, we prioritized clarity by incorporating sufficient comments that elucidate the functionality of each code segment. This documentation will assist others in grasping the methodology and rationale behind our solutions, promoting knowledge sharing and educational growth within the computational biology community.

### Libraries and Their Purposes:

#### 1) R.utils

```r
library(R.utils) #use BiocManager::install("Biostrings") if it is not already installed in your Rstudio
```
Purpose: Provides utility functions that enhance base R capabilities, including file handling, string manipulation, and data compression (Bengtsson, 2023). It helps streamline various operations within the R environment.

#### 2) Biostrings:

```r
library(Biostrings)
```
Purpose: Provides efficient methods for manipulating biological sequences, particularly DNA, RNA, and protein sequences (Pagès et al., 2023). It allows for reading, writing, and performing operations on sequences with high performance.

#### 3) seqinr:

```r
library(seqinr)
```
Purpose: Facilitates biological sequence analysis, allowing for functions to compute frequency distributions of k-mers (short subsequences). It is useful for various biological computations involving nucleotide and amino acid sequences (Charif and Lobry, 2007).

#### 4) dplyr:

```r
library(dplyr)
```
Purpose: A powerful package for data manipulation and transformation. It provides a set of functions for filtering, arranging, and summarizing data frames, making data analysis more efficient efficient (Wickham et al., 2023).

#### 5) ggplot2:

```r
library(ggplot2)
```
Purpose: A widely-used data visualization package that implements the Grammar of Graphics. It enables the creation of complex and customizable plots, allowing for clear data representation (Wickham and Wickham, 2016).

#### 6) readr:

```r
library(readr)
```
Purpose: Provides functions for reading and writing data in various formats (e.g., CSV, TSV). It is optimized for speed and ease of use, making data import and export straightforward (Wickham et al., 2018).

#### 7) tidyr:

```r
library(tidyr)
```
Purpose: Focuses on tidying data, ensuring datasets are in a structured format suitable for analysis. It includes functions for reshaping and organizing data frames, enhancing the clarity of data preparation steps (Wickham and Henry, 2020).

#### 8) knitr:

```r
library(knitr)
```
Purpose: Allows dynamic report generation, allowing users to integrate R code with text to create high-quality documents in various formats like HTML, PDF, and Word.

#### Inputs:

Sequence Files: FASTA or other formats containing DNA or protein sequences.
Data Parameters: Configurable settings for analysis, such as k-mer lengths, frequency thresholds, and visualization preferences.

#### Outputs:

Frequency Tables: Data frames summarizing the counts or proportions of codons or k-mers.
Visualizations: Plots generated using ggplot2, showcasing the distribution of codon usage or other metrics.
Cleaned Data Frames: Manipulated datasets ready for further analysis or export.

### Images format

```r
knitr::opts_chunk$set(fig.width=7, fig.height=5, out.width="80%")
```
Purpose: this code sets global options for all subsequent R code chunks in the document.

## Part 1. Task 1

In this task, we aim to analyze gene expression data by reading a file containing gene identifiers and their expression values. We'll display the first six genes to understand the dataset, calculate the mean expression for each gene, identify the top 10 genes with the highest mean expression, count how many genes have a mean expression below 10, and visualize the distribution of low-expression genes using a histogram plot. This analysis will help us understand gene expression patterns and identify potential candidates for further investigation.


### Download data and read it

The following code downloads a gene expression dataset from a specified GitHub repository and then reads that dataset into R as a data frame. This allows for further analysis and processing of the gene expression data in subsequent steps.

```r
#Download the data from the github link provided
URL = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv"
download.file(URL, destfile = "gene_expression.tsv")

# Read the downloaded TSV file into R
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
# View the first few rows of the data
gene_expression_show <- head(n=6, gene_expression)
#Give format to the table
kable (gene_expression_show, align = 'c', caption = "First 6 genes of gene_expression file")
```

#### Inputs:
gene_expression: This is the data frame that was previously created by reading the gene expression data from a TSV file. It contains gene expression information where rows correspond to genes and columns correspond to different conditions or samples.

n=6: This argument specifies the number of rows to be displayed, in this case, the first six rows.

#### Outputs:
Console Output: The variable gene_expression_show will hold the first six rows of the gene_expression data frame, which will be printed to the console when called.

Formatted Table: The kable() function formats the gene_expression_show data frame into a neatly structured table with the following attributes:

  -Alignment: Centered columns (align = 'c').
  -Caption: "First 6 genes of gene_expression file", providing context to the displayed data.
  

### Calculate Mean Expression

This code calculates the mean expression level for each gene across all samples and add it as a new column to the gene_expression data frame.

```r
# Calculate the mean across the samples and add as a new column
gene_expression <- gene_expression %>%
  mutate(mean_expression = rowMeans(select(., everything())))
# Show a table of values for the first six genes including the mean
gene_expression_show <-head(n=6, gene_expression)
# Give format to a table
kable (gene_expression_show, align = 'c', caption = "First 6 genes of gene_expression file")
```

#### Inputs:
gene_expression: The data frame containing gene expression data, with rows for genes and columns for different expression values across samples.

select(., everything()): This function selects all columns in the data frame to compute the row means.

#### Outputs:
Updated gene_expression Data Frame: The original data frame is modified to include a new column named mean_expression, which contains the mean expression value for each gene calculated across all samples.

Console Output: The variable gene_expression_show holds the first six rows of the updated gene_expression data frame, which will be printed to the console when called.

Formatted Table: The kable() function creates a neatly formatted table displaying the first six genes along with their mean expression values. This table includes:

  -Alignment: Centered columns (align = 'c').
  -Caption: "First 6 genes of gene_expression file", providing context to the displayed data.
  
  
###Identify Top 10 Genes

The top 10 genes can be identified and list using the highest mean expression values from the gene_expression data frame, showed in the following code.

```r
# List the 10 genes with the highest mean expression
top_genes <- gene_expression %>%
  arrange(desc(mean_expression)) %>%
  head (n=10) 
# Print the top genes in a table
kable (top_genes, align = 'c', caption = "Top 10 expressed genes")
```

#### Inputs:
gene_expression: The data frame containing gene expression data, including the mean_expression column.

#### Outputs:
Updated top_genes Data Frame: A data frame containing the top ten genes sorted by their mean expression values in descending order.

Formatted Table: The kable() function generates a neatly formatted table displaying the top ten expressed genes, including:

  -Alignment: Centered columns (align = 'c').
  -Caption: "Top 10 expressed genes", providing context for the data presented.
  

### Count and Plot the Genes with Low Expression (Mean < 10)

This code determines the number of genes with a mean expression value less than 10, prints this count, and creates a histogram to visualize the distribution of these low mean expression values. It also saves the histogram plot to a file for inclusion in a report.

```r
# Determine the number of genes with a mean < 10
num_genes_below_10 <- sum(gene_expression$mean_expression < 10)

# Print the number of genes
cat(format(num_genes_below_10), "\n")

# Make a histogram plot of the mean values
filtered_data <- gene_expression[gene_expression$mean_expression < 10, ]
ggplot(filtered_data, aes(x = mean_expression)) +
  geom_histogram( binwidth = 1, fill = "orange", color = "black") +
  labs(title = "Histogram of Mean Low Gene Expression", x = "Mean Expression", y = "Frequency") +
  theme_minimal()

# Save the plot to your report
ggsave("histogram_mean_low_gene_expression.png")
```

#### Inputs:
gene_expression: The data frame containing gene expression data, including the mean_expression column.

#### Outputs:
num_genes_below_10: An integer value that counts the number of genes with a mean expression less than 10.

Console Output: The count of genes below the threshold is printed to the console in a formatted manner.

Histogram Plot: A histogram that visualizes the frequency of mean expression values for genes with mean expression less than 10.

Saved Plot: The histogram plot is saved as a PNG file named "histogram_mean_low_gene_expression.png" in the current working directory.

## Part 1. Task 2

In this task, we aim to analyze tree circumference measurements collected from two sites planted 20 years ago. We will import the data from a CSV file and examine the column names to understand the dataset structure. Next, we'll calculate the mean and standard deviation of tree circumference at both sites at the start and end of the study. We will visualize the data using a box plot to compare tree circumferences at both time points. Additionally, we'll compute the mean growth over the last 10 years for each site and use a t-test to assess whether there is a statistically significant difference in growth between the two sites. This analysis will provide insights into the effects of the treatment on tree growth.

### Read Data to perfrorm a Tree Circumference Data Analysis

The code below was used to import tree circumference growth data from a CSV file in a Github repository. It also displays the column names for further analysis.

```r
# Read in the growth data
#Download the data from the github link provided
URL = "https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv"
download.file(URL, destfile = "growth_data.csv")

# Read the downloaded TSV file into R
growth_data <- read.csv("growth_data.csv")

#print the table with the 6 first rows
growth_data_show <- head(growth_data)
kable (growth_data_show, align = 'c', caption = "First 6 rows of growth_data file")

# Show column names
cat("The column names are:")
print(colnames(growth_data))
```

#### Inputs:
URL: A string containing the link to the CSV file on GitHub.

destfile: A string specifying the name of the file where the downloaded data will be saved locally ("growth_data.csv").

#### Outputs:
growth_data: A data frame containing the growth data read from the CSV file.

growth_data_show: A table displaying the first six rows of the growth_data data frame, formatted using kable() with a caption "First 6 rows of growth_data file".

Console Output: The message "The column names are:" followed by the names of the columns in the growth_data data frame printed to the console.


### Calculate the Mean and Standard Deviation Statistics

The purpose of this code is to calculate and summarize the mean and standard deviation of tree circumference measurements at two different sites (southwest and northeast) for the years 2005 and 2020.

```r
# Calculate mean and standard deviation for tree circumference
summary_stats <- growth_data %>%
  summarise(mean_start_2005_southwest = mean(Circumf_2005_cm[Site == "southwest"]),
    sd_start_2005_southwest = sd(Circumf_2005_cm[Site == "southwest"]),
    mean_start_2005_northeast = mean(Circumf_2005_cm[Site == "northeast"]),
    sd_start_2005_northeast = sd(Circumf_2005_cm[Site == "northeast"]),
    mean_end_2020_southwest = mean(Circumf_2020_cm[Site == "southwest"]),
    sd_end_2020_southwest = sd(Circumf_2020_cm[Site == "southwest"]),
    mean_end_2020_northeast = mean(Circumf_2020_cm[Site == "northeast"]),
    sd_end_2020_northeast = sd(Circumf_2020_cm[Site == "northeast"])
  )

# Print summary statistics
cat(paste(capture.output(print(summary_stats)), collapse = "\n"), "\n")
```

#### Inputs:
growth_data: A data frame containing tree circumference measurements across different sites and years.

#### Outputs:
summary_stats: A data frame containing the calculated mean and standard deviation for tree circumference at both sites for the years 2005 and 2020.

Console Output: The summary statistics printed to the console, formatted as a table.


### Plot the Tree Circumferences Comparisson between 2005 and 2020 

The purpose of this code is to create a box plot that visually compares the tree circumference measurements for the years 2005 and 2020 across two different sites.

```r
# Reshape data from wide to long format for Circumf_2005_cm and Circumf_2020_cm
long_data <- growth_data %>%
  select(Site, TreeID, Circumf_2005_cm, Circumf_2020_cm) %>%
  pivot_longer(cols = starts_with("Circumf"), 
               names_to = "Year", 
               values_to = "Circumference")

# Filter for only the start and end years
long_data <- long_data %>%
  filter(Year %in% c("Circumf_2005_cm", "Circumf_2020_cm"))

# Create a box plot
ggplot(long_data, aes(x = Year, y = Circumference, fill = Site)) +
  geom_boxplot() +
  labs(title = "Tree Circumference at Start (2005) and End (2020) of Study", 
       x = "Year", 
       y = "Circumference (cm)") +
  scale_fill_manual(values = c("northeast" = "red", "southwest" = "pink")) +
  theme_minimal()

# Save the box plot to your report
ggsave("boxplot_tree_circumference.png")
```

#### Inputs:
growth_data: A data frame containing tree circumference measurements, site information, and tree identifiers.

#### Outputs:
long_data: A reshaped data frame in long format containing site, year, and circumference values.

Box Plot: A box plot visualizing the distribution of tree circumferences for 2005 and 2020, differentiated by site.

Saved Plot: The box plot is saved as "boxplot_tree_circumference.png".

### Mean Tree Growth Calculation in 10 years

The following code calculates the growth of the trees over the last 10 years for each site and to compute the mean growth for each site.

```r
# Calculate growth over the last 10 years for each tree
growth_data <- growth_data %>%
  mutate(Growth_10_years = Circumf_2020_cm - Circumf_2010_cm)

# Calculate mean growth at each site
mean_growth <- growth_data %>%
  group_by(Site) %>%
  summarise(mean_growth = mean(Growth_10_years, na.rm = TRUE), .groups = 'drop')

# Show the mean growth
kable(mean_growth, align='c', caption = "Mean Growth Calculation")
```

#### Inputs:
growth_data: A data frame containing tree circumference measurements for the years 2010 and 2020, along with site information.

#### Outputs:
Updated growth_data: The original data frame augmented with a new column, Growth_10_years, representing the growth of each tree over the last 10 years.

mean_growth: A data frame summarizing the mean growth for each site.

Formatted Table: The mean growth data is displayed in a table format using kable(), providing a clear and organized view of the results.


### T-Test for Growth Difference

The purpose of this code is to perform a t-test to compare the growth of trees between the two different sites.
The t-test is a statistical method used to determine if there is a significant difference between the means of two groups. In this case, used to determine if the difference in mean growth between the two sites (control and treatment) is statistically significant.

```r
# Perform t-test to compare growth between sites
t_test_result <- t.test(Growth_10_years ~ Site, data = growth_data)

# Print t-test results
t_test_result

```

#### Inputs:
growth_data: A data frame containing tree growth data, including the Growth_10_years variable and Site classification.

#### Outputs:
t_test_result: An object containing the results of the t-test, including the t-statistic, degrees of freedom, p-value, and confidence interval.

Console Output: The t-test results are printed to the console, providing insights into whether there is a statistically significant difference in growth between the sites.



## Part 2: Examining Biological Sequence Diversity of *Streptacidiphilus jiangxiensis* and comparing it with *Escherichia Coli*

In this task, we aim to analyze and compare the sequence features of *Streptacidiphilus jiangxiensis* with those of *Escherichia coli*. We will begin by downloading and examining the coding DNA sequences of both organisms to determine the number of coding sequences and their total length, presenting our findings in tables. We will also analyze the length of coding sequences using boxplots and calculate the mean and median lengths for each organism. Further, we will assess the frequency of DNA bases and amino acids, visualizing these frequencies through bar plots. Additionally, we will create a codon usage table to quantify codon usage bias and identify k-mers in the protein sequences, comparing their representation between the two organisms. Through this comprehensive analysis, we will highlight the differences in sequence features, providing insights into the evolutionary and functional implications of these variations.


### Sequences

This code is used to download and read the coding DNA sequences for *Escherichia coli* and Streptacidiphilus jiangxiensis from specified URLs, decompress the files, and load the sequences into R for further analysis.

```r
# URLs for the coding DNA sequences
URL_Ecoli <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_117_collection/escherichia_coli_110957_gca_000485615/cds/Escherichia_coli_110957_gca_000485615.ASM48561v1.cds.all.fa.gz"
URL_Streptacidiphilus <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_57_collection/streptacidiphilus_jiangxiensis_gca_900109465/cds/Streptacidiphilus_jiangxiensis_gca_900109465.IMG-taxon_2675903135_annotated_assembly.cds.all.fa.gz"

# Downloading the sequences
download.file(URL_Ecoli, destfile = "e_coli_cds.fa.gz")
download.file(URL_Streptacidiphilus, destfile = "streptacidiphilus_cds.fa.gz")

#Decompress the files

gunzip("e_coli_cds.fa.gz")
gunzip("streptacidiphilus_cds.fa.gz")

# Reading the sequences
ecoli_seqs <- seqinr::read.fasta ("e_coli_cds.fa")
streptacidiphilus_seqs <- seqinr::read.fasta ("streptacidiphilus_cds.fa")
```

####: Inputs
URL_Ecoli: URL for the *Escherichia coli* coding sequences (FASTA format, gzipped).

URL_Streptacidiphilus: URL for the *Streptacidiphilus jiangxiensis* coding sequences (FASTA format, gzipped).

#### Outputs:
ecoli_seqs: A list of coding DNA sequences for *Escherichia coli*.

streptacidiphilus_seqs: A list of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

Decompressed files: e_coli_cds.fa and streptacidiphilus_cds.fa in the working directory.


### CDS count for both organisms

The purpose of this code is to count the number of coding sequences in the *Escherichia coli* and *Streptacidiphilus jiangxiensis* datasets and create a summary table of these counts.

```r
# Count coding sequences
ecoli_count <- length (ecoli_seqs)
streptacidiphilus_count <- length (streptacidiphilus_seqs)

# Creating a summary table
coding_counts <- data.frame(
  Organism = c("Escherichia coli", "Streptacidiphilus jiangxiensis"),
  Coding_Sequences = c(ecoli_count, streptacidiphilus_count)
)
kable (coding_counts, align = 'c', caption = "Coding Sequences of each organism")
```

#### Inputs:
ecoli_seqs: A list of coding DNA sequences for *Escherichia coli*.

streptacidiphilus_seqs: A list of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
ecoli_count: The total number of coding sequences in *Escherichia coli*.

streptacidiphilus_count: The total number of coding sequences in *Streptacidiphilus jiangxiensis*.

coding_counts: A data frame summarizing the counts of coding sequences for both organisms.

Table Output: A formatted table displayed using kable, showing the number of coding sequences for each organism with a caption "Coding Sequences of each organism."


### Total Coding DNA Length of both organisms

The code below calculates the total length of coding DNA sequences for *Escherichia coli* and *Streptacidiphilus jiangxiensis*, and to create a summary table displaying these lengths.

```r
# Calculate total coding DNA length
ecoli_length <- as.numeric(summary(ecoli_seqs)[,1])
streptacidiphilus_length <- as.numeric(summary(streptacidiphilus_seqs)[,1])

# Creating a summary table
total_lengths <- data.frame(
  Organism = c("Escherichia coli", "Streptacidiphilus jiangxiensis"),
  Total_Length = c(sum(ecoli_length), sum(streptacidiphilus_length))
)
# Give format
kable (total_lengths, align = 'c', caption = "Total length of DNA sequences from each organism")
```

#### Inputs:
ecoli_seqs: A list of coding DNA sequences for *Escherichia coli*.

streptacidiphilus_seqs: A list of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
ecoli_length: A numeric vector containing the lengths of coding sequences for *Escherichia coli*.

streptacidiphilus_length: A numeric vector containing the lengths of coding sequences for *Streptacidiphilus jiangxiensis*.

total_lengths: A data frame summarizing the total lengths of DNA sequences for both organisms, with columns for the organism names and their respective total lengths.

Table Output: A formatted table displayed using kable, showing the total length of DNA sequences for each organism with a caption "Total length of DNA sequences from each organism."


### Coding Sequence Length Distribution

The purpose of this code is to visualize the distribution of coding sequence lengths for *Escherichia coli* and *Streptacidiphilus jiangxiensis* using a boxplot.

```r
# Boxplot of coding sequence lengths

boxplot(list(`Escherichia coli` = ecoli_length, 
             `Streptacidiphilus jiangxiensis` = streptacidiphilus_length),
        col = "yellow",
        xlab = "Organism",
        ylab = "Length (bp)",
        main = "Coding Sequence Length Distribution")
```

#### Inputs:
ecoli_length: A numeric vector containing the lengths of coding sequences for *Escherichia coli*.

streptacidiphilus_length: A numeric vector containing the lengths of coding sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
A boxplot displaying the distribution of coding sequence lengths for both organisms, with labeled axes and a title.

### Mean and Median Coding Sequence Length

The purpose of this code is to calculate and summarize the mean and median coding sequence lengths for *Escherichia coli* and *Streptacidiphilus jiangxiensis*.

```r
mean_median <- data.frame(
  Organism = c("Escherichia coli", "Streptacidiphilus jiangxiensis"),
  Mean_Length = c(mean(ecoli_length), mean(streptacidiphilus_length)),
  Median_Length = c(median(ecoli_length), median(streptacidiphilus_length))
)

kable (mean_median, align = 'c', caption = "Mean and Median of the CDS length each organism")
```

#### Inputs:
ecoli_length: A numeric vector containing the lengths of coding sequences for *Escherichia coli*.

streptacidiphilus_length: A numeric vector containing the lengths of coding sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
A data frame named mean_median that includes:

  -Organism: Names of the two organisms.

  -Mean_Length: The mean coding sequence length for each organism.

  -Median_Length: The median coding sequence length for each organism.
  
Table Output: A formatted table displayed using kable, showing the mean and median lengths of coding sequences for each organism with a caption "Mean and Median of the CDS length each organism."


### Frequency of DNA Bases

The frecuency of DNA bases in the coding sequences of *Escherichia coli* and *Streptacidiphilus jiangxiensis*, is calculated and visualized using the following code.

```r
# Calculate base frequencies

dna_ecoli <- unlist(ecoli_seqs)
ecoli_dna_df <- data.frame(base = dna_ecoli) 
ecoli_dna_composition <- ecoli_dna_df %>% count(base)  

dna_streptacidiphilus <- unlist(streptacidiphilus_seqs)
streptacidiphilus_dna_df <- data.frame(base = dna_streptacidiphilus)
streptacidiphilus_dna_composition <- streptacidiphilus_dna_df %>% count(base)

# Bar plots

barplot(height=ecoli_dna_composition$n,
        names.arg=ecoli_dna_composition$base,
        col = "lightgreen",
        xlab="nucleotides",
        ylab="frequency", 
        main="E. coli CDS composition")
grid()
barplot(height=streptacidiphilus_dna_composition$n,
        names.arg = streptacidiphilus_dna_composition$base,
        col="lightblue",
        xlab="nucleotides",
        ylab="frequency", 
        main="S. jiangxiensis CDS composition")
grid()


```

#### Inputs:
ecoli_seqs: A list of coding DNA sequences for *Escherichia coli*.

streptacidiphilus_seqs: A list of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
ecoli_dna_composition: A vector containing the frequency of each nucleotide in *Escherichia coli* coding sequences.

streptacidiphilus_dna_composition: A vector containing the frequency of each nucleotide in *Streptacidiphilus jiangxiensis* coding sequences.

Bar Plots: Visual representations of nucleotide frequencies for both organisms, labeled appropriately.


### Amino Acid Frequency

The purpose of this code is to convert coding DNA sequences into protein sequences and then calculate and visualize the frequency of amino acids in *Escherichia coli* and *Streptacidiphilus jiangxiensis*.

```r
# Convert coding DNA to protein sequences

ecoli_seqs <- seqinr::read.fasta ("e_coli_cds.fa")
streptacidiphilus_seqs <- seqinr::read.fasta ("streptacidiphilus_cds.fa")

ecoli_prot <- lapply(ecoli_seqs, translate)
streptacidiphilus_prot <- lapply(streptacidiphilus_seqs, translate)


# Calculate amino acid frequencies
ecoli_proteins <- unlist (ecoli_prot)
streptacidiphilus_proteins <-unlist (streptacidiphilus_prot)

aa_ecoli <-unique(ecoli_proteins)
aa_ecoli <- aa_ecoli[aa_ecoli != "*"]
aa_streptacidiphilus <- unique(streptacidiphilus_proteins)
aa_streptacidiphilus <- aa_streptacidiphilus[aa_streptacidiphilus != "*"]

ecoli_proteins_df <- data.frame(aa = ecoli_proteins, stringsAsFactors = FALSE)
streptacidiphilus_proteins_df <- data.frame(aa = streptacidiphilus_proteins, stringsAsFactors = FALSE)

ecoli_aa_freq <- ecoli_proteins_df %>%
  count(aa)

streptacidiphilus_aa_freq <- streptacidiphilus_proteins_df %>%
  count(aa)

# Bar plots
barplot(height=ecoli_aa_freq$n,
        names.arg = ecoli_aa_freq$aa,
        col = "lightgreen",
        xlab="Aminoacids",
        ylab="Frequency", 
        main="Amino Acid Frequency in Escherichia coli")
grid()
barplot(height=streptacidiphilus_aa_freq$n, 
        names.arg = streptacidiphilus_aa_freq$aa,
        col="lightblue",
        xlab="Aminoacids",
        ylab="Frequency", 
        main="Amino Acid Frequency in Streptacidiphilus jiangxiensis")
grid()
```

#### Inputs:
ecoli_seqs: A list of coding DNA sequences for *Escherichia coli*.

streptacidiphilus_seqs: A list of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

#### Outputs:
ecoli_aa_freq: A vector containing the frequency of each amino acid in *Escherichia coli* protein sequences.

streptacidiphilus_aa_freq: A vector containing the frequency of each amino acid in *Streptacidiphilus jiangxiensis* protein sequences.

Bar Plots: Visual representations of amino acid frequencies for both organisms, labeled appropriately.


### Codon Usage Bias

This code calculates and visualizes the codon usage bias in *Escherichia coli* and *Streptacidiphilus jiangxiensis* by creating relative synonymous codon usage (RSCU) tables and generating bar plots for comparison.

```r
# Create codon usage tables
RSCU_ecoli <- uco(dna_ecoli,index="rscu",as.data.frame=TRUE)
RSCU_streptacidiphilus <- uco(dna_streptacidiphilus,index="rscu",as.data.frame=TRUE)

# Create a combined table
RSCU_combined <- cbind(RSCU_ecoli, "S. jiangxiensis RSCU"=RSCU_streptacidiphilus$RSCU)
colnames(RSCU_combined)[5] <- "E. coli RSCU"
#Give format to the table
kable (head(RSCU_combined), align = 'c', caption = "Codon Usage Bias from each from each organism") 

# Bar plots:
# E. coli

barplot(height=RSCU_ecoli$RSCU, 
        names.arg=RSCU_ecoli$codon,
        col = "lightgreen",
        xlab= "Codons",
        ylab="Relative Synonymous Codon Usage",
        main="Codon Usage Bias in Escherichia coli",
        las = 2,  # Rotate labels
        space = 0.5,  # Reduce space between bars
        width = 0.4,
        ylim = c(0, 4))
grid()

# Streptacidiphilus jiangxiensis
 
barplot(height=RSCU_streptacidiphilus$RSCU, 
        names.arg=RSCU_streptacidiphilus$codon,
        col = "lightblue",
        xlab= "Codons",
        ylab="Relative Synonymous Codon Usage", 
        main="Codon Usage in Streptacidiphilus jiangxiensis",
        las = 2,  # Rotate labels
        space = 0.5,  # Reduce space between bars
        width = 0.4,
        ylim = c(0, 4))
grid()
```

### Inputs:
dna_ecoli: A vector of coding DNA sequences for *Escherichia coli*.

dna_streptacidiphilus: A vector of coding DNA sequences for *Streptacidiphilus jiangxiensis*.

### Outputs:
RSCU_ecoli: A data frame containing the codon usage data for *Escherichia coli*.

RSCU_streptacidiphilus: A data frame containing the codon usage data for *Streptacidiphilus jiangxiensis*.

RSCU_combined: A combined table of RSCU values for both organisms.

Formatted Table Output: A table generated using kable, showing the first few rows of the combined RSCU data with the caption "Codon Usage Bias from each organism."

Bar Plots: Visual representations of codon usage bias for both *Escherichia coli* and *Streptacidiphilus jiangxiensis*, labeled appropriately.


### Overlay the two barplots to compare the RSCU between the two organisms.

The purpose of this code is to create an overlaid bar plot to visually compare the relative synonymous codon usage (RSCU) between *Escherichia coli* and *Streptacidiphilus jiangxiensis*.

```r
# Comparing two barplots

bar_heights_ecoli <- barplot(height = RSCU_ecoli$RSCU, 
                              names.arg = RSCU_ecoli$codon,
                              col = rgb(0, 1, 0, 0.5),  # Transparent green
                              xlab = "Codons",
                              ylab = "Relative Synonymous Codon Usage",
                              main = "Codon Usage Bias",
                              las = 2,  # Rotate labels
                              space = 0.5,  # Reduce space between bars
                              width = 0.4,
                              ylim = c(0, 3.8))  # Set proper limits for y-axis


bar_heights_streptacidiphilus <- barplot(height = RSCU_streptacidiphilus$RSCU, 
                                          col = rgb(0, 0, 1, 0.5),  # Transparent blue
                                          add = TRUE,  # Overlay on the existing plot
                                          las = 2,  # Rotate labels
                                         space = 0.5,  # Same space for consistency
                                          width = 0.4)

# Legend 
legend("topright", 
       legend = c("Escherichia coli", "Streptacidiphilus jiangxiensis"), 
       fill = c(rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)))

grid()
```

#### Inputs:
RSCU_ecoli: A data frame containing RSCU values for *Escherichia coli codons*.

RSCU_streptacidiphilus: A data frame containing RSCU values for *Streptacidiphilus jiangxiensis* codons.

#### Outputs:
Overlaid Bar Plot: A single plot displaying the codon usage for both organisms with transparency to distinguish between the two datasets.

Legend: A legend indicating which color corresponds to each organism.


### K-mer Analysis

The purpose of this code is to analyze the frequency of k-mers (3-mers, 4-mers, and 5-mers) in the protein sequences of *Escherichia coli* and *Streptacidiphilus jiangxiensis*, and to visualize the top and bottom k-mers for both organisms.

```r
# Function to calculate k-mers

# 3-mers Escherichia coli
ecoli_prot_freq_3 <- seqinr::count(ecoli_proteins, wordsize=3, alphabet=aa_ecoli,freq=TRUE)
ecoli_prot_freq_3 <- as.data.frame (ecoli_prot_freq_3)
colnames(ecoli_prot_freq_3)[1] <- "3-mer"
#If needed to confirm: head(ecoli_prot_freq_3)

# 3-mers Streptacidiphilus jiangxiensis
streptacidiphilus_prot_freq_3 <- seqinr::count(streptacidiphilus_proteins, wordsize=3, alphabet=aa_streptacidiphilus,freq=TRUE)
streptacidiphilus_prot_freq_3 <- as.data.frame(streptacidiphilus_prot_freq_3)
colnames(streptacidiphilus_prot_freq_3)[1] <- "3-mer"
#If needed to confirm: head(streptacidiphilus_prot_freq_3)

# 4-mers Escherichia coli
ecoli_prot_freq_4 <- seqinr::count(ecoli_proteins, wordsize=4, alphabet=aa_ecoli, freq=TRUE)
ecoli_prot_freq_4 <- as.data.frame (ecoli_prot_freq_4)
colnames(ecoli_prot_freq_4)[1] <- "4-mer"
#If needed to confirm: head(ecoli_prot_freq_4)

# 4-mers Streptacidiphilus jiangxiensis
streptacidiphilus_prot_freq_4 <- seqinr::count(streptacidiphilus_proteins, wordsize=4, alphabet=aa_streptacidiphilus,freq=TRUE)
streptacidiphilus_prot_freq_4 <- as.data.frame(streptacidiphilus_prot_freq_4)
colnames(streptacidiphilus_prot_freq_4)[1] <- "4-mer"
#If needed to confirm: head(streptacidiphilus_prot_freq_4)

# 5-mers Escherichia coli
ecoli_prot_freq_5 <- seqinr::count(ecoli_proteins, wordsize=5, alphabet=aa_ecoli, freq=TRUE)
ecoli_prot_freq_5 <- as.data.frame (ecoli_prot_freq_5)
colnames(ecoli_prot_freq_5)[1] <- "5-mer"
#If needed to confirm: head(ecoli_prot_freq_5)

# 5-mers Streptacidiphilus jiangxiensis
streptacidiphilus_prot_freq_5 <- seqinr::count(streptacidiphilus_proteins, wordsize=5, alphabet=aa_streptacidiphilus,freq=TRUE)
streptacidiphilus_prot_freq_5 <- as.data.frame(streptacidiphilus_prot_freq_5)
colnames(streptacidiphilus_prot_freq_5)[1] <- "5-mer"
#If needed to confirm: head(streptacidiphilus_prot_freq_5)

# Create a table for each organism with the k-mers and their frequencies

#E. coli
# Combine k-mer data into one column
ecoli_kmers <- data.frame(
  Kmer = c(
    as.character(ecoli_prot_freq_3$`3-mer`),
    as.character(ecoli_prot_freq_4$`4-mer`),
    as.character(ecoli_prot_freq_5$`5-mer`)
  ),
  Frequency = c(
    as.numeric(ecoli_prot_freq_3$Freq),
    as.numeric(ecoli_prot_freq_4$Freq),
    as.numeric(ecoli_prot_freq_5$Freq)
  )
)
# Ordering the data from most frequent to least frequent
ecoli_kmers <- ecoli_kmers[order(ecoli_kmers$Frequency, decreasing = TRUE), ]
#If needed to confirm: head(ecoli_kmers)


#Filtering E.coli K-mers with >0 frequency
ecoli_kmers_filtered <- ecoli_kmers[ecoli_kmers$Frequency > 0, ]


#Selecting the 10 most frequent protein sequence k-mers in E.coli
ecoli_top_10 <- ecoli_kmers_filtered[order(ecoli_kmers_filtered$Frequency, decreasing = TRUE), ][1:10, ]
kable(head(ecoli_top_10), align = 'c', caption = "Top 10 most frequent k-mers in E.coli")

#Selecting the 10 least frequent protein sequence k-mers in E.coli
ecoli_bottom_10 <- ecoli_kmers_filtered[order(ecoli_kmers_filtered$Frequency), ][1:10, ]
kable (head(ecoli_bottom_10), align = 'c', caption = "Least 10 frequent k-mers in S.jiangxiensis")

#Selecting and countig E.coli K-mers with 0 frequency
ecoli_kmers_null <- ecoli_kmers[ecoli_kmers$Frequency == 0, ]
ecoli_num_kmers_null <- nrow(ecoli_kmers_null)
cat("Number of E. coli K-mers with 0 frequency:")
print (ecoli_num_kmers_null)

#S. jiangxiensis
# Combine k-mer data into one column
streptacidiphilus_kmers <- data.frame(
  Kmer = c(
    as.character(streptacidiphilus_prot_freq_3$`3-mer`),
    as.character(streptacidiphilus_prot_freq_4$`4-mer`),
    as.character(streptacidiphilus_prot_freq_5$`5-mer`)
  ),
  Frequency = c(
    as.numeric(streptacidiphilus_prot_freq_3$Freq),
    as.numeric(streptacidiphilus_prot_freq_4$Freq),
    as.numeric(streptacidiphilus_prot_freq_5$Freq)
  )
)
# Ordering the data from most frequent to least frequent
streptacidiphilus_kmers <- streptacidiphilus_kmers[order(streptacidiphilus_kmers$Frequency, decreasing = TRUE), ]
#If needed to confirm: head(streptacidiphilus_kmers)

#Filtering S.jiangxiensis K-mers with >0 frequency
streptacidiphilus_kmers_filtered <- streptacidiphilus_kmers[streptacidiphilus_kmers$Frequency > 0, ]

#Selecting the 10 most frequent protein sequence k-mers in S.jiangxiensis
streptacidiphilus_top_10 <- streptacidiphilus_kmers_filtered[order(streptacidiphilus_kmers_filtered$Frequency, decreasing = TRUE), ][1:10, ]
kable (head(streptacidiphilus_top_10), align = 'c', caption = "Top 10 most frequent k-mers in S.jiangxiensis")

#Selecting the 10 least frequent protein sequence k-mers in S.jiangxiensis
streptacidiphilus_bottom_10 <- streptacidiphilus_kmers_filtered[order(streptacidiphilus_kmers_filtered$Frequency), ][1:10, ]
kable (head(streptacidiphilus_bottom_10), align = 'c', caption = "Least 10 frequent k-mers in S.jiangxiensis")


#Selecting and countig S.jiangxiensis K-mers with 0 frequency
streptacidiphilus_kmers_null <- streptacidiphilus_kmers[streptacidiphilus_kmers$Frequency == 0, ]
streptacidiphilus_num_kmers_null <- nrow(streptacidiphilus_kmers_null)
cat("Number of S. jiangxiensis K-mers with 0 frequency:")
print (streptacidiphilus_num_kmers_null)

#Plot of top 10 E.coli K-mers
barplot(height=ecoli_top_10$Frequency, 
        names.arg=ecoli_top_10$Kmer,
        col = "lightgreen",
        xlab= "K-mers",
        ylab="Frequency",
        main="Top 10 E.coli K-mers",
        las = 2,  # Rotate labels
        space = 0.5,  
        width = 0.4,
        ylim = c(0, 0.002),
        cex.axis = 0.6)

grid()

#Plot of top 10 S.jiangxiensis K-mers
barplot(height=streptacidiphilus_top_10$Frequency, 
        names.arg=streptacidiphilus_top_10$Kmer,
        col = "lightblue",
        xlab= "K-mers",
        ylab="Frequency",
        main="Top 10 S.jiangxiensis K-mers",
        las = 2,  # Rotate labels
        space = 0.5,  
        width = 0.4,
        ylim = c(0, 0.005),
        cex.axis = 0.6)
grid()

#Plot of bottom 10 E.coli K-mers
barplot(height=ecoli_bottom_10$Frequency, 
        names.arg=ecoli_bottom_10$Kmer,
        col = "lightgreen",
        xlab= "K-mers",
        ylab="Frequency",
        main="Bottom 10 E.coli K-mers",
        las = 2,  # Rotate labels
        space = 0.5,  
        width = 0.4,
        cex.axis = 0.6
       )
grid()

#Plot of bottom 10 S.jiangxiensis K-mers
barplot(height=streptacidiphilus_bottom_10$Frequency, 
        names.arg=streptacidiphilus_bottom_10$Kmer,
        col = "lightblue",
        xlab= "K-mers",
        ylab="Frequency",
        main="Bottom 10 S.jiangxiensis K-mers",
        las = 2,  # Rotate labels
        space = 0.5,  
        width = 0.4,
        cex.axis = 0.6
        )
grid()

```

#### Inputs:
ecoli_proteins: Protein sequence data for Escherichia coli.

streptacidiphilus_proteins: Protein sequence data for Streptacidiphilus jiangxiensis.

aa_ecoli: Alphabet used for Escherichia coli, likely a set of amino acids.

aa_streptacidiphilus: Alphabet used for Streptacidiphilus jiangxiensis.

#### Outputs:
K-mer Frequency Tables:
  -ecoli_kmers: A data frame containing k-mers of lengths 3, 4, and 5 for Escherichia coli along with their frequencies.
  -streptacidiphilus_kmers: A data frame for Streptacidiphilus jiangxiensis containing the same information.
  
Formatted Table Outputs:
  -A table displaying the top 10 most frequent k-mers for both organisms with the caption "Top 10 most frequent k-mers in E. coli" and "Top 10 most frequent k-mers in S. jiangxiensis."
  -A table for the least frequent k-mers with the respective captions.

Count of K-mers with Zero Frequency:
  -A printed statement indicating the number of k-mers with a frequency of zero for both organisms.

Bar Plots:
  -Bar plots visualizing the top and bottom 10 k-mers for both Escherichia coli and Streptacidiphilus jiangxiensis, showing the frequency of each k-mer.
  

### Session Info

```r
sessionInfo()
```

give a comprehensive overview of the working R environment. This is particularly helpful for:

-Debugging and reproducibility.
-Sharing your session details with others when reporting issues.
-Keeping track of package versions and R settings.

#### Outputs
R Version: The version of R that is currently running.

Platform: Information about the operating system (e.g., Windows, macOS, Linux).

Locale: The current locale settings, which affect how R handles text encoding and date formats.

Attached Base Packages: A list of base R packages that are loaded by default.

Other Loaded Packages: A list of any additional packages that have been loaded during the session, along with their versions.


## Conclusion

In this assignment, we explored various aspects of gene expression and biological sequence diversity. Initially, we analyzed RNA-seq count data from the file "gene_expression.tsv," calculating the mean expression for each gene and identifying the top 10 genes with the highest mean values, along with determining how many genes had a mean expression below 10. The histogram of mean expression values provided a visual representation of the distribution, highlighting patterns in gene activity across the samples.

Next, we examined tree growth data from "growth_data.csv," comparing tree circumference measurements taken at a control site and a treatment site over a 20-year period. We calculated the mean and standard deviation of tree circumference at both the start and end of the study, visualized the data through box plots, and assessed growth over the last decade using t-tests to evaluate significant differences between the sites.

Part 2 of the assignment focused on comparing the sequence features of *Streptacidiphilus jiangxiensis* to *E. coli*. We downloaded and analyzed coding DNA sequences, tabulating the number of coding sequences and their total length, and identified differences in coding sequence lengths between the organisms. We also calculated the frequency of DNA bases and amino acids, generating bar plots to illustrate these frequencies. Codon usage bias was quantified and compared, providing insights into the evolutionary strategies of the two organisms. Finally, we identified and compared k-mers, revealing significant differences in sequence representation between the two genomes.

Through this comprehensive analysis, we not only gained a deeper understanding of gene expression and tree growth patterns but also explored the genetic diversity and evolutionary implications of coding sequences in two distinct bacterial species. This integrative approach highlighted the importance of bioinformatics in elucidating the complexities of biological systems.


## References:

BENGTSSON, H. 2023. Various Programming Utilities [R package R. utils version 2.12. 3].

CHARIF, D. & LOBRY, J. R. 2007. SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. Structural approaches to sequence evolution: Molecules, networks, populations. Springer.

PAGÈS, H., ABOYOUN, P., GENTLEMAN, R. & DEBROY, S. 2023. Biostrings: efficient manipulation of biological strings. Biostrings, R package version 2.70. 1.

WICKHAM, H., FRANÇOIS, R., HENRY, L., MÜLLER, K. & VAUGHAN, D. 2023. dplyr: a grammar of data manipulation. R package version 1.1. 2. Computer software.

WICKHAM, H. & HENRY, L. 2020. tidyr: Tidy Messy Data. R package version 1.1. 2. CRAN. R-project. org/package= tidyr.

WICKHAM, H., HESTER, J. & FRANCOIS, R. 2018. Readr: Read rectangular text data.

WICKHAM, H. & WICKHAM, H. 2016. Data analysis, Springer.

