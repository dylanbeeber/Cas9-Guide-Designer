# Cas9-Guide-Designer
Software used to design guide RNA sequences for CRISPR/Cas9 genome editing

This software aims to provide all scientifically pertinent information when designing guide RNA sequences for Cas9 genome editing. When provided a target DNA sequence for editing, a genome to check for off-targets in, and a genome annotation file (.gtf) to provide addition information about off-target matches it will out put information for two separate data tables. The first table contains all information on the generated sgRNA themselves (sgRNA sequence, PAM, Direction, Start, End, GC content, Presence of Homopolymers, Effciency Score (Doench 2014), and Genomic Matches). The second table contains all information on the found off-target sequences (Original sgRNA Sequence, Chromosome, Start, End, Number of Mismatches, Direction, Matched Sequence, Gene ID, Gene Name, Sequence Type, and Exon Number)

# Required files:
StandaloneFindsgRNAfunction_Doench2014.R - The main script that contains all the code needed to design sgRNA. This is the only R file most users will need.
Doench_Model_Weights_Singleonly.csv and Doench_Model_Weights_Doubleonly.csv - Two data tables used to assist with efficiency scoring. These must be put in the working directory when using the sgRNA_design function.

# Optional files:
RunShiny.R - A script that contains code for a user interface. This UI requires installation of several addition packages and is currently only available for the human and yeast genomes.
FindsgRNAfunction_Doench2014.R - A script designed to be used with the user optional Shiny user interface
Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz - an example of a gene annotation file (.gtf) that needs to be used with the the sgRNA_design function

# Instructions for StandaloneFindsgRNAfunction_Doench2014.R
The working directory must be set to a folder that contains this script, a .gtf file from a species of your choice, and the files: "Doench_Model_Weights_Singleonly.csv" and "Doench_Model_Weights_Doubleonly.csv"

Example: setwd("C://Users//User//Desktop//Folder")

if packages need to be installed input the following code:
install.packages("stringr", repos='http://cran.us.r-project.org')
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("BSgenome")

Your organism's genome must also be obtained from BSgenome
To obtain a list of available genomes type:
available.genomes()
When a genome has been selected use the following code:
biocLite("your.genenome")
Example: biocLite("BSgenome.Hsapiens.UCSC.hg19")

Finally, a genome annotation file (.gtf) must be obtained for your organism and placed in the working directory
These can be found at: https://useast.ensembl.org/info/data/ftp/index.html

Simply source this file and run the sgRNA_desgin function withthe target sequence, genome, and gtf as arguments.
Example:
source("StandaloneFindsgRNAfunction_Doench2014.R")
alldata <- sgRNA_design("ATTCGAGGAGACTATAGAGCAGGATTAGGACAGAGACCATGTGACAGAA", Scerevisiae, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")

Important Note: When designing sgRNA for large genomes (billions of base pairs), use short query DNA sequences (under 250 bp). Depending on your hardware checking for off-targets can be quite computationally intensive and may take several hours if not limited to smaller query sequences.
