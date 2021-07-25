Outdated project: Please see crispRdesignR.

# Cas9-Guide-Designer
Software used to design guide RNA sequences for CRISPR/Cas9 genome editing

This directory contains files used to develop the crispRdesignR package and a User Interface (in Shiny) to run the crispRdesignR tool through a browser.
For the R package crispRdesignR, see https://github.com/dylanbeeber/crispRdesignR

This software aims to provide all scientifically pertinent information when designing guide RNA sequences for Cas9 genome editing. When provided a target DNA sequence for editing, a genome to check for off-targets in, and a genome annotation file (.gtf) to provide addition information about off-target matches it will out put information for two separate data tables. The first table contains all information on the generated sgRNA themselves (sgRNA sequence, PAM, Direction, Start, End, GC content, Presence of Homopolymers, Self-Complementary Sequences, Effciency Score (Doench 2016), and Genomic Matches). The second table contains all information on the found off-target sequences (Original sgRNA Sequence, Chromosome, Start, End, Number of Mismatches, Direction, CFD Scores, Matched Sequence, Gene ID, Gene Name, Sequence Type, and Exon Number)

# Required files:
RunShiny.R - A script that contains code for a user interface. This UI requires installation of several addition packages and is currently only available for the human and yeast genomes.

FindsgRNAfunction.R - A script designed to be used from within the Shiny user interface

Rule_Set_2_Model.rds - A gradient boosted regression model trained on data from the Doench 2016 paper.

CFD_Scoring.csv - A data table that contains the information used to calculate the off-target effects of off-target sequences.

# Optional files
Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz - an example of a gene annotation file (.gtf) that needs to be used with the the sgRNA_design function. In order to run the program with default settings, a .gtf file for the target organism must be provided.

StandaloneFindsgRNAfunction_Doench2014.R (outdated) - The script that contains all the code needed to design sgRNA without the Shiny UI. This is useful for debugging and testing.

FindsgRNAfunction_Doench2014.R(outdated) - An older version of the FindsgRNA script that uses the Doench 2014 rule set to predict the efficiency of sgRNA. Requires the files: Doench_Model_Weights_Singleonly.csv, Doench_Model_Weights_Doubleonly.csv.

# Instructions for the Shiny UI
Download RunShiny.R, FindsgRNAfunction_Doench2014.R, CFD_Scoring.csv, and Rule_Set_2_Model.rds. Open RunShiny.R and set your working directory to a location that contains all of the previously downloaded files for this program. Install packages as necesary.

A list of compatible genomes to check for off-targets in may be located by using the command `available.genomes()` in the R console. These genomes may then be installed using the following command:
`biocLite("your.genome")`
Example: `biocLite("BSgenome.Hsapiens.UCSC.hg19")`

A a genome annotation file (.gtf) specfic to your genome of choice is required to provide detailed information on possible off-target sequences. These can be found at: https://useast.ensembl.org/info/data/ftp/index.html. This file can be saved anywhere as the UI will ask for it to be provided manually. A sample .gtf file has been provided in this directory (Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz).

The user interface can then be run using the command: `runApp("RunShiny.R")`

In the UI, provide either a DNA sequence or a .fasta file to design sgRNA for. Then choose your desired genome from the drop-down genome list. Provide a .gtf file to annatote off-target sequences. Options can be set to skip the annotation of off-target sequences or prevent the program from calling off-target sequences. Finally click "Find sgRNA" and the program will begin.

Important Note: When designing sgRNA for large genomes (billions of base pairs), use short query DNA sequences (under 250 bp). Depending on your hardware checking for off-targets can be quite computationally intensive and may take several hours if not limited to smaller query sequences.

# Instructions for StandaloneFindsgRNAfunction_Doench2014.R
The working directory must be set to a folder that contains this script, a .gtf file from a species of your choice, and the files: "Doench_Model_Weights_Singleonly.csv" and "Doench_Model_Weights_Doubleonly.csv"


Example: `setwd("C://Users//User//Desktop//Folder")`

if packages need to be installed input the following code:

`install.packages("stringr", repos='http://cran.us.r-project.org')`

`source("https://bioconductor.org/biocLite.R")`

`biocLite("Biostrings")`

`biocLite("BSgenome")`


Your organism's genome must also be obtained from BSgenome

To obtain a list of available genomes type:

`available.genomes()`

When a genome has been selected use the following code:

`biocLite("your.genome")`

Example: `biocLite("BSgenome.Hsapiens.UCSC.hg19")`


Finally, a genome annotation file (.gtf) must be obtained for your organism and placed in the working directory

These can be found at: https://useast.ensembl.org/info/data/ftp/index.html

Simply source this file and run the sgRNA_desgin function withthe target sequence, genome, and gtf as arguments.

Example:

`source("StandaloneFindsgRNAfunction_Doench2014.R")`

`alldata <- sgRNA_design("ATTCGAGGAGACTATAGAGCAGGATTAGGACAGAGACCATGTGACAGAA", Scerevisiae, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")`


Important Note: When designing sgRNA for large genomes (billions of base pairs), use short query DNA sequences (under 250 bp). Depending on your hardware checking for off-targets can be quite computationally intensive and may take several hours if not limited to smaller query sequences.
