## The working directory must be set to a folder that contains this
## script, a gtf file from a species of your choice, and the files:
## "Doench_Model_Weights_Singleonly.csv" and "Doench_Model_Weights_Doubleonly.csv"
##
## Example: setwd("C://Users//User//Desktop//Folder")
##
## if packages need to be installed input the following code:
## install.packages("stringr", repos='http://cran.us.r-project.org')
## source("https://bioconductor.org/biocLite.R")
## biocLite("Biostrings")
## biocLite("BSgenome")
##
## Your organism's genome must also be obtained from BSgenome
## To obtain a list of available genomes type:
## available.genomes()
## When a genome has been selected use the following code:
## biocLite("your.genenome")
## library(your.genome)
## Example: biocLite("BSgenome.Hsapiens.UCSC.hg19")
## library(BSgenome.Hsapiens.UCSC.hg19)
##
## Finally, a genome annotation file (.gtf) must be obtained for your organism and placed in the working directory
## These can be found at: https://useast.ensembl.org/info/data/ftp/index.html
##
## Simply source this file and run the sgRNA_desgin function with
## the target sequence, genome, and gtf as arguments.
## Example:
## source("StandaloneFindsgRNAfunction_Doench2014.R")
## alldata <- sgRNA_design("GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGGCCTATGTGAGTTTTTGAAGGGGGTTCAAAGCGCCTTGTAA", Scerevisiae, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")
##
## Important Note: When designing sgRNA for large genomes (billions of base pairs),
## use short query DNA sequences (under 250 bp). Depending on your hardware
## checking for off-targets can be quite computationally intensive and may
## take several hours if not limited to smaller query sequences.
##
library(Biostrings)
library(stringr)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

sgRNA_design <- function(usersequence, genomename, gtfname){
  sequence <- paste(usersequence, collapse = "")
  sequence <- str_replace_all(sequence, fixed(" "), "")
  Biostrings_sequence <- DNAString(sequence)
  ## Creates a character string that contains the
  ## complementary sequence (Both in the reverse
  ## direction and with substituted nucleotides)
  Biostrings_rev_seq <- reverseComplement(Biostrings_sequence)
  rev_seq <- as.character(Biostrings_rev_seq)
  ## Create an empty list for the forward sgRNA to go
  sgRNA_list_f <- c()
  ## Create an empty list for the reverse sgRNA to go
  sgRNA_list_r <- c()
  ## Create four empty lists for forward and reverse start and end positions
  ## Check to make the start and end numbers are correct
  sgRNA_f_start <- c()
  sgRNA_f_end <- c()
  sgRNA_r_start <- c()
  sgRNA_r_end <- c()
  ## Sets number that helps determine when to stop looking
  ## for possible sgRNA (This prevents it from choosing
  ## an incomplete sgRNA at the very end of the sequence)
  num_char_in_seq <- nchar(sequence) - 29
  ## Sets the PAM sequence and determines what the program
  ## will describe as a possible sgRNA. Even though most sgRNA is
  ## only 20 nucleotides long, nucleotides surrounding the sgRNA
  ## are used for study-based scoring
  PAM <- ".........................GG..."
  ## Sets n to zero, which is used to incrementally increase
  ## and search for all possible sgRNA in the sequence
  n <- 0
  ## Searches all 23 nt streches in the sequence for
  ## possible matches to the PAM, then puts entire 30 nt
  ## matches into a list (including the PAM)
  for (x in 0:num_char_in_seq){
    poss_sgRNA <- substr(sequence, 1+n, 30+n)
    if (str_detect(poss_sgRNA, PAM) == TRUE){
      sgRNA_list_f[[length(sgRNA_list_f)+1]] <- poss_sgRNA
      sgRNA_f_start[[length(sgRNA_f_start)+1]] <- n+5
      sgRNA_f_end[[length(sgRNA_f_end)+1]] <- n+27
    }
    n <- n+1
  }
  ## Same as above but with the reverse sequence
  n <- 0
  for (x in 0:num_char_in_seq){
    poss_sgRNA <- substr(rev_seq, 1+n, 30+n)
    if (str_detect(poss_sgRNA, PAM) == TRUE){
      sgRNA_list_r[[length(sgRNA_list_r)+1]] <- poss_sgRNA
      sgRNA_r_start[[length(sgRNA_r_start)+1]] <- nchar(rev_seq)-n+5
      sgRNA_r_end[[length(sgRNA_r_end)+1]] <- nchar(rev_seq)-n+27
    }
    n <- n+1
  }
  ## Removes any sgRNA that contain degerate bases
  sgRNA_list_f <- sgRNA_list_f[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_f_start <- sgRNA_f_start[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_f_end <- sgRNA_f_end[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_list_r <- sgRNA_list_r[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  sgRNA_r_start <- sgRNA_r_start[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  sgRNA_r_end <- sgRNA_r_end[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  ## Creates list for all sgRNA to go
  sgRNA_list <- c(sgRNA_list_f, sgRNA_list_r)
  if (is.null(sgRNA_list) == FALSE) {
    sgRNA_start <- c(sgRNA_f_start, sgRNA_r_start)
    sgRNA_end <- c(sgRNA_f_end, sgRNA_r_end)
    ## Creates a list with only the sgRNA sequences (no PAMs or
    ## flanking sequence)
    breakseq <- function(seqlist){
      str_sub(seqlist, 5, 24)
    }
    sgRNA_seq <- sapply(sgRNA_list, breakseq)
    ## Creates a list with only the PAM sequences
    breakPAM <- function(seqlist){
      str_sub(seqlist, 25, 27)
    }
    sgRNA_PAM <- sapply(sgRNA_list, breakPAM)
    ## Makes a list of all of the sgRNA sequences with their PAM
    sgRNA_with_PAM <- paste(sgRNA_seq, sgRNA_PAM, sep = "")
    ## Makes a list of whether sgRNA are forward or reverse
    sgRNA_fow <- rep("+", each = length(sgRNA_list_f))
    sgRNA_rev <- rep("-", each = length(sgRNA_list_r))
    sgRNA_fow_or_rev <- c(sgRNA_fow, sgRNA_rev)
    ## Find GC percentage for each sgRNA and puts that data into
    ## a list called "GCinstance"
    FindGC <- function(seqlist){
      ((str_count(seqlist, "G") + str_count(seqlist, "C")) / 20)
    }
    GCinstance <- sapply(sgRNA_seq, FindGC)
    ## Creates a list that determines if the GC percentage
    ## is within the user-defined (WIP) threshold
    EvalGC <- function(GC){
      isTRUE(30<=GC&GC<= 70)
    }
    GClist <- sapply(GCinstance, EvalGC)
    ## Find homopolmers
    Findhomopolymer <- function(seqlist){
      str_detect(seqlist, "TTTT|AAAA|GGGG|CCCC")
    }
    Homopolymerdetect <- sapply(sgRNA_seq, Findhomopolymer)
    Homopolymerdetect
    ## Find TTTT homopolymers
    FindTTTThomopolymer <- function(seqlist){
      str_detect(seqlist, "TTTT")
    }
    TTTTHomopolymerdetect <- sapply(sgRNA_seq, FindTTTThomopolymer)
    ## Detect if self complementary (not currently implemented)
    ## Assign a study-based efficiency score
    ## ***Need to add GC Penalty
    ## Following two lines retrieve the penalty constants (one for single nucleotides, the other for paired nucleotides)
    Doench_model_weights_singleonly <- read.csv("Doench_Model_Weights_Singleonly.csv", header = FALSE)
    Doench_model_weights_doubleonly <- read.csv("Doench_Model_Weights_Doubleonly.csv", header = FALSE)
    ## Creates an empty list for Doench study-based scores to go into
    Doench_Score <- c()
    ## "g" allows this script to go through the sgRNA_list, one by one
    g <- 0
    for (h in 1:length(sgRNA_list)){
      ## Splits the sgRNA into individual nucleotides
      split_sgRNA <- str_split(sgRNA_list[1+g], "", simplify = TRUE)
      ## Creates the sgRNA_model_weight list and adds the intercept to it
      sgRNA_model_weights <- c(0.597636154)
      ## Finds the model weights for each individual nucleotide in the sgRNA and adds them to sgRNA_model_weights
      n <- 0
      for (t in 1:39){
        if ((split_sgRNA[,Doench_model_weights_singleonly[1+n, 2]] == Doench_model_weights_singleonly[1+n, 1]) == TRUE){
          sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- Doench_model_weights_singleonly[1+n, 3]
        }
        n <- n+1
      }
      ## Splits the sgRNA into pieces two nucleotides long
      ## Serves the same purpose as: "split_sgRNA <- str_split(sgRNA_list[1+g], "", simplify = TRUE)"
      double <- ".."
      split_double_sgRNA <- c()
      n <- 0
      for (xp in 1:29){
        double_sgRNA <- substr(sgRNA_list[1+g], 1+n, 2+n)
        split_double_sgRNA[[length(split_double_sgRNA)+1]] <- double_sgRNA
        n <- n+1
      }
      ## Finds the model weights for each double nucleotide in the sgRNA and adds them to sgRNA_model_weights
      n <- 0
      for (t in 1:31){
        if ((split_double_sgRNA[Doench_model_weights_doubleonly[1+n, 2]] == Doench_model_weights_doubleonly[1+n, 1]) == TRUE){
          sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- Doench_model_weights_doubleonly[1+n, 3]
        }
        n <- n+1
      }
      GC <- (((str_count(sgRNA_seq[1+g], "G") + str_count(sgRNA_seq[1+g], "C")) / 20) * 100)
      if (isTRUE(GC>70 & 80>=GC) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.202625894258
      }
      if (isTRUE(GC<30 & GC>=20) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.202625894258
      }
      if (isTRUE(GC<20 | GC>80) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.166587751983
      }
      ## Completes the calculation described in Doench et al. (2014) and adds it to the Doench Score list
      Doench_Score[[length(Doench_Score)+1]] <- round(1/(1+exp(-(sum(sgRNA_model_weights)))), digits = 3)
      g <- g+1
    }
    ## Check for off-targets in the genome
    ## Creates Function that converts all sgRNAs into a format readable
    ## by Biostrings
    multiple_DNAString <- function(seqlist){
      DNAString(seqlist)
    }
    Biostrings_sgRNA <- lapply(sgRNA_with_PAM, multiple_DNAString)
    ## Define genome
    usegenome <- genomename
    seqnames <- seqnames(usegenome)
    ## Creates a series of lists to store the incoming mismatch information
    mm0_list <- c()
    mm1_list <- c()
    mm2_list <- c()
    mm3_list <- c()
    mm4_list <- c()
    off_start <- c()
    off_end <- c()
    off_direction <- c()
    off_sgRNAseq <- c()
    off_offseq <- c()
    off_chr <- c()
    off_mismatch <- c()
    ## Decides whether to run a function for a sequence that contains multiple sgRNA or only a single sgRNA
    if (isTRUE(class(Biostrings_sgRNA) == "list")){
      ## Creates lists that will store information about off-target info on the individual sgRNA
      for (pattern in Biostrings_sgRNA) {
        ini_mm0_list <- c()
        ini_mm1_list <- c()
        ini_mm2_list <- c()
        ini_mm3_list <- c()  
        ini_mm4_list <- c()
        ## Creates separate lists for the reverse direction
        rev_ini_mm0_list <- c()
        rev_ini_mm1_list <- c()
        rev_ini_mm2_list <- c()
        rev_ini_mm3_list <- c()
        rev_ini_mm4_list <- c()
        ## Finds the mismatch (MM) info using Biostrings
        for (seqname in seqnames) {
          subject <- usegenome[[seqname]]
          off_info <- matchPattern(pattern, subject, max.mismatch = 3, min.mismatch = 0)
          mis_info <- mismatch(pattern, off_info)
          rev_pattern <- reverseComplement(pattern)
          rev_off_info <- matchPattern(rev_pattern, subject, max.mismatch = 3, min.mismatch = 0)
          rev_mis_info <- mismatch(rev_pattern, rev_off_info)
          ## Puts all forward MM info into lists to later be put into data frame
          if (length(off_info) > 0) {
            for (f in 1:length(off_info)) {
              off_start[[length(off_start)+1]] <- start(off_info)[f]
              off_end[[length(off_end)+1]] <- end(off_info)[f]
              off_direction[[length(off_direction)+1]] <- "+"
              off_chr[[length(off_chr)+1]] <- seqname
              off_mismatch[[length(off_mismatch)+1]] <- length(mis_info[[f]])
              off_sgRNAseq[[length(off_sgRNAseq)+1]] <- as.character(pattern)
              off_offseq[[length(off_offseq)+1]] <- as.character(off_info[[f]])
            }
          }
          ## Adds forward MM counts to lists
          if (length(mis_info) > 0) {
            for (f in 1:length(mis_info)) {
              ini_mm0_list[[length(ini_mm0_list)+1]] <- sum(length(mis_info[[f]]) == 0)
              ini_mm1_list[[length(ini_mm1_list)+1]] <- sum(length(mis_info[[f]]) == 1)
              ini_mm2_list[[length(ini_mm2_list)+1]] <- sum(length(mis_info[[f]]) == 2)
              ini_mm3_list[[length(ini_mm3_list)+1]] <- sum(length(mis_info[[f]]) == 3)
            }
            ## Adds zeroes to all lists if there is no MM info
          } else {
            ini_mm0_list[[length(ini_mm0_list)+1]] <- 0
            ini_mm1_list[[length(ini_mm1_list)+1]] <- 0
            ini_mm2_list[[length(ini_mm2_list)+1]] <- 0
            ini_mm3_list[[length(ini_mm3_list)+1]] <- 0
          }
          ## Puts all reverse MM info into lists to later be put into data frame
          if (length(rev_off_info) > 0) {
            for (f in 1:length(rev_off_info)) {
              off_start[[length(off_start)+1]] <- start(rev_off_info)[f]
              off_end[[length(off_end)+1]] <- end(rev_off_info)[f]
              off_direction[[length(off_direction)+1]] <- "-"
              off_chr[[length(off_chr)+1]] <- seqname
              off_mismatch[[length(off_mismatch)+1]] <- length(rev_mis_info[[f]])
              off_sgRNAseq[[length(off_sgRNAseq)+1]] <- as.character(pattern)
              off_offseq[[length(off_offseq)+1]] <- as.character(rev_off_info[[f]])
            }
          }
          ## Adds reverse MM counts to lists
          if (length(rev_mis_info) > 0) {
            for (f in 1:length(rev_mis_info)) {
              rev_ini_mm0_list[[length(rev_ini_mm0_list)+1]] <- sum(length(rev_mis_info[[f]]) == 0)
              rev_ini_mm1_list[[length(rev_ini_mm1_list)+1]] <- sum(length(rev_mis_info[[f]]) == 1)
              rev_ini_mm2_list[[length(rev_ini_mm2_list)+1]] <- sum(length(rev_mis_info[[f]]) == 2)
              rev_ini_mm3_list[[length(rev_ini_mm3_list)+1]] <- sum(length(rev_mis_info[[f]]) == 3)
            }
            ## Adds zeroes to all lists if there is no MM info
          } else {
            rev_ini_mm0_list[[length(rev_ini_mm0_list)+1]] <- 0
            rev_ini_mm1_list[[length(rev_ini_mm1_list)+1]] <- 0
            rev_ini_mm2_list[[length(rev_ini_mm2_list)+1]] <- 0
            rev_ini_mm3_list[[length(rev_ini_mm3_list)+1]] <- 0
          }
        }
        ## Compiles list of all MM counts for one sgRNA
        mm0_list[[length(mm0_list)+1]] <- (sum(ini_mm0_list) + sum(rev_ini_mm0_list))
        mm1_list[[length(mm1_list)+1]] <- (sum(ini_mm1_list) + sum(rev_ini_mm1_list))
        mm2_list[[length(mm2_list)+1]] <- (sum(ini_mm2_list) + sum(rev_ini_mm2_list))
        mm3_list[[length(mm3_list)+1]] <- (sum(ini_mm3_list) + sum(rev_ini_mm3_list))
      }
      ## Finds the mismatch (MM) info using Biostrings (only used if there is a single sgRNA)
    } else {
      ini_mm0_list <- c()
      ini_mm1_list <- c()
      ini_mm2_list <- c()
      ini_mm3_list <- c()
      rev_ini_mm0_list <- c()
      rev_ini_mm1_list <- c()
      rev_ini_mm2_list <- c()
      rev_ini_mm3_list <- c()
      for (seqname in seqnames) {
        subject <- usegenome[[seqname]]
        off_info <- matchPattern(Biostrings_sgRNA, subject, max.mismatch = 4, min.mismatch = 0)
        mis_info <- mismatch(Biostrings_sgRNA, off_info)
        rev_Biostrings_sgRNA <- reverseComplement(Biostrings_sgRNA)
        rev_off_info <- matchPattern(rev_Biostrings_sgRNA, subject, max.mismatch = 4, min.mismatch = 0)
        rev_mis_info <- mismatch(rev_Biostrings_sgRNA, rev_off_info)
        if (length(mis_info) > 0) {
          for (f in 1:length(mis_info)) {
            ini_mm0_list[[length(ini_mm0_list)+1]] <- sum(length(mis_info[[f]]) == 0)
            ini_mm1_list[[length(ini_mm1_list)+1]] <- sum(length(mis_info[[f]]) == 1)
            ini_mm2_list[[length(ini_mm2_list)+1]] <- sum(length(mis_info[[f]]) == 2)
            ini_mm3_list[[length(ini_mm3_list)+1]] <- sum(length(mis_info[[f]]) == 3)
          }
        } else {
          ini_mm0_list[[length(ini_mm0_list)+1]] <- 0
          ini_mm1_list[[length(ini_mm1_list)+1]] <- 0
          ini_mm2_list[[length(ini_mm2_list)+1]] <- 0
          ini_mm3_list[[length(ini_mm3_list)+1]] <- 0
        }
        if (length(rev_mis_info) > 0) {
          for (f in 1:length(rev_mis_info)) {
            rev_ini_mm0_list[[length(rev_ini_mm0_list)+1]] <- sum(length(rev_mis_info[[f]]) == 0)
            rev_ini_mm1_list[[length(rev_ini_mm1_list)+1]] <- sum(length(rev_mis_info[[f]]) == 1)
            rev_ini_mm2_list[[length(rev_ini_mm2_list)+1]] <- sum(length(rev_mis_info[[f]]) == 2)
            rev_ini_mm3_list[[length(rev_ini_mm3_list)+1]] <- sum(length(rev_mis_info[[f]]) == 3)
          }
        } else {
          rev_ini_mm0_list[[length(rev_ini_mm0_list)+1]] <- 0
          rev_ini_mm1_list[[length(rev_ini_mm1_list)+1]] <- 0
          rev_ini_mm2_list[[length(rev_ini_mm2_list)+1]] <- 0
          rev_ini_mm3_list[[length(rev_ini_mm3_list)+1]] <- 0
        }
      }
      ## Creates the final MM list that includes all sgRNA sequences
      mm0_list[[length(mm0_list)+1]] <- (sum(ini_mm0_list) + sum(rev_ini_mm0_list))
      mm1_list[[length(mm1_list)+1]] <- (sum(ini_mm1_list) + sum(rev_ini_mm1_list))
      mm2_list[[length(mm2_list)+1]] <- (sum(ini_mm2_list) + sum(rev_ini_mm2_list))
      mm3_list[[length(mm3_list)+1]] <- (sum(ini_mm3_list) + sum(rev_ini_mm3_list))
    }
    if ((sum(mm0_list) + sum(mm1_list) + sum(mm2_list) + sum(mm3_list)) == 0) {
      ## Put lists in data frame
      sgRNA_data <- data.frame(sgRNA_seq, sgRNA_PAM, sgRNA_fow_or_rev, sgRNA_start, sgRNA_end, GCinstance, TTTTHomopolymerdetect, Homopolymerdetect, Doench_Score, mm0_list, mm1_list, mm2_list, mm3_list)
      ## Set the names of each column
      colnames(sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content", "TTTT Homopolymer", "Homopolymer", "Doench Score", "MM0", "MM1", "MM2", "MM3")
      sgRNA_data <- sgRNA_data[order(-sgRNA_data$`Doench Score`),]
      sgRNA_data
    } else {
      ## Creates a function that annotates the off-targets called above
      annotate_genome <- function(ochr, ostart, oend, odir, gtfname) {  
        gtf <- import(gtfname)
        seqlevelsStyle(gtf) <- "UCSC"
        seqer <- unlist(ochr)
        starter <- as.numeric(ostart)
        ender <- as.numeric(unlist(oend))
        strander <- unlist(odir)
        off_ranges <- GRanges(seqer, IRanges(starter, ender), strander)
        olaps <- findOverlaps(off_ranges, gtf)
        geneid <- c()
        geneidlist <- c()
        genename <- c()
        genenamelist <- c()
        sequencetype <- c()
        sequencetypelist <- c()
        exonnumber <- c()
        exonnumberlist <- c()
        mcols(off_ranges)$gene_id <- c()
        for (p in 1:length(off_ranges)) {
          if (p %in% queryHits(olaps)) {
            geneid <- mcols(gtf)$gene_id[subjectHits(olaps[which(p == queryHits(olaps))])]
            geneid <- unique(geneid)
            geneidlist[[length(geneidlist)+1]] <- paste(geneid, collapse = ", ")
            genename <- mcols(gtf)$gene_name[subjectHits(olaps[which(p == queryHits(olaps))])]
            genename <- unique(genename)
            genenamelist[[length(genenamelist)+1]] <- paste(genename, collapse = ", ")
            sequencetype <- mcols(gtf)$type[subjectHits(olaps[which(p == queryHits(olaps))])]
            sequencetype <- unique(sequencetype)
            sequencetypelist[[length(sequencetypelist)+1]] <- paste(sequencetype, collapse = ", ")
            exonnumber <- mcols(gtf)$exon_number[subjectHits(olaps[which(p == queryHits(olaps))])]
            exonnumber <- unique(exonnumber)
            exonnumberlist[[length(exonnumberlist)+1]] <- paste(exonnumber, collapse = ", ")
          } else {
            geneidlist[[length(geneidlist)+1]] <- "NA"
            genenamelist[[length(genenamelist)+1]] <- "NA"
            sequencetypelist[[length(sequencetypelist)+1]] <- "NA"
            exonnumberlist[[length(exonnumberlist)+1]] <- "NA"
          }
        }
        mcols(off_ranges)$gene_id <- geneidlist
        more_off_info <- data.frame(geneidlist, genenamelist, sequencetypelist, exonnumberlist)
        more_off_info
      }
      ## Compiles data frame of all off-target annotations
      more_off_info <- annotate_genome(off_chr, off_start, off_end, off_direction, gtfname)
      ## Complies all extra sgRNA info into a separate data frame
      all_offtarget_info <- data.frame(off_sgRNAseq, off_chr, off_start, off_end, off_mismatch, off_direction, off_offseq, more_off_info$geneidlist, more_off_info$genenamelist, more_off_info$sequencetypelist, more_off_info$exonnumberlist)
      colnames(all_offtarget_info) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
      ## Put lists in data frame
      sgRNA_data <- data.frame(sgRNA_seq, sgRNA_PAM, sgRNA_fow_or_rev, sgRNA_start, sgRNA_end, GCinstance, TTTTHomopolymerdetect, Homopolymerdetect, Doench_Score, mm0_list, mm1_list, mm2_list, mm3_list)
      ## Set the names of each column
      colnames(sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content", "TTTT Homopolymer", "Homopolymer", "Doench Score", "MM0", "MM1", "MM2", "MM3")
      sgRNA_data <- sgRNA_data[order(-sgRNA_data$`Doench Score`),]
      data_list <- c("sgRNA_data" = sgRNA_data, "all_offtarget_info" = all_offtarget_info)
      data_list
    }
  } else {
    data_list <- data.frame()
  }
}
