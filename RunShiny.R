## This script is designed to be used with "FindsgRNAfunction_Doench2014.R",
## "Genome_annotation.R", "Doench_Model_Weights_Singleonly.csv",
## and "Doench_Model_Weights_Doubleonly.csv"
## The above files must be in the working directory
##
## example: setwd("C://Users//Dylan//Desktop//SP")
##
##
## To run this program, a computer must have offline Java installed
##
## if packages need to be installed
## install.packages("shiny")
## install.packages("stringr", repos='http://cran.us.r-project.org')
## install.packages("stringi")
## install.packages("mailR", dep = TRUE)
## install.packages("rJava")
## source("https://bioconductor.org/biocLite.R")
## biocLite("Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "AnnotationHub")
## biocLite("BSgenome.Scerevisiae.UCSC.sacCer2") 
##
## In order to access shinyapps server
## install.packages('rsconnect')
## Additionally, to access the shinyapps server your desktop must
## be authorized
##
## Activate packages
##library(rsconnect)
library(shiny)
library(stringr)
library(rJava)
library(mailR)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
##library(BSgenome.Mmusculus.UCSC.mm10)
##library(BSgenome.Dmelanogaster.UCSC.dm6)
##library(BSgenome.Ptroglodytes.UCSC.panTro5)
##library(BSgenome.Ecoli.NCBI.20080805)
##library(BSgenome.Celegans.UCSC.ce11)
##library(BSgenome.Rnorvegicus.UCSC.rn6)
##library(BSgenome.Athaliana.TAIR.04232008)

ui <- fluidPage(
  navbarPage("Cas9 Guide Designer",
    tabPanel("sgRNA Designer",
      titlePanel("sgRNA Designer"),
      
      sidebarLayout(
        sidebarPanel(
          textInput("sequence", "Sequence", placeholder = "Paste sequence here"),
          actionButton("run", "Find sgRNA"),
          selectInput("genome_select", "Select Genome",
                      c("Homo sapiens (UCSC.hg19)" = "Hsapiens",
                        "Saccharomyces cerevisiae (UCSC.sacCer2)" = "Scerevisiae",
                        "Mus musculus (UCSC.mm10)" = "Mmusculus",
                        "Drosphila melanogaster (UCSC.dm6)" = "Dmelanogaster",
                        "Pan troglodytes (UCSC.panTro5)" = "Ptroglodytes",
                        "Escheria coli (NCBI.20080805)" = "Ecoli",
                        "Caenorhabditis elegans (UCSC.ce11)" = "Celegans",
                        "Ratus norvegicus (UCSC.rn6)" = "Rnorvegicus",
                        "Arabidopsis thaliana (TAIR.04232008)" = "Athaliana")),
          selectInput("scoring_select", "Select Scoring Method (WIP)",
                      c("Doench Rule Set 1 (2014)" = "Doench_2014",
                        "Doench Rule Set 2 (2016)" = "Doench_2016"),
                      selected = "Doench_2014"),
          checkboxInput("email", "Send me an email with the results", value = FALSE),
          tags$div(id = "placeholder2")
        ),
        mainPanel(  
          tags$div(id = "placeholder3"),
          dataTableOutput("sgRNA_data"),
          tags$div(id = "placeholder4"),
          dataTableOutput("offtarget_data")
        )
      )
    ),
    tabPanel("Options",
      titlePanel("Options")
    ),
    tabPanel("About",
      titlePanel("About"),
      column(8, "The Cas9 Guide Finder designs guide RNA sequences (sgRNA) for Cas9 DNA editing.
             To begin, enter a sequence into the sequence box, select a genome to search for
             Off-Targets, and click find sgRNA."),
      column(8, "Note about Off-target calling in large genomes: When using a large genome like
             Homo sapiens, we reccomend using sequences under 500 base pairs. The time it can take
             to search these genomes can be multiple hours if too many sgRNA are generated.")
    )
  )
)

server <- function(input, output) {
  ## Sources the file that contains the function for sgRNA design
  ## and the file for off-target searches
  source("FindsgRNAfunction_Doench2014.R")
  source("Offtargetsearchfunction.R")
  
  ## Creates a list of reactive values that allows the program to
  ## update only when the action button is pressed
  maindf <- reactiveValues(data = NULL)
  offtargetdf <- reactiveValues(data = NULL)
  
  ## Runs the sgRNA_design function when the action button is pressed
  observeEvent(input$run, {
    # Check to see if input is valid
    sequence <- paste(input$'sequence', collapse = "")
    sequence <- str_replace_all(sequence, fixed(" "), "")
    if (isTRUE(try(class(DNAString(sequence)) == "DNAString"))) {
      # Create a Progress object
      designprogress <- shiny::Progress$new()
      designprogress$set(message = "Finding sgRNA", value = 0, detail = "This may take a while")
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(designprogress$close())
      all_data <- sgRNA_design(usersequence = input$'sequence', genomename = input$'genome_select', designprogress)
      if ((length(all_data) == 0) == FALSE) {
        int_sgRNA_data <- data.frame(all_data[1:13])
        colnames(int_sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content",
                                    "TTTT Homopolymer", "Homopolymer", "Doench Score", "MM0", "MM1", "MM2", "MM3")
        if (input$run == 1) {
          insertUI(
            selector = "#placeholder3",
            where = "afterEnd",
            ui = tags$div(id = 'sgRNAdftext',
                          titlePanel("sgRNA Table"),
                          downloadButton("Download_sgRNA", "Download sgRNA")
            )
          )
        }
        maindf$sgRNA_data <- int_sgRNA_data
        int_offtarget_data <- data.frame(all_data[14:24])
        colnames(int_offtarget_data) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction",
                                          "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
        if (input$run == 1) {
          insertUI(
            selector = "#placeholder4",
            where = "afterEnd",
            ui = tags$div(id = 'sgRNAofftext',
                          titlePanel("Additional Off-target Information"),
                          downloadButton("Download_off", "Download Off-Targets")
            )
          )
        }
        offtargetdf$data <- int_offtarget_data
        if (input$email == TRUE) {
          recipient <- input$recipientbox
          project <- input$project_title
          write.csv(maindf$sgRNA_data, file = "sgRNA data.csv")
          write.csv(offtargetdf$data, file = "Off-Target data.csv")
          ## Sends an email to an address entered into the UI
          send.mail(from = "<uml.sgRNA.design@gmail.com>",
                    to = paste("<", recipient, ">", sep = ""),
                    subject = paste("sgRNA Design data for", project),
                    body = "Your sgRNA data is attached",
                    smtp = list(host.name = "smtp.gmail.com", port = 465,
                                user.name = "uml.sgRNA.design@gmail.com",
                                passwd = "tg5pvm19zq", ssl = TRUE),
                    authenticate = TRUE,
                    send = TRUE,
                    attach.files = "sgRNA_data.csv", "Off-Target data.csv")
        }
      } else {
        showModal(modalDialog(
          title = "Error",
          "Error! No sgRNA were generated from sequence"
        ))
      }
    } else {
      showModal(modalDialog(
        title = "Error",
        "Error! Sequence may contain unsupported characters"
      ))
    }
  })
  
  output$Download_sgRNA <- downloadHandler(
    filename = function(){"sgRNA.csv"},
    content = function(file) {
      write.csv(maindf$sgRNA_data, file, row.names = TRUE)
    }
  )
  
  output$Download_off <- downloadHandler(
    filename = function(){"Offtarget.csv"},
    content = function(file) {
      write.csv(offtargetdf$data, file, row.names = TRUE)
    }
  )
  
  ## Reactively outputs an sgRNA table when the function is complete
  output$sgRNA_data <- renderDataTable(maindf$sgRNA_data)
  output$offtarget_data <- renderDataTable(offtargetdf$data)
  
  ## Add email input to the UI
  observeEvent(input$email, {
    if (input$email == TRUE) {
      insertUI(
        selector = "#placeholder2",
        where = "afterEnd",
        ui = tags$div(id = 'emailinfo',
                      textInput("recipientbox", NULL, placeholder = "Enter your email here"),
                      textInput("project_title", NULL, placeholder = "Enter a name for your project")
        )
      )
    } else {
      removeUI(
        selector = 'div#emailinfo',
        multiple = TRUE
      )
    }
  })
}

shinyApp(ui=ui, server=server)