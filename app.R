#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyLP)
library(data.table)
source('src/data_visualization.R')
load('data/INFIMA-shiny-data.RData')
# prior, infima, gwas, input_queries_small

# Define UI for application that draws a histogram
ui <- navbarPage(
   'INFIMA',
   theme = shinytheme("sandstone"),
   tabPanel(
     'Home', icon = icon('home'),
     jumbotron(
       'Welcome to INFIMA!',
       'INFIMA is an R package for the Integrative Fine-Mapping with Model Organism Multi-Omics Data. INFIMA utilizes the diversity outbred (DO) mice population as the model organsim. The major usage of the INFIMA package is to fine-map the eQTL markers in DO studies (DO-eQTL). INFIMA implements an empirical Bayes model which quantifies how well each non-coding SNP explains the observed DO allelic patern through consistency of founder mice RNA-seq data, founder mice ATAC-seq data (including the existence and consistency of a footprint) with the observed allelic pattern.',
       button = F),
     fluidRow(
       column(
         width = 8,
         panel_div(
           class_type = 'primary',
           panel_title = 'What do we obtain from INFIMA?',
           HTML('Given a library of mouse SNPs, the corresponding local ATAC-seq signals, DO-eQTL data as well as the gene expression in founder mice, INFIMA provides the following functionalities: <br> (1) Fine-mapping DO-eQTLs and estimate SNP-level posterior probabilities. <br> (2) Linking SNPs or the local ATAC-seq peaks to effector genes.')
         ),
         
         
         panel_div(
           class_type = 'success',
           panel_title = 'How do we utilize INFIMA predictions for human GWAS?',
           content = 'INFIMA results from DO studies can be mapped to human orthologs using peak-based lift-over strategies, which provides putative effector genes of human GWAS SNPs. See the "GWAS Effector Genes" page for the results validated by promoter capture Hi-C data.'
         ),
         
         panel_div(
           class_type = 'default',
           panel_title = 'About',
           HTML("Email us: <a href='mailto:cdong@stat.wisc.edu'>Chenyang Dong</a>, <a href='mailto:keles@stat.wisc.edu'>Sunduz Keles</a><br><br>Copyright (c) Chenyang Dong and Sunduz Keles")
         )
         
       ),
       column(
         width = 4,
         h3('INFIMA model overview'),
         imageOutput('infima_overview')         
       )
     )
   ),
   
   
   tabPanel(
     'GWAS Effector Genes', icon = icon('table'),
     
     fluidPage(
       titlePanel('Effector genes of islet human GWAS SNPs validated by pcHi-C'),
       fluidRow(
         column(
           width = 6,
           selectInput('select_trait', label = h4('Pancreatic islet traits'), 
                       choices = c('All', sort(unique(gwas$trait))), 
                       selected = 'All')           
         )
       ),
       
       mainPanel(
         DT::dataTableOutput('gwas_data')
       )
       
     )
   ),
   
   
   ### The panel showing the infima data
   tabPanel(
     'INFIMA Predictions', icon = icon('search'),
     
     fluidPage(
       titlePanel('INFIMA Predictions'),
       
       fluidRow(
         column(
           width = 12,
           h4('DO mouse QTL markers were fine-mapped to local-ATAC-QTLs by using INFIMA.'),
           h4('Click any row of the following table to see the relevant multi-omics data in the "Data Visualization" page.')
         )
       ),
       
       mainPanel(
         DT::dataTableOutput('infima_data')
       )
     )
   ),
   
   
   ### The panel showing the infima data
   tabPanel(
     'Data Visualization', icon = icon('bar-chart-o'),
     
     fluidPage(
       titlePanel('Data Visualization'),
       sidebarLayout(
         
         sidebarPanel(
           selectInput('plot_which', 
                       label = h4('Select data to plot'), 
                       choices = list('DO allele effect' = 1, 
                                      'Local ATAC-seq signal' = 2, 
                                      'Founder gene expression' = 3,
                                      'Founder allele effect & edit distance' = 4), 
                       selected = 1),
           h4('Four priors:'),
           h5('cor(A,E): correlation between ATAC-seq and founder allele effects;'),
           h5('cor(A,B): correlation between ATAC-seq and founder gene expression;'),
           h5('footprint: in-silico mutation motif and footprint indicator;'),
           h5('distance: distance score.'),
           plotOutput('plot_prior')
         ),
         
         
         mainPanel(
           plotOutput('plot_data', width = "700px"),
           h5('Selected INFIMA prediction'),
           tableOutput('infima_data_selected'),
           h5('Genotype of the local-ATAC-QTL'),
           tableOutput('genotype'),
           h5('Posterior probability of the local-ATAC-QTL'),
           tableOutput('pprob')
         )
         
         
       )
     )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  output$infima_overview <- renderImage({
    width <- session$clientData$output_infima_overview_width
    height <- session$clientData$output_infima_overview_height
    len <- min(width, height)*1.2
    list(src = 'figure/infima.png',
         width = len,
         height = len/40*45,
         alt = 'INFIMA overview.')
  }, deleteFile = F)
  
  output$gwas_data <- DT::renderDataTable({
    if(input$select_trait != 'All'){
      DT::datatable(gwas[gwas$trait == input$select_trait,], 
                    selection = 'single',
                    options = list(lengthMenu = c(20, 50, 100)))
    }
    else{
      DT::datatable(gwas, 
                    selection = 'single',
                    options = list(lengthMenu = c(20, 50, 100)))
    }
  })
  
  output$infima_data <- DT::renderDataTable({
    DT::datatable(infima,
                  selection = 'single',
                  options = list(lengthMenu = c(20, 50, 100)))
  })
  
  observe({
    req(input$infima_data_rows_selected)
    infima_selected <- infima[input$infima_data_rows_selected,]
    prior_selected <- prior[input$infima_data_rows_selected,]
    input_query_selected <- input_queries_small[input$infima_data_rows_selected][[1]]
    
    output$infima_data_selected <- renderTable({
      infima_selected
    })
    
    output$genotype <- renderTable({
      genotype_data <- as.data.frame(t(as.integer(input_query_selected$genotype)))
      colnames(genotype_data) <- c('129', 'AJ', 'B6', 'CAST',
                                   'NOD', 'NZO', 'PWK', 'WSB')
      genotype_data
    })
    
    output$pprob <- renderTable({
      pprob_data <- prior_selected[, c('p', 'k', 'Z', 'Z.rs')]
      pprob_data$p <- as.integer(pprob_data$p)
      pprob_data$k <- as.integer(pprob_data$k)
      colnames(pprob_data) <- c('Total # of candidates',
                                'Candidates in the credible set',
                                'Posterior probability',
                                'Posterior probability (rank score)')
      pprob_data
    })
    
    
    output$plot_prior <- renderPlot(
      plot_prior(prior_selected)
    )
    
    
    output$plot_data <- renderPlot(
      plot_input(input_query_selected, option = input$plot_which)
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

