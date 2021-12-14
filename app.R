library(shiny)
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))

load('app/sample_meta.Rdata')


sva_counts <-  read_tsv('app/sva_counts.tsv.gz')

sample_meta_D <- sample_meta %>% filter(Sample %in% colnames(sva_counts)) %>%
  dplyr::select(Sample:Section, Layout:Fusion) %>%
  mutate(S2 = case_when(Section == 'OF' ~ 'OF', TRUE ~ 'OC')) %>%
  unique()

box_maker <- function(table, genes, section = c('OF','OC'), type = 'temporal'){
  if ('A1CF' %in% row.names(table)){
    table <- table %>%
      as_tibble(rownames = 'Gene')
  }
  if (type == 'temporal'){
    plot <- table %>%
      pivot_longer(-Gene, names_to = 'Sample', values_to = 'Expression') %>%
      mutate(Sample = gsub('_.*|.CEL.*','',Sample)) %>%
      left_join(sample_meta_D) %>%
      mutate(S2 = case_when(Section == 'OF' ~ 'OF',TRUE ~ 'OC')) %>%
      filter(Gene %in% genes, S2 %in% section) %>%
      #filter(Gene %in% row.names(top.table_OF_AD %>% head(10))) %>%
      mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
      mutate(OrgTech = paste(Organism, Technology, sep = '_')) %>%
      ggplot(aes(x=Fusion, y=Expression, color = Organism, shape = Technology)) +
      # geom_boxplot(aes(group = Fusion), color = 'Black', outlier.colour = NA) +
      # geom_boxplot(aes(group = interaction(Organism,Fusion)), outlier.colour = NA) +
      ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7) +
      cowplot::theme_cowplot() +
      facet_grid(~Gene + S2, scales = 'free_y') +
      ggsci::scale_color_aaas() +
      ylab('log2 (corrected counts)') +
      stat_summary(fun=mean, geom = 'line', aes(group = OrgTech, color = Organism)) }
  else {
    plot <- table %>%
      pivot_longer(-Gene, names_to = 'Sample', values_to = 'Expression') %>%
      mutate(Sample = gsub('_.*|.CEL.*','',Sample)) %>%
      left_join(sample_meta_D) %>%
      mutate(S2 = case_when(Section == 'OF' ~ 'OF',TRUE ~ 'OC')) %>%
      filter(Gene %in% genes, S2 %in% section) %>%
      mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
      filter(Fusion == 'During') %>%
      mutate(OrgTech = paste(Organism, Technology, sep = '_')) %>%
      ggplot(aes(x=S2, y=Expression, color = Organism, shape = Technology)) +
      # geom_boxplot(aes(group = Fusion), color = 'Black', outlier.colour = NA) +
      # geom_boxplot(aes(group = interaction(Organism,Fusion)), outlier.colour = NA) +
      ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7) +
      cowplot::theme_cowplot() +
      ggsci::scale_color_aaas() +
      ylab('log2 (corrected counts)') +
      xlab('Section') +
      stat_summary(fun=mean, geom = 'line', aes(group = OrgTech, color = Organism)) + facet_wrap(~Gene)
  }
  plot
}


ui <- fluidPage(

  selectizeInput(inputId = "Gene",
                 "Gene Name",
                 choices = NULL,
                 selected = 'NTN1'),
  plotOutput(outputId = "plotter")
)


server <- function(input, output, session) {

  updateSelectizeInput(session, 'Gene',
                       choices = sva_counts$Gene,
                       selected = 'MITF',
                       server = TRUE)



  output$plotter <- renderPlot({

    req(input$Gene)

    plot_temporal <- box_maker(sva_counts,
                               genes = input$Gene,
                               section = c('OF','OC'), type = 'temporal')


    plotOC_vs_OF <- box_maker(sva_counts,
                              genes = input$Gene,
                              section = c('OF','OC'), type = 'stage')
    cowplot::plot_grid(plot_temporal + ggtitle('Stage') + theme(legend.position = "none"), plotOC_vs_OF + ggtitle('During Fusion'), nrow = 1, rel_widths = c(1.5,1))

  }

  )
}

shinyApp(ui = ui, server = server)
