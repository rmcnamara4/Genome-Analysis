# ---------------------------------------------------------------------------------------------------------------------

# Import libraries

library(magrittr)
library(readxl)
library(xlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(devtools)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(ggthemes)
library(scales)
library(ggrepel)
library(ggforce)
library(shiny)

load_var <- c("aedes_aegypti", "human", "dengue_1_isolate", "dengue_2_isolate", "dengue_3_isolate", "dengue_4_isolate", "dengue_consensus", "dengue_2_replicon_strain",
              "all_virus")

species <- c("Aedes aegypti", "Human", "CHIKV", "DV1", "HHV1", "HIV1", "FLUAV", "GU280", "ZIKV", "Dengue 1 Isolate", "Dengue 2 Isolate", "Dengue 3 Isolate", 
             "Dengue 4 Isolate", "Dengue 1 Consensus", "Dengue 2 Consensus", "Dengue 3 Consensus", "Dengue 4 Consensus", "Dengue 2 Replicon", 
             "Dengue 2 - Strain TH36")

data_names <- c("aedes_aegypti", "human", "CHIKV", "DV1", "HHV1", "HIV1", "FLUAV", "GU280", "ZIKV", "Dengue_1_isolate", "Dengue_2_isolate", "Dengue_3_isolate", "Dengue_4_isolate",
                "Dengue_1", "Dengue_2", "Dengue_3", "Dengue_4", "Dengue_2_replicon", "Dengue_2_strain_TH36")

names(data_names) <- species

data_type <- c("Nucleotides", "Codons", "Amino_Acids", "RSCU")

nuc_vars <- c("A", "G", "C", "T", "A1", "G1", "C1", "T1", "A2", "G2", "C2", "T2", "A3", "G3", "C3", "T3", "AT", "GC", "AT3", "GC3")

synonymous_codons <- list(
  Phe = c("TTT", "TTC"),
  Leu = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"), 
  Ile = c("ATT", "ATC", "ATA"), 
  Val = c("GTT", "GTC", "GTA", "GTG"), 
  Ser = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"), 
  Pro = c("CCT", "CCC", "CCA", "CCG"), 
  Thr = c("ACT", "ACC", "ACA", "ACG"), 
  Ala = c("GCT", "GCC", "GCA", "GCG"), 
  Tyr = c("TAT", "TAC"), 
  His = c("CAT", "CAC"), 
  Gln = c("CAA", "CAG"), 
  Asn = c("AAT", "AAC"), 
  Lys = c("AAA", "AAG"), 
  Asp = c("GAT", "GAC"), 
  Glu = c("GAA", "GAG"), 
  Cys = c("TGT", "TGC"), 
  Arg = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"), 
  Gly = c("GGT", "GGC", "GGA", "GGG"),
  Met = c("ATG"), 
  Trp = c("TGG"), 
  Stp = c("TAG", "TAA", "TGA")
)

nucleotides_all <- data.frame()
codons_all <- data.frame()
amino_acids_all <- data.frame()
rscu_all <- data.frame()
for (d in data_type) {
  for (s in load_var) {
    assign(paste0(s, "_", tolower(d)), 
           read_excel(paste0("/n/projects/rm2498/Data_Files/Full_Species_Data/", s, "_data.xls"), 
                      sheet = d)) 
    assign(paste0(tolower(d), "_all"),
           bind_rows(get(paste0(s, "_", tolower(d))), get(paste0(tolower(d), "_all"))))
    
  }
}

nucleotides_all <- select(nucleotides_all, -c(gene_ID, coding)) %>%
  mutate(
    AT = .[["A"]] + .[["T"]], 
    GC = .[["G"]] + .[["C"]],
    AT3 = .[["A3"]] + .[["T3"]],
    GC3 = .[["G3"]] + .[["C3"]]
  )
codons_all <- select(codons_all, -c(species, gene_ID, coding))
amino_acids_all <- select(amino_acids_all, -c(species, gene_ID, coding))
rscu_all <- select(rscu_all, -c(species, gene_ID, coding))

names(rscu_all) <- paste0(names(rscu_all), "_rscu")

all_data <- cbind(nucleotides_all, codons_all, amino_acids_all, rscu_all)
all_data[is.na(all_data)] <- 0

all_data <- all_data %>%
  group_by(species) %>%
  summarise_all(median)

codon_optimality <- read_excel("/n/projects/rm2498/Data_Files/codon_optimality.xls") %>%
  select(1, 3:5) %>%
  mutate(
    score = (.[[2]] + .[[3]] + .[[4]]) / 3
  ) %>%
  select(1, 5)
codon_optimality <- codon_optimality[order(codon_optimality[[2]]), ]
codon_optimality <- rbind(codon_optimality, c("TAA", NA), c("TGA", NA), c("TAG", NA))
names(codon_optimality)[1] <- "codon"

ordered_codons <- codon_optimality[[1]]

ordered_aa <- read_excel("/n/projects/rm2498/Data_Files/human_order_of_aa.xls")
ordered_aa <- ordered_aa[[2]]

ui <- fluidPage(
  
  titlePanel(
    
    strong("Species/Virus Analysis")
    
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      strong("Plot Options"), 
      
      br(),
      br(), 
      
      selectInput(
        inputId = "plot_data", 
        label = "Choose which data you want to plot:", 
        choices = data_type,
        multiple = FALSE
      ),
      
      br(), 
      br(), 
      
      checkboxGroupInput(
        inputId = "species_choice", 
        label = "Choose which species to plot:", 
        choices = data_names,
        inline = TRUE
      ),
      
      checkboxInput(
        inputId = "fold_change", 
        label = "Plot as Fold Change", 
        value = FALSE
      ),
      
      conditionalPanel(
        condition = "input.fold_change == 1",
        radioButtons(
          inputId = "plot_type", 
          label = "Choose a type of plot:", 
          choices = c("Bar Plot", "Heat Map"), 
          inline = TRUE
        )
      ),
      
      conditionalPanel(
        condition = "input.plot_data == 'Nucleotides' && input.plot_type == 'Bar Plot'", 
        checkboxGroupInput(
          inputId = "nucleotide_variables", 
          label = "Choose which variables to plot:", 
          choices = nuc_vars,
          inline = TRUE
        )
      ),
      
      conditionalPanel(
        condition = "input.plot_data == 'Codons' && input.plot_type == 'Bar Plot'",
        selectInput(
          inputId = "codon_variables", 
          label = "Choose which variables to plot:", 
          choices = ordered_codons,
          multiple = TRUE
        )
      ), 
      
      conditionalPanel(
        condition = "input.plot_data == 'Amino_Acids' && input.plot_type == 'Bar Plot'",
        selectInput(
          inputId = "aa_variables", 
          label = "Choose which variables to plot:",
          choices = ordered_aa_human, 
          multiple = TRUE
        )
      ),
      
      conditionalPanel(
        condition = "input.plot_data == 'RSCU' && input.plot_type == 'Bar Plot'", 
        selectInput(
          inputId = "rscu_variables", 
          label = "Choose which variable to plot:",
          choices = ordered_aa, 
          multiple = FALSE
        )
      ),
      
      conditionalPanel(
        condition = "input.fold_change == 1", 
        uiOutput("chosen_species")
      ), 
      
      actionButton("submit", "PLOT")
      
    ),
    
    mainPanel(
      
      plotOutput(
        outputId = "main_plot"
      )
      
    )
  )
)

server <- function(input, output, session) {
  
  output$chosen_species <- renderUI({
    
    selectInput(
      inputId = "fold_change_variable",
      label = "Choose which species to base fold change on:", 
      choices = setdiff(data_names, input$species_choice), 
      multiple = FALSE
    )
    
  })
  
  output$main_plot <- renderPlot({
    
    if (input$submit == 0) {
      
      return()
      
    }
    
    isolate({
      
      if (input$fold_change == TRUE) {
        
        fold_var <- which(all_data$species %in% input$fold_change_variable)
        plot_vals <- apply(all_data[-1], 2, function(x) x / x[fold_var])
        plot_vals <- apply(plot_vals, 2, function(x) log2(x)) %>%
          as.data.frame(.) %>%
          add_column(species = all_data$species, .before = 1)
        
        text <- "Median Fold Change and log2"
        text2 <- paste0(" with Respect to ", str_to_title(input$fold_change_variable))
        
        # ordered_aa_places <- plot_vals[fold_var, ] %>%
        #   select(names(synonymous_codons)) %>%
        #   order()
        # ordered_aa <- names(plot_vals[fold_var, ] %>% select(names(synonymous_codons)))[ordered_aa_places]
        
      } else {
        
        plot_vals <- all_data
        
        text <- "Median"
        text2 <- ""
        
        #ordered_aa <- ordered_aa_human
        
      }
      
      selected_data <- plot_vals %>%
        filter(species %in% input$species_choice)
      
      if (input$plot_data == "Nucleotides") {
        
        if (input$plot_type == "Bar Plot") {
          
          selected_data <- selected_data %>%
            select(species, input$nucleotide_variables)
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = variable, y = value, fill = species)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(
              x = "Nucleotides", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Nucleotide Composition", text2),
              fill = "Species"
            )
          
        } else {
          
          selected_data <- selected_data %>%
            select(species, all_of(nuc_vars))
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = variable, y = species, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradientn(colors = c("plum4", "plum4", "plum4", "white", "springgreen4", "springgreen4"), 
                                 values = rescale(c(-4, -2, -1.5, 0, 1, 2), to = c(0, 1))) +
            theme(aspect.ratio = 1) +
            theme(axis.text.x = element_blank(),      
                  axis.ticks.x = element_blank(),     
                  strip.background = element_blank()) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            facet_grid(~ variable, scales = "free", space = "free", switch = "y") +
            labs(
              x = "Nucleotides", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Nucleotide Composition", text2), 
              fill = "Species"
            )
          
        }
        
      }
      
      if (input$plot_data == "Codons") {
        
        if (input$plot_type == "Bar Plot") {
          
          selected_data <- selected_data %>%
            select(species, input$codon_variables)
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_codons), y = value, fill = species)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(
              x = "Codons", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Codon Composition", text2),
              fill = "Species"
            )
            
          
        } else {
          
          selected_data <- selected_data %>%
            select(species, all_of(ordered_codons))
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_codons), y = species, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradientn(colors = c("plum4", "plum4", "plum4", "white", "springgreen4", "springgreen4"), 
                                 values = rescale(c(-4, -2, -1.5, 0, 1, 2), to = c(0, 1))) +
            theme(aspect.ratio = 1) +
            theme(axis.text.x = element_blank(),      
                  axis.ticks.x = element_blank(),     
                  strip.background = element_blank()) +
            facet_grid(~ factor(variable, levels = ordered_codons), scales = "free", space = "free", switch = "y") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(
              x = "Codons", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Codon Composition", text2), 
              fill = "Species"
            )
          
        }
        
      }
      
      if (input$plot_data == "Amino_Acids") {
        
        if (input$plot_type == "Bar Plot") {
          
          selected_data <- selected_data %>%
            select(species, input$aa_variables)
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_aa), y = value, fill = species)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(
              x = "Amino Acids", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Amino Acid Composition", text2),
              fill = "Species"
            )
          
        } else {
          
          selected_data <- selected_data %>%
            select(species, all_of(ordered_aa))
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_aa), y = species, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradientn(colors = c("plum4", "plum4", "plum4", "white", "springgreen4", "springgreen4"), 
                                 values = rescale(c(-4, -2, -1.5, 0, 1, 2), to = c(0, 1))) +
            theme(aspect.ratio = 1) + 
            theme(axis.text.x = element_blank(),      
                  axis.ticks.x = element_blank(),     
                  strip.background = element_blank()) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            facet_grid(~ factor(variable, levels = ordered_aa), scales = "free", space = "free", switch = "y") +
            labs(
              x = "Amino Acids", 
              y = paste0(text, " of Percent Composition"), 
              title = paste0(text, " of Amino Acid Composition", text2), 
              fill = "Species"
            )
          
        }
        
      }
      
      if (input$plot_data == "RSCU") {
        
        if (input$plot_type == "Bar Plot") {
          
          rscu_choices <- vector()
          for (a in input$rscu_variables) {
            rscu_choices <- append(rscu_choices, synonymous_codons[[a]])
          }
          
          selected_data <- selected_data %>%
            select(species, paste0(rscu_choices, "_rscu"))
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          melted_data$variable <- str_replace_all(melted_data$variable, "_rscu", "")
          
          aa <- vector()
          optimality <- vector()
          for (c in melted_data$variable) {
            aa <- append(aa, names(synonymous_codons)[which(grepl(c, synonymous_codons))])
            optimality <- append(optimality, codon_optimality[[which(codon_optimality[[1]] == c), 2]])
          }
          
          melted_data$aa <- aa
          melted_data$optimality <- optimality
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_codons), y = value, fill = species)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(
              x = "Codons", 
              y = paste0(text, " of RSCU of ", aa[1]), 
              title = paste0(text, " of RSCU of ", aa[1], text2),
              fill = "Species"
            )
          
        } else {
          
          selected_data <- selected_data %>%
            select(species, paste0(ordered_codons, "_rscu"))
          
          melted_data <- reshape2::melt(selected_data, id = "species")
          
          melted_data$variable <- str_replace_all(melted_data$variable, "_rscu", "")
          
          aa <- vector()
          for (c in melted_data$variable) {
            aa <- append(aa, names(synonymous_codons)[which(grepl(c, synonymous_codons))])
          }
          
          melted_data$aa <- aa
          
          figure <- melted_data %>%
            ggplot(aes(x = factor(variable, levels = ordered_codons), y = species, fill = value, label = aa)) +
            geom_tile(color = "white") +
            scale_fill_gradientn(colors = c("plum4", "plum4", "plum4", "white", "springgreen4", "springgreen4"), 
                                 values = rescale(c(-4, -2, -1.5, 0, 1, 2), to = c(0, 1))) +
            theme(aspect.ratio = 1) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            facet_grid(~ factor(aa, levels = ordered_aa), scales = "free", space = "free") +
            labs(
              x = "Codons", 
              y = paste0(text, " of RSCU"), 
              title = paste0(text, " of RSCU", text2), 
              fill = "Species"
            )
          
        }
        
      }
      
    })
  
    figure
    
  })
  
}

shinyApp(ui = ui, server = server) 







