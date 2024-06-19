
# VSC user day 130624 ----

# packrat: Helps manage package dependencies in R projects
#install.packages("packrat")
#packrat::init()
#library(packrat)
#.libPaths() # check the library paths


# renv: A modern replacement for packrat
install.packages("renv")
renv::init()

# Use tools like lintr to analyze your code for potential issues and to 
# identify dependencies.
install.packages("lintr")
lintr::lint("data_assembling.R") # for aesthetic reason only

# these 3 only help on packages management or aesthetics.
# should look at R profillers instead.

# Load packages ----
library(readxl)
library(tidyverse)
library(reshape)
library(stringr)
library(ggpubr)
library(randomcoloR)
library(pracma) # for calculating AUC trapezoidal

## Set seed for reproducible research
set.seed(12345)

# Extracting data ----

## define a function to extract the data ----

setwd("./Datasets")

extract_data <- function(filename) {
  
  temp_name <- str_split(file_path_sans_ext(filename), "_")[[1]][2]
  
  # Second chunk: Capitalize the string and extract the first word
  capitalized_name <- str_to_title(temp_name)
  NAME <- str_split(capitalized_name, " ")[[1]][1] # extract only the drug name, with the first letter capitalized
  
  bio_rep = extract_numeric(word(filename, 4)) # extract the biological replicate number
  
  output_name <- str_to_title(str_split(file_path_sans_ext(filename), "_")[[1]][2]) # extract only the part after "_" and without ".xlsx""
  #################    slice datasets by celline and ATUX concentration    ####
  # WT cellline
  data1 <- read_excel(path = filename, sheet = "WT", range = "A2:J5"  , na = "") #ATUX=6 µM
  data2 <- read_excel(path = filename, sheet = "WT", range = "A7:J10"  , na = "") #ATUX=4,5 µM
  data3 <- read_excel(path = filename, sheet = "WT", range = "A14:J17"  , na = "") #ATUX=3 µM
  data4 <- read_excel(path = filename, sheet = "WT", range = "A19:J22"  , na = "") #ATUX=0 µM
 
  # R183W cellline
  data5 <- read_excel(path = filename, sheet = "R183W", range = "A2:J5"  , na = "") #ATUX=6 µM
  data6 <- read_excel(path = filename, sheet = "R183W", range = "A7:J10"  , na = "") #ATUX=4,5 µM
  data7 <- read_excel(path = filename, sheet = "R183W", range = "A14:J17"  , na = "") #ATUX=3 µM
  data8 <- read_excel(path = filename, sheet = "R183W", range = "A19:J22"  , na = "") #ATUX=0 µM
  
  #################    assemble WT cell_line data   ############################
  data1 <- data1  |> 
    gather(key = "DRUG_CONC", value = "DV", -"ATUX (6 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (6 µM)`),
           ADDITIVE_CONC = 6,
           CELLINE = 1,
           CELLINE_NAME = "WT",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO", 0, as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (6 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data2 <- data2 |> 
    gather(key = "DRUG_CONC", value = "DV", -"ATUX (4,5 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (4,5 µM)`),                                
           ADDITIVE_CONC = 4.5,
           CELLINE = 1,
           CELLINE_NAME = "WT",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (4,5 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data3 <- data3 |> 
    gather(key = "DRUG_CONC", value = "DV", -"ATUX (3 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (3 µM)`),                                
           ADDITIVE_CONC = 3,
           CELLINE = 1,
           CELLINE_NAME = "WT",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (3 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data4 <- data4 |> 
    gather(key = "DRUG_CONC", value = "DV", -"DMSO...1") |>  
    mutate(REP = extract_numeric(`DMSO...1`),                                
           ADDITIVE_CONC = 0,
           CELLINE = 1,
           CELLINE_NAME = "WT",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO...10",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`DMSO...1`) |>  arrange(REP,DRUG_CONC)
  
  ##################    assemble R183W cell_line data   ############################
  data5 <- data5  |> 
    gather(key = "DRUG_CONC", value = "DV", -"ATUX (6 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (6 µM)`),
           ADDITIVE_CONC = 6,
           CELLINE = 2,
           CELLINE_NAME = "R183W",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (6 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data6 <- data6  |> 
    gather(key = "DRUG_CONC",value = "DV", -"ATUX (4,5 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (4,5 µM)`),                                
           ADDITIVE_CONC = 4.5,
           CELLINE = 2,
           CELLINE_NAME = "R183W",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (4,5 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data7 <- data7  |> 
    gather(key = "DRUG_CONC",value = "DV", -"ATUX (3 µM)") |>  
    mutate(REP = extract_numeric(`ATUX (3 µM)`),                                
           ADDITIVE_CONC = 3,
           CELLINE = 2,
           CELLINE_NAME = "R183W",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`ATUX (3 µM)`) |>  arrange(REP,DRUG_CONC)
  
  data8 <- data8  |> 
    gather(key = "DRUG_CONC",value = "DV", -"DMSO...1") |>  
    mutate(REP = extract_numeric(`DMSO...1`),                                
           ADDITIVE_CONC = 0,
           CELLINE = 2,
           CELLINE_NAME = "R183W",
           DRUG_NAME = NAME,
           REPEAT = bio_rep,
           DRUG_CONC = ifelse(DRUG_CONC == "DMSO...10",0,as.numeric(DRUG_CONC))) |>  
    select(CELLINE,CELLINE_NAME,DRUG_NAME,ADDITIVE_CONC,REPEAT,REP,DRUG_CONC,DV,-`DMSO...1`) |>  arrange(REP,DRUG_CONC)
  
  ##################   bind data 1-8 together   ##################################
  data <-  bind_rows(data1, data2, data3, data4, data5, data6, data7, data8)
  write.csv(data, paste0(output_name,".csv"), quote = F, row.names = F)
}

## generate the extracted datasets
path  <- getwd()
file_info <- data.frame(filename = list.files(getwd(), pattern = "\\.xlsx$"))

purrr::pmap(list(file_info$filename), extract_data)

## combine all datasets

# List files with .csv extension in the specified path
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# Initialize an empty data frame
combined_data <- data.frame()

# Loop through each file and bind rows
for (i in files) {
  data <- read.csv(file = i)
  combined_data <- bind_rows(combined_data, data)
}

## generate final dataset for modeling ----------------------------

onco_combine_therapy <- combined_data |>
  filter(!is.na(DV)) |>  
  mutate(ADDITIVE_NAME = "ATUX", 
         ADDITIVE = 1,
         DV = DV*0.01,
         CELLINE_DRUG_ADDITIVE_NAME = paste0(CELLINE_NAME,"_",DRUG_NAME,"_",ADDITIVE_NAME)) |>  
  group_by(DRUG_NAME) |>  mutate(DRUG = cur_group_id()) |>  ungroup() |>  
  arrange(CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, REPEAT,REP)

## Define ID in different ways
onco_combine_therapy <- onco_combine_therapy |>
  group_by(CELLINE, DRUG, ADDITIVE) |>
  mutate(CELLINE_DRUG_ADDITIVE = cur_group_id()) |>
  ungroup() # ID1

onco_combine_therapy <- onco_combine_therapy |>  
  group_by(CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC) |>
  mutate(CELLINE_DRUG_ADDITIVE_ADDITIVECONC = cur_group_id()) |>
  ungroup() # ID2

onco_combine_therapy <- onco_combine_therapy |>
  group_by(CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC,REPEAT) |>
  mutate(CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT = cur_group_id()) |>
  ungroup() # ID3

## rearrange the data and export the dataset
onco_combine_therapy <- onco_combine_therapy |>  
  select(CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, REPEAT, REP, DRUG_CONC, DV, 
         CELLINE_DRUG_ADDITIVE, CELLINE_DRUG_ADDITIVE_ADDITIVECONC, 
         CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT,CELLINE_NAME,
         DRUG_NAME,ADDITIVE_NAME,
         CELLINE_DRUG_ADDITIVE_NAME) 
write.csv(onco_combine_therapy, "onco_combine_therapy.csv", quote = F, row.names = F)

# return to the original directory
Path <- getwd()
setwd(dirname(Path))

# Some EDAs ----

## unique  conc per drug ----

Dactolisib <- onco_combine_therapy |>
  dplyr::filter(DRUG == 1) |>
  dplyr::select(DRUG_CONC) |>
  unique() |>
  dplyr::pull(DRUG_CONC) |>
  round(5) |>
  paste(collapse = ", ")
Dactolisib

Erlotinib <- onco_combine_therapy |>
  dplyr::filter(DRUG == 2) |>
  dplyr::select(DRUG_CONC) |>
  unique() |>
  dplyr::pull(DRUG_CONC) |>
  round(5) |>
  paste(collapse = ", ")
Erlotinib

## TV bottom vs ADDITIVE CONC ----

TVbottom_erlotinib <- onco_combine_therapy_out |>
dplyr::filter(DRUG==2) |>
ggplot(aes(x = ADDITIVECONC, y = TVBOTTOM)) +
geom_line() +
facet_wrap(~CELLINE)

ggsave("./plots/TVbottom_erlotinib.png", width = 10, height = 10)
# this is in line with what observed in dDSS0

# Data visualization ----

## Cell viability plot ----

onco_combine_therapy <- read.csv(file = './Datasets/onco_combine_therapy.csv')

data_plot <- onco_combine_therapy |>  
  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC),
                           'ATUX concentration (µM)' = as.factor(ADDITIVE_CONC)) |>  
  mutate(CELLINE_DRUG_ADDITIVE_NAME_ADDITIVE_CONC_DRUG_CONC = paste0(CELLINE_DRUG_ADDITIVE_NAME,"_", ADDITIVE_CONC,"_", DRUG_CONC)) |> 
  group_by(CELLINE_DRUG_ADDITIVE_NAME, ADDITIVE_CONC, DRUG_CONC) |>
  mutate(meanDV = mean(DV), medianDV = median(DV)) |>
  ungroup()

# set seed for reproducibility
set.seed(321)
palette_cell_viability <- unname(distinctColorPalette(4)) # color palette for cell viability plot
# "#B76DD5" "#B6D876" "#9FD2D3" "#D89BA5"

cell_viability_dactolisib <- data_plot |> 
  dplyr::filter(DRUG_NAME == "Dactolisib") |>
  ggplot(aes(x = DRUG_CONC, colour = `ATUX concentration (µM)`)) +
  geom_point(aes(y = DV), alpha = 0.5) +
  geom_line(aes(y = meanDV), size = 1) +
  facet_grid(CELLINE_NAME ~ DRUG_NAME, scales = "free") +
  scale_x_continuous(breaks = seq(0, 1.7, by = 0.1), limits = c(0, 1.7)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.2), limits = c(0, 2)) +
  scale_color_manual(values = palette_cell_viability) +
  xlab("Drug concentration (µM)") +
  ylab("Cell viability") +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_text(size = 10, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        strip.text = element_text(size = 10, family = "sans"))
cell_viability_dactolisib

cell_viability_erlotinib <- data_plot |> 
  dplyr::filter(DRUG_NAME == "Erlotinib") |>
  ggplot(aes(x = DRUG_CONC, colour = `ATUX concentration (µM)`)) +
  geom_point(aes(y = DV), alpha = 0.5) +
  geom_line(aes(y = meanDV), size = 1) +
  facet_grid(CELLINE_NAME ~ DRUG_NAME, scales = "free") +
  scale_x_continuous(breaks = seq(0, 60, by = 5), limits = c(0, 60)) +
  scale_y_continuous(breaks = seq(0, 2.4, by = 0.2), limits = c(0, 2.4)) +
  scale_color_manual(values = palette_cell_viability) +
  xlab("Drug concentration (µM)") +
  ylab("Cell viability") +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_text(size = 10, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        strip.text = element_text(size = 10, family = "sans"))
cell_viability_erlotinib

# export the plot
setwd("./plots")
ggsave('cell_viability_dactolisib.svg', cell_viability_dactolisib , width = 19, height = 19, unit = "cm", dpi = 300)
ggsave('cell_viability_erlotinib.svg', cell_viability_erlotinib , width = 19, height = 19, unit = "cm", dpi = 300)
# return to the original directory
Path <- getwd()
setwd(dirname(Path))


## Goodness-of-fit ----

onco_combine_therapy_out <- read.table(file = './NONMEM_modelling/onco_combine_therapy_out_003.csv', skip = 1, header = T)
onco_combine_therapy_out$CELLINE <- ordered(onco_combine_therapy_out$CELLINE, levels = c(1, 2),
        labels = c("WT", "R183W"))

#setwd("./plots/")

onco_IPRED_DV <- ggplot(data = onco_combine_therapy_out, aes(x = (IPRED), y = (DV))) + 
  #choose:
  facet_grid(~CELLINE) +
  geom_point(colour = "dodgerblue", shape = 1, size = 2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1.7), breaks = seq(0, 1.7, by = 0.4), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by = 0.5), expand = c(0,0.01)) +
  geom_smooth(se = F, aes(colour = "red"), method = "loess", size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18),
        panel.spacing = unit(1, "lines")) +
  labs(x = "Predicted cell viability", y = "Observed cell viability")
onco_IPRED_DV
ggsave('./plots/onco_IPRED_DV.png', onco_IPRED_DV , device = "png", width = 16, height = 4, dpi = 300)

onco_IPRED_CWRES <- ggplot(data = onco_combine_therapy_out, aes(x = IPRED, y = CWRES)) + 
  #choose:
  facet_grid(~CELLINE) +
  geom_point(colour = "dodgerblue", shape = 1, size = 2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(limits = c(0, 1.7), breaks = seq(0, 1.7, by = 0.2), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-100, 999, by = 1), expand = c(0,0.01)) +
  geom_smooth(se = F, aes(colour = "red"), method = "loess", size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18),
        panel.spacing = unit(1, "lines")) +
  geom_abline(intercept = -2, slope = 0, linetype = "dashed") +
  geom_abline(intercept =  2, slope = 0, linetype = "dashed") +
  labs(x = "Predicted cell viability", y = "CWRES")
onco_IPRED_CWRES
ggsave('./plots/onco_IPRED_CWRES.png', onco_IPRED_CWRES, device = "png", width = 16, height = 4, dpi = 300)


## Fitted cell viability ----

## Create more values for drug concentration

# Calculate the number of observations, min and max DRUGCONC for each ID
id_summary <- onco_combine_therapy_out |>
  group_by(ID) |>
  summarise(
    n = n(),
    min_drugconc = min(DRUGCONC),
    max_drugconc = max(DRUGCONC)
  )

# Generate the DRUG_CONC values for each ID and add them to the dataset

ids_group1 <- c(1:4, 9:12)
ids_group2 <- setdiff(unique(onco_combine_therapy_out$ID), ids_group1)


onco_combine_therapy_out <- onco_combine_therapy_out |>
  mutate(
    DRUG_CONC = case_when(
      ID %in% ids_group1 ~ rep(seq(0, 1.67, length.out = 108), length.out = n()),
      ID %in% ids_group2 ~ rep(seq(0, 58.3, length.out = 81), length.out = n())
    )
  )

# add predicted_cell_viability to the dataset
onco_combine_therapy_out <- onco_combine_therapy_out |>
  mutate(predicted_cell_viability = TVBOTTOM + (TVTOP - TVBOTTOM) / (1 + (DRUG_CONC / TVIC50)^TVGAMMA) + 0.00001)

# add DRUG_NAME to onco_combine_therapy_out, equal "Dacotlisib' if DRUG == 1, 
# and "Erlotinib" otherwise
onco_combine_therapy_out <- onco_combine_therapy_out |>
  mutate(DRUG_NAME = ifelse(DRUG == 1, "Dactolisib", "Erlotinib"))

# Plot the fitted cell viability

# "palette_cell_viability" taken from above

predicted_cell_viability_dactolisib <- ggplot(data = onco_combine_therapy_out |>
         filter(DRUG == 1)) +
  geom_point(aes(x = DRUGCONC, y = DV, colour = as.factor(ADDITIVECONC)), alpha = 0.5) +
  geom_line(aes(x = DRUG_CONC, y = predicted_cell_viability, colour = as.factor(ADDITIVECONC)), size = 1) +
  geom_segment(aes(x = TVIC50, y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), linetype = "dashed") +
  geom_segment(aes(x = 0, y = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, xend = TVIC50, yend = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), linetype = "dashed") +
  geom_point(aes(x = TVIC50, y = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), shape = 8, stroke = 1, size = 1) +
  facet_grid(CELLINE ~ DRUG_NAME) +
  scale_x_continuous(breaks = seq(0, 1.7, by = 0.1), limits = c(0, 1.7)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.2), limits = c(0, 2)) +
  labs(x = expression(paste("Dactolisib concentration (µM)")), col = expression(paste("ATUX concentration (µM)")), y = "Cell viability") +
  scale_color_manual(values = palette_cell_viability) +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_text(size = 10, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        strip.text = element_text(size = 10, family = "sans"))
predicted_cell_viability_dactolisib

predicted_cell_viability_erlotinib <- ggplot(data = onco_combine_therapy_out |>
                                               filter(DRUG == 2)) +
  geom_point(aes(x = DRUGCONC, y = DV, colour = as.factor(ADDITIVECONC)), alpha = 0.5) +
  geom_line(aes(x = DRUG_CONC, y = predicted_cell_viability, colour = as.factor(ADDITIVECONC)), size = 1) +
  geom_segment(aes(x = TVIC50, y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), linetype = "dashed") +
  geom_segment(aes(x = 0, y = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, xend = TVIC50, yend = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), linetype = "dashed") +
  geom_point(aes(x = TVIC50, y = TVBOTTOM + (TVTOP - TVBOTTOM) / 2, colour = as.factor(ADDITIVECONC)), shape = 8, stroke = 1, size = 1) +
  facet_grid(CELLINE ~ DRUG_NAME) +
  scale_x_continuous(breaks = seq(0, 60, by = 5), limits = c(0, 60), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2.4, by = 0.2), limits = c(0, 2.4), expand = c(0, 0)) +
  labs(x = expression(paste("Erlotinib concentration (µM)")), col = expression(paste("ATUX concentration (µM)")), y = "Cell viability") +
  scale_color_manual(values = palette_cell_viability) +
  theme_bw() +
  theme(legend.position = "top", 
        axis.title = element_text(size = 10, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        strip.text = element_text(size = 10, family = "sans"))
predicted_cell_viability_erlotinib

# export the plot
setwd("./plots")
ggsave('predicted_cell_viability_dactolisib.svg', predicted_cell_viability_dactolisib , width = 19, height = 19, unit = "cm", dpi = 300)
ggsave('predicted_cell_viability_erlotinib.svg', predicted_cell_viability_erlotinib , width = 19, height = 19, unit = "cm", dpi = 300)
# return to the original directory
Path <- getwd()
setwd(dirname(Path))


## AUC - drug sensitivity score (DSS) ---------------------------

# Function to create a model function given a dataset
create_model <- function(dat) {
  function(x) {
    dat$TVBOTTOM + ((dat$TVTOP - dat$TVBOTTOM) / (1 + ((x / dat$TVIC50)^dat$TVGAMMA)))
  }
}

# Initialize an empty list to store models
models <- list()

# Loop through IDs 1 to 16, subset the data, and create models
for (id in 1:16) {
  dat <- onco_combine_therapy_out |> filter(ID == id) |> distinct(ID, TVBOTTOM, TVTOP, TVIC50, TVGAMMA)
  models[[paste0("model_id", id)]] <- create_model(dat)
}

# upper concentration varies by drug
onco_upper_drugconc <- onco_combine_therapy_out |>  
  group_by(ID) |>  
  slice(n()) |>  
  ungroup() |>   
  select(ID, DRUGCONC)

# calculate effect AUC from 0 to upper drug concentration
onco_DSS <- list()  # Initialize onco_DSS as a list to store integration results
for (i in 1:16) {
  result <- integrate(models[[i]], lower = 0, upper = onco_upper_drugconc$DRUGCONC[i])  # onco_DSS = integrand
  onco_DSS[[i]] <- result$value  # Store the value of the integration result in the onco_DSS list
  print(onco_DSS[[i]])
}
onco_DSS <- as.data.frame(onco_DSS)
onco_DSS <- onco_DSS |> 
  t() |>  
  as.data.frame() |>
  rownames_to_column() |>
  select(MODEL = rowname, AUC = V1)

# Since BOTTOM was not always reached at upper drug concentration --> 
# Calculate the area from 0 to upper drug concentration below Y = function(X = onco_upper_drugconc$DRUGCONC)
onco_DSS_bttm <- data.frame(
  models[[1]](x = onco_upper_drugconc$DRUGCONC[1]),
  models[[2]](x = onco_upper_drugconc$DRUGCONC[2]),
  models[[3]](x = onco_upper_drugconc$DRUGCONC[3]),
  models[[4]](x = onco_upper_drugconc$DRUGCONC[4]),
  models[[5]](x = onco_upper_drugconc$DRUGCONC[5]),
  models[[6]](x = onco_upper_drugconc$DRUGCONC[6]),
  models[[7]](x = onco_upper_drugconc$DRUGCONC[7]),
  models[[8]](x = onco_upper_drugconc$DRUGCONC[8]),
  models[[9]](x = onco_upper_drugconc$DRUGCONC[9]),
  models[[10]](x = onco_upper_drugconc$DRUGCONC[10]),
  models[[11]](x = onco_upper_drugconc$DRUGCONC[11]),
  models[[12]](x = onco_upper_drugconc$DRUGCONC[12]),
  models[[13]](x = onco_upper_drugconc$DRUGCONC[13]),
  models[[14]](x = onco_upper_drugconc$DRUGCONC[14]),
  models[[15]](x = onco_upper_drugconc$DRUGCONC[15]),
  models[[16]](x = onco_upper_drugconc$DRUGCONC[16])
  )
  
onco_DSS_bttm <- onco_DSS_bttm |> 
  t() |>
  as.data.frame() |>
  rownames_to_column()

onco_DSS_bttm <- onco_DSS_bttm |> 
  mutate(rowname = sub(".x...*.", "", rowname),
         upper_DRUGCONC = onco_upper_drugconc$DRUGCONC) |>  
  mutate(BTTM_AUC = V1 * upper_DRUGCONC) |>
  select(MODEL = rowname, BTTM = V1, BTTM_AUC,upper_DRUGCONC)

onco_DSS <- onco_DSS |>
  mutate(ID = row_number())
onco_DSS_bttm <- onco_DSS_bttm |>
  mutate(ID = row_number())

onco_combine_therapy_out <- full_join(onco_combine_therapy_out, onco_DSS
                                      |>  select(ID, AUC), by = c("ID"))

onco_combine_therapy_out <- full_join(onco_combine_therapy_out, onco_DSS_bttm |>
                                        select(ID, BTTM, BTTM_AUC, 
                                               upper_DRUGCONC), by = c("ID"))

onco_combine_therapy_out <- onco_combine_therapy_out |>
  mutate(DSS0  = (AUC - (TVBOTTOM * upper_DRUGCONC)), # BOTTOM is the model EGFRameter estimate. BTTM = the value of Y where X=45. BTTM_AUC = BTTM*upper_DRUGCONC = the AUC below BTTM from 0 to 45.
                      DSS0i = (AUC - (BTTM_AUC)), # Bottom was not always reached at Sora/Regora = 45 --> therefore this.
                      DSS1  = (AUC - (TVBOTTOM * upper_DRUGCONC)) / ((TVTOP - TVBOTTOM) * upper_DRUGCONC),
                      DSS1i = (AUC - (BTTM_AUC   )) / ((TVTOP - BTTM ) * upper_DRUGCONC),
                      DSS2  = DSS1/log10(TVTOP),
                      DSS2i = DSS1i / log10(TVTOP),
                      DSS3  = DSS2   * ((upper_DRUGCONC)/(upper_DRUGCONC)),
                      DSS3i = DSS2i  * ((upper_DRUGCONC)/(upper_DRUGCONC))) 

onco_DSS_full <- full_join(onco_combine_therapy |>
                      select(ID = CELLINE_DRUG_ADDITIVE_ADDITIVECONC, CELLINE_DRUG_ADDITIVE,CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, CELLINE_NAME, DRUG_NAME, ADDITIVE_NAME, CELLINE_DRUG_ADDITIVE_NAME) |>  unique(),
                 onco_combine_therapy_out |>  select(ID, TVTOP, TVBOTTOM, TVIC50, TVGAMMA, AUC,upper_DRUGCONC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i),
                 by = c("ID"))


### Difference in DSS as comEGFRed to no ATUX (ATUX = 0) ----

onco_DSS_reference_0 <- onco_DSS_full |>  
  filter(ADDITIVE_CONC == 0) |>  
  distinct(CELLINE_DRUG_ADDITIVE, AUC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i) |>  
  select(CELLINE_DRUG_ADDITIVE, 
         ref_AUC = AUC,
         ref_BTTM = BTTM,
         ref_BTTM_AUC = BTTM_AUC,
         ref_DSS0 = DSS0,
         ref_DSS1 = DSS1,
         ref_DSS2 = DSS2,
         ref_DSS3 = DSS3,
         ref_DSS0i = DSS0i,
         ref_DSS1i = DSS1i,
         ref_DSS2i = DSS2i,
         ref_DSS3i = DSS3i)

onco_DSS <- left_join(onco_DSS_full, onco_DSS_reference_0, by = "CELLINE_DRUG_ADDITIVE")
onco_DSS <- onco_DSS |>  
  mutate(dDSS0   = DSS0   - ref_DSS0,
         dDSS1   = DSS1   - ref_DSS1,
         dDSS2   = DSS2   - ref_DSS2,
         dDSS3   = DSS3   - ref_DSS3,
         dDSS0i  = DSS0i  - ref_DSS0i,
         dDSS1i  = DSS1i  - ref_DSS1i,
         dDSS2i  = DSS2i  - ref_DSS2i,
         dDSS3i  = DSS3i  - ref_DSS3i) 

### dDSS dsitribution ----

group.colors <- c(Dactolisib  = "#d39200", 
                  Erlotinib = "#93aa00")
ggplot(data = onco_DSS, mapping = aes(x = ADDITIVE_CONC, y = dDSS0,  colour = as.factor(DRUG_NAME))) +
  facet_wrap(~ CELLINE_NAME, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey") +
  geom_point(size = 3) +
  geom_line(size = 0.75, linetype = "dashed") +
  labs(x = expression(paste("ATUX concentration (µM)")), y = "dDSS0", col = "Compound") +
  scale_x_continuous(breaks = c(0, 3, 4.5, 6)) +
  #scale_y_continuous(breaks = seq(-100, 100, 1), limits = c(-6, 0.5)) +
  scale_colour_manual(values = group.colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0), 
        panel.grid.minor = element_blank(),
        legend.position = "top", 
        axis.title = element_text(size = 10, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        strip.text = element_text(size = 10, family = "sans"))
ggsave('./plots/dDSS0.svg',
       width = 19, height = 19, unit = "cm", dpi = 300)

### AUC trapezoidal rule ----

# calculate AUC based on trapezoidal rule 

onco_AUC_trapezoid <- cbind(onco_DSS, onco_combine_therapy |> 
                              select(DRUG_CONC,DV,CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT,REP))

onco_AUC_trapezoid <- onco_AUC_trapezoid |> 
  group_by(CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT,REP)  |>  
  mutate(ID_trapez = cur_group_id()) |>  # ID_trapez based on CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT_REP
  ungroup() |>
  mutate(AUC_trapez_i = NA)


IDs = unique(onco_AUC_trapezoid$ID_trapez)
for (i in IDs) {
  data <- onco_AUC_trapezoid |> 
    filter(ID_trapez == i)|> 
    select(ID_trapez, DRUG_CONC,DV)
  AUC_value <- trapz(data$DRUG_CONC, data$DV)
  onco_AUC_trapezoid$AUC_trapez_i[onco_AUC_trapezoid$ID_trapez == i] = AUC_value
}

onco_AUC_trapezoid <- onco_AUC_trapezoid |> 
  group_by(ID) |> 
  mutate(AUC_trapez = median(AUC_trapez_i)) |>
  ungroup() 

group.colors <- c(WT = "#d39200", 
                  R183W = "#93aa00")

ggplot(data = onco_AUC_trapezoid , aes(x = AUC_trapez, y = AUC, colour = CELLINE_NAME)) + 
  geom_point( size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5), expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, by = 5), expand = c(0,0.01)) +
  geom_smooth(se = F,colour = "red", method = "lm", size = 1.2, alpha=0.1) +
  scale_colour_manual(values=group.colors)+
  theme_bw()+
  theme(text = element_text(size = 12),
        legend.position = "top") +
  labs(x = "Predicted AUC", y = "Observed AUC", col = "Cell line") 
ggsave('./plots/Observed AUC vs Predicted AUC.png', width = 7.5, height = 6, dpi = 300)


## Simulations ----

onco_combine_therapy_sim_001 <- read.csv(file = './Datasets/onco_combine_therapy_sim_001.csv', skip = 1, header = TRUE, sep = "")

onco_combine_therapy_sim_001 <- onco_combine_therapy_sim_001 |>   
  filter(ID != "TABLE" & ID != "ID" )          # remove two rows for each simulation

# sim dataset for dactolisib
onco_combine_therapy_sim_001_dactolisib <- onco_combine_therapy_sim_001 |>
  mutate_all(function(x) as.numeric(as.character(x))) |> # convert dataframe to numeric
  filter(ID %in% c(1:4, 9:12)) |>
  filter(DRUGCONC == 0 | DRUGCONC == 1.666700) |> 
  mutate(SIM = ceiling(row_number()/192)) # label each time of simulation

# sim dataset for erlotinib
onco_combine_therapy_sim_001_erlotinib <- onco_combine_therapy_sim_001 |>
  mutate_all(function(x) as.numeric(as.character(x))) |> # convert dataframe to numeric
  filter(ID %in% c(5:8, 13:16)) |>
  filter(DRUGCONC == 0 | DRUGCONC == 58.33300) |> 
  mutate(SIM = ceiling(row_number()/144)) # label each time of simulation

## first, for dactolisib analysis

onco_combine_therapy_sim_001_dactolisib <- onco_combine_therapy_sim_001_dactolisib |>  
  arrange(ID, DAY, REP, SIM) |> 
  group_by(ID, DAY, REP)     |>  
  mutate(ID_DAY_REP = cur_group_id()) |>  
  ungroup() |>  
  group_by(ID, DAY, REP, SIM) |>  
  mutate(ID_DAY_REP_SIM = cur_group_id()) |>  
  ungroup() |>  
  select(ID_DAY_REP, ID_DAY_REP_SIM, SIM, everything())

onco_combine_therapy_sim_001_dactolisib$CELLINE <- ordered(
  onco_combine_therapy_sim_001_dactolisib$CELLINE, 
  levels = c(1, 2),
  labels = c("WT", "R183W")
)

upper_DRUGCONC_dactolisib <- onco_combine_therapy_sim_001_dactolisib |>  
  group_by(ID_DAY_REP_SIM) |>  
  slice(n()) |>  
  ungroup() |>  
  select(ID_DAY_REP_SIM, DRUGCONC) 

compute_DSS <- function(dat) {
  equation <- function(x) { dat$BOTTOM + ((dat$TOP - dat$BOTTOM) / (1 + ((x / dat$IC50) ^ dat$GAMMA))) }
  DSS <- integrate(equation, lower = 0, upper = 1.6667)$value
  DSS_bttm <- equation(1) * 1
  list(DSS = DSS, DSS_bttm = DSS_bttm)
}

DSS_results <- onco_combine_therapy_sim_001_dactolisib |> 
  distinct(ID_DAY_REP_SIM, CELLINE, ADDITIVECONC, TOP, BOTTOM, IC50, GAMMA) |> 
  nest(data = c(CELLINE, ADDITIVECONC, TOP, BOTTOM, IC50, GAMMA)) |> 
  mutate(results = map(data, compute_DSS)) |> 
  unnest_wider(results) |> 
  select(ID_DAY_REP_SIM, DSS, DSS_bttm)

DSS_dactolisib <- DSS_results$DSS
DSS_bttm_dactolisib <- DSS_results$DSS_bttm

DSS_dactolisib <- as.data.frame(DSS_dactolisib) |>  
  rownames_to_column() |>  
  select(MODEL = rowname, AUC = DSS_dactolisib) |>  
  mutate(ID_DAY_REP_SIM = row_number())

DSS_bttm_dactolisib <- as.data.frame(DSS_bttm_dactolisib) |>  
  mutate(ID_DAY_REP_SIM = row_number(), 
         BTTM = DSS_bttm_dactolisib, 
         upper_DRUGCONC = 1.6667) |>  
  select(ID_DAY_REP_SIM, 
         BTTM_AUC = DSS_bttm_dactolisib, 
         BTTM, 
         upper_DRUGCONC)

onco_combine_therapy_sim_001_dactolisib <- full_join(onco_combine_therapy_sim_001_dactolisib, DSS_dactolisib |>  
                                                       select(ID_DAY_REP_SIM, AUC), 
                                                     by = c("ID_DAY_REP_SIM"))

onco_combine_therapy_sim_001_dactolisib <- full_join(onco_combine_therapy_sim_001_dactolisib, 
                                                     DSS_bttm_dactolisib,  
                                                     by = c("ID_DAY_REP_SIM"))

onco_combine_therapy_sim_001_dactolisib <- onco_combine_therapy_sim_001_dactolisib |>  
  distinct(ID_DAY_REP_SIM, .keep_all = T) |>  
  group_by(ID, SIM) |>  
  mutate(AUC_median = median(AUC), 
         BTTM_AUC_median = median(BTTM_AUC), 
         BOTTOM_median = median(BOTTOM), 
         TOP_median = median(TOP), 
         BTTM_median = median(BTTM))

onco_combine_therapy_sim_001_dactolisib <- onco_combine_therapy_sim_001_dactolisib |>  
  mutate(DSS0 = (AUC_median - (BOTTOM_median * upper_DRUGCONC)), 
         DSS0i = (AUC_median - (BTTM_AUC_median)), 
         DSS1 = (AUC_median - (BOTTOM_median * upper_DRUGCONC)) / ((TOP_median - BOTTOM_median) * upper_DRUGCONC), 
         DSS1i = (AUC_median - (BTTM_AUC_median)) / ((TOP_median - BTTM_median) * upper_DRUGCONC), 
         DSS2 = DSS1 / log10(TOP_median), 
         DSS2i = DSS1i / log10(TOP_median), 
         DSS3 = DSS2 * ((upper_DRUGCONC) / (upper_DRUGCONC)), 
         DSS3i = DSS2i * ((upper_DRUGCONC) / (upper_DRUGCONC))) |>  
  distinct(ID, SIM, ADDITIVECONC, .keep_all = T)

# Calculate the difference in DSS as compared to no OA (OA = 0)
DSS_reference_0_dactolisib <- onco_combine_therapy_sim_001_dactolisib |>  
  filter(ADDITIVECONC == 0) |>  
  ungroup() |>  
  select(SIM, CELLINE, ref_AUC = AUC_median, ref_BTTM = BTTM_median, 
         ref_BTTM_AUC = BTTM_AUC_median, ref_DSS0 = DSS0, ref_DSS1 = DSS1, 
         ref_DSS2 = DSS2, ref_DSS3 = DSS3, ref_DSS0i = DSS0i, ref_DSS1i = DSS1i, 
         ref_DSS2i = DSS2i, ref_DSS3i = DSS3i)

DSS_dactolisib <- left_join(onco_combine_therapy_sim_001_dactolisib, 
                            DSS_reference_0_dactolisib, 
                            by = c("CELLINE", "SIM"))
DSS_dactolisib <- DSS_dactolisib |>  
  mutate(dDSS0 = DSS0 - ref_DSS0, dDSS1 = DSS1 - ref_DSS1, 
         dDSS2 = DSS2 - ref_DSS2, dDSS3 = DSS3 - ref_DSS3, 
         dDSS0i = DSS0i - ref_DSS0i, dDSS1i = DSS1i - ref_DSS1i, 
         dDSS2i = DSS2i - ref_DSS2i, dDSS3i = DSS3i - ref_DSS3i)

# Plot the dDSS distribution
DSS_ordered_dactolisib <- DSS_dactolisib |>  
  filter(ADDITIVECONC != 0) 
DSS_ordered_dactolisib$CELLINE <- factor(DSS_ordered_dactolisib$CELLINE, 
                                         levels = c("WT", "R183W"))
mycomparison <- list(c("WT", "R183W"))

p1_dactolisib <- ggboxplot(DSS_ordered_dactolisib, x = "CELLINE", y = "dDSS0") +
  stat_compare_means(comparisons = mycomparison, 
                     method = "wilcox.test", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  labs(x = "Cell line", y = "dDSS0 (ATUX = 6 µM)") +
  facet_wrap(~c("Dactolisib")) +
  #scale_y_continuous(breaks = seq(-0.3, 0, 0.1), limits = c(-0.3, 0)) +
  theme_bw()
p1_dactolisib

ggsave('./plots/dDSS0_dactolisib_1000_simulation.png', 
       p1_dactolisib, 
       width = 6, 
       height = 5.5, 
       dpi = 300)

## second, for erlotinib analysis

# Process the dataset
onco_combine_therapy_sim_001_erlotinib <- onco_combine_therapy_sim_001_erlotinib |>  
  arrange(ID, DAY, REP, SIM) |> 
  group_by(ID, DAY, REP)     |>  
  mutate(ID_DAY_REP = cur_group_id()) |>  
  ungroup() |>  
  group_by(ID, DAY, REP, SIM) |>  
  mutate(ID_DAY_REP_SIM = cur_group_id()) |>  
  ungroup() |>  
  select(ID_DAY_REP, ID_DAY_REP_SIM, SIM, everything())

onco_combine_therapy_sim_001_erlotinib$CELLINE <- ordered(
  onco_combine_therapy_sim_001_erlotinib$CELLINE, 
  levels = c(1, 2),
  labels = c("WT", "R183W")
)

upper_DRUGCONC_erlotinib <- onco_combine_therapy_sim_001_erlotinib |>  
  group_by(ID_DAY_REP_SIM) |>  
  slice(n()) |>  
  ungroup() |>  
  select(ID_DAY_REP_SIM, DRUGCONC) 

compute_DSS_erlotinib <- function(dat) {
  equation <- function(x) { dat$BOTTOM + ((dat$TOP - dat$BOTTOM) / (1 + ((x / dat$IC50) ^ dat$GAMMA))) }
  DSS <- integrate(equation, lower = 0, upper = 58.33300)$value
  DSS_bttm <- equation(1) * 1
  list(DSS = DSS, DSS_bttm = DSS_bttm)
}

DSS_results_erlotinib <- onco_combine_therapy_sim_001_erlotinib |> 
  distinct(ID_DAY_REP_SIM, CELLINE, ADDITIVECONC, TOP, BOTTOM, IC50, GAMMA) |> 
  nest(data = c(CELLINE, ADDITIVECONC, TOP, BOTTOM, IC50, GAMMA)) |> 
  mutate(results = map(data, compute_DSS_erlotinib)) |> 
  unnest_wider(results) |> 
  select(ID_DAY_REP_SIM, DSS, DSS_bttm)

DSS_erlotinib <- DSS_results_erlotinib$DSS
DSS_bttm_erlotinib <- DSS_results_erlotinib$DSS_bttm

DSS_erlotinib <- as.data.frame(DSS_erlotinib) |>  
  rownames_to_column() |>  
  select(MODEL = rowname, AUC = DSS_erlotinib) |>  
  mutate(ID_DAY_REP_SIM = row_number())

DSS_bttm_erlotinib <- as.data.frame(DSS_bttm_erlotinib) |>  
  mutate(ID_DAY_REP_SIM = row_number(), 
         BTTM = DSS_bttm_erlotinib, 
         upper_DRUGCONC = 58.33300) |>  
  select(ID_DAY_REP_SIM, 
         BTTM_AUC = DSS_bttm_erlotinib, 
         BTTM, 
         upper_DRUGCONC)

onco_combine_therapy_sim_001_erlotinib <- full_join(onco_combine_therapy_sim_001_erlotinib, DSS_erlotinib |>  
                                                      select(ID_DAY_REP_SIM, AUC), 
                                                    by = c("ID_DAY_REP_SIM"))

onco_combine_therapy_sim_001_erlotinib <- full_join(onco_combine_therapy_sim_001_erlotinib, 
                                                    DSS_bttm_erlotinib,  
                                                    by = c("ID_DAY_REP_SIM"))

onco_combine_therapy_sim_001_erlotinib <- onco_combine_therapy_sim_001_erlotinib |>  
  distinct(ID_DAY_REP_SIM, .keep_all = T) |>  
  group_by(ID, SIM) |>  
  mutate(AUC_median = median(AUC), 
         BTTM_AUC_median = median(BTTM_AUC), 
         BOTTOM_median = median(BOTTOM), 
         TOP_median = median(TOP), 
         BTTM_median = median(BTTM))

onco_combine_therapy_sim_001_erlotinib <- onco_combine_therapy_sim_001_erlotinib |>  
  mutate(DSS0 = (AUC_median - (BOTTOM_median * upper_DRUGCONC)), 
         DSS0i = (AUC_median - (BTTM_AUC_median)), 
         DSS1 = (AUC_median - (BOTTOM_median * upper_DRUGCONC)) / ((TOP_median - BOTTOM_median) * upper_DRUGCONC), 
         DSS1i = (AUC_median - (BTTM_AUC_median)) / ((TOP_median - BTTM_median) * upper_DRUGCONC), 
         DSS2 = DSS1 / log10(TOP_median), 
         DSS2i = DSS1i / log10(TOP_median), 
         DSS3 = DSS2 * ((upper_DRUGCONC) / (upper_DRUGCONC)), 
         DSS3i = DSS2i * ((upper_DRUGCONC) / (upper_DRUGCONC))) |>  
  distinct(ID, SIM, ADDITIVECONC, .keep_all = T)

# Calculate the difference in DSS as compared to no OA (OA = 0)
DSS_reference_0_erlotinib <- onco_combine_therapy_sim_001_erlotinib |>  
  filter(ADDITIVECONC == 0) |>  
  ungroup() |>  
  select(SIM, CELLINE, ref_AUC = AUC_median, ref_BTTM = BTTM_median, 
         ref_BTTM_AUC = BTTM_AUC_median, ref_DSS0 = DSS0, ref_DSS1 = DSS1, 
         ref_DSS2 = DSS2, ref_DSS3 = DSS3, ref_DSS0i = DSS0i, ref_DSS1i = DSS1i, 
         ref_DSS2i = DSS2i, ref_DSS3i = DSS3i)

DSS_erlotinib <- left_join(onco_combine_therapy_sim_001_erlotinib, 
                           DSS_reference_0_erlotinib, 
                           by = c("CELLINE", "SIM"))
DSS_erlotinib <- DSS_erlotinib |>  
  mutate(dDSS0 = DSS0 - ref_DSS0, dDSS1 = DSS1 - ref_DSS1, 
         dDSS2 = DSS2 - ref_DSS2, dDSS3 = DSS3 - ref_DSS3, 
         dDSS0i = DSS0i - ref_DSS0i, dDSS1i = DSS1i - ref_DSS1i, 
         dDSS2i = DSS2i - ref_DSS2i, dDSS3i = DSS3i - ref_DSS3i)

# Plot the dDSS distribution
DSS_ordered_erlotinib <- DSS_erlotinib |>  
  filter(ADDITIVECONC != 0) 
DSS_ordered_erlotinib$CELLINE <- factor(DSS_ordered_erlotinib$CELLINE, 
                                        levels = c("WT", "R183W"))
mycomparison <- list(c("WT", "R183W"))

p1_erlotinib <- ggboxplot(DSS_ordered_erlotinib, x = "CELLINE", y = "dDSS0") +
  stat_compare_means(comparisons = mycomparison, 
                     method = "wilcox.test", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  facet_wrap(~c("Erlotinib")) +
  labs(x = "Cell line", y = "dDSS0 (ATUX = 6 µM)") +
  theme_bw()
p1_erlotinib

ggsave('./plots/dDSS0_erlotinib_1000_simulation.png', 
       p1_erlotinib, 
       width = 6, 
       height = 5.5, 
       dpi = 300)

