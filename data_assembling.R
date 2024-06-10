# Load packages ----
library(readxl)
library(tidyverse)
library(reshape)
library(stringr)
library(ggpubr)


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
###################generate final dataset for modeling   #######################
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

#####################  END OF DATA MANAGEMENT  #################################
################################################################################

##########################  DATA VISUALISATION  #################################
data_plot <- onco_combine_therapy |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC),
                           'OA concentration (nM)'=as.factor(ADDITIVE_CONC)) |>  
  mutate(CELLINE_DRUG_ADDITIVE_NAME_ADDITIVE_CONC_DRUG_CONC=paste0(CELLINE_DRUG_ADDITIVE_NAME,"_", ADDITIVE_CONC,"_", DRUG_CONC)) |> 
  group_by(CELLINE_DRUG_ADDITIVE_NAME, ADDITIVE_CONC, DRUG_CONC) |>  mutate(meanDV = mean(DV),
                                                                      medianDV= median(DV)) |>  ungroup()
p <- ggplot(data_plot,aes(x= DRUG_CONC, colour=`OA concentration (nM)`)) +
  geom_point(aes(y= DV), alpha = 0.2)+
  geom_line(aes(y= meanDV))+
  facet_grid(CELLINE_NAME~., scales = "free_x")+
  # scale_x_log10(breaks = c(0,0.01, 0.025,  0.05, 0.1, 0.2,0.4,1),
  #               labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  )) +
  #scale_y_continuous(breaks = seq(0, 150, by = 0.5), limits = c(0, 1.8)) +
  xlab("Clofarabine concentration (µM)")+
  ylab("Cell viability")+
  theme_bw()+
  theme(legend.position = "top")
p
ggsave('Cell viability plot.png', p , device = "png",width = 20, height = 9, dpi = 500)
################################################################################

############################# Goodness-of-fit ##################################
dat_2023 <- read.table (file='~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023//modeling/run006_out.csv', skip=1, header=T)
dat_2023$CELLINE <- ordered(dat_2023$CELLINE, levels=c(1, 2, 3, 4, 5),
        labels=c("EGFP", "WT", "R183W", "P179R", "S256F"))

setwd("~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023/plots/")

p1 <- ggplot (data = dat_2023, aes(x=   (IPRED), y=   (DV))) + 
  #choose:
  facet_grid(~CELLINE)+
  geom_point(colour = "dodgerblue", shape = 1, size = 2, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 2.15), breaks = seq(0, 999, by = 0.5)) +
  scale_y_continuous(limits = c(0, 2.15), breaks = seq(0, 999, by = 0.5)) +
  geom_smooth(se = F, aes(colour = "red"), method = "loess", size = 1.2) +
  theme_bw()+
  theme(legend.position="none") +
  theme(text = element_text(size = 18)) +
  labs(x = "Predicted cell viability", y = "Observed cell viability")
p1
ggsave('IPRED vs DV.png', p1 , device = "png",width = 16, height = 4, dpi = 500)

p2 <- ggplot (data = dat_2023, aes(x= PRED, y=CWRES)) + # LnDV run_9999vi
  #choose:
  facet_grid(~CELLINE)+
  geom_point(colour = "dodgerblue", shape = 1, size = 2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 0) +
  scale_x_continuous(limits = c(0.001, 1.2), breaks = seq(-100, 999, by = 0.2)) +
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-100, 999, by = 5)) +
  geom_smooth(se = F, aes(colour = "red"), method = "loess", size = 1.2) +
  theme_bw()+
  theme(legend.position="none") +
  theme(text = element_text(size = 18)) +
  geom_abline(intercept = -5, slope = 0, linetype = "dashed") +
  geom_abline(intercept =  5, slope = 0, linetype = "dashed") +
  labs(x = "Predicted cell viability", y = "CWRES")
p2
ggsave('CWRES vs PRED.png', p2 , device = "png",width = 16, height = 4, dpi = 500)
################################################################################

#### define the equations for effect curve ####################################
dat1 <- dat_2023 |>  filter(ID == 1) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id1 <- function(x) { dat1$TVBOTTOM + ( (dat1$TVTOP-dat1$TVBOTTOM) / ( 1+ ((x/dat1$TVIC50)^dat1$TVGAMMA) ) )}
dat2 <- dat_2023 |>  filter(ID == 2) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id2 <- function(x) { dat2$TVBOTTOM + ( (dat2$TVTOP-dat2$TVBOTTOM) / ( 1+ ((x/dat2$TVIC50)^dat2$TVGAMMA) ) )}
dat3 <- dat_2023 |>  filter(ID == 3) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id3 <- function(x) { dat3$TVBOTTOM + ( (dat3$TVTOP-dat3$TVBOTTOM) / ( 1+ ((x/dat3$TVIC50)^dat3$TVGAMMA) ) )}
dat4 <- dat_2023 |>  filter(ID == 4) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id4 <- function(x) { dat4$TVBOTTOM + ( (dat4$TVTOP-dat4$TVBOTTOM) / ( 1+ ((x/dat4$TVIC50)^dat4$TVGAMMA) ) )}
dat5 <- dat_2023 |>  filter(ID == 5) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id5 <- function(x) { dat5$TVBOTTOM + ( (dat5$TVTOP-dat5$TVBOTTOM) / ( 1+ ((x/dat5$TVIC50)^dat5$TVGAMMA) ) )}
dat6 <- dat_2023 |>  filter(ID == 6) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id6 <- function(x) { dat6$TVBOTTOM + ( (dat6$TVTOP-dat6$TVBOTTOM) / ( 1+ ((x/dat6$TVIC50)^dat6$TVGAMMA) ) )}
dat7 <- dat_2023 |>  filter(ID == 7) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id7 <- function(x) { dat7$TVBOTTOM + ( (dat7$TVTOP-dat7$TVBOTTOM) / ( 1+ ((x/dat7$TVIC50)^dat7$TVGAMMA) ) )}
dat8 <- dat_2023 |>  filter(ID == 8) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id8 <- function(x) { dat8$TVBOTTOM + ( (dat8$TVTOP-dat8$TVBOTTOM) / ( 1+ ((x/dat8$TVIC50)^dat8$TVGAMMA) ) )}
dat9 <- dat_2023 |>  filter(ID == 9) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id9 <- function(x) { dat9$TVBOTTOM + ( (dat9$TVTOP-dat9$TVBOTTOM) / ( 1+ ((x/dat9$TVIC50)^dat9$TVGAMMA) ) )}

# ID 10-19
dat10 <- dat_2023 |>  filter(ID == 10) |>  distinct(ID,CELLINE,ADDITIVECONC,TVTOP,TVBOTTOM,TVIC50,TVGAMMA)
model_id10 <- function(x) { dat10$TVBOTTOM + ( (dat10$TVTOP-dat10$TVBOTTOM) / ( 1+ ((x/dat10$TVIC50)^dat10$TVGAMMA) ) )}

TV_parameters <- bind_rows(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10) |>  
gather(Parameter, Value, -ID, -CELLINE,-ADDITIVECONC) |>  mutate(comb=paste0(CELLINE," ,OA=",ADDITIVECONC))

paramater <- ggplot(TV_parameters |>  filter(Parameter!="TVGAMMA"),aes(x=comb,y=Value,color=Parameter))+ 
  geom_point() +
  scale_y_continuous("Parameters",breaks = seq(0,1.5,by=0.2))+
  theme_bw()
ggsave('parameters.png', paramater  , device = "png",width = 15, height = 7, dpi = 500)
################################################################################

#### plot the predicted curve with observed data ###############################
# 1 Clofarabine - OA ############
p1 <- ggplot() +
  geom_point(data = onco_combine_therapy |>  filter(CELLINE_NAME == "EGFP" & DRUG_NAME == "Clofarabine" & ADDITIVE_NAME == "OA") |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC)),
             mapping = aes(x = DRUG_CONC, y = DV, colour = as.factor(ADDITIVE_CONC)), alpha = 0.2) +  # ADDITIVECONC is a covariate
  # scale_x_log10(breaks = c(0, 0.001,0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  ),
  #               labels = scales::number_format(accuracy = 0.001)) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  )) +
  scale_y_continuous(breaks = seq(0, 99, by = 0.5), limits = c(0, 1.8), expand = c(0, 0)) +
  facet_grid(.~CELLINE_NAME)+
  labs(x = expression(paste("Clofarabine concentration (µM)")), col=expression(paste("OA concentration (nM)")), y = "Cell viability") +
  scale_color_manual(values=c("red","blue", "purple", "orange", "green"))+
  stat_function(fun = model_id1, color = "red", size = 1, linetype = "solid") +
  stat_function(fun = model_id2, color = "blue", size = 1, linetype = "solid") +
  geom_segment(data=dat1, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_segment(data=dat1, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_point(data=dat1, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1, colour = "red", size = 2.5)+
  geom_segment(data=dat2, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_segment(data=dat2, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_point(data=dat2, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1, colour = "blue", size = 2.5)+
  theme_bw()
p1

p2 <- ggplot() +
  geom_point(data = onco_combine_therapy |>  filter(CELLINE_NAME == "WT" & DRUG_NAME == "Clofarabine" & ADDITIVE_NAME == "OA") |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC)),
             mapping = aes(x = DRUG_CONC, y = DV, colour = as.factor(ADDITIVE_CONC)), alpha = 0.2) +  # ADDITIVECONC is a covariate
  # scale_x_log10(breaks = c(0, 0.001,0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  ),
  #               labels = scales::number_format(accuracy = 0.001)) +
  scale_x_continuous(breaks = c(0,0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  )) +
  scale_y_continuous(breaks = seq(0, 99, by = 0.5), limits = c(0, 1.8), expand = c(0, 0)) +
  facet_grid(.~CELLINE_NAME)+
  labs(x = expression(paste("Clofarabine concentration (µM)")), col=expression(paste("OA concentration (nM)")), y = "Cell viability") +
  scale_color_manual(values=c("red","blue", "purple", "orange", "green"))+
  stat_function(fun = model_id3, color = "red", size = 1, linetype = "solid") +
  stat_function(fun = model_id4, color = "blue", size = 1, linetype = "solid") +
  geom_segment(data=dat3, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_segment(data=dat3, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_point(data=dat3, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "red", size = 2.5)+
  geom_segment(data=dat4, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_segment(data=dat4, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_point(data=dat4, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "blue", size = 2.5)+
  theme_bw()
p2

p3 <- ggplot() +
  geom_point(data = onco_combine_therapy |>  filter(CELLINE_NAME == "R183W" & DRUG_NAME == "Clofarabine" & ADDITIVE_NAME == "OA") |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC)),
             mapping = aes(x = DRUG_CONC, y = DV, colour = as.factor(ADDITIVE_CONC)), alpha = 0.2) +  # ADDITIVECONC is a covariate
  # scale_x_log10(breaks = c(0, 0.001,0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  ),
  #               labels = scales::number_format(accuracy = 0.001)) +
  scale_x_continuous(breaks = c(0,0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  )) +
  scale_y_continuous(breaks = seq(0, 99, by = 0.5), limits = c(0, 1.8), expand = c(0, 0)) +
  facet_grid(.~CELLINE_NAME)+
  labs(x = expression(paste("Clofarabine concentration (µM)")), col=expression(paste("OA concentration (nM)")), y = "Cell viability") +
  scale_color_manual(values=c("red","blue", "purple", "orange", "green"))+
  stat_function(fun = model_id5, color = "red", size = 1, linetype = "solid") +
  stat_function(fun = model_id6, color = "blue", size = 1, linetype = "solid") +
  geom_segment(data=dat5, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_segment(data=dat5, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_point(data=dat5, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "red", size = 2.5)+
  geom_segment(data=dat6, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_segment(data=dat6, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_point(data=dat6, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "blue", size = 2.5)+
  theme_bw()
p3

p4 <- ggplot() +
  geom_point(data = onco_combine_therapy |>  filter(CELLINE_NAME == "P179R" & DRUG_NAME == "Clofarabine" & ADDITIVE_NAME == "OA") |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC)),
             mapping = aes(x = DRUG_CONC, y = DV, colour = as.factor(ADDITIVE_CONC)), alpha = 0.2) +  # ADDITIVECONC is a covariate
  # scale_x_log10(breaks = c(0, 0.001,0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  ),
  #               labels = scales::number_format(accuracy = 0.001)) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  )) +
  scale_y_continuous(breaks = seq(0, 99, by = 0.5), limits = c(0, 1.8), expand = c(0, 0)) +
  facet_grid(.~CELLINE_NAME)+
  labs(x = expression(paste("Clofarabine concentration (µM)")), col=expression(paste("OA concentration (nM)")), y = "Cell viability") +
  scale_color_manual(values=c("red","blue", "purple", "orange", "green"))+
  stat_function(fun = model_id7, color = "red", size = 1, linetype = "solid") +
  stat_function(fun = model_id8, color = "blue", size = 1, linetype = "solid") +
  geom_segment(data=dat7, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_segment(data=dat7, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_point(data=dat7, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "red", size = 2.5)+
  geom_segment(data=dat8, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_segment(data=dat8, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_point(data=dat8, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "blue", size = 2.5)+
  theme_bw()
p4

p5 <- ggplot() +
  geom_point(data = onco_combine_therapy |>  filter(CELLINE_NAME == "S256F" & DRUG_NAME == "Clofarabine" & ADDITIVE_NAME == "OA") |>  mutate(DRUG_CONC = ifelse(DRUG_CONC == 0, 0.0001, DRUG_CONC)),
             mapping = aes(x = DRUG_CONC, y = DV, colour = as.factor(ADDITIVE_CONC)), alpha = 0.2) +  # ADDITIVECONC is a covariate
  # scale_x_log10(breaks = c(0, 0.001,0.025, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1  ),
  #               labels = scales::number_format(accuracy = 0.001)) +
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(breaks = seq(0, 99, by = 0.5), limits = c(0, 1.8), expand = c(0, 0)) +
  facet_grid(.~CELLINE_NAME)+
  labs(x = expression(paste("Clofarabine concentration (µM)")), col=expression(paste("OA concentration (nM)")), y = "Cell viability") +
  scale_color_manual(values=c("red","blue", "purple", "orange", "green"))+
  stat_function(fun = model_id9, color = "red", size = 1, linetype = "solid") +
  stat_function(fun = model_id10, color = "blue", size = 1, linetype = "solid") +
  geom_segment(data=dat9, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_segment(data=dat9, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "red", alpha = 0.5)+
  geom_point(data=dat9, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "red", size = 2.5)+
  geom_segment(data=dat10, aes(x = TVIC50 , y = 0, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_segment(data=dat10, aes(x = 0 , y = TVBOTTOM + (TVTOP-TVBOTTOM)/2, xend = TVIC50, yend = TVBOTTOM + (TVTOP-TVBOTTOM)/2), linetype = "dashed", colour = "blue", alpha = 0.5)+
  geom_point(data=dat10, aes(x = TVIC50,  y= TVBOTTOM + (TVTOP-TVBOTTOM)/2), shape=8, stroke = 1,colour = "blue", size = 2.5)+
  theme_bw()
p5

figure <- ggarrange(p1,p2,p3,p4,p5, 
                    ncol = 1, nrow = 5,
                    common.legend = TRUE, legend = "top",
                    align = "hv", 
                    font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "bottom"))
figure
ggsave('PRED vs OBS Clofarabine.png', figure , device = "png",width = 15, height = 20, dpi = 500)

################################################################################

### Calculating AUC = drug sensitivity score (DSS) = Integrate #################
models <- c(paste0("model_id",1:10))

# upper concentration varies by drug
upper_DRUGCONC <- dat_2023 |>  group_by(ID) |>  slice(n()) |>  ungroup() |>   select(ID, DRUGCONC)

# calculate effect AUC from 0 to upper drug concentration
DSS <- 0
for(i in 1:10){
  DSS[i] <- integrate(models[i], lower = 0, upper = upper_DRUGCONC$DRUGCONC[i]) # DSS = integrand
  print(DSS[i])
}
DSS <- as.data.frame(DSS)
DSS <- DSS |>  t() |>  as.data.frame() |>  rownames_to_column() |>  select(MODEL = rowname, AUC = V1)

# Since BOTTOM was not always reached at upper drug concentration --> Calculate the area from 0 to upper drug concentration below Y=function(X=upper_DRUGCONC$DRUGCONC)
DSS_bttm <- data.frame(
  model_id1(x=upper_DRUGCONC$DRUGCONC[1]),
  model_id2(x=upper_DRUGCONC$DRUGCONC[2]),
  model_id3(x=upper_DRUGCONC$DRUGCONC[3]),
  model_id4(x=upper_DRUGCONC$DRUGCONC[4]),
  model_id5(x=upper_DRUGCONC$DRUGCONC[5]),
  model_id6(x=upper_DRUGCONC$DRUGCONC[6]),
  model_id7(x=upper_DRUGCONC$DRUGCONC[7]),
  model_id8(x=upper_DRUGCONC$DRUGCONC[8]),
  model_id9(x=upper_DRUGCONC$DRUGCONC[9]),
  model_id10(x=upper_DRUGCONC$DRUGCONC[10])
  )
  
DSS_bttm <- DSS_bttm |>  t() |>  as.data.frame() |>  rownames_to_column()
DSS_bttm <- DSS_bttm |>   mutate(rowname = sub(".x...*.", "", rowname),
                                 upper_DRUGCONC=upper_DRUGCONC$DRUGCONC) |>  mutate(BTTM_AUC = V1*upper_DRUGCONC) |>  select(MODEL = rowname, BTTM = V1, BTTM_AUC,upper_DRUGCONC)
DSS      <- DSS      |>  mutate(ID = row_number())
DSS_bttm <- DSS_bttm |>  mutate(ID = row_number())

dat_2023 <- full_join(dat_2023, DSS      |>  select(ID, AUC           ), by = c("ID"))
dat_2023 <- full_join(dat_2023, DSS_bttm |>  select(ID, BTTM, BTTM_AUC, upper_DRUGCONC), by = c("ID"))

dat_2023 <- dat_2023 |>  mutate(DSS0  = (AUC-(TVBOTTOM*upper_DRUGCONC)), # BOTTOM is the model EGFRameter estimate. BTTM = the value of Y where X=45. BTTM_AUC = BTTM*upper_DRUGCONC = the AUC below BTTM from 0 to 45.
                      DSS0i = (AUC-(BTTM_AUC   )), # Bottom was not always reached at Sora/Regora = 45 --> therefore this.
                      DSS1  = (AUC-(TVBOTTOM*upper_DRUGCONC)) / ((TVTOP-TVBOTTOM)*upper_DRUGCONC),
                      DSS1i = (AUC-(BTTM_AUC   )) / ((TVTOP-BTTM )*upper_DRUGCONC),
                      DSS2  = DSS1/log10(TVTOP),
                      DSS2i = DSS1i / log10(TVTOP),
                      DSS3  = DSS2  *((upper_DRUGCONC)/(upper_DRUGCONC)),
                      DSS3i = DSS2i *((upper_DRUGCONC)/(upper_DRUGCONC))) 

DS_dat <- full_join(onco_combine_therapy |>   select(ID=CELLINE_DRUG_ADDITIVE_ADDITIVECONC, CELLINE_DRUG_ADDITIVE,CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, CELLINE_NAME, DRUG_NAME, ADDITIVE_NAME, CELLINE_DRUG_ADDITIVE_NAME) |>  unique(),
                 dat_2023 |>  select(ID, TVTOP, TVBOTTOM, TVIC50, TVGAMMA, AUC,upper_DRUGCONC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i),
                 by = c("ID"))


################################################################################

### calculate difference in DSS as comEGFRed to no OA (OA=0)#################
DSS_reference_0 <- DS_dat |>  filter(ADDITIVE_CONC==0) |>  
  distinct(CELLINE_DRUG_ADDITIVE, AUC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i) |>  
  select(CELLINE_DRUG_ADDITIVE, 
         ref_AUC=AUC,
         ref_BTTM=BTTM,
         ref_BTTM_AUC=BTTM_AUC,
         ref_DSS0=DSS0,
         ref_DSS1=DSS1,
         ref_DSS2=DSS2,
         ref_DSS3=DSS3,
         ref_DSS0i=DSS0i,
         ref_DSS1i=DSS1i,
         ref_DSS2i=DSS2i,
         ref_DSS3i=DSS3i)

DSS <- left_join(DS_dat, DSS_reference_0,by="CELLINE_DRUG_ADDITIVE")
DSS <- DSS |>  mutate(dDSS0   = DSS0   - ref_DSS0,
                      dDSS1   = DSS1   - ref_DSS1,
                      dDSS2   = DSS2   - ref_DSS2,
                      dDSS3   = DSS3   - ref_DSS3,
                      dDSS0i  = DSS0i  - ref_DSS0i,
                      dDSS1i  = DSS1i  - ref_DSS1i,
                      dDSS2i  = DSS2i  - ref_DSS2i,
                      dDSS3i  = DSS3i  - ref_DSS3i) 


################################################################################

### plot the dDSS dsitribution #################################################
# Plots
group.colors <- c(EGFP  ="#f8766d",
                  WT   ="#d39200", 
                  R183W ="#93aa00",
                  P179R  ="#00ba38",
                  S256F   ="#FFDB6D")
p1 <- ggplot(data = DSS, mapping = aes(x = ADDITIVE_CONC, y = dDSS0,  colour = CELLINE_NAME)) +
  #facet_wrap(~ CELLINE_NAME, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey") +
  geom_point(size = 3) +
  geom_line(size = 0.75, linetype = "dashed") +
  labs(x = expression(paste("OA concentration (nM)")), y = "dDSS0", col = "Cell line") +
  scale_x_continuous(breaks = c(0, 10)) +
  #scale_y_continuous(breaks = seq(-100, 100, 1), limits = c(-6, 0.5)) +
  scale_colour_manual(values=group.colors)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0), legend.position = "right", panel.grid.minor = element_blank())
p1
ggsave('dDSS0.png', p1 , device = "png",width = 7, height = 5.5, dpi = 500)



### Plot predicted AUC vs AUC based on trapezoidal rule ########################
### calculate AUC based on trapezoidal rule 
dat_trapezoid <- cbind(DSS, onco_combine_therapy |>  select(DRUG_CONC,DV,CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT,REP))
dat_trapezoid <- dat_trapezoid |>  
  group_by(CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT,REP) |>  
  mutate(ID_trapez = cur_group_id()) |>   # ID_trapez based on CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT_REP
  ungroup() |> 
  mutate(AUC_trapez_i=NA)

library(pracma)
IDs= unique(dat_trapezoid$ID_trapez)
for (i in IDs) {
  data <- dat_trapezoid |>  filter(ID_trapez==i)|>  select(ID_trapez, DRUG_CONC,DV)
  AUC_value <- trapz(data$DRUG_CONC, data$DV)
  dat_trapezoid$AUC_trapez_i[dat_trapezoid$ID_trapez==i] = AUC_value
}
dat_trapezoid <- dat_trapezoid |>  
  group_by(ID) |>  
  #mutate(AUC_trapez=mean(AUC_trapez_i)) |> 
  mutate(AUC_trapez=median(AUC_trapez_i)) |> 
  ungroup() 

group.colors <- c(EGFP  ="#f8766d",
                  WT   ="#d39200", 
                  R183W ="#93aa00",
                  P179R  ="#00ba38",
                  S256F   ="#FFDB6D")
p1 <- ggplot (data = dat_trapezoid , aes(x= AUC_trapez, y= AUC, colour=CELLINE_NAME)) + 
  geom_point( size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_x_continuous(limits=c(0,0.5)) +
  scale_y_continuous(limits=c(0,0.5))+
  geom_smooth(se = F,colour = "red", method = "lm", size = 1.2, alpha=0.1) +
  scale_colour_manual(values=group.colors)+
  theme_bw()+
  theme(text = element_text(size = 12)) +
  labs(x = "Predicted AUC", y = "Observed AUC", col = "Cell line") 
p1
ggsave('Observed AUC vs Predicted AUC.png', p1 , device = "png",width = 7.5, height = 6, dpi = 500)


# ############################# use ID+REPEAT as new ID ##################################
# dat_2023 <- read.table (file='~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023//modeling/run003_out.csv', skip=1, header=T)
# dat_2023$CELLINE <- ordered(dat_2023$CELLINE, levels=c(1, 2, 3, 4, 5),
#                             labels=c("EGFP", "WT", "R183W", "P179R", "S256F"))
# dat_2023 <- dat_2023 |>  group_by(ID,DAY)   |>  mutate(ID_DAY   = cur_group_id()) |>  ungroup()
# 
# setwd("~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023/plots/")
# 
# #### define the equations for effect curve ####################################
# # 1-10
# dat1 <- dat_2023 |>  filter(ID_DAY == 1) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id1 <- function(x) { dat1$BOTTOM + ( (dat1$TOP-dat1$BOTTOM) / ( 1+ ((x/dat1$IC50)^dat1$GAMMA) ) )}
# dat2 <- dat_2023 |>  filter(ID_DAY == 2) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id2 <- function(x) { dat2$BOTTOM + ( (dat2$TOP-dat2$BOTTOM) / ( 1+ ((x/dat2$IC50)^dat2$GAMMA) ) )}
# dat3 <- dat_2023 |>  filter(ID_DAY == 3) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id3 <- function(x) { dat3$BOTTOM + ( (dat3$TOP-dat3$BOTTOM) / ( 1+ ((x/dat3$IC50)^dat3$GAMMA) ) )}
# dat4 <- dat_2023 |>  filter(ID_DAY == 4) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id4 <- function(x) { dat4$BOTTOM + ( (dat4$TOP-dat4$BOTTOM) / ( 1+ ((x/dat4$IC50)^dat4$GAMMA) ) )}
# dat5 <- dat_2023 |>  filter(ID_DAY == 5) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id5 <- function(x) { dat5$BOTTOM + ( (dat5$TOP-dat5$BOTTOM) / ( 1+ ((x/dat5$IC50)^dat5$GAMMA) ) )}
# dat6 <- dat_2023 |>  filter(ID_DAY == 6) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id6 <- function(x) { dat6$BOTTOM + ( (dat6$TOP-dat6$BOTTOM) / ( 1+ ((x/dat6$IC50)^dat6$GAMMA) ) )}
# dat7 <- dat_2023 |>  filter(ID_DAY == 7) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id7 <- function(x) { dat7$BOTTOM + ( (dat7$TOP-dat7$BOTTOM) / ( 1+ ((x/dat7$IC50)^dat7$GAMMA) ) )}
# dat8 <- dat_2023 |>  filter(ID_DAY == 8) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id8 <- function(x) { dat8$BOTTOM + ( (dat8$TOP-dat8$BOTTOM) / ( 1+ ((x/dat8$IC50)^dat8$GAMMA) ) )}
# dat9 <- dat_2023 |>  filter(ID_DAY == 9) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id9 <- function(x) { dat9$BOTTOM + ( (dat9$TOP-dat9$BOTTOM) / ( 1+ ((x/dat9$IC50)^dat9$GAMMA) ) )}
# 
# # ID_DAY 10-19
# dat10 <- dat_2023 |>  filter(ID_DAY == 10) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id10 <- function(x) { dat10$BOTTOM + ( (dat10$TOP-dat10$BOTTOM) / ( 1+ ((x/dat10$IC50)^dat10$GAMMA) ) )}
# dat11 <- dat_2023 |>  filter(ID_DAY == 11) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id11 <- function(x) { dat11$BOTTOM + ( (dat11$TOP-dat11$BOTTOM) / ( 1+ ((x/dat11$IC50)^dat11$GAMMA) ) )}
# dat12 <- dat_2023 |>  filter(ID_DAY == 12) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id12 <- function(x) { dat12$BOTTOM + ( (dat12$TOP-dat12$BOTTOM) / ( 1+ ((x/dat12$IC50)^dat12$GAMMA) ) )}
# dat13 <- dat_2023 |>  filter(ID_DAY == 13) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id13 <- function(x) { dat13$BOTTOM + ( (dat13$TOP-dat13$BOTTOM) / ( 1+ ((x/dat13$IC50)^dat13$GAMMA) ) )}
# dat14 <- dat_2023 |>  filter(ID_DAY == 14) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id14 <- function(x) { dat14$BOTTOM + ( (dat14$TOP-dat14$BOTTOM) / ( 1+ ((x/dat14$IC50)^dat14$GAMMA) ) )}
# dat15 <- dat_2023 |>  filter(ID_DAY == 15) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id15 <- function(x) { dat15$BOTTOM + ( (dat15$TOP-dat15$BOTTOM) / ( 1+ ((x/dat15$IC50)^dat15$GAMMA) ) )}
# dat16 <- dat_2023 |>  filter(ID_DAY == 16) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id16 <- function(x) { dat16$BOTTOM + ( (dat16$TOP-dat16$BOTTOM) / ( 1+ ((x/dat16$IC50)^dat16$GAMMA) ) )}
# dat17 <- dat_2023 |>  filter(ID_DAY == 17) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id17 <- function(x) { dat17$BOTTOM + ( (dat17$TOP-dat17$BOTTOM) / ( 1+ ((x/dat17$IC50)^dat17$GAMMA) ) )}
# dat18 <- dat_2023 |>  filter(ID_DAY == 18) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id18 <- function(x) { dat18$BOTTOM + ( (dat18$TOP-dat18$BOTTOM) / ( 1+ ((x/dat18$IC50)^dat18$GAMMA) ) )}
# dat19 <- dat_2023 |>  filter(ID_DAY == 19) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id19 <- function(x) { dat19$BOTTOM + ( (dat19$TOP-dat19$BOTTOM) / ( 1+ ((x/dat19$IC50)^dat19$GAMMA) ) )}
# 
# # ID_DAY 20-29
# dat20 <- dat_2023 |>  filter(ID_DAY == 20) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id20 <- function(x) { dat20$BOTTOM + ( (dat20$TOP-dat20$BOTTOM) / ( 1+ ((x/dat20$IC50)^dat20$GAMMA) ) )}
# dat21 <- dat_2023 |>  filter(ID_DAY == 21) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id21 <- function(x) { dat21$BOTTOM + ( (dat21$TOP-dat21$BOTTOM) / ( 1+ ((x/dat21$IC50)^dat21$GAMMA) ) )}
# dat22 <- dat_2023 |>  filter(ID_DAY == 22) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id22 <- function(x) { dat22$BOTTOM + ( (dat22$TOP-dat22$BOTTOM) / ( 1+ ((x/dat22$IC50)^dat22$GAMMA) ) )}
# dat23 <- dat_2023 |>  filter(ID_DAY == 23) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id23 <- function(x) { dat23$BOTTOM + ( (dat23$TOP-dat23$BOTTOM) / ( 1+ ((x/dat23$IC50)^dat23$GAMMA) ) )}
# dat24 <- dat_2023 |>  filter(ID_DAY == 24) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id24 <- function(x) { dat24$BOTTOM + ( (dat24$TOP-dat24$BOTTOM) / ( 1+ ((x/dat24$IC50)^dat24$GAMMA) ) )}
# dat25 <- dat_2023 |>  filter(ID_DAY == 25) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id25 <- function(x) { dat25$BOTTOM + ( (dat25$TOP-dat25$BOTTOM) / ( 1+ ((x/dat25$IC50)^dat25$GAMMA) ) )}
# dat26 <- dat_2023 |>  filter(ID_DAY == 26) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id26 <- function(x) { dat26$BOTTOM + ( (dat26$TOP-dat26$BOTTOM) / ( 1+ ((x/dat26$IC50)^dat26$GAMMA) ) )}
# dat27 <- dat_2023 |>  filter(ID_DAY == 27) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id27 <- function(x) { dat27$BOTTOM + ( (dat27$TOP-dat27$BOTTOM) / ( 1+ ((x/dat27$IC50)^dat27$GAMMA) ) )}
# dat28 <- dat_2023 |>  filter(ID_DAY == 28) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id28 <- function(x) { dat28$BOTTOM + ( (dat28$TOP-dat28$BOTTOM) / ( 1+ ((x/dat28$IC50)^dat28$GAMMA) ) )}
# dat29 <- dat_2023 |>  filter(ID_DAY == 29) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id29 <- function(x) { dat29$BOTTOM + ( (dat29$TOP-dat29$BOTTOM) / ( 1+ ((x/dat29$IC50)^dat29$GAMMA) ) )}
# 
# # ID_DAY 30-39
# dat30 <- dat_2023 |>  filter(ID_DAY == 30) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id30 <- function(x) { dat30$BOTTOM + ( (dat30$TOP-dat30$BOTTOM) / ( 1+ ((x/dat30$IC50)^dat30$GAMMA) ) )}
# dat31 <- dat_2023 |>  filter(ID_DAY == 31) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id31 <- function(x) { dat31$BOTTOM + ( (dat31$TOP-dat31$BOTTOM) / ( 1+ ((x/dat31$IC50)^dat31$GAMMA) ) )}
# dat32 <- dat_2023 |>  filter(ID_DAY == 32) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id32 <- function(x) { dat32$BOTTOM + ( (dat32$TOP-dat32$BOTTOM) / ( 1+ ((x/dat32$IC50)^dat32$GAMMA) ) )}
# dat33 <- dat_2023 |>  filter(ID_DAY == 33) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id33 <- function(x) { dat33$BOTTOM + ( (dat33$TOP-dat33$BOTTOM) / ( 1+ ((x/dat33$IC50)^dat33$GAMMA) ) )}
# dat34 <- dat_2023 |>  filter(ID_DAY == 34) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id34 <- function(x) { dat34$BOTTOM + ( (dat34$TOP-dat34$BOTTOM) / ( 1+ ((x/dat34$IC50)^dat34$GAMMA) ) )}
# dat35 <- dat_2023 |>  filter(ID_DAY == 35) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id35 <- function(x) { dat35$BOTTOM + ( (dat35$TOP-dat35$BOTTOM) / ( 1+ ((x/dat35$IC50)^dat35$GAMMA) ) )}
# dat36 <- dat_2023 |>  filter(ID_DAY == 36) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id36 <- function(x) { dat36$BOTTOM + ( (dat36$TOP-dat36$BOTTOM) / ( 1+ ((x/dat36$IC50)^dat36$GAMMA) ) )}
# dat37 <- dat_2023 |>  filter(ID_DAY == 37) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id37 <- function(x) { dat37$BOTTOM + ( (dat37$TOP-dat37$BOTTOM) / ( 1+ ((x/dat37$IC50)^dat37$GAMMA) ) )}
# dat38 <- dat_2023 |>  filter(ID_DAY == 38) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id38 <- function(x) { dat38$BOTTOM + ( (dat38$TOP-dat38$BOTTOM) / ( 1+ ((x/dat38$IC50)^dat38$GAMMA) ) )}
# dat39 <- dat_2023 |>  filter(ID_DAY == 39) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id39 <- function(x) { dat39$BOTTOM + ( (dat39$TOP-dat39$BOTTOM) / ( 1+ ((x/dat39$IC50)^dat39$GAMMA) ) )}
# 
# # ID_DAY 40-49
# dat40 <- dat_2023 |>  filter(ID_DAY == 40) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id40 <- function(x) { dat40$BOTTOM + ( (dat40$TOP-dat40$BOTTOM) / ( 1+ ((x/dat40$IC50)^dat40$GAMMA) ) )}
# dat41 <- dat_2023 |>  filter(ID_DAY == 41) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id41 <- function(x) { dat41$BOTTOM + ( (dat41$TOP-dat41$BOTTOM) / ( 1+ ((x/dat41$IC50)^dat41$GAMMA) ) )}
# dat42 <- dat_2023 |>  filter(ID_DAY == 42) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id42 <- function(x) { dat42$BOTTOM + ( (dat42$TOP-dat42$BOTTOM) / ( 1+ ((x/dat42$IC50)^dat42$GAMMA) ) )}
# dat43 <- dat_2023 |>  filter(ID_DAY == 43) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id43 <- function(x) { dat43$BOTTOM + ( (dat43$TOP-dat43$BOTTOM) / ( 1+ ((x/dat43$IC50)^dat43$GAMMA) ) )}
# dat44 <- dat_2023 |>  filter(ID_DAY == 44) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id44 <- function(x) { dat44$BOTTOM + ( (dat44$TOP-dat44$BOTTOM) / ( 1+ ((x/dat44$IC50)^dat44$GAMMA) ) )}
# dat45 <- dat_2023 |>  filter(ID_DAY == 45) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id45 <- function(x) { dat45$BOTTOM + ( (dat45$TOP-dat45$BOTTOM) / ( 1+ ((x/dat45$IC50)^dat45$GAMMA) ) )}
# dat46 <- dat_2023 |>  filter(ID_DAY == 46) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id46 <- function(x) { dat46$BOTTOM + ( (dat46$TOP-dat46$BOTTOM) / ( 1+ ((x/dat46$IC50)^dat46$GAMMA) ) )}
# dat47 <- dat_2023 |>  filter(ID_DAY == 47) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id47 <- function(x) { dat47$BOTTOM + ( (dat47$TOP-dat47$BOTTOM) / ( 1+ ((x/dat47$IC50)^dat47$GAMMA) ) )}
# dat48 <- dat_2023 |>  filter(ID_DAY == 48) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id48 <- function(x) { dat48$BOTTOM + ( (dat48$TOP-dat48$BOTTOM) / ( 1+ ((x/dat48$IC50)^dat48$GAMMA) ) )}
# dat49 <- dat_2023 |>  filter(ID_DAY == 49) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id49 <- function(x) { dat49$BOTTOM + ( (dat49$TOP-dat49$BOTTOM) / ( 1+ ((x/dat49$IC50)^dat49$GAMMA) ) )}
# 
# # ID_DAY 50-59
# dat50 <- dat_2023 |>  filter(ID_DAY == 50) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id50 <- function(x) { dat50$BOTTOM + ( (dat50$TOP-dat50$BOTTOM) / ( 1+ ((x/dat50$IC50)^dat50$GAMMA) ) )}
# dat51 <- dat_2023 |>  filter(ID_DAY == 51) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id51 <- function(x) { dat51$BOTTOM + ( (dat51$TOP-dat51$BOTTOM) / ( 1+ ((x/dat51$IC50)^dat51$GAMMA) ) )}
# dat52 <- dat_2023 |>  filter(ID_DAY == 52) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id52 <- function(x) { dat52$BOTTOM + ( (dat52$TOP-dat52$BOTTOM) / ( 1+ ((x/dat52$IC50)^dat52$GAMMA) ) )}
# dat53 <- dat_2023 |>  filter(ID_DAY == 53) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id53 <- function(x) { dat53$BOTTOM + ( (dat53$TOP-dat53$BOTTOM) / ( 1+ ((x/dat53$IC50)^dat53$GAMMA) ) )}
# dat54 <- dat_2023 |>  filter(ID_DAY == 54) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id54 <- function(x) { dat54$BOTTOM + ( (dat54$TOP-dat54$BOTTOM) / ( 1+ ((x/dat54$IC50)^dat54$GAMMA) ) )}
# dat55 <- dat_2023 |>  filter(ID_DAY == 55) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id55 <- function(x) { dat55$BOTTOM + ( (dat55$TOP-dat55$BOTTOM) / ( 1+ ((x/dat55$IC50)^dat55$GAMMA) ) )}
# dat56 <- dat_2023 |>  filter(ID_DAY == 56) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id56 <- function(x) { dat56$BOTTOM + ( (dat56$TOP-dat56$BOTTOM) / ( 1+ ((x/dat56$IC50)^dat56$GAMMA) ) )}
# dat57 <- dat_2023 |>  filter(ID_DAY == 57) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id57 <- function(x) { dat57$BOTTOM + ( (dat57$TOP-dat57$BOTTOM) / ( 1+ ((x/dat57$IC50)^dat57$GAMMA) ) )}
# dat58 <- dat_2023 |>  filter(ID_DAY == 58) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id58 <- function(x) { dat58$BOTTOM + ( (dat58$TOP-dat58$BOTTOM) / ( 1+ ((x/dat58$IC50)^dat58$GAMMA) ) )}
# dat59 <- dat_2023 |>  filter(ID_DAY == 59) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id59 <- function(x) { dat59$BOTTOM + ( (dat59$TOP-dat59$BOTTOM) / ( 1+ ((x/dat59$IC50)^dat59$GAMMA) ) )}
# 
# # ID_DAY 60-66
# dat60 <- dat_2023 |>  filter(ID_DAY == 60) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id60 <- function(x) { dat60$BOTTOM + ( (dat60$TOP-dat60$BOTTOM) / ( 1+ ((x/dat60$IC50)^dat60$GAMMA) ) )}
# dat61 <- dat_2023 |>  filter(ID_DAY == 61) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id61 <- function(x) { dat61$BOTTOM + ( (dat61$TOP-dat61$BOTTOM) / ( 1+ ((x/dat61$IC50)^dat61$GAMMA) ) )}
# dat62 <- dat_2023 |>  filter(ID_DAY == 62) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id62 <- function(x) { dat62$BOTTOM + ( (dat62$TOP-dat62$BOTTOM) / ( 1+ ((x/dat62$IC50)^dat62$GAMMA) ) )}
# dat63 <- dat_2023 |>  filter(ID_DAY == 63) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id63 <- function(x) { dat63$BOTTOM + ( (dat63$TOP-dat63$BOTTOM) / ( 1+ ((x/dat63$IC50)^dat63$GAMMA) ) )}
# dat64 <- dat_2023 |>  filter(ID_DAY == 64) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id64 <- function(x) { dat64$BOTTOM + ( (dat64$TOP-dat64$BOTTOM) / ( 1+ ((x/dat64$IC50)^dat64$GAMMA) ) )}
# dat65 <- dat_2023 |>  filter(ID_DAY == 65) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id65 <- function(x) { dat65$BOTTOM + ( (dat65$TOP-dat65$BOTTOM) / ( 1+ ((x/dat65$IC50)^dat65$GAMMA) ) )}
# dat66 <- dat_2023 |>  filter(ID_DAY == 66) |>  distinct(ID_DAY,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)
# model_id66 <- function(x) { dat66$BOTTOM + ( (dat66$TOP-dat66$BOTTOM) / ( 1+ ((x/dat66$IC50)^dat66$GAMMA) ) )}
# 
# 
# ### Calculating AUC = drug sensitivity score (DSS) = Integrate #################
# models <- c(paste0("model_id",1:66))
# 
# # upper concentration varies by drug
# upper_DRUGCONC <- dat_2023 |>  group_by(ID_DAY) |>  slice(n()) |>  ungroup() |>   select(ID_DAY, DRUGCONC)
# 
# # calculate effect AUC from 0 to upper drug concentration
# DSS <- 0
# for(i in 1:66){
#   DSS[i] <- integrate(models[i], lower = 0, upper = upper_DRUGCONC$DRUGCONC[i]) # DSS = integrand
#   print(DSS[i])
# }
# DSS <- as.data.frame(DSS)
# DSS <- DSS |>  t() |>  as.data.frame() |>  rownames_to_column() |>  select(MODEL = rowname, AUC = V1)
# 
# # Since BOTTOM was not always reached at upper drug concentration --> Calculate the area from 0 to upper drug concentration below Y=function(X=upper_DRUGCONC$DRUGCONC)
# DSS_bttm <- data.frame(
#   model_id1(x=upper_DRUGCONC$DRUGCONC[1]),
#   model_id2(x=upper_DRUGCONC$DRUGCONC[2]),
#   model_id3(x=upper_DRUGCONC$DRUGCONC[3]),
#   model_id4(x=upper_DRUGCONC$DRUGCONC[4]),
#   model_id5(x=upper_DRUGCONC$DRUGCONC[5]),
#   model_id6(x=upper_DRUGCONC$DRUGCONC[6]),
#   model_id7(x=upper_DRUGCONC$DRUGCONC[7]),
#   model_id8(x=upper_DRUGCONC$DRUGCONC[8]),
#   model_id9(x=upper_DRUGCONC$DRUGCONC[9]),
#   model_id10(x=upper_DRUGCONC$DRUGCONC[10]),
#   model_id11(x=upper_DRUGCONC$DRUGCONC[11]),
#   model_id12(x=upper_DRUGCONC$DRUGCONC[12]),
#   model_id13(x=upper_DRUGCONC$DRUGCONC[13]),
#   model_id14(x=upper_DRUGCONC$DRUGCONC[14]),
#   model_id15(x=upper_DRUGCONC$DRUGCONC[15]),
#   model_id16(x=upper_DRUGCONC$DRUGCONC[16]),
#   model_id17(x=upper_DRUGCONC$DRUGCONC[17]),
#   model_id18(x=upper_DRUGCONC$DRUGCONC[18]),
#   model_id19(x=upper_DRUGCONC$DRUGCONC[19]),
#   model_id20(x=upper_DRUGCONC$DRUGCONC[20]),
#   model_id21(x=upper_DRUGCONC$DRUGCONC[21]),
#   model_id22(x=upper_DRUGCONC$DRUGCONC[22]),
#   model_id23(x=upper_DRUGCONC$DRUGCONC[23]),
#   model_id24(x=upper_DRUGCONC$DRUGCONC[24]),
#   model_id25(x=upper_DRUGCONC$DRUGCONC[25]),
#   model_id26(x=upper_DRUGCONC$DRUGCONC[26]),
#   model_id27(x=upper_DRUGCONC$DRUGCONC[27]),
#   model_id28(x=upper_DRUGCONC$DRUGCONC[28]),
#   model_id29(x=upper_DRUGCONC$DRUGCONC[29]),
#   model_id30(x=upper_DRUGCONC$DRUGCONC[30]),
#   model_id31(x=upper_DRUGCONC$DRUGCONC[31]),
#   model_id32(x=upper_DRUGCONC$DRUGCONC[32]),
#   model_id33(x=upper_DRUGCONC$DRUGCONC[33]),
#   model_id34(x=upper_DRUGCONC$DRUGCONC[34]),
#   model_id35(x=upper_DRUGCONC$DRUGCONC[35]),
#   model_id36(x=upper_DRUGCONC$DRUGCONC[36]),
#   model_id37(x=upper_DRUGCONC$DRUGCONC[37]),
#   model_id38(x=upper_DRUGCONC$DRUGCONC[38]),
#   model_id39(x=upper_DRUGCONC$DRUGCONC[39]),
#   model_id40(x=upper_DRUGCONC$DRUGCONC[40]),
#   model_id41(x=upper_DRUGCONC$DRUGCONC[41]),
#   model_id42(x=upper_DRUGCONC$DRUGCONC[42]),
#   model_id43(x=upper_DRUGCONC$DRUGCONC[43]),
#   model_id44(x=upper_DRUGCONC$DRUGCONC[44]),
#   model_id45(x=upper_DRUGCONC$DRUGCONC[45]),
#   model_id46(x=upper_DRUGCONC$DRUGCONC[46]),
#   model_id47(x=upper_DRUGCONC$DRUGCONC[47]),
#   model_id48(x=upper_DRUGCONC$DRUGCONC[48]),
#   model_id49(x=upper_DRUGCONC$DRUGCONC[49]),
#   model_id50(x=upper_DRUGCONC$DRUGCONC[50]),
#   model_id51(x=upper_DRUGCONC$DRUGCONC[51]),
#   model_id52(x=upper_DRUGCONC$DRUGCONC[52]),
#   model_id53(x=upper_DRUGCONC$DRUGCONC[53]),
#   model_id54(x=upper_DRUGCONC$DRUGCONC[54]),
#   model_id55(x=upper_DRUGCONC$DRUGCONC[55]),
#   model_id56(x=upper_DRUGCONC$DRUGCONC[56]),
#   model_id57(x=upper_DRUGCONC$DRUGCONC[57]),
#   model_id58(x=upper_DRUGCONC$DRUGCONC[58]),
#   model_id59(x=upper_DRUGCONC$DRUGCONC[59]),
#   model_id60(x=upper_DRUGCONC$DRUGCONC[60]),
#   model_id61(x=upper_DRUGCONC$DRUGCONC[61]),
#   model_id62(x=upper_DRUGCONC$DRUGCONC[62]),
#   model_id63(x=upper_DRUGCONC$DRUGCONC[63]),
#   model_id64(x=upper_DRUGCONC$DRUGCONC[64]),
#   model_id65(x=upper_DRUGCONC$DRUGCONC[65]),
#   model_id66(x=upper_DRUGCONC$DRUGCONC[66]))
# 
# DSS_bttm <- DSS_bttm |>  t() |>  as.data.frame() |>  rownames_to_column()
# DSS_bttm <- DSS_bttm |>   mutate(rowname = sub(".x...*.", "", rowname),
#                                  upper_DRUGCONC=upper_DRUGCONC$DRUGCONC) |>  mutate(BTTM_AUC = V1*upper_DRUGCONC) |>  select(MODEL = rowname, BTTM = V1, BTTM_AUC,upper_DRUGCONC)
# DSS      <- DSS      |>  mutate(ID_DAY = row_number())
# DSS_bttm <- DSS_bttm |>  mutate(ID_DAY = row_number())
# 
# dat_2023 <- full_join(dat_2023, DSS      |>  select(ID_DAY, AUC           ), by = c("ID_DAY"))
# dat_2023 <- full_join(dat_2023, DSS_bttm |>  select(ID_DAY, BTTM, BTTM_AUC, upper_DRUGCONC), by = c("ID_DAY"))
# 
# dat_2023 <- dat_2023 |>  mutate(DSS0  = (AUC-(BOTTOM*upper_DRUGCONC)), # BOTTOM is the model EGFRameter estimate. BTTM = the value of Y where X=45. BTTM_AUC = BTTM*upper_DRUGCONC = the AUC below BTTM from 0 to 45.
#                                 DSS0i = (AUC-(BTTM_AUC   )), # Bottom was not always reached at Sora/Regora = 45 --> therefore this.
#                                 DSS1  = (AUC-(BOTTOM*upper_DRUGCONC)) / ((TOP-BOTTOM)*upper_DRUGCONC),
#                                 DSS1i = (AUC-(BTTM_AUC   )) / ((TOP-BTTM )*upper_DRUGCONC),
#                                 DSS2  = DSS1/log10(TOP),
#                                 DSS2i = DSS1i / log10(TOP),
#                                 DSS3  = DSS2  *((upper_DRUGCONC)/(upper_DRUGCONC)),
#                                 DSS3i = DSS2i *((upper_DRUGCONC)/(upper_DRUGCONC))) 
# 
# DS_dat <- full_join(onco_combine_therapy |>   select(ID_DAY=CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT, ID = CELLINE_DRUG_ADDITIVE_ADDITIVECONC, CELLINE_DRUG_ADDITIVE,CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, DAY=REPEAT, CELLINE_NAME, DRUG_NAME, ADDITIVE_NAME, CELLINE_DRUG_ADDITIVE_NAME) |>  unique() |> 
#                       group_by(CELLINE_DRUG_ADDITIVE, DAY)   |>  mutate(CELLINE_DRUG_ADDITIVE_DAY   = cur_group_id()) |>  ungroup(), 
#                     dat_2023 |>  select(ID_DAY, TOP, BOTTOM, IC50, GAMMA, AUC,upper_DRUGCONC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i),
#                     by = c("ID_DAY"))
# 
# 
# ################################################################################
# 
# ### calculate difference in DSS as compared to no OA (OA=0)#################
# DSS_reference_0 <- DS_dat |>  filter(ADDITIVE_CONC==0) |>  
#   distinct(CELLINE_DRUG_ADDITIVE_DAY, ID,AUC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i) |>   # group by CELLINE_DRUG_ADDITIVE_DAY
#   select(CELLINE_DRUG_ADDITIVE_DAY, 
#          ref_AUC=AUC,
#          ref_BTTM=BTTM,
#          ref_BTTM_AUC=BTTM_AUC,
#          ref_DSS0=DSS0,
#          ref_DSS1=DSS1,
#          ref_DSS2=DSS2,
#          ref_DSS3=DSS3,
#          ref_DSS0i=DSS0i,
#          ref_DSS1i=DSS1i,
#          ref_DSS2i=DSS2i,
#          ref_DSS3i=DSS3i)
# 
# DSS <- left_join(DS_dat, DSS_reference_0,by="CELLINE_DRUG_ADDITIVE_DAY")
# DSS <- DSS |>  mutate(dDSS0   = DSS0   - ref_DSS0,
#                       dDSS1   = DSS1   - ref_DSS1,
#                       dDSS2   = DSS2   - ref_DSS2,
#                       dDSS3   = DSS3   - ref_DSS3,
#                       dDSS0i  = DSS0i  - ref_DSS0i,
#                       dDSS1i  = DSS1i  - ref_DSS1i,
#                       dDSS2i  = DSS2i  - ref_DSS2i,
#                       dDSS3i  = DSS3i  - ref_DSS3i) 
# 
# ### plot the dDSS dsitribution (with distribution bar) #################################################
# # Plots
# group.colors <- c(EGFP  ="#f8766d",
#                   WT   ="#d39200", 
#                   R183W ="#93aa00",
#                   P179R  ="#00ba38",
#                   S256F   ="#FFDB6D")
# DSS_ordered <- DSS |>  filter(ADDITIVE_CONC!=0) |>  distinct(ID_DAY,.keep_all = T)
# DSS_ordered$CELLINE_NAME <- factor(DSS_ordered$CELLINE_NAME, 
#                                    levels = c( "WT", "S256F","EGFP", "P179R","R183W" ))
# 
# p1 <- ggplot(data = DSS_ordered, mapping = aes(x = CELLINE_NAME, y = dDSS0)) +
#   #facet_wrap(~ CELLINE_NAME, scales = "free_x") +
#   geom_hline(yintercept = 0, linetype = "solid", colour = "grey") +
#   geom_boxplot()+
#   geom_point(size = 1.5) +
#   #geom_line(size = 0.75, linetype = "dashed") +
#   labs(x = expression(paste("Cell line")), y = "dDSS0 (OA= 10 nM)", col = "Cell line") +
#   #scale_x_continuous(breaks = c(0, 10)) +
#   #scale_y_continuous(breaks = seq(-100, 100, 1), limits = c(-6, 0.5)) +
#   scale_colour_manual(values=group.colors)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 0), legend.position = "right", panel.grid.minor = element_blank())
# p1
# ggsave('dDSS0_Group by DAY.png', p1 , device = "png",width = 7, height = 5.5, dpi = 500)
# ################################################################################






############################# simulations ##################################
dat_sim <- read.csv(file='~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023//modeling/run006_sim_out.csv', skip = 1, header = TRUE, sep = "")

dat_sim <-dat_sim |>   filter(ID!="TABLE" & ID!="ID" ) |>           # remove two rows for each simulation 
  mutate_all(function(x) as.numeric(as.character(x))) |>            # convert dataframe to numeric
  filter(DRUGCONC==0|DRUGCONC==1) |> 
  mutate(SIM=ceiling(row_number()/396))                             # label each time of simulation

dat_sim <- dat_sim |>  
  arrange(ID,DAY,REP,SIM) |> 
  group_by(ID,DAY,REP)     |>  mutate(ID_DAY_REP   = cur_group_id())     |>  ungroup() |>  # For each replicate, there is a 
  group_by(ID,DAY,REP,SIM) |>  mutate(ID_DAY_REP_SIM   = cur_group_id()) |>  ungroup() |>  
  select( ID_DAY_REP,ID_DAY_REP_SIM,SIM, everything())
  
dat_sim$CELLINE <- ordered(dat_sim$CELLINE, levels=c(1, 2, 3, 4, 5),
                            labels=c("EGFP", "WT", "R183W", "P179R", "S256F"))

setwd("~/Library/CloudStorage/OneDrive-KULeuven/PhD Garin/new Michiel_combinational_therapy_analysis_2023/plots/")

# upper concentration varies by drug
upper_DRUGCONC <- dat_sim |>  group_by(ID_DAY_REP_SIM) |>  slice(n()) |>  ungroup() |>   select(ID_DAY_REP_SIM, DRUGCONC) 

#### define the equations for effect curve ####################################
DSS <- 0
DSS_bttm <- 0
for(i in 1:198000){                                                                                                        # in total 198000 ID_DAY_REP_SIMs
  dat <- dat_sim |>  filter(ID_DAY_REP_SIM == i) |>  distinct(ID_DAY_REP_SIM,CELLINE,ADDITIVECONC,TOP,BOTTOM,IC50,GAMMA)   # get the parameter estimates of each ID
  equation = function(x) { dat$BOTTOM + ( (dat$TOP-dat$BOTTOM) / ( 1+ ((x/dat$IC50)^dat$GAMMA) ) )}                        # define the equation for cell viability calculation
  DSS[i] <- integrate(equation, lower = 0, upper = 1)                                                                      # calculate the AUC from the lowest drug concentration to the highest drug concentration
  DSS_bttm[i] = equation(1) *1                                                                                             # calculate the area from 0-1 x cell viability at the highest drug concentration
  print(i)
}
DSS <- as.data.frame(DSS)
DSS <- DSS |>  t() |>  as.data.frame() |>  rownames_to_column() |>  select(MODEL = rowname, AUC = V1) |>  mutate(ID_DAY_REP_SIM = row_number())

DSS_bttm <- as.data.frame(DSS_bttm)
DSS_bttm <- DSS_bttm  |>  mutate(ID_DAY_REP_SIM = row_number(),
                                 BTTM=DSS_bttm,         # upper drug concentration is 1, thereby the minimum cell viability = BTTM_AUC
                                 upper_DRUGCONC=1) |>  select(ID_DAY_REP_SIM, BTTM_AUC = DSS_bttm,BTTM,upper_DRUGCONC)

dat_sim <- full_join(dat_sim, DSS      |>  select(ID_DAY_REP_SIM, AUC), by = c("ID_DAY_REP_SIM"))
dat_sim <- full_join(dat_sim, DSS_bttm,  by = c("ID_DAY_REP_SIM"))
dat_sim <- dat_sim |>  distinct(ID_DAY_REP_SIM, .keep_all=T) |>  group_by(ID,SIM) |> 
  mutate(AUC_median      = median(AUC),
         BTTM_AUC_median = median(BTTM_AUC),
         BOTTOM_median   = median(BOTTOM),
         TOP_median      = median(TOP),
         BTTM_median     = median(BTTM))
  

dat_sim <- dat_sim |>    mutate(DSS0  = (AUC_median -(BOTTOM_median*upper_DRUGCONC)), # BOTTOM is the model parameter estimate. BTTM = the value of Y where X=1. BTTM_AUC = BTTM*upper_DRUGCONC = the AUC below BTTM from 0 to 1.
                                DSS0i = (AUC_median-(BTTM_AUC_median   )), # Bottom was not always reached at Clofarabine=1 --> therefore this.
                                DSS1  = (AUC_median-(BOTTOM_median*upper_DRUGCONC)) / ((TOP_median-BOTTOM_median)*upper_DRUGCONC),
                                DSS1i = (AUC_median-(BTTM_AUC_median   )) / ((TOP_median-BTTM_median )*upper_DRUGCONC),
                                DSS2  = DSS1/log10(TOP_median),
                                DSS2i = DSS1i / log10(TOP_median),
                                DSS3  = DSS2  *((upper_DRUGCONC)/(upper_DRUGCONC)),
                                DSS3i = DSS2i *((upper_DRUGCONC)/(upper_DRUGCONC))) |> 
  distinct(ID,SIM,ADDITIVECONC,.keep_all = T) 
  
  

# DS_dat <- full_join(onco_combine_therapy |>   select(ID_DAY=CELLINE_DRUG_ADDITIVE_ADDITIVECONC_REPEAT, ID = CELLINE_DRUG_ADDITIVE_ADDITIVECONC, CELLINE_DRUG_ADDITIVE,CELLINE, DRUG, ADDITIVE, ADDITIVE_CONC, DAY=REPEAT, CELLINE_NAME, DRUG_NAME, ADDITIVE_NAME, CELLINE_DRUG_ADDITIVE_NAME) |>  unique() |> 
#                       group_by(CELLINE_DRUG_ADDITIVE, DAY)   |>  mutate(CELLINE_DRUG_ADDITIVE_DAY   = cur_group_id()) |>  ungroup(), 
#                     dat_2023 |>  select(ID_DAY, TOP, BOTTOM, IC50, GAMMA, AUC,upper_DRUGCONC, BTTM, BTTM_AUC, DSS0, DSS1, DSS2, DSS3, DSS0i, DSS1i, DSS2i, DSS3i),
#                     by = c("ID_DAY"))

################################################################################

### calculate difference in DSS as compared to no OA (OA=0)#################
DSS_reference_0 <- dat_sim |>  filter(ADDITIVECONC==0) |>  ungroup() |> 
  select(SIM,
         CELLINE,
         ref_AUC=AUC_median,
         ref_BTTM=BTTM_median,
         ref_BTTM_AUC=BTTM_AUC_median,
         ref_DSS0=DSS0,
         ref_DSS1=DSS1,
         ref_DSS2=DSS2,
         ref_DSS3=DSS3,
         ref_DSS0i=DSS0i,
         ref_DSS1i=DSS1i,
         ref_DSS2i=DSS2i,
         ref_DSS3i=DSS3i)

DSS <- left_join(dat_sim, DSS_reference_0,by=c("CELLINE","SIM"))
DSS <- DSS |>  mutate(dDSS0   = DSS0   - ref_DSS0,
                      dDSS1   = DSS1   - ref_DSS1,
                      dDSS2   = DSS2   - ref_DSS2,
                      dDSS3   = DSS3   - ref_DSS3,
                      dDSS0i  = DSS0i  - ref_DSS0i,
                      dDSS1i  = DSS1i  - ref_DSS1i,
                      dDSS2i  = DSS2i  - ref_DSS2i,
                      dDSS3i  = DSS3i  - ref_DSS3i) 

### plot the dDSS dsitribution (with distribution bar) #################################################
# Plots
group.colors <- c(EGFP  ="#f8766d",
                  WT   ="#d39200", 
                  R183W ="#93aa00",
                  P179R  ="#00ba38",
                  S256F   ="#FFDB6D")
DSS_ordered <- DSS |>  filter(ADDITIVECONC!=0) 
DSS_ordered$CELLINE <- factor(DSS_ordered$CELLINE, 
                                   levels = c( "R183W","WT", "S256F", "P179R","EGFP"))
mycomparison <- list(c("R183W", "WT"),
                     c("WT"   , "S256F"),
                     c("S256F", "P179R"),
                     c("P179R", "EGFP"))
p1 <- ggboxplot(DSS_ordered, x = "CELLINE", y = "dDSS0")+
  stat_compare_means(comparisons = mycomparison,
                     method="wilcox.test", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")),
                     label.y = c(-0.09, -0.07, -0.05, -0.03,-0.01), label.x = 1)+
  labs(x = "Cell line", y = "dDSS0 (OA= 10 nM)") +
  scale_y_continuous(breaks = seq(-0.3, 0, 0.1), limits = c(-0.3, 0)) +
  theme_bw()
p1
ggsave('dDSS0_1000 simulation.png', p1 , device = "png",width = 6, height = 5.5, dpi = 500)

p2 <- ggboxplot(DSS_ordered, x = "CELLINE", y = "dDSS0i")+
  stat_compare_means(comparisons = mycomparison,
                     method="wilcox.test", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))+
  labs(x = "Cell line", y = "dDSS0i (OA= 10 nM)") +
  scale_y_continuous(breaks = seq(-0.3, 0, 0.1), limits = c(-0.3, 0.05)) +
  theme_bw()
p2
