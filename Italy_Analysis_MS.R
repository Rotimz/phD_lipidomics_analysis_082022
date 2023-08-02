#############################################
# References #
# title = {bspr-workshop-2018},
# author = {Alistair Bailey},
#year = {2019},
#note = {Unknown},
#url = {https://github.com/ab604/bspr-workshop-2018}
#############################################
## 31 August 2022 ##
## Replicate Code from r_Script_Generic_Packages_v2 ##


sessionInfo()
R.version
# click "About RStudio for RStudio version"

# Software Version Index ----------------------------------------------------------------
# R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
# RStudio version 2022.07.1+554 "Spotted Wakerobin" 
####???????# Bioconductor version 3.15 (BiocManager 1.30.18), R 4.2.1 (2022-06-23)
# Software Version Index ----------------------------------------------------------------

# Setting Up ----------------------------------------------------------------
### Cleaning Workspace, Setting Directory and Installing Necessary Packages
rm(list = ls())
getwd()
setwd("/Users/rotimi/Library/CloudStorage/OneDrive-Personal/OneDrive QMUL LIDo PhD - All 4 Years/Year 4/Quarter 14_15_Jan 2023/Data/Full_lipidomics")
getwd()

Bar_italy_data_ordered <- read.csv("bar_extract_Evaluations4PCApmol.csv")
type_sum(Bar_italy_data_ordered)

### need to group data by columns then rows ##
## need to sum columns together by lipid class ##
# Method 1 - didn't work
?lapply
base <- c("Cer.NDS.", "Cer.NP.", "Cer.NH.")
cols <- lapply(base, grep, names(Italy_data_ordered), fixed = TRUE)
for (i in seq_along(base)) {
  Italy_data_ordered[, base[i]] <- rowSums(Italy_data_ordered[, cols[[i]]])
}
# Method 2
install.packages("tidyverse")
library(tidyverse)

a<-as.data.frame(c(2000:2005))
a$Col1<-c(1:6)
a$Col2<-seq(2,12,2)

colnames(a)<-c("year","Col1","Col2")

a %>%
  mutate(Total = select(., Col1:Col2) %>% rowSums(na.rm = TRUE))

condensed_italy_data <- Italy_data_ordered %>%
  mutate(Cer.NDS = select(., Cer.NDS.32.0:Cer.NDS.48.0) %>% rowSums(na.rm = TRUE))
condensed_italy_data <- Italy_data_ordered %>%
  mutate(Cer.NP = select(., Cer.NP.36.0:Cer.NP.48.0) %>% rowSums(na.rm = TRUE))
condensed_italy_data <- Italy_data_ordered %>%
  mutate(Cer.NH = select(., Cer.NH.34.1:Cer.NH.48.1) %>% rowSums(na.rm = TRUE))

# Method 3
multiple.func <- function(x) {
  c(mutate(Cer.NDS = select(., Cer.NDS.32.0:Cer.NDS.48.0) %>% rowSums(na.rm = TRUE)),
    mutate(Cer.NP = select(., Cer.NP.36.0:Cer.NP.48.0) %>% rowSums(na.rm = TRUE)),
    mutate(Cer.NH = select(., Cer.NH.34.1:Cer.NH.48.1) %>% rowSums(na.rm = TRUE)))
}

condensed_italy_data <- sapply(Italy_data_ordered, multiple.func)

# Method 4
condensed_italy_data_NDS <- Bar_italy_data_ordered %>%
  mutate(Cer.NDS = select(., Cer.NDS.32.0:Cer.NDS.48.0) %>% rowSums(na.rm = TRUE))
condensed_italy_data_NP <- Bar_italy_data_ordered %>%
  mutate(Cer.NP = select(., Cer.NP.36.0:Cer.NP.48.0) %>% rowSums(na.rm = TRUE))
condensed_italy_data_NH <- Bar_italy_data_ordered %>%
  mutate(Cer.NH = select(., Cer.NH.34.1:Cer.NH.48.1) %>% rowSums(na.rm = TRUE))
         
Italy_data_ordered_2 <-bind_cols(condensed_italy_data_NDS,
                                 Cer.NP = condensed_italy_data_NP[,50],
                                 Cer.NH = condensed_italy_data_NH[,50])

Italy_shortened <- Italy_data_ordered_2 %>%
  select(c(Name, Cer.NDS, Cer.NP, Cer.NH)) %>%
  slice(4:21)

## need to sum rows together by condition ## --> actually don't need to do this
# Method 1

#Italy <- Italy_shortened %>%
 # group_by(Name) %>%
 # summarise(mean_Cer.68.2 = mean(Cer.68.2, na.rm = T))
###

### Graphing
plot_NDS <- ggplot(Italy_shortened) +
  geom_col(aes(x = Name, y = Cer.NDS, fill = Name)) +
             labs(x = "Condition")

plot_NP <- ggplot(Italy_shortened) +
  geom_col(aes(x = Name, y = Cer.NP, fill = Name)) +
             labs(x = "Condition")

plot_NH <- ggplot(Italy_shortened) +
  geom_col(aes(x = Name, y = Cer.NH, fill = Name)) +
             labs(x = "Condition")

plot_NDS;plot_NH;plot_NP

install.packages('patchwork')
library(patchwork)

plot_NDS + plot_NH + plot_NP +
  plot_layout(ncol = 3) + plot_layout(guides = 'collect')

### Need to Transpose Data for -omics analysis
omics_italy_data_ordered <- read.csv("omics_extract_Evaluations4PCApmol.csv")

omics_condensed_italy_data_NDS <- omics_italy_data_ordered %>%
  mutate(Cer.NDS = select(., Cer.NDS.32.0:Cer.NDS.48.0) %>% rowSums(na.rm = TRUE))
omics_condensed_italy_data_NP <- omics_italy_data_ordered %>%
  mutate(Cer.NP = select(., Cer.NP.36.0:Cer.NP.48.0) %>% rowSums(na.rm = TRUE))
omics_condensed_italy_data_NH <- omics_italy_data_ordered %>%
  mutate(Cer.NH = select(., Cer.NH.34.1:Cer.NH.48.1) %>% rowSums(na.rm = TRUE))

omics_italy_data_ordered_2 <-bind_cols(omics_condensed_italy_data_NDS,
                                 Cer.NP = omics_condensed_italy_data_NP[,50],
                                 Cer.NH = omics_condensed_italy_data_NH[,50])

omics_italy_shortened <- omics_italy_data_ordered_2 %>%
  select(c(Name, Cer.NDS, Cer.NP, Cer.NH)) %>%
  slice(4:21)

omics_transpose_italy <- data.frame(t(omics_italy_shortened[-1]))
colnames(omics_transpose_italy) <- omics_italy_shortened[, 1]



library('dplyr')
#omics_transpose_italy <- omics_transpose_italy %>%
  #mutate(id = row_number())
#data.frame(id = 0, omics_transpose_italy)
# decided to transpose in excel as I can't find a way to 


#omics_transpose_italy_ordered <- omics_transpose_italy #%>%
  #slice(match(row_order, Name))
## can't transform a data frame with duplicate names

t_test <- function(dt,grp1,grp2){
  # Subset control group and convert to numeric
  x <- dt[grp1] %>% unlist %>% as.numeric()
  # Subset treatment group and convert to numeric
  y <- dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
}

# omics_transpose_italy column 1 is not retained, will just upload an excel version
write.csv(omics_transpose_italy, "omics_transpose_italy.csv", row.names=TRUE)
omics_transpose_italy <- read.csv("omics_transpose_italy.csv")

t.test(as.numeric(omics_transpose_italy[1,11:13]),
       as.numeric(omics_transpose_italy[1,14:16]))$p.value

Italy_t.testpos_NTC_vs_sh1 <- plyr::adply(omics_transpose_italy,.margins = 1, .fun = t_test,
                                          grp1 = c(11:13), grp2 = c(14:16)) %>% as_tibble()
typeof(Italy_t.testpos_NTC_vs_sh1)
is.data.frame(Italy_t.testpos_NTC_vs_sh1)

# Plot histogram
Italy_t.testpos_NTC_vs_sh1 %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.003, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()

Italy_t.testpos_NTC_vs_sh2 <- plyr::adply(omics_transpose_italy,.margins = 1, .fun = t_test,
                                          grp1 = c(11:13), grp2 = c(17:19)) %>% as_tibble()
# Plot histogram
Italy_t.testpos_NTC_vs_sh2 %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.02, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()

# Select columns and log data
colnames(Italy_t.testpos_NTC_vs_sh1)
logtransform_Italy_pos_NTC_vs_sh1 <- Italy_t.testpos_NTC_vs_sh1 %>% 
  select(-c(Lipid_Class,neg_NTC.DOX.1,neg_NTC.DOX.2,neg_NTC.DOX.3,neg_Sh1.DOX.1,neg_Sh1.DOX.2,
            neg_Sh1.DOX.3,neg_Sh2.DOX.1,neg_Sh2.DOX.2,neg_Sh2.DOX.3,pos_Sh2.DOX.1,pos_Sh2.DOX.2,
            pos_Sh2.DOX.3,p_val)) %>% 
  log2()

logtransform_Italy_pos_NTC_vs_sh2 <- Italy_t.testpos_NTC_vs_sh2 %>% 
  select(-c(Lipid_Class,neg_NTC.DOX.1,neg_NTC.DOX.2,neg_NTC.DOX.3,neg_Sh1.DOX.1,neg_Sh1.DOX.2,
            neg_Sh1.DOX.3,neg_Sh2.DOX.1,neg_Sh2.DOX.2,neg_Sh2.DOX.3,pos_Sh1.DOX.1,pos_Sh1.DOX.2,
            pos_Sh1.DOX.3,p_val)) %>% 
  log2()

# Bind columns to create transformed data frame
combine_Italy_pos_NTC_vs_sh1 <- bind_cols(Italy_t.testpos_NTC_vs_sh1[,c(1)], logtransform_Italy_pos_NTC_vs_sh1, Italy_t.testpos_NTC_vs_sh1[,20])
combine_Italy_pos_NTC_vs_sh2 <- bind_cols(Italy_t.testpos_NTC_vs_sh2[,c(1)], logtransform_Italy_pos_NTC_vs_sh2, Italy_t.testpos_NTC_vs_sh1[,20])

# Fold Change
colnames(combine_Italy_pos_NTC_vs_sh1)
fc_Italy_pos_NTC_vs_sh1 <- combine_Italy_pos_NTC_vs_sh1 %>% 
  group_by(Lipid_Class) %>% 
  mutate(mean_NTC = mean(c(pos_NTC.DOX.1,
                           pos_NTC.DOX.2,
                           pos_NTC.DOX.3)),
         mean_knockdown= mean(c(pos_Sh1.DOX.1,
                                pos_Sh1.DOX.2,
                                pos_Sh1.DOX.3)),
         log_fc = mean_knockdown - mean_NTC,
         log_pval = -1*log10(p_val))

fc_Italy_pos_NTC_vs_sh2 <- combine_Italy_pos_NTC_vs_sh2 %>% 
  group_by(Lipid_Class) %>% 
  mutate(mean_NTC = mean(c(pos_NTC.DOX.1,
                           pos_NTC.DOX.2,
                           pos_NTC.DOX.3)),
         mean_knockdown= mean(c(pos_Sh2.DOX.1,
                                pos_Sh2.DOX.2,
                                pos_Sh2.DOX.3)),
         log_fc = mean_knockdown - mean_NTC,
         log_pval = -1*log10(p_val))


condensed_fc_Italy_pos_NTC_vs_sh1 <- fc_Italy_pos_NTC_vs_sh1 %>%
  select(Lipid_Class,log_fc, log_pval)

condensed_fc_Italy_pos_NTC_vs_sh2 <- fc_Italy_pos_NTC_vs_sh2 %>%
  select(Lipid_Class,log_fc, log_pval)

# Plot a histogram - values far above or below zero suggests knockdown effect
condensed_fc_Italy_pos_NTC_vs_sh1 %>%
  ggplot(aes(log_fc)) + 
  geom_histogram(binwidth = 0.5,
                 boundary = 0.5,
                 fill = "darkblue",
                 colour = "white") +
  xlab("log2 fold change") +
  ylab("Frequency") +
  theme_minimal()

condensed_fc_Italy_pos_NTC_vs_sh2 %>%
  ggplot(aes(log_fc)) + 
  geom_histogram(binwidth = 0.5,
                 boundary = 0.5,
                 fill = "darkblue",
                 colour = "white") +
  xlab("log2 fold change") +
  ylab("Frequency") +
  theme_minimal()

## Volcano Plot

condensed_fc_Italy_pos_NTC_vs_sh1 %>% ggplot(aes(log_fc,log_pval)) + geom_point()
condensed_fc_Italy_pos_NTC_vs_sh2 %>% ggplot(aes(log_fc,log_pval)) + geom_point()

## Thresholding
# p-value on a log10 scale means that a p-value of 0.05 or below is transformed
# to 1.3 or above and a p-value of 0.01 is equal to 2

condensed_fc_Italy_pos_NTC_vs_sh1 %>%
  # Add a threshold for significant observations
  mutate(threshold = if_else(log_fc >= 1.0 & log_pval >= 1.3 |
                               log_fc <= -1.0 & log_pval >= 1.3,"A", "B")) %>%
  # Plot with points coloured according to the threshold
  ggplot(aes(log_fc,log_pval, colour = threshold)) +
  geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 1.0, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.0, linetype = 2, alpha = 0.5) +
  # Set the colour of the points
  scale_colour_manual(values = c("A"= "red", "B"= "black")) +
  xlab("log2 fold change") + ylab("-log10 p-value") + # Relabel the axes
  theme_minimal() + # Set the theme
  theme(legend.position="none") # Hide the legend

condensed_fc_Italy_pos_NTC_vs_sh2 %>%
  # Add a threshold for significant observations
  mutate(threshold = if_else(log_fc >= 1.0 & log_pval >= 1.3 |
                               log_fc <= -1.0 & log_pval >= 1.3,"A", "B")) %>%
  # Plot with points coloured according to the threshold
  ggplot(aes(log_fc,log_pval, colour = threshold)) +
  geom_point(alpha = 0.5) + # Alpha sets the transparency of the points
  # Add dotted lines to indicate the threshold, semi-transparent
  geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
  geom_vline(xintercept = 1.0, linetype = 2, alpha = 0.5) +
  geom_vline(xintercept = -1.0, linetype = 2, alpha = 0.5) +
  # Set the colour of the points
  scale_colour_manual(values = c("A"= "red", "B"= "black")) +
  xlab("log2 fold change") + ylab("-log10 p-value") + # Relabel the axes
  theme_minimal() + # Set the theme
  theme(legend.position="none") # Hide the legend
