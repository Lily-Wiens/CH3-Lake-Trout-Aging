# Find and set the working directory to the histology folder on the desktop


# computer working directory
setwd("C:/Users/user/Desktop/CH3")


#### Checking for loaded objects ####
# list objects
ls()

# remove list of objects
rm(list = ls())

# done
#### next get the data ####

library(tidyverse)
library(magrittr)
#### Packages used ####

library(FSA)
library(nlstools)
library(plotrix)

# An R script to analyze a panel of metabolites using linear models

library(sjstats)
library(performance)
library(forestmodel)
library(emmeans)
library(pwr)

library(ggpubr)
library(rstatix)
library(FSAdata) # for data
library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(car)     # for Boot()
library(dplyr)   # for filter(), mutate()
library(ggplot2) 


SBB_OD1 <- read.csv("C:/Users/user/Desktop/Histological pics/SBB_heart_OD1.csv")
View(SBB_OD1)
# Put csv on an object to be read by R
OD_sheet <- read.csv("C:/Users/user/Desktop/CH3/OD_sheet.csv")
View(OD_sheet)
names(OD_sheet)
summary(OD_sheet)
 # lake           sample_date             id           
# Length:66          Length:66          Length:66         
# Class :character   Class :character   Class :character  
# Mode  :character   Mode  :character   Mode  :character  
# 
# 
# mass_post_resp_g  fork_length     total_length  
# Min.   : 347.0   Min.   :332.0   Min.   :369.0  
# 1st Qu.: 647.2   1st Qu.:389.5   1st Qu.:428.8  
# Median : 979.0   Median :450.0   Median :491.5  
# Mean   :1252.2   Mean   :468.4   Mean   :514.2  
# 3rd Qu.:1798.5   3rd Qu.:551.0   3rd Qu.:599.5  
# Max.   :3081.0   Max.   :659.0   Max.   :719.0  
# 
# standard_mass        sex              maturity        
# Min.   : 314.0   Length:66          Length:66         
# 1st Qu.: 578.2   Class :character   Class :character  
# Median : 860.5   Mode  :character   Mode  :character  
# Mean   :1078.8                                        
# 3rd Qu.:1526.8                                        
# Max.   :2763.0                                        
# 
# gonad_mass      ventricle_mass    liver_mass    
# Min.   :  0.407   Min.   :0.271   Min.   : 2.308  
# 1st Qu.: 16.755   1st Qu.:0.501   1st Qu.: 5.909  
# Median : 29.514   Median :0.814   Median :11.411  
# Mean   : 77.022   Mean   :1.092   Mean   :15.326  
# 3rd Qu.:105.328   3rd Qu.:1.621   3rd Qu.:20.831  
# Max.   :380.364   Max.   :3.171   Max.   :40.702  
# 
# Age           OD_heart         OD_liver     
# Min.   : 4.00   Min.   :0.7840   Min.   :0.4000  
# 1st Qu.: 7.00   1st Qu.:0.9563   1st Qu.:0.8735  
# Median :10.00   Median :0.9954   Median :0.9410  
# Mean   :11.62   Mean   :1.0147   Mean   :0.9185  
# 3rd Qu.:15.00   3rd Qu.:1.0976   3rd Qu.:0.9995  
# Max.   :26.00   Max.   :1.2930   Max.   :1.1480  
# NA's   :3                        NA's   :17      
# Calculate standard deviation for each lake, excluding NA values
lake_sds <- apply(OD_sheet, 2, sd, na.rm = TRUE)

# Print results
print(lake_sds)
# lake       sample_date                id 
# NA                NA                NA 
# mass_post_resp_g       fork_length      total_length 
# 753.8868762        91.3083980        99.0045195 
# standard_mass               sex          maturity 
# 636.3845937                NA                NA 
# gonad_mass    ventricle_mass        liver_mass 
# 90.7186834         0.7288732        10.6573775 
# Age          OD_heart          OD_liver 
# 5.5019896         0.1036847         0.1517693 
# predicted_hrtmass 
# 0.7077048 

# Check for non-numeric values in each column
sapply(OD_sheet, function(x) sum(!is.numeric(x)))

# Replace non-numeric values with NA (if needed)
OD_sheet[] <- lapply(OD_sheet, function(x) as.numeric(as.character(x)))

# Print results
print(lake_means)
print(lake_sds)

# Ensure the Lake column is treated as a factor
OD_sheet$lake <- as.factor(OD_sheet$lake)

# Calculate mean and standard deviation for mass within each lake
summary_stats_mass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_mass = mean(mass_post_resp_g, na.rm = TRUE),
    sd_mass = sd(mass_post_resp_g, na.rm = TRUE)
  )

# Print results
print(summary_stats_mass)

# Load dplyr package
library(dplyr)

# Calculate mean and standard deviation for mass within each lake, and round to 2 decimal places
summary_stats_mass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_mass = round(mean(mass_post_resp_g, na.rm = TRUE), 2),
    sd_mass = round(sd(mass_post_resp_g, na.rm = TRUE), 2)
  )

# Print results
print(summary_stats_mass)

# Load dplyr package
library(dplyr)

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$mass_post_resp_g <- as.numeric(OD_sheet$mass_post_resp_g)

# Calculate mean and standard deviation for mass within each lake, and round to 2 decimal places
summary_stats_mass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_mass = round(mean(mass_post_resp_g, na.rm = TRUE), 2),
    sd_mass = round(sd(mass_post_resp_g, na.rm = TRUE), 2)
  )

# Print results
print(summary_stats_mass)

head(OD_sheet)

# Load dplyr package
library(dplyr)

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$mass_post_resp_g <- as.numeric(OD_sheet$mass_post_resp_g)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_mass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_mass = formatC(mean(mass_post_resp_g, na.rm = TRUE), format = "f", digits = 2),
    sd_mass = formatC(sd(mass_post_resp_g, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_mass)

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$fork_length <- as.numeric(OD_sheet$fork_length)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_fl <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_fl = formatC(mean(fork_length, na.rm = TRUE), format = "f", digits = 2),
    sd_fl = formatC(sd(fork_length, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_fl)
# A tibble: 4 × 3
# lake    mean_fl sd_fl
# <fct>   <chr>   <chr>
#   1 Hogan   556.00  67.12
# 2 LOTR    402.76  35.33
# 3 Opeongo 531.12  54.11
# 4 Shirley 382.25  37.02

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$Age <- as.numeric(OD_sheet$Age)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_age <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_fl = formatC(mean(Age, na.rm = TRUE), format = "f", digits = 2),
    sd_fl = formatC(sd(Age, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_age)
# A tibble: 4 × 3
# lake    mean_fl sd_fl
# <fct>   <chr>   <chr>
#   1 Hogan   13.82   6.00 
# 2 LOTR    7.47    3.62 
# 3 Opeongo 12.07   5.76 
# 4 Shirley 12.75   4.37 

# Check the column names
names(OD_sheet)
# Check the structure of the dataframe
str(OD_sheet)

# Load dplyr package
library(dplyr)

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$Lake <- as.factor(OD_sheet$Lake)
OD_sheet$Mass <- as.numeric(OD_sheet$Mass)

# Calculate the number of individuals per lake
summary_stats_sample_size <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    sample_fl = n()
  )

# Print results
print(summary_stats_sample_size)



#### checking the column sex ####
library(dplyr)

# Remove leading and trailing spaces
OD_sheet$sex <- trimws(OD_sheet$sex)

# Convert to lowercase for consistency
OD_sheet$sex <- tolower(OD_sheet$sex)

# Check unique values
unique_values <- unique(OD_sheet$sex)
print(unique_values)

# Check the length of each entry
lengths <- unique(nchar(OD_sheet$sex))
print(lengths)

# Assuming your data frame is named df and you want to calculate the standard deviation for the variable 'var1'

# Calculate standard deviation for 'var1'
std_dev_mass_post <- sd(OD_sheet$mass_post_resp_g, na.rm = TRUE)

# Print the standard deviation
print(std_dev_mass_post)
#753.8869
# Calculate standard deviation for 'var1'
std_dev_fork_length <- sd(OD_sheet$fork_length, na.rm = TRUE)

# Print the standard deviation
print(std_dev_fork_length)
#91.3084

# Calculate standard deviation for 'var1'
std_dev_age <- sd(OD_sheet$Age, na.rm = TRUE)

# Print the standard deviation
print(std_dev_age)
# [1] 5.50199

## gonad mass

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$gonad_mass <- as.numeric(OD_sheet$gonad_mass)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_gnmass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_fl = formatC(mean(gonad_mass, na.rm = TRUE), format = "f", digits = 2),
    sd_fl = formatC(sd(gonad_mass, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_gnmass)

# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$Lake <- as.factor(OD_sheet$Lake)
OD_sheet$Mass <- as.numeric(OD_sheet$Mass)

# Calculate the number of individuals per lake
summary_stats_sample_size <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    sample_fl = n()
  )

# Print results
print(summary_stats_sample_size)

# liver mass
# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$liver_mass <- as.numeric(OD_sheet$liver_mass)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_lvmass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_fl = formatC(mean(liver_mass, na.rm = TRUE), format = "f", digits = 2),
    sd_fl = formatC(sd(liver_mass, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_lvmass)


# ventricle mass
# Ensure the Lake column is treated as a factor and Mass as numeric
OD_sheet$lake <- as.factor(OD_sheet$lake)
OD_sheet$ventricle_mass <- as.numeric(OD_sheet$ventricle_mass)

# Calculate mean and standard deviation for mass within each lake, and format to 2 decimal places
summary_stats_vmass <- OD_sheet %>%
  group_by(lake) %>%
  summarise(
    mean_fl = formatC(mean(ventricle_mass, na.rm = TRUE), format = "f", digits = 2),
    sd_fl = formatC(sd(ventricle_mass, na.rm = TRUE), format = "f", digits = 2)
  )

# Print results
print(summary_stats_vmass)

#### correlations across morfological data ####
# to test the correlation
#### Age vs FL ####
# Perform correlation test
correlation_age_forklength <- cor.test(OD_sheet$Age, OD_sheet$fork_length)

# Print the test results
print(correlation_age_forklength)
#Pearson's product-moment correlation
# 
# data:  OD_sheet$Age and OD_sheet$fork_length
# t = 4.4003, df = 61, p-value = 4.42e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.2767490 0.6585328
# sample estimates:
#   cor 
# 0.4908566 


#### age vs mass post ####
# Perform correlation test
correlation_age_masspost <- cor.test(OD_sheet$Age, 
                                     OD_sheet$mass_post_resp_g)

# Print the test results
print(correlation_age_masspost)
# Pearson's product-moment correlation
# 
# data:  OD_sheet$Age and OD_sheet$mass_post_resp_g
# t = 4.7442, df = 61, p-value = 1.303e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3114596 0.6795199
# sample estimates:
#       cor 
# 0.5191616 


#### FL vs mass post correlation ####
# Perform correlation test
correlation_fl_mass <- cor.test(OD_sheet$fork_length, 
                                OD_sheet$mass_post_resp_g)

# Print the test results
print(correlation_fl_mass)

# Pearson's product-moment correlation
# 
# data:  OD_sheet$fork_length and OD_sheet$mass_post_resp_g
# t = 32.655, df = 64, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9533625 0.9823733
# sample estimates:
#       cor 
# 0.9712778 


#### cor FL vs standard mass  ####
# Perform correlation test
correlation_fl_stmass <- cor.test(OD_sheet$fork_length, 
                                  OD_sheet$standard_mass)

# Print the test results
print(correlation_fl_stmass)
# Pearson's product-moment correlation
# 
# data:  OD_sheet$fork_length and OD_sheet$standard_mass
# t = 31.402, df = 64, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9497763 0.9809961
# sample estimates:
#       cor 
# 0.9690474 

#### FL vs mass per lake separate ####
library(dplyr)
library(purrr)
library(broom)


# Function to perform correlation test
cor_test <- function(df) {
  cor.test(df$fork_length, df$standard_mass)
}

# Apply correlation test to each lake
correlation_results <- OD_sheet %>%
  group_by(lake) %>%
  nest() %>%
  mutate(correlation = map(data, cor_test))

# Extract and print results
correlation_results %>%
  mutate(correlation_summary = map(correlation, broom::tidy)) %>%
  select(lake, correlation_summary) %>%
  unnest(correlation_summary)
# A tibble: 4 × 9
# Groups:   lake [4]
# lake    estimate statistic  p.value parameter conf.low conf.high
# <fct>      <dbl>     <dbl>    <dbl>     <int>    <dbl>     <dbl>
#   1 Opeongo    0.969     14.6  7.35e-10        14    0.910     0.989
# 2 LOTR       0.647      3.29 4.95e- 3        15    0.242     0.860
# 3 Shirley    0.978     17.4  6.88e-11        14    0.935     0.992
# 4 Hogan      0.965     14.3  3.93e-10        15    0.904     0.988
# # ℹ 2 more variables: method <chr>, alternative <chr>


#### graficos para morfologia ####
#### boxplots ####
# Create the line graph with regression lines
# Create the box plot 
ggplot(OD_sheet, aes(x = lake, y = fork_length)) + geom_boxplot()

ggplot(OD_sheet, aes(x = lake, y = fork_length, fill = lake)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("A" = "red", "B" = "blue", "C" = "green")) + 
  theme_minimal(base_size = 15) + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
# Check the unique levels in your lake variable
unique(OD_sheet$lake)
ggplot(OD_sheet, aes(x = lake, y = fork_length, fill = lake)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Lake1" = "red", "Lake2" = "blue", "Lake3" = "green")) + 
  theme_minimal(base_size = 15) + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

# Check the unique levels in your lake variable
levels(OD_sheet$lake)

# Ensure the lake variable is a factor
OD_sheet$lake <- as.factor(OD_sheet$lake)


# Create the box plot with custom colors and a white background
ggplot(OD_sheet, aes(x = lake, y = fork_length, fill = lake)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("Lake1" = "red", "Lake2" = "blue", "Lake3" = "green")) + 
  theme_minimal(base_size = 15) + 
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )
#### graphs lakes vs gonad, ventrcle, liver mass ####
ggplot(data = OD_sheet, mapping = aes(x = lake, y = liver_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Liver Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$liver_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$liver_mass) + 1, 
           label = "*", size = 10, color = "red")

liver_mass_bp <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = liver_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Liver Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$liver_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$liver_mass) + 1, 
           label = "*", size = 10, color = "red")

ggsave("liver_mass_bp.png", 
       plot = liver_mass_bp, width = 10, height = 8, dpi = 300)
print(liver_mass_bp)

ggplot(data = OD_sheet, mapping = aes(x = lake, y = gonad_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Gonad Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$gonad_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$gonad_mass) + 1, 
           label = "*", size = 10, color = "red")

gonad_mass_bp <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = gonad_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Gonad Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$gonad_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$gonad_mass) + 1, 
           label = "*", size = 10, color = "red")

ggsave("gonad_mass_bp.png", 
       plot = gonad_mass_bp, width = 10, height = 8, dpi = 300)
print(gonad_mass_bp)

ggplot(data = OD_sheet, mapping = aes(x = lake, y = ventricle_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Ventricle Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$ventricle_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$ventricle_mass) + 1, 
           label = "*", size = 10, color = "red")

ventricle_mass_bp <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = ventricle_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Ventricle Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$ventricle_mass) + 1, 
           label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$ventricle_mass) + 1, 
           label = "*", size = 10, color = "red")


ggsave("ventricle_mass_bp.png", 
       plot = ventricle_mass_bp, width = 10, height = 8, dpi = 300)
print(ventricle_mass_bp)


# y.labs <- expression("Optical Density SBB")
# x.labs <- expression("Sex")
ggplot(data = OD_sheet, mapping = aes(x = lake, y = fork_length)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Fork Length (mm)")
#### Fork Length vs Lake boxplot ####
ggplot(data = OD_sheet, mapping = aes(x = lake, y = fork_length)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Fork Length (mm)") +
  annotate("text", x = 1, y = max(OD_sheet$fork_length) + 1, label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$fork_length) + 1, label = "*", size = 10, color = "red")

lake_vs_fork_length_boxplot <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = fork_length)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Fork Length (mm)") +
  annotate("text", x = 1, y = max(OD_sheet$fork_length) + 1, label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$fork_length) + 1, label = "*", size = 10, color = "red")


ggsave("lake_vs_fork_length_boxplot.png", 
       plot = lake_vs_fork_length_boxplot, width = 10, height = 8, dpi = 300)
print(lake_vs_fork_length_boxplot)

ggplot(data = OD_sheet, mapping = aes(x = lake, 
                                      y = standard_mass)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Standard Mass (g)")

ggplot(data = OD_sheet, mapping = aes(x = lake, y = standard_mass)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Standard Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$standard_mass) + 1, label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$standard_mass) + 1, label = "*", size = 10, color = "red")

lake_vs_standar_mass_boxplot <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = standard_mass)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Standard Mass (g)") +
  annotate("text", x = 1, y = max(OD_sheet$standard_mass) + 1, label = "*", size = 10, color = "red") +
  annotate("text", x = 3, y = max(OD_sheet$standard_mass) + 1, label = "*", size = 10, color = "red")

ggsave("lake_vs_standar_mass_boxplot.png", 
       plot = lake_vs_fork_length_boxplot, width = 10, height = 8, dpi = 300)
print(lake_vs_standar_mass_boxplot)

#### age vs lake boxplot ####
ggplot(data = OD_sheet, mapping = aes(x = lake, 
                                      y = Age)) + 
  geom_boxplot(aes(fill=lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Age (years)")

ggplot(data = OD_sheet, mapping = aes(x = lake, y = Age )) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Age (years)") + annotate("text", x = 2, 
                                                 y = max(OD_sheet$Age) + 1, 
                                                 label = "*", size = 10, 
                                                 color = "red")
  
library(ggplot2)

# Ensure OD_sheet is correctly loaded with lake and Age columns
ggplot(data = OD_sheet, mapping = aes(x = lake, y = Age)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Age (years)") + 
  annotate("text", x = 2, y = max(OD_sheet$Age) + 1, label = "*", size = 10, color = "red")

# cleaning the data of nas
# Remove rows with non-finite or missing values
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(Age) & !is.na(Age) & !is.na(lake))

library(ggplot2)
library(dplyr)

# Clean the dataset
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(Age) & !is.na(Age) & !is.na(lake))

# Create the plot
ggplot(data = OD_sheet_clean, mapping = aes(x = lake, y = Age)) + 
  geom_boxplot(aes(fill = lake)) +
  theme_classic() +
  labs(x = "Lake", y = "Age (years)") + 
  annotate("text", x = 2, y = max(OD_sheet_clean$Age) + 1, label = "*", size = 10, color = "red")

#### mass vs fl per lake graphs ####

# Create the plot with no formulas
ggplot(OD_sheet, aes(x = fork_length, y = mass_post_resp_g, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Mass vs Length for Each Lake",
       x = "Fork Length (mm)",
       y = "Mass (g)",
       color = "Lake")

#### liver mass vs fl bylake ####

# Fit linear models for each lake and extract formulas
formulaslm <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(liver_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("liver_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
plm <- ggplot(OD_sheet, aes(x = fork_length, y = liver_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Liver Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Liver Mass (g)",
       color = "Lake")
plm
# Add formulas to each facet
plm + geom_text(data = formulaslm, aes(x = Inf, y = Inf, label = formulaslm), 
                hjust = 1.1, vjust = 2, size = 3, color = "black")
#### fork length age bylake ####

# Fit linear models for each lake and extract formulas
formula <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(fork_length ~ Age, data = .)) %>%
  mutate(formula = paste0("fork_length = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
alm <- ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Fork Legth vs Age",
       x = "Age (y)",
       y = "Fork Length (mm)",
       color = "Lake")
alm
# Add formulas to each facet
alm + geom_text(data = formula, aes(x = Inf, y = Inf, label = formula), 
                hjust = 1.5, vjust = 2, size = 3, color = "black")


#### trying ln for fork length ####
library(ggplot2)
library(dplyr)
library(purrr)
library(broom)


# Fit linear models for each lake and extract formulas
formulas <- OD_sheet %>%
  group_by(lake) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log(fork_length) ~ Age, data = .x))) %>%
  mutate(formula = map_chr(model, ~ paste0("fork_length = ", round(coef(.x)[1], 2), " + ", round(coef(.x)[2], 2), " * log(Age)")))

# Create the base plot
lalm <- ggplot(OD_sheet, aes(x = Age, y = log(fork_length), color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Fork Length vs Age",
       x = "Age (y)",
       y = "ln Fork Length (mm)",
       color = "Lake")

# Add formulas to each facet
lalm + geom_text(data = formulas, aes(x = Inf, y = Inf, label = formula), 
                 hjust = 1.5, vjust = 2, size = 3, color = "black")

#### organ mass vs standard mass by lake logged ####
# liver
# Fit linear models for each lake and extract formulas
formulasll <- OD_sheet %>%
  group_by(lake) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log(liver_mass) ~ log(standard_mass), data = .x))) %>%
  mutate(formula = map_chr(model, ~ paste0("liver_mass = ", round(coef(.x)[1], 2), " + ", round(coef(.x)[2], 2), " * log(standard_mass)")))

# Create the base plot
llalm <- ggplot(OD_sheet, aes(x = log(standard_mass), y = log(liver_mass), color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "ln Liver Mass vs ln Standard Mass",
       x = "ln Standard Mass (g)",
       y = "ln Liver Mass (g)",
       color = "Lake")

# Add formulas to each facet
llalm + geom_text(data = formulasll, aes(x = Inf, y = Inf, label = formula), 
                  hjust = 1.3, vjust = 2, size = 3, color = "black")

# ventricle
# Fit linear models for each lake and extract formulas
formulasvl <- OD_sheet %>%
  group_by(lake) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log(ventricle_mass) ~ log(standard_mass), data = .x))) %>%
  mutate(formula = map_chr(model, ~ paste0("ventricle_mass = ", round(coef(.x)[1], 2), " + ", round(coef(.x)[2], 2), " * log(standard_mass)")))

# Create the base plot
vlalm <- ggplot(OD_sheet, aes(x = log(standard_mass), y = log(ventricle_mass), color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "ln Ventricle Mass vs ln Standard Mass",
       x = "ln Standard Mass (g)",
       y = "ln Ventricle Mass (g)",
       color = "Lake")

# Add formulas to each facet
vlalm + geom_text(data = formulasvl, aes(x = Inf, y = Inf, label = formula), 
                  hjust = 1.3, vjust = 2, size = 3, color = "black")


# gonad
# Fit linear models for each lake and extract formulas
formulasgl <- OD_sheet %>%
  group_by(lake) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(log(gonad_mass) ~ log(standard_mass), data = .x))) %>%
  mutate(formula = map_chr(model, ~ paste0("gonad_mass = ", round(coef(.x)[1], 2), " + ", round(coef(.x)[2], 2), " * log(standard_mass)")))

# Create the base plot
glalm <- ggplot(OD_sheet, aes(x = log(standard_mass), y = log(gonad_mass), color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "ln Gonad Mass vs ln Standard Mass",
       x = "ln Standard Mass (g)",
       y = "ln Gonad Mass (g)",
       color = "Lake")

# Add formulas to each facet
glalm + geom_text(data = formulasvl, aes(x = Inf, y = Inf, label = formula), 
                  hjust = 1.3, vjust = 2, size = 3, color = "black")

#### standard mass vs fl bylake ####

# Fit linear models for each lake and extract formulas
formula <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(standard_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("standard_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
stmlm <- ggplot(OD_sheet, aes(x = fork_length, y = standard_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Standard Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Standard Mass (g)",
       color = "Lake")
stmlm
# Add formulas to each facet
stmlm + geom_text(data = formula, aes(x = Inf, y = Inf, label = formula), 
                  hjust = 1.1, vjust = 2, size = 3, color = "black")


#### mass post resp vs fl and lake ####

# Fit linear models for each lake and extract formulas
formulas <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(mass_post_resp_g ~ fork_length, data = .)) %>%
  mutate(formula = paste0("mass_post_resp_g = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
p <- ggplot(OD_sheet, aes(x = fork_length, y = mass_post_resp_g, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Mass (g)",
       color = "Lake")

# Add formulas to each facet and assign the plot to an object
final_plot <- p + geom_text(data = formulas, aes(x = Inf, y = Inf, label = formula), 
                            hjust = 1.1, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("mass_vs_fork_length_plot.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(final_plot)

#### graph log mass vs fork length ####


# Transform variables to natural log
OD_sheet <- OD_sheet %>%
  mutate(log_mass = log(mass_post_resp_g),
         log_length = log(fork_length))

# Fit linear models for each lake and extract formulas
formulas <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(log_mass ~ log_length, data = .)) %>%
  mutate(formula = paste0("log(mass) = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * log(length)"))

# Create the base plot
p <- ggplot(OD_sheet, aes(x = log_length, y = log_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Log(Mass) vs Log(Fork Length)",
       x = "Log(Fork Length) (log(mm))",
       y = "Log(Mass) (log(g))",
       color = "Lake")

# Add formulas to each facet
final_plot <- p + geom_text(data = formulas, aes(x = Inf, y = Inf, label = formula), 
                            hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("log_mass_vs_log_length_plot.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(final_plot)

#### standard mass vs fl and lake ####

# Fit linear models for each lake and extract formulas
formulasstm <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(standard_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("standard_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
pstm <- ggplot(OD_sheet, aes(x = fork_length, y = standard_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Standard Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Mass (g)",
       color = "Lake")

# Add formulas to each facet and assign the plot to an object
final_plotst <- pstm + geom_text(data = formulasstm, aes(x = Inf, y = Inf, label = formula), 
                                 hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("standard_mass_vs_fork_length_plot.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(final_plotst)

#### graph log mass vs fork length ####


# Transform variables to natural log
OD_sheet <- OD_sheet %>%
  mutate(log_standard_mass = log(standard_mass),
         log_length = log(fork_length))

# Fit linear models for each lake and extract formulas
formulaslst <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(log_standard_mass ~ log_length, data = .)) %>%
  mutate(formula = paste0("log(mass) = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * log(length)"))

# Create the base plot
plst <- ggplot(OD_sheet, aes(x = log_length, y = log_standard_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Log(Mass) vs Log(Fork Length)",
       x = "Log(Fork Length) (log(mm))",
       y = "LogStandard (Mass) (log(g))",
       color = "Lake")

# Add formulas to each facet
final_plotlst <- plst + geom_text(data = formulas, aes(x = Inf, y = Inf, label = formula), 
                                  hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("log_standard_mass_vs_log_length_plot.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(final_plotlst)

#### ventricle, liver, gonad graphs ####
#### liver mass vs fl and lake ####

# Fit linear models for each lake and extract formulas
formulaliver <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(liver_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("liver_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
pl <- ggplot(OD_sheet, aes(x = fork_length, y = liver_mass, 
                           color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Liver Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Liver Mass (g)",
       color = "Lake")

# Add formulas to each facet and assign the plot to an object
final_plotpl <- pl + geom_text(data = formulaliver, aes(x = Inf, y = Inf, label = formula), 
                               hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("liver_mass_vs_fork_length_plot.png", plot = final_plotpl, width = 10, height = 8, dpi = 300)
print(final_plotpl)

#### ventricle mass vs fl and lake ####

# Fit linear models for each lake and extract formulas
formulaventricle <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(ventricle_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("ventricle_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
pv <- ggplot(OD_sheet, aes(x = fork_length, y = ventricle_mass, 
                           color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Ventricle Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Ventricle Mass (g)",
       color = "Lake")

# Add formulas to each facet and assign the plot to an object
final_plotpv <- pv + geom_text(data = formulaventricle, aes(x = Inf, y = Inf, label = formula), 
                               hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("ventricle_mass_vs_fork_length_plot.png", plot = final_plotpv, width = 10, height = 8, dpi = 300)
print(final_plotpv)

#### gonad mass vs fl and lake ####

# Fit linear models for each lake and extract formulas
formulagonad <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(gonad_mass ~ fork_length, data = .)) %>%
  mutate(formula = paste0("gonad_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * fork_length"))

# Create the base plot
pg <- ggplot(OD_sheet, aes(x = fork_length, y = gonad_mass, 
                           color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Gonad Mass vs Fork Length",
       x = "Fork Length (mm)",
       y = "Gonad Mass (g)",
       color = "Lake")

# Add formulas to each facet and assign the plot to an object
final_plotpg <- pg + geom_text(data = formulagonad, aes(x = Inf, y = Inf, label = formula), 
                               hjust = 1.3, vjust = 2, size = 3, color = "black")

# Save the plot
ggsave("gonad_mass_vs_fork_length_plot.png", plot = final_plotpg, width = 10, height = 8, dpi = 300)
print(final_plotpg)


#### data morfologica ####

#### st mass =age + fl  linear models ####

# mass vs fl
# Fit the linear model
model_mass <- lm(standard_mass ~ fork_length, data = OD_sheet)
summary(model_mass)
# Call:
#   lm(formula = standard_mass ~ fork_length, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -659.08  -80.09  -12.81   94.26  407.13 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2084.5603   102.6058  -20.32   <2e-16 ***
#   fork_length     6.7539     0.2151   31.40   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 158.3 on 64 degrees of freedom
# Multiple R-squared:  0.9391,	Adjusted R-squared:  0.9381 
# F-statistic: 986.1 on 1 and 64 DF,  p-value: < 2.2e-16
check_model(model_mass) #
# not too bad, no collinearity
forest_model(model_mass)

# Tidy the model output
tidy_model_mass <- tidy(model_mass)
library(kableExtra)
# Create a table
tidy_model_mass %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)



# Fit the linear model
model_mass2 <- lm(standard_mass ~ fork_length + lake, data = OD_sheet)

# Print the summary of the model
summary(model_mass2)
# Call:
#   lm(formula = standard_mass ~ fork_length + lake, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -640.36  -76.87  -18.07   98.98  422.65 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2201.517    224.535  -9.805 3.77e-14 ***
#   fork_length     6.892      0.398  17.319  < 2e-16 ***
#   lakeLOTR       31.225     81.456   0.383    0.703    
# lakeOpeongo    88.657     55.727   1.591    0.117    
# lakeShirley    93.651     88.253   1.061    0.293    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 157.4 on 61 degrees of freedom
# Multiple R-squared:  0.9426,	Adjusted R-squared:  0.9388 
# F-statistic: 250.2 on 4 and 61 DF,  p-value: < 2.2e-16


# F-statistic:   509 on 2 and 58 DF,  p-value: < 2.2e-16
check_model(model_mass2) #
# not too bad, no collinearity
forest_model(model_mass2)

# Tidy the model output
tidy_model_mass2 <- tidy(model_mass2)

# Create a table
tidy_model_mass2 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Fit the linear model
model_mass3 <- lm(standard_mass ~ fork_length + Age + lake, 
                  data = OD_sheet)

# Print the summary of the model
summary(model_mass3)
# Call:
#   lm(formula = standard_mass ~ fork_length + Age + lake, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -602.94  -85.89  -19.46   97.57  377.11 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2095.7752   248.8909  -8.420 1.38e-11 ***
#   fork_length     6.5139     0.5054  12.889  < 2e-16 ***
#   Age             7.5592     5.0647   1.493   0.1411    
# lakeLOTR       26.1014    86.2055   0.303   0.7632    
# lakeOpeongo    98.6557    57.1914   1.725   0.0899 .  
# lakeShirley    36.0666   101.0741   0.357   0.7225    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 158.4 on 57 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9444,	Adjusted R-squared:  0.9395 
# F-statistic: 193.5 on 5 and 57 DF,  p-value: < 2.2e-16
check_model(model_mass3) #
# not too bad, no collinearity
forest_model(model_mass3)
# Tidy the model output
tidy_model_mass3 <- tidy(model_mass3)

# Create a table
tidy_model_mass3 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-2095.775202	248.8909292	-8.4204563	0.0000000
# fork_length	6.513941	0.5053723	12.8893895	0.0000000
# Age	7.559203	5.0647079	1.4925249	0.1410762
# lakeLOTR	26.101433	86.2055347	0.3027814	0.7631587
# lakeOpeongo	98.655742	57.1914334	1.7250091	0.0899447
# lakeShirley	36.066603	101.0740876	0.3568333	0.7225342


# Fit the linear model
model_mass4 <- lm(standard_mass ~ fork_length + Age, 
                  data = OD_sheet)

# Print the summary of the model
summary(model_mass4)
# Call:
#   lm(formula = standard_mass ~ fork_length + Age, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -611.72  -84.03   -8.68   76.55  409.05 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2094.0860   106.9003 -19.589   <2e-16 ***
#   fork_length     6.5803     0.2489  26.438   <2e-16 ***
#   Age             7.7299     4.2419   1.822   0.0736 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 159.8 on 58 degrees of freedom
# (5 observations deleted due to missingness)
# Multiple R-squared:  0.9417,	Adjusted R-squared:  0.9397 
# F-statistic: 468.8 on 2 and 58 DF,  p-value: < 2.2e-16
check_model(model_mass4)
forest_model(model_mass4)


# Tidy the model output
tidy_model_mass4 <- tidy(model_mass4)

# Create a table
tidy_model_mass4 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)



# Fit the linear model
model_mass5 <- lm(standard_mass ~ 
                    fork_length + Age + lake + sex, data = OD_sheet)

# Print the summary of the model
summary(model_mass5)
# Call:
#   lm(formula = standard_mass ~ fork_length + Age + lake + sex, 
#      data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -585.76  -87.18    0.00   90.72  346.62 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2043.3413   275.7789  -7.409 9.90e-10 ***
#   fork_length     6.3873     0.5458  11.703 2.59e-16 ***
#   Age             8.4794     5.0693   1.673    0.100    
# lakeLOTR      -26.2480   101.4508  -0.259    0.797    
# lakeOpeongo    82.2246    62.2955   1.320    0.193    
# lakeShirley    -8.6654   118.9003  -0.073    0.942    
# sexmale        38.6184    47.0920   0.820    0.416    
# sexunknown    -26.7202   177.6196  -0.150    0.881    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 161 on 53 degrees of freedom
# (5 observations deleted due to missingness)
# Multiple R-squared:  0.946,	Adjusted R-squared:  0.9389 
# F-statistic: 132.6 on 7 and 53 DF,  p-value: < 2.2e-16
check_model(model_mass5) 
forest_model(model_mass5)
# Tidy the model output
tidy_model_mass5 <- tidy(model_mass5)

# Create a table
tidy_model_mass5 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

#### checking the column sex ####
library(dplyr)

# Remove leading and trailing spaces
OD_sheet$sex <- trimws(OD_sheet$sex)

# Convert to lowercase for consistency
OD_sheet$sex <- tolower(OD_sheet$sex)

# Check unique values
unique_values <- unique(OD_sheet$sex)
print(unique_values)

# Check the length of each entry
lengths <- unique(nchar(OD_sheet$sex))
print(lengths)


# Fit the linear model
model_growth4 <- lm(log(mass_post_resp_g) ~ 
                      log(fork_length), data = OD_sheet)

# Print the summary of the model
summary(model_growth4)
# lm(formula = log(mass_post_resp_g) ~ log(fork_length), data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.74162 -0.04455  0.00641  0.07507  0.17206 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -11.89200    0.48452  -24.54   <2e-16 ***
#   log(fork_length)   3.07404    0.07899   38.92   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1229 on 64 degrees of freedom
# Multiple R-squared:  0.9595,	Adjusted R-squared:  0.9588 
# F-statistic:  1514 on 1 and 64 DF,  p-value: < 2.2e-16
# growth by age

check_model(model_growth4) # good

# Tidy the model output
tidy_model_growth4 <- tidy(model_growth4)

# Create a table
tidy_model_growth4 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Fit the linear model
model_growth5 <- lm(fork_length ~ 
                      Age, data = OD_sheet)

summary(model_growth5)
# Call:
#   lm(formula = fork_length ~ Age, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -179.04  -55.85   21.58   63.38  147.95 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   372.44      24.01   15.51  < 2e-16 ***
#   Age             8.23       1.87    4.40 4.42e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 81.03 on 61 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.2409,	Adjusted R-squared:  0.2285 
# F-statistic: 19.36 on 1 and 61 DF,  p-value: 4.42e-05
check_model(model_growth5)

# Tidy the model output
tidy_model_growth5 <- tidy(model_growth5)

# Create a table
tidy_model_growth5 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	372.438795	24.009783	15.511960	0.00e+00
# Age	8.229994	1.870329	4.400292	4.42e-05

library(kableExtra)
# Fit the linear model
model_growth6 <- lm(fork_length ~ 
                      log(Age), data = OD_sheet)

summary(model_growth6)
# 
# Call:
#   lm(formula = fork_length ~ log(Age), data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -169.09  -58.92   28.46   62.87  145.65 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   256.75      51.94   4.943 6.33e-06 ***
#   log(Age)       90.24      21.74   4.152 0.000104 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 82.12 on 61 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.2203,	Adjusted R-squared:  0.2075 
# F-statistic: 17.24 on 1 and 61 DF,  p-value: 0.0001042

check_model(model_growth6)

# Tidy the model output
tidy_model_growth6 <- tidy(model_growth6)

# Create a table
tidy_model_growth6 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	256.7522	51.94049	4.943199	0.0000063
# log(Age)	90.2411	21.73681	4.151535	0.0001042

#### FL vs mass age ####
model_growth7 <- lm(log(fork_length)~ Age + log(mass_post_resp_g) + lake, data = OD_sheet)

# Print the summary of the model
summary(model_growth7)
# Call:
#   lm(formula = log(fork_length) ~ Age + log(mass_post_resp_g) + 
#        lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.046486 -0.020078 -0.002878  0.011036  0.231822 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            4.067496   0.134454  30.252   <2e-16 ***
#   Age                    0.002122   0.001142   1.858   0.0683 .  
# log(mass_post_resp_g)  0.295669   0.018881  15.660   <2e-16 ***
#   lakeLOTR              -0.013304   0.021455  -0.620   0.5377    
# lakeOpeongo           -0.027288   0.013604  -2.006   0.0496 *  
#   lakeShirley           -0.028683   0.024982  -1.148   0.2557    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03809 on 57 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9649,	Adjusted R-squared:  0.9619 
# F-statistic: 313.7 on 5 and 57 DF,  p-value: < 2.2e-16
check_model(model_growth7)
forest_model(model_growth7)
# Tidy the model output
tidy_model_growth7 <- tidy(model_growth7)

# Create a table
tidy_model_growth7 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# OD_sheet and contains the columns 'fork_length' and 'age'

library(ggplot2)

ggplot(data = OD_sheet, aes(x = Age, y = fork_length)) + 
  geom_point(color = "blue", size = 3) + 
  theme_classic() +
  labs(title = "Scatterplot of Fork Length vs. Age",
       x = "Age (years)",
       y = "Fork Length (mm)")


                              
                              
#### mass post =age + fl ####
# Fit the linear model
model_growth <- lm(mass_post_resp_g ~ fork_length, data = OD_sheet)
# Print the summary of the model
summary(model_growth)
# Call:
#   lm(formula = mass_post_resp_g ~ fork_length, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -797.49 -107.95   21.59  119.09  368.07 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2503.8945   117.1559  -21.37   <2e-16 ***
#   fork_length     8.0193     0.2456   32.66   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 180.8 on 64 degrees of freedom
# Multiple R-squared:  0.9434,	Adjusted R-squared:  0.9425 
# F-statistic:  1066 on 1 and 64 DF,  p-value: < 2.2e-16 
check_model(model_growth) #
# not too bad, no collinearity
forest_model(model_growth)
                            
# Tidy the model output
tidy_model_growth <- tidy(model_growth)
                            
# Create a table
tidy_model_growth %>%
      kable("html", caption = "Linear Regression Results") %>%
      kable_styling(full_width = FALSE)
                            
# ln fl vs mass
lnmodel_growth <- lm(mass_post_resp_g ~ fork_length, data = OD_sheet)
# Print the summary of the model
summary(lnmodel_growth)
# Call:
#   lm(formula = mass_post_resp_g ~ fork_length, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -797.49 -107.95   21.59  119.09  368.07 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2503.8945   117.1559  -21.37   <2e-16 ***
#   fork_length     8.0193     0.2456   32.66   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 180.8 on 64 degrees of freedom
# Multiple R-squared:  0.9434,	Adjusted R-squared:  0.9425 
# F-statistic:  1066 on 1 and 64 DF,  p-value: < 2.2e-16

 
# F-statistic:  1066 on 1 and 64 DF,  p-value: < 2.2e-16 
check_model(lnmodel_growth) #
# not too bad, no collinearity
forest_model(lnmodel_growth)

# Tidy the model output
tidy_model_growthln <- tidy(lnmodel_growth)

# Create a table
tidy_model_growthln %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# 
# 
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-2503.894467	117.1559051	-21.37233	0
# fork_length	8.019346	0.2455772	32.65510	0
                         
#### mass vs fl and lake ####                            
# Fit the linear model
model_growth_mfl <- lm(mass_post_resp_g ~ fork_length + lake, data = OD_sheet)
                            
        # Print the summary of the model
summary(model_growth_mfl)
# Call:
#   lm(formula = mass_post_resp_g ~ fork_length + lake, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -779.76 -110.46   27.48  111.86  338.20 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2634.1251   258.8020 -10.178 9.04e-15 ***
#   fork_length     8.1909     0.4587  17.857  < 2e-16 ***
#   lakeLOTR       29.3072    93.8874   0.312    0.756    
# lakeOpeongo    75.7480    64.2322   1.179    0.243    
# lakeShirley    98.9141   101.7215   0.972    0.335    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 181.5 on 61 degrees of freedom
# Multiple R-squared:  0.9456,	Adjusted R-squared:  0.9421 
# F-statistic: 265.2 on 4 and 61 DF,  p-value: < 2.2e-16  

check_model(model_growth_mfl)
    # not too bad, no collinearity
forest_model(model_growth_mfl)
                            
    # Tidy the model output
tidy_model_growth_mfl <- tidy(model_growth_mfl)
                            
  # Create a table
tidy_model_growth_mfl %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)

# 
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-2634.125108	258.8019863	-10.1781487	0.0000000
# fork_length	8.190872	0.4586903	17.8570883	0.0000000
# lakeLOTR	29.307227	93.8874393	0.3121528	0.7559889
# lakeOpeongo	75.747953	64.2322211	1.1792828	0.2428647
# lakeShirley	98.914096	101.7215451	0.9724006	0.3346915

# Create a dataframe with the model predictions
predictions <- OD_sheet %>%
  mutate(predicted_mass = predict(model_growth_mfl, newdata = .))

# Plot the scatterplot with regression lines
ggplot(OD_sheet, aes(x = fork_length, y = mass_post_resp_g, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = fork_length, y = predicted_mass, color = lake)) +
  labs(title = "ANCOVA: Mass vs Length by Lake",
       x = "Fork Length (mm)",
       y = "Mass (g)") +
  theme_minimal()



#### mass vs fl and lake ln ####
lnmodel_growth_mfl <- lm(log(mass_post_resp_g) ~ log(fork_length) + lake, data = OD_sheet)

# Print the summary of the model
summary(lnmodel_growth_mfl)
# 
# Call:
#   lm(formula = log(mass_post_resp_g) ~ log(fork_length) + lake, 
#      data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.66860 -0.06010  0.02938  0.06854  0.14254 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -9.95213    0.90180 -11.036 3.59e-16 ***
#   log(fork_length)  2.76382    0.14276  19.359  < 2e-16 ***
#   lakeLOTR         -0.09440    0.06038  -1.563    0.123    
# lakeOpeongo       0.06913    0.04076   1.696    0.095 .  
# lakeShirley      -0.12524    0.06664  -1.879    0.065 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1157 on 61 degrees of freedom
# Multiple R-squared:  0.9658,	Adjusted R-squared:  0.9635 
# F-statistic: 430.5 on 4 and 61 DF,  p-value: < 2.2e-16
 

check_model(lnmodel_growth_mfl)
# not too bad, no collinearity
forest_model(lnmodel_growth_mfl)

# Tidy the model output
tidy_lnmodel_growth_mfl <- tidy(lnmodel_growth_mfl)

# Create a table
tidy_lnmodel_growth_mfl %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-9.9521263	0.9017954	-11.035903	0.0000000
# log(fork_length)	2.7638160	0.1427632	19.359447	0.0000000
# lakeLOTR	-0.0943953	0.0603759	-1.563461	0.1231176
# lakeOpeongo	0.0691283	0.0407601	1.695980	0.0949883
# lakeShirley	-0.1252400	0.0666448	-1.879215	0.0649944

# Create a dataframe with the model predictions
predictionsln <- OD_sheet %>%
  mutate(predicted_lnmass = predict(lnmodel_growth_mfl, newdata = .))
# plot
ggplot(OD_sheet, aes(x = fork_length, y = mass_post_resp_g, color = lake)) +
  geom_point() +
  geom_line(data = predictionsln, aes(x = fork_length, y = predicted_lnmass, color = lake)) +
  labs(title = "ANCOVA: Mass vs Length by Lake",
       x = "Fork Length (mm)",
       y = "Mass (g)") +
  theme_minimal()

#### mass vs fl and lake and sex ln ####
lnmodel_growth_mfls <- lm(log(mass_post_resp_g) ~ log(fork_length) + 
                            lake + sex, data = OD_sheet)

# Print the summary of the model
summary(lnmodel_growth_mfls)
# Call:
#   lm(formula = log(mass_post_resp_g) ~ log(fork_length) + lake + 
#        sex, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.64803 -0.05354  0.01801  0.07208  0.15282 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -9.85631    0.96389 -10.226 1.11e-14 ***
#   log(fork_length)  2.75088    0.15234  18.058  < 2e-16 ***
#   lakeLOTR         -0.08679    0.06754  -1.285   0.2038    
# lakeOpeongo       0.07641    0.04313   1.772   0.0816 .  
# lakeShirley      -0.11120    0.07600  -1.463   0.1487    
# sexmale          -0.04402    0.03157  -1.394   0.1685    
# sexunknown       -0.10872    0.12585  -0.864   0.3911    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.115 on 59 degrees of freedom
# Multiple R-squared:  0.9673,	Adjusted R-squared:  0.964 
# F-statistic: 290.7 on 6 and 59 DF,  p-value: < 2.2e-16

check_model(lnmodel_growth_mfls)
# not too bad, no collinearity
forest_model(lnmodel_growth_mfls)

# Tidy the model output
tidy_lnmodel_growth_mfls <- tidy(lnmodel_growth_mfls)

# Create a table
tidy_lnmodel_growth_mfls %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)



# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-9.8563115	0.9638949	-10.2255048	0.0000000
# log(fork_length)	2.7508835	0.1523359	18.0580059	0.0000000
# lakeLOTR	-0.0867882	0.0675438	-1.2849170	0.2038429
# lakeOpeongo	0.0764114	0.0431320	1.7715712	0.0816300
# lakeShirley	-0.1111985	0.0759956	-1.4632232	0.1487136
# sexmale	-0.0440188	0.0315718	-1.3942443	0.1684723
# sexunknown	-0.1087168	0.1258455	-0.8638915	0.3911461

predictions_sex <- OD_sheet %>%
  mutate(predicted_lnmass = predict(lnmodel_growth_mfls, newdata = .),
         predicted_mass = exp(predicted_lnmass))  # Back-transform to original scale

ggplot(OD_sheet, aes(x = fork_length, y = mass_post_resp_g, color = lake, shape = sex)) +
  geom_point() +
  geom_line(data = predictions_sex, aes(x = fork_length, y = predicted_mass, color = lake, linetype = sex)) +
  labs(title = "ANCOVA: Mass vs Length by Lake and Sex",
       x = "Fork Length (mm)",
       y = "Mass (g)") +
  theme_minimal()


#### mass vs fl and lake and sex and age ln ####
lnmodel_growth_mflsa <- lm(log(mass_post_resp_g) ~ log(fork_length) + 
                            lake + sex + Age, data = OD_sheet)

# Print the summary of the model
summary(lnmodel_growth_mflsa)
# Call:
#   lm(formula = log(mass_post_resp_g) ~ log(fork_length) + lake + 
#        sex + Age, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.63121 -0.05170  0.02330  0.06946  0.14764 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -9.6276398  1.1562951  -8.326 2.57e-11 ***
#   log(fork_length)  2.7137366  0.1868999  14.520  < 2e-16 ***
#   lakeLOTR         -0.1193507  0.0733170  -1.628   0.1093    
# lakeOpeongo       0.0747834  0.0445040   1.680   0.0986 .  
# lakeShirley      -0.1338490  0.0855755  -1.564   0.1235    
# sexmale          -0.0287496  0.0334677  -0.859   0.3941    
# sexunknown       -0.1146959  0.1279424  -0.896   0.3739    
# Age               0.0002545  0.0036014   0.071   0.9439    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1165 on 55 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9681,	Adjusted R-squared:  0.964 
# F-statistic: 238.4 on 7 and 55 DF,  p-value: < 2.2e-16


check_model(lnmodel_growth_mflsa)
# not too bad, no collinearity
forest_model(lnmodel_growth_mflsa)

# Tidy the model output
tidy_lnmodel_growth_mflsa <- tidy(lnmodel_growth_mflsa)

# Create a table
tidy_lnmodel_growth_mflsa %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-9.6276398	1.1562951	-8.3262825	0.0000000
# log(fork_length)	2.7137366	0.1868999	14.5197325	0.0000000
# lakeLOTR	-0.1193507	0.0733170	-1.6278712	0.1092680
# lakeOpeongo	0.0747834	0.0445040	1.6803742	0.0985563
# lakeShirley	-0.1338490	0.0855755	-1.5641029	0.1235295
# sexmale	-0.0287496	0.0334677	-0.8590239	0.3940559
# sexunknown	-0.1146959	0.1279424	-0.8964656	0.3739113
# Age	0.0002545	0.0036014	0.0706777	0.9439106



#### Model of age w/ morfological ####

# Fit the linear model
model_growth7 <- lm(fork_length~ log(Age) + 
                      mass_post_resp_g + lake, data = OD_sheet)
                            
# Print the summary of the model
summary(model_growth7)
# Call:
#   lm(formula = fork_length ~ log(Age) + mass_post_resp_g + lake, 
#      data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -26.268 -13.416  -3.255   7.075  93.988 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      347.740532  16.826219  20.667  < 2e-16 ***
#   log(Age)           8.382300   7.679126   1.092 0.279615    
# mass_post_resp_g   0.097381   0.007129  13.660  < 2e-16 ***
#   lakeLOTR         -29.007474  10.293167  -2.818 0.006629 ** 
#   lakeOpeongo      -12.718486   7.297851  -1.743 0.086766 .  
# lakeShirley      -44.385067  11.644432  -3.812 0.000341 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 20.39 on 57 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.9551,	Adjusted R-squared:  0.9512 
# F-statistic: 242.5 on 5 and 57 DF,  p-value: < 2.2e-16                            
                            

check_model(model_growth7)
forest_model(model_growth7)
                            # Tidy the model output
tidy_model_growth7 <- tidy(model_growth7)
                            
                            # Create a table
tidy_model_growth7 %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)

#### fl vs age by lake ####
# Fit the linear model
model_growth8 <- lm(fork_length~ log(Age) * lake, data = OD_sheet)

# Print the summary of the model
summary(model_growth8)

# Call:
#   lm(formula = fork_length ~ log(Age) * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -115.343  -20.115   -5.607   23.119   87.552 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)          290.4843    58.1269   4.997
# log(Age)             104.5444    22.5846   4.629
# lakeLOTR              83.0203    76.3669   1.087
# lakeOpeongo           -0.1653    76.8316  -0.002
# lakeShirley           36.0379    93.6729   0.385
# log(Age):lakeLOTR    -91.1812    33.8583  -2.693
# log(Age):lakeOpeongo  -4.1984    30.6465  -0.137
# log(Age):lakeShirley -82.1725    36.9390  -2.225
# Pr(>|t|)    
# (Intercept)          6.25e-06 ***
#   log(Age)             2.28e-05 ***
#   lakeLOTR              0.28172    
# lakeOpeongo           0.99829    
# lakeShirley           0.70193    
# log(Age):lakeLOTR     0.00937 ** 
#   log(Age):lakeOpeongo  0.89154    
# log(Age):lakeShirley  0.03023 *  
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 38.82 on 55 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.8429,	Adjusted R-squared:  0.8229 
# F-statistic: 42.15 on 7 and 55 DF,  p-value: < 2.2e-16

                            
#### scatterplots of age vs length and stuff by lake ####
                            
library(ggplot2)
                            
ggplot(data = OD_sheet, aes(x = Age, y = fork_length)) + 
        geom_point(color = "blue", size = 3) + 
        theme_classic() +
        labs(title = "Scatterplot of Fork Length vs. Age",
        x = "Age (years)",
        y = "Fork Length (mm)")
                            
                            
ggplot(data = OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age by Lake",
          x = "Age (years)",
          y = "Fork Length (mm)")
                            
                            
                            
ggplot(data = OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
        geom_point(size = 3) + 
        geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "black") +  # Logarithmic curve
        theme_classic() +
        labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Curve",
        x = "Age (years)",
        y = "Fork Length (mm)",
        color = "Locality")

# Fit the model
model <- lm(fork_length ~ log(Age), data = OD_sheet)

# Extract the coefficients
coefficients <- coef(model)
intercept <- coefficients[1]
slope <- coefficients[2]

# Display the formula
cat("Fork Length (mm) =", intercept, "+", slope, "* log(Age)")



ggplot(data = OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "black") +  # Logarithmic curve
  theme_classic() +
  labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Curve",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Locality") +
  annotate("text", x = max(OD_sheet$Age), y = max(OD_sheet$fork_length), 
           label = paste("Fork Length (mm) =", round(intercept, 2), "+", round(slope, 2), "* log(Age)"), 
           hjust = 1, vjust = 1, color = "black")

ggplot(data = OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "black") +  # Logarithmic curve
  theme_classic() +
  labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Curve",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Locality") +
  annotate("text", x = max(OD_sheet$Age), y = max(OD_sheet$fork_length), 
           label = paste("Fork Length (mm) =", round(intercept, 2), "+", round(slope, 2), "* log(Age)"), 
           hjust = 1, vjust = 1, color = "black") +
 # Display the formula
cat("Fork Length (mm) =", intercept, "+", slope, "* log(Age)")

library(dplyr)

# Remove rows with non-finite or missing values
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(Age) & !is.na(Age) & is.finite(fork_length) & !is.na(fork_length))

# Fit the model
lnmodel_growth_mfls <- lm(log(mass_post_resp_g) ~ log(fork_length) + lake + sex, data = OD_sheet_clean)

# Extract the coefficients
coefficients <- coef(lnmodel_growth_mfls)
intercept <- coefficients[1]
slope <- coefficients[2]

# Display the formula
cat("Fork Length (mm) =", intercept, "+", slope, "* log(Age)")

predictions_sex <- OD_sheet_clean %>%
  mutate(predicted_lnmass = predict(lnmodel_growth_mfls, newdata = .),
         predicted_mass = exp(predicted_lnmass))  # Back-transform to original scale
library(ggplot2)

ggplot(OD_sheet_clean, aes(x = Age, y = fork_length, color = lake, shape = sex)) +
  geom_point(size = 3) +
  geom_line(data = predictions_sex, aes(x = Age, y = predicted_mass, color = lake, linetype = sex)) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "black") +
  theme_classic() +
  labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Curve",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Locality") +
  annotate("text", x = max(OD_sheet_clean$Age), y = max(OD_sheet_clean$fork_length), 
           label = paste("Fork Length (mm) =", round(intercept, 2), "+", round(slope, 2), "* log(Age)"), 
           hjust = 1, vjust = 1, color = "black")

                            
# Example data for fitting (replace with your actual data)
fish_data <- OD_sheet
                            
# Fit the von Bertalanffy growth function
vbgf_model <- nls(fork_length ~ L_inf * (1 - exp(-k * (Age - t_0))), 
    data = fish_data, 
    start = list(L_inf = 7000, k = 0.2, t_0 = -1))
                            
# Add predicted values to the data frame
fish_data$predicted_fork_length <- predict(vbgf_model)
                            
# Plot with the von Bertalanffy growth curve
ggplot(fish_data, aes(x = age, y = fork_length, color = locality)) + 
      geom_point(size = 3) + 
      geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
      theme_classic() +
      labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
      x = "Age (years)",
      y = "Fork Length (mm)",
      color = "Locality")
                            
                            
# another try
                            
# Example data for fitting (replace with your actual data)
fish_data <- OD_sheet
                            
# Initial plotting to understand the data
plot(fish_data$Age, fish_data$fork_length)
                            
# Fit the von Bertalanffy growth function with better starting values
vbgf_model <- nls(fork_length ~ L_inf * (1 - exp(-k * (Age - t_0))), 
              data = fish_data, 
              start = list(L_inf = max(fish_data$fork_length), k = 0.1, t_0 = 0))
                            
# Add predicted values to the data frame
fish_data$predicted_fork_length <- predict(vbgf_model)
                            
# Plot with the von Bertalanffy growth curve
ggplot(fish_data, aes(x = age, y = fork_length, color = locality)) + 
      geom_point(size = 3) + 
      geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
      theme_classic() +
      labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
      x = "Age (years)",
      y = "Fork Length (mm)",
      color = "Locality")
                            

                            
# Scale your data
fish_data <- OD_sheet
fish_data$scaled_age <- scale(fish_data$Age)
fish_data$scaled_fork_length <- scale(fish_data$fork_length)
                            
# Fit the von Bertalanffy growth function with scaled data and higher evaluation limits
vbgf_model <- nls(scaled_fork_length ~ L_inf * (1 - exp(-k * (scaled_age - t_0))), 
                  data = fish_data, 
                  start = list(L_inf = 1, k = 0.1, t_0 = 0),
                  algorithm = "port",
                  control = nls.control(maxiter = 500, maxeval = 1000)
# Add predicted values to the data frame and scale back
fish_data$predicted_scaled_fork_length <- predict(vbgf_model)
fish_data$predicted_fork_length <- (fish_data$predicted_scaled_fork_length * sd(fish_data$fork_length)) + mean(fish_data$fork_length)
                            
 # Plot with the von Bertalanffy growth curve
 ggplot(fish_data, aes(x = age, y = fork_length, color = locality)) + 
                              geom_point(size = 3) + 
                              geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
                              theme_classic() +
                              labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Locality")
                            
## one more time
# Fit the von Bertalanffy growth function with increased iterations
vbgf_model <- nls(scaled_fork_length ~ L_inf * (1 - exp(-k * (scaled_age - t_0))), 
                                              data = fish_data, 
                                              start = list(L_inf = 1, k = 0.1, t_0 = 0),
                                              algorithm = "port",
                                              control = nls.control(maxiter = 500))
                            
# Add predicted values to the data frame and scale back
fish_data$predicted_scaled_fork_length <- predict(vbgf_model)
fish_data$predicted_fork_length <- (fish_data$predicted_scaled_fork_length * sd(fish_data$fork_length)) + mean(fish_data$fork_length)
                            
# Plot with the von Bertalanffy growth curve
ggplot(fish_data, aes(x = Age, y = fork_length, color = locality)) + 
                              geom_point(size = 3) + 
                              geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
                              theme_classic() +
                              labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Locality")
                            
                            
                            
#### trying another package on VB ####
                            
install.packages("minpack.lm")
library(minpack.lm)
# try again :)
                            
# Fit the von Bertalanffy growth function using nlsLM from minpack.lm
vbgf_model <- nlsLM(scaled_fork_length ~ L_inf * (1 - exp(-k * (scaled_age - t_0))), 
                  data = fish_data, 
                  start = list(L_inf = 1, k = 0.1, t_0 = 0),
                  control = nls.lm.control(maxiter = 500))
                            
# Add predicted values to the data frame and scale back
fish_data$predicted_scaled_fork_length <- predict(vbgf_model)
fish_data$predicted_fork_length <- (fish_data$predicted_scaled_fork_length * sd(fish_data$fork_length)) + mean(fish_data$fork_length)
                            
# Plot with the von Bertalanffy growth curve
ggplot(fish_data, aes(x = age, y = fork_length, color = locality)) + 
                              geom_point(size = 3) + 
                              geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
                              theme_classic() +
                              labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Locality")
                            
#### getting rid of missing values for VB ####
# Identify rows with missing values
# Identify rows with missing values
missing_values <- fish_data[is.na(fish_data$fork_length) | is.na(fish_data$age), ]
                            print(missing_values)
                            
#omit them
# Remove rows with missing values
fish_data <- na.omit(fish_data)
# fit again or try
# Fit the von Bertalanffy growth function using nlsLM from minpack.lm
vbgf_model <- nlsLM(scaled_fork_length ~ L_inf * (1 - exp(-k * (scaled_age - t_0))), 
                                                data = fish_data, 
                                                start = list(L_inf = 1, k = 0.1, t_0 = 0),
                                                control = nls.lm.control(maxiter = 500))
                            
# Add predicted values to the data frame and scale back
fish_data$predicted_scaled_fork_length <- predict(vbgf_model)
fish_data$predicted_fork_length <- (fish_data$predicted_scaled_fork_length * sd(fish_data$fork_length)) + mean(fish_data$fork_length)
                            
                            
                            
ggplot(fish_data, aes(x = Age, y = fork_length, color = lake)) + 
                              geom_point(size = 3) + 
                              geom_line(aes(y = predicted_fork_length), color = "black") +  # VBGF curve
                              theme_classic() +
                              labs(title = "Scatterplot of Fork Length vs. Age with VBGF Curve",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Locality")
                            
                            
#### fork length as y just trying simple model ####
#### fl=age+lake ####
# Fit a linear model
lm_model_fl_age_lake <- lm(fork_length ~ Age + lake, data = fish_data)
# Fit a linear model
lm_model_fl_age_lake <- lm(fork_length ~ Age + lake, data = OD_sheet)
                            
# Print the summary of the model
summary(lm_model_fl_age_lake)
# Call:
#   lm(formula = fork_length ~ Age + lake, data = fish_data)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -95.388 -19.770  -3.321  24.511  76.820 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  452.612     21.631  20.924  < 2e-16 ***
#   Age            6.736      1.251   5.383 3.05e-06 ***
#   lakeLOTR    -114.964     18.912  -6.079 3.06e-07 ***
#   lakeOpeongo   -3.859     17.017  -0.227    0.822    
# lakeShirley -152.642     16.207  -9.418 6.51e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 39.95 on 42 degrees of freedom
# Multiple R-squared:  0.8431,	Adjusted R-squared:  0.8281 
# F-statistic: 56.41 on 4 and 42 DF,  p-value: 2.405e-16
                            
check_model(lm_model_fl_age_lake)
forest_model(lm_model_fl_age_lake)

report(lm_model_fl_age_lake)
# We fitted a linear model (estimated using OLS) to predict
# fork_length with Age and lake (formula: fork_length ~ Age +
#                                  lake). The model explains a statistically significant and
# substantial proportion of variance (R2 = 0.81, F(4, 58) = 63.39,
#                                     p < .001, adj. R2 = 0.80). The model's intercept, corresponding
# to Age = 0 and lake = Hogan, is at 473.66 (95% CI [438.21,
# 509.11], t(58) = 26.75, p < .001). Within this model:
# 
#   - The effect of Age is statistically significant and positive
# (beta = 5.96, 95% CI [3.84, 8.07], t(58) = 5.63, p < .001; Std.
# beta = 0.36, 95% CI [0.23, 0.48])
#   - The effect of lake [LOTR] is statistically significant and
# negative (beta = -118.94, 95% CI [-151.07, -86.80], t(58) =
# -7.41, p < .001; Std. beta = -1.29, 95% CI [-1.64, -0.94])
#   - The effect of lake [Opeongo] is statistically non-significant
# and negative (beta = -16.74, 95% CI [-46.15, 12.68], t(58) =
# -1.14, p = 0.259; Std. beta = -0.18, 95% CI [-0.50, 0.14])
#   - The effect of lake [Shirley] is statistically significant and
# negative (beta = -167.36, 95% CI [-196.14, -138.57], t(58) =
# -11.64, p < .001; Std. beta = -1.81, 95% CI [-2.13, -1.50])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.

                            
# Tidy the model output
tidy_model_growthfla <- tidy(lm_model_fl_age_lake)
                            
# Create a table
tidy_model_growthfla %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	473.662392	17.709618	26.746053	0.0000000
# Age	5.956338	1.058278	5.628330	0.0000006
# lakeLOTR	-118.936379	16.055445	-7.407853	0.0000000
# lakeOpeongo	-16.735532	14.696162	-1.138769	0.2594813
# lakeShirley	-167.355696	14.378968	-11.638922	0.0000000


predictions <- OD_sheet %>%
  mutate(predicted_fork_length = predict(lm_model_fl_age_lake, newdata = .))
library(ggplot2)

ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = Age, y = predicted_fork_length, color = lake)) +
  labs(title = "Fork Length vs Age by Lake",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Lake") +
  theme_minimal()

                            
#### fork length age + lake logged as y just trying simple model ####
#### no interaction ####
# Fit a linear model
lm_model_fl_age_lakeln <- lm(log(fork_length) ~ Age + lake, data = Fish_data)

lm_model_fl_age_lakeln <- lm(log(fork_length) ~ Age + lake, data = OD_sheet)

# Print the summary of the model
summary(lm_model_fl_age_lakeln)
# Call:
#   lm(formula = log(fork_length) ~ Age + lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.19805 -0.05135 -0.01018  0.05139  0.21462 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  6.157280   0.037421 164.542  < 2e-16 ***
#   Age          0.011313   0.002236   5.059 4.53e-06 ***
#   lakeLOTR    -0.255627   0.033925  -7.535 3.70e-10 ***
#   lakeOpeongo -0.028241   0.031053  -0.909    0.367    
# lakeShirley -0.359747   0.030383 -11.840  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08696 on 58 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.8141,	Adjusted R-squared:  0.8012 
# F-statistic: 63.48 on 4 and 58 DF,  p-value: < 2.2e-16

check_model(lm_model_fl_age_lakeln)
forest_model(lm_model_fl_age_lakeln)
report(lm_model_fl_age_lakeln)
# We fitted a linear model (estimated using OLS) to predict fork_length with Age
# and lake (formula: log(fork_length) ~ Age + lake). The model explains a
# statistically significant and substantial proportion of variance (R2 = 0.81, F(4,
#                                                                                58) = 63.48, p < .001, adj. R2 = 0.80). The model's intercept, corresponding to
# Age = 0 and lake = Hogan, is at 6.16 (95% CI [6.08, 6.23], t(58) = 164.54, p <
# .001). Within this model:
# 
#   - The effect of Age is statistically significant and positive (beta = 0.01, 95%
# CI [6.84e-03, 0.02], t(58) = 5.06, p < .001; Std. beta = 0.12, 95% CI [0.06,
# 0.17])
#   - The effect of lake [LOTR] is statistically significant and negative (beta =
# -0.26, 95% CI [-0.32, -0.19], t(58) = -7.53, p < .001; Std. beta = -0.54, 95% CI
# [-0.70, -0.39])
#   - The effect of lake [Opeongo] is statistically non-significant and negative
# (beta = -0.03, 95% CI [-0.09, 0.03], t(58) = -0.91, p = 0.367; Std. beta = -0.05,
# 95% CI [-0.19, 0.09])
#   - The effect of lake [Shirley] is statistically significant and negative (beta =
# -0.36, 95% CI [-0.42, -0.30], t(58) = -11.84, p < .001; Std. beta = -0.78, 95% CI
# [-0.92, -0.65])
# 
# Standardized parameters were obtained by fitting the model on a standardized
# version of the dataset. 95% Confidence Intervals (CIs) and p-values were computed
#using a Wald t-distribution approximation.
# Tidy the model output
tidy_model_fl_age_lakeln <- tidy(lm_model_fl_age_lakeln)

# Create a table
tidy_model_fl_age_lakeln %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	6.1572804	0.0374207	164.5422540	0.0000000
# Age	0.0113130	0.0022362	5.0591380	0.0000045
# lakeLOTR	-0.2556267	0.0339254	-7.5349705	0.0000000
# lakeOpeongo	-0.0282408	0.0310532	-0.9094317	0.3668861
# lakeShirley	-0.3597466	0.0303830	-11.8404083	0.0000000

lnmodel_growth <- lm(log(fork_length) ~ Age + lake, data = fish_data)
predictionsln <- OD_sheet %>%
  mutate(predicted_lnfork_length = predict(lm_model_fl_age_lakeln, newdata = .),
         predicted_fork_length = exp(predicted_lnfork_length))  # Back-transform to original scale
library(ggplot2)

ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) +
  geom_point() +
  geom_line(data = predictionsln, aes(x = Age, y = predicted_fork_length, color = lake)) +
  labs(title = "Fork Length vs Age by Lake",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Lake") +
  theme_minimal()

#### Fork length vs age by lake with interaction ####


# Fit a linear model
lm_model_fl_age_lakelni <- lm(log(fork_length) ~ Age * lake, data = OD_sheet)

# Print the summary of the model
summary(lm_model_fl_age_lakelni)
# Call:
#   lm(formula = log(fork_length) ~ Age * lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.201167 -0.045921 -0.007094  0.041759  0.204162 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      6.126061   0.053077 115.418  < 2e-16 ***
#   Age              0.013571   0.003538   3.835 0.000325 ***
#   lakeLOTR        -0.171183   0.074085  -2.311 0.024630 *  
#   lakeOpeongo     -0.056319   0.074558  -0.755 0.453251    
# lakeShirley     -0.227375   0.085771  -2.651 0.010463 *  
#   Age:lakeLOTR    -0.009387   0.007198  -1.304 0.197616    
# Age:lakeOpeongo  0.002656   0.005296   0.501 0.618040    
# Age:lakeShirley -0.010192   0.006138  -1.661 0.102493    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08496 on 55 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.8317,	Adjusted R-squared:  0.8103 
# F-statistic: 38.82 on 7 and 55 DF,  p-value: < 2.2e-16

check_model(lm_model_fl_age_lakelni) # some collinearities
forest_model(lm_model_fl_age_lakelni)

# Tidy the model output
tidy_model_fl_age_lakelni <- tidy(lm_model_fl_age_lakelni)

# Create a table
tidy_model_fl_age_lakelni %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# 
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	6.1260613	0.0530774	115.4176083	0.0000000
# Age	0.0135714	0.0035385	3.8353706	0.0003251
# lakeLOTR	-0.1711827	0.0740849	-2.3106278	0.0246295
# lakeOpeongo	-0.0563185	0.0745578	-0.7553672	0.4532515
# lakeShirley	-0.2273752	0.0857713	-2.6509489	0.0104635
# Age:lakeLOTR	-0.0093867	0.0071976	-1.3041443	0.1976160
# Age:lakeOpeongo	0.0026557	0.0052958	0.5014721	0.6180403
# Age:lakeShirley	-0.0101919	0.0061376	-1.6605612	0.1024928



predictionslnai <- OD_sheet %>%
  mutate(predicted_lnai = predict(lm_model_fl_age_lakelni, newdata = .),
         predicted_lnai = exp(predicted_lnai))  
# Back-transform to original scale

ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
  geom_point() +
  geom_line(data = predictionslnai, aes(x = Age, y = predicted_lnai, color = lake)) +
  labs(title = "Fork Length vs Age by Lake",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Lake") +
  theme_minimal()
                            
# trying stargazer
                            
# Install and load the stargazer package
library(stargazer)
                            
# Create a summary table
stargazer(lm_model_fl_age_lake, type = "text",
                                      title = "Linear Regression Results",
                                      dep.var.labels = "Fork Length",
                                      covariate.labels = c("Age"),
                                      out = "linear_model_results.txt")
                            
# Linear Regression Results
# ===============================================
#   Dependent variable:    
#   ---------------------------
#   Fork Length        
# -----------------------------------------------
#   Age                          6.374***          
#   (1.289)          
# 
# lakeLOTR                    -118.337***        
#   (20.240)          
# 
# lakeOpeongo                    2.129           
# (17.440)          
# 
# lakeShirley                 -153.839***        
#   (16.894)          
# 
# Constant                    457.981***         
#   (22.357)          
# 
# -----------------------------------------------
#   Observations                    46             
# R2                             0.831           
# Adjusted R2                    0.814           
# Residual Std. Error      41.678 (df = 41)      
# F Statistic           50.382*** (df = 4; 41)   
# ===============================================
#   Note:               *p<0.1; **p<0.05; ***p<0.01
                            
                            
                            
                            
                            
                            
# Fit a linear model
lm_model <- lm(fork_length ~ Age, data = OD_sheet)
# Print the summary of the model
summary(lm_model)
# Call:
#   lm(formula = fork_length ~ Age, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -179.04  -55.85   21.58   63.38  147.95 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   372.44      24.01   15.51  < 2e-16 ***
#   Age             8.23       1.87    4.40 4.42e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 81.03 on 61 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.2409,	Adjusted R-squared:  0.2285 
# F-statistic: 19.36 on 1 and 61 DF,  p-value: 4.42e-05                         
ggplot(fish_data, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear regression line
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age with Linear Regression",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Locality")
                            
                            
# trying to add formula
# Fit the linear model
lm_model <- lm(log(fork_length) ~ Age, data = OD_sheet)
summary(lm_model)
# Call:
#   lm(formula = log(fork_length) ~ Age, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38524 -0.12205  0.04704  0.14365  0.30475 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 5.941235   0.051802 114.690  < 2e-16
# Age         0.016227   0.004035   4.021 0.000162
# 
# (Intercept) ***
#   Age         ***
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1748 on 61 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.2095,	Adjusted R-squared:  0.1966 
# F-statistic: 16.17 on 1 and 61 DF,  p-value: 0.0001619

# Extract the coefficients
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]
# Create the formula string
formula_text <- paste("fork_length =", round(intercept, 2), "+", round(slope, 2), "* Age")
# Plot with annotation
ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
theme_classic() +
labs(title = "Scatterplot of ln Fork Length vs. Age with Linear Regression",
                                   x = "Age (years)",
                                   y = "ln Fork Length (mm)",
                                   color = "Lake") +
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 1.5, vjust = 2, size = 4, color = "black")
                            
myplot_fl_agelm<- ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
theme_classic() +
labs(title = "Scatterplot of Natural Log Fork Length vs. Age with Linear Regression",
                                   x = "Age (years)",
                                   y = "ln Fork Length (mm)",
                                   color = "Lake") +
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 1.5, vjust = 2, size = 4, color = "black")
                            
                            
# Save the specific plot object with a proper file extension
ggsave("myplot_fl_agelm.png", plot = myplot_fl_agelm, width = 6, height = 4, dpi = 300)
                            
                            
#### fork length vs age ####
# Fit the linear model
lm_model <- lm(fork_length ~ Age, data = OD_sheet)
                            
# Extract the coefficients
intercept <- coef(lm_model)[1]
slope <- coef(lm_model)[2]
                            
# Create the formula string
formula_text <- paste("fork_length =", round(intercept, 2), "+", round(slope, 2), "* Age")
                            
# Plot with annotation
ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age with Linear Regression",
x = "Age (years)",
y = "Fork Length (mm)",
color = "Lake") +
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 1.5, vjust = 2, size = 4, color = "black")
                            
myplot_fl_agelm<- ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age with Linear Regression",
x = "Age (years)",
y = "Fork Length (mm)",
color = "Lake") +
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 1.5, vjust = 2, size = 4, color = "black")
                            
                            
# Save the specific plot object with a proper file extension
ggsave("myplot_fl_agelm.png", plot = myplot_fl_agelm, width = 6, height = 4, dpi = 300)
                            
                            
#### fitting a log curve ####
                            
# Fit a logarithmic model
log_model <- lm(fork_length ~ log(Age), data = OD_sheet)
                            
# Print the summary of the model
summary(log_model)
# Call:
#   lm(formula = fork_length ~ log(Age), data = fish_data)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -179.21  -51.61   17.25   61.98  163.95 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   200.91      66.50   3.021 0.004184 ** 
#   log(Age)      112.26      27.52   4.079 0.000187 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 83.35 on 44 degrees of freedom
# Multiple R-squared:  0.2744,	Adjusted R-squared:  0.2579 
# F-statistic: 16.64 on 1 and 44 DF,  p-value: 0.0001868
# Extract the coefficients
intercept <- coef(log_model)[1]
slope <- coef(log_model)[2]
# Create the formula string for annotation
formula_text <- paste("fork_length =", round(intercept, 2), "+", round(slope, 2), "* log(Age)")
                            
#### Plot with logarithmic fit ####
ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "green") +  # Logarithmic curve
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Fit",
x = "Age (years)",
y = "Fork Length (mm)",
color = "Lake") +
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 1.1, vjust = 2, size = 5, color = "black")
                            
# Plot with logarithmic fit and adjusted text position
ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
geom_point(size = 3) + 
geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "black") +  # Logarithmic curve
theme_classic() +
labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Fit",
x = "Age (years)",
y = "Fork Length (mm)",
color = "Lake") +
annotate("text", x = Inf, y = Inf, 
label = formula_text, hjust = 1.5, 
vjust = 2, size = 4, color = "black")
# adjusting the hlust
                            
annotate("text", x = Inf, y = Inf, label = formula_text, hjust = 0.8, vjust = 2, size = 5, color = "black")
                            
 # Save the specific plot object 
 myplot_fl_age <- ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) + 
 geom_point(size = 3) + 
  geom_smooth(method = "lm", formula = y ~ log(x), se = FALSE, color = "green") +  # Logarithmic curve
                   theme_classic() +
                  labs(title = "Scatterplot of Fork Length vs. Age with Logarithmic Fit",
                                   x = "Age (years)",
                                   y = "Fork Length (mm)",
                                   color = "Lake") +
                              annotate("text", x = Inf, y = Inf, 
                                       label = formula_text, hjust = 1.5, 
                                       vjust = 2, size = 4, color = "black")
                            myplot_fl_age
                            ggsave("myplot_fl_age", plot = myplot_fl_age, 
                                   width = 6, height = 4, dpi = 300)
                            
                            
   # Save the specific plot object with a proper file extension
  ggsave("myplot_fl_age.png", 
                 plot = myplot_fl_age, width = 6, height = 4, dpi = 300)
                            
                            
#### anova length ####
anova_length <- aov(fork_length ~ lake, data = OD_sheet)
summary(anova_length)
# Df Sum Sq Mean Sq F value Pr(>F)    
# lake         3 385390  128463   50.88 <2e-16 ***
#   Residuals   62 156530    2525                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
library(forestmodel)
forest_model(anova_length)
                            
TukeyHSD(anova_length)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = fork_length ~ lake, data = OD_sheet)
# 
# $lake
# diff        lwr        upr     p adj
# LOTR-Hogan      -153.23529 -198.73557 -107.73501 0.0000000
# Opeongo-Hogan    -24.87500  -71.08075   21.33075 0.4910795
# Shirley-Hogan   -173.75000 -219.95575 -127.54425 0.0000000
# Opeongo-LOTR     128.36029   82.15454  174.56605 0.0000000
# Shirley-LOTR     -20.51471  -66.72046   25.69105 0.6465079
# Shirley-Opeongo -148.87500 -195.77562 -101.97438 0.0000000
                            
anova_mass <- aov(mass_post_resp_g ~ lake, data = OD_sheet)
summary(anova_mass)
# Df   Sum Sq Mean Sq F value   Pr(>F)    
# lake         3 24431871 8143957   40.36 1.38e-14 ***
#   Residuals   62 12510581  201784                     
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_mass)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = mass_post_resp_g ~ lake, data = OD_sheet)
# 
# $lake
# diff        lwr       upr
# LOTR-Hogan      -1225.82353 -1632.5986 -819.0485
# Opeongo-Hogan    -128.00000  -541.0820  285.0820
# Shirley-Hogan   -1324.25000 -1737.3320 -911.1680
# Opeongo-LOTR     1097.82353   684.7415 1510.9055
# Shirley-LOTR      -98.42647  -511.5085  314.6555
# Shirley-Opeongo -1196.25000 -1615.5441 -776.9559
# p adj
# LOTR-Hogan      0.0000000
# Opeongo-Hogan   0.8457039
# Shirley-Hogan   0.0000000
# Opeongo-LOTR    0.0000000
# Shirley-LOTR    0.9223144
# Shirley-Opeongo 0.0000000
                            
anova_age <- aov(Age ~ lake, data = OD_sheet)
summary(anova_age)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# lake         3  364.7  121.57   4.744 0.00496 **
#   Residuals   59 1512.1   25.63                   
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 3 observations deleted due to missingness
TukeyHSD(anova_age)
# Tukey mu ltiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Age ~ lake, data = OD_sheet)
# 
# $lake
# diff        lwr       upr
# LOTR-Hogan      -6.3568627 -11.098212 -1.615513
# Opeongo-Hogan   -1.7568627  -6.498212  2.984487
# Shirley-Hogan   -1.0735294  -5.735500  3.588442
# Opeongo-LOTR     4.6000000  -0.287271  9.487271
# Shirley-LOTR     5.2833333   0.473032 10.093635
# Shirley-Opeongo  0.6833333  -4.126968  5.493635
# p adj
# LOTR-Hogan      0.0042130
# Opeongo-Hogan   0.7615317
# Shirley-Hogan   0.9288924
# Opeongo-LOTR    0.0721897
# Shirley-LOTR    0.0259264
# Shirley-Opeongo 0.9817670
                            
anova_lmass <- aov(liver_mass ~ lake, data = OD_sheet)
summary(anova_lmass)
  # Df Sum Sq Mean Sq F value   Pr(>F)    
# lake         3   4538  1512.6   32.97 7.29e-13 ***
#   Residuals   62   2845    45.9                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_lmass)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = liver_mass ~ lake, data = OD_sheet)
# 
# $lake
# diff        lwr        upr     p adj
# LOTR-Hogan      -16.658647 -22.792626 -10.524668 0.0000000
# Opeongo-Hogan    -1.478713  -7.707799   4.750372 0.9230915
# Shirley-Hogan   -17.869338 -24.098424 -11.640253 0.0000000
# Opeongo-LOTR     15.179934   8.950848  21.409019 0.0000001
# Shirley-LOTR     -1.210691  -7.439777   5.018394 0.9556576
# Shirley-Opeongo -16.390625 -22.713386 -10.067864 0.0000000
                            
anova_gmass <- aov(gonad_mass ~ lake, data = OD_sheet)
summary(anova_gmass)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# lake         3 164058   54686   9.142 4.27e-05 ***
#   Residuals   62 370884    5982                     
# ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_gmass)
# Tukey multiple comparisons of means
 # 95% family-wise confidence level
# 
 # Fit: aov(formula = gonad_mass ~ lake, data = OD_sheet)
  # 
 # $lake
# diff         lwr       upr     p adj
# LOTR-Hogan      -107.068529 -177.106696 -37.03036 0.0008556
# Opeongo-Hogan    -30.573360 -101.697455  40.55073 0.6695630
# Shirley-Hogan   -117.057485 -188.181580 -45.93339 0.0003009
# Opeongo-LOTR      76.495169    5.371075 147.61926 0.0302810
# Shirley-LOTR      -9.988956  -81.113050  61.13514 0.9824385
# Shirley-Opeongo  -86.484125 -158.677815 -14.29044 0.0126147
anova_vmass <- aov(ventricle_mass ~ lake, data = OD_sheet)
summary(anova_vmass)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# lake         3  20.96   6.987   31.93 1.33e-12 ***
#   Residuals   62  13.57   0.219                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(anova_vmass)
# Tukey multiple comparisons of means
 # 95% family-wise confidence level
# 
# Fit: aov(formula = ventricle_mass ~ lake, data = OD_sheet)
# 
 # $lake
 # diff        lwr        upr     p adj
  # LOTR-Hogan      -1.01294118 -1.4365785 -0.5893038 0.0000002
# Opeongo-Hogan    0.17029412 -0.2599116  0.6004999 0.7237411
# Shirley-Hogan   -1.06370588 -1.4939116 -0.6335001 0.0000001
 # Opeongo-LOTR     1.18323529  0.7530295  1.6134411 0.0000000
# Shirley-LOTR    -0.05076471 -0.4809705  0.3794411 0.9894189
# Shirley-Opeongo -1.23400000 -1.6706754 -0.7973246 0.0000000
                            
library(broom)
lm_length <- aov(fork_length ~ lake, data = OD_sheet)
summary(lm_length)
 # Df Sum Sq Mean Sq F value Pr(>F)    
# lake         3 385390  128463   50.88 <2e-16 ***
#   Residuals   62 156530    2525                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Generate a tidy summary of the linear model
forest_model(lm_length)
plot(lm_length)
anova(lm_length)
# Analysis of Variance Table
# 
# Response: fork_length
 # Df Sum Sq Mean Sq F value    Pr(>F)    
# lake       3 385390  128463  50.883 < 2.2e-16 ***
#   Residuals 62 156530    2525                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            
                            
anova_mass_post <- aov(mass_post_resp_g ~ lake, data = OD_sheet)
summary(anova_mass_post)
# Df   Sum Sq Mean Sq F value   Pr(>F)    
 # lake         3 24431871 8143957   40.36 1.38e-14 ***
#   Residuals   62 12510581  201784                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            
TukeyHSD(anova_mass_post)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
 # Fit: aov(formula = mass_post_resp_g ~ lake, data = OD_sheet)
# 
# $lake
# diff        lwr       upr     p adj
# LOTR-Hogan      -1225.82353 -1632.5986 -819.0485 0.0000000
 # Opeongo-Hogan    -128.00000  -541.0820  285.0820 0.8457039
# Shirley-Hogan   -1324.25000 -1737.3320 -911.1680 0.0000000
# Opeongo-LOTR     1097.82353   684.7415 1510.9055 0.0000000
# Shirley-LOTR      -98.42647  -511.5085  314.6555 0.9223144
# Shirley-Opeongo -1196.25000 -1615.5441 -776.9559 0.0000000
                            
                            
                            
 # Perform correlation test
result_mass <- cor.test(OD_sheet$standard_mass, OD_sheet$mass_post_resp_g)
                            
 # Print the result
 print(result_mass)
 # Pearson's product-moment correlation
# 
# data:  OD_sheet$standard_mass and OD_sheet$mass_post_resp_g
# t = 87.914, df = 64, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
 # 95 percent confidence interval:
#  0.9932662 0.9974869
# sample estimates:
#       cor 
 # 0.9958852 
                            
#### linear models per lake ####
library(dplyr)
library(purrr)
library(broom)
                            
#### lms mass port vs fl per separate lake ####
# Fit linear models for each lake
models <- OD_sheet %>%group_by(lake) %>% nest() %>% mutate(model = map(data, ~ lm(mass_post_resp_g ~ fork_length, data = .x)))
 # Extract model summaries
model_summaries <- models %>%
mutate(summary = map(model, glance)) %>% unnest(summary)
                            
# Print model summaries
                            
print(model_summaries)
# A tibble: 4 × 15
# Groups:   lake [4]
# lake    data     model  r.squared adj.r.squared sigma statistic
# <fct>   <list>   <list>     <dbl>         <dbl> <dbl>     <dbl>
#   1 Opeongo <tibble> <lm>       0.940         0.936 139.     219.  
# 2 LOTR    <tibble> <lm>       0.399         0.359 122.       9.96
# 3 Shirley <tibble> <lm>       0.950         0.947  43.4    267.  
# 4 Hogan   <tibble> <lm>       0.943         0.939 165.     247.  
# ℹ 8 more variables: p.value <dbl>, df <dbl>, logLik <dbl>,
#   AIC <dbl>, BIC <dbl>, deviance <dbl>, df.residual <int>,
#   nobs <int>                           
                            
# Extract detailed model summaries 
model_summaries <- models %>% mutate(tidy_summary = map(model, tidy), 
                                                        glance_summary = map(model, glance)) %>% select(lake, tidy_summary, glance_summary) %>% unnest(c(tidy_summary, glance_summary)) 
# Print model summaries 
print(model_summaries)
# A tibble: 4 × 15
# Groups:   lake [4]
# lake    data     model  r.squared adj.r.squared sigma statistic
# <fct>   <list>   <list>     <dbl>         <dbl> <dbl>     <dbl>
#   1 Opeongo <tibble> <lm>       0.940         0.936 139.     219.  
# 2 LOTR    <tibble> <lm>       0.399         0.359 122.       9.96
# 3 Shirley <tibble> <lm>       0.950         0.947  43.4    267.  
# 4 Hogan   <tibble> <lm>       0.943         0.939 165.     247.  
# ℹ 8 more variables: p.value <dbl>, df <dbl>, logLik <dbl>,
#   AIC <dbl>, BIC <dbl>, deviance <dbl>, df.residual <int>,
#   nobs <int>
                            
# one more time
                            
# Fit linear models for each lake
models <- OD_sheet %>%
      group_by(lake) %>%
      nest() %>%
mutate(model = map(data, ~ lm(mass_post_resp_g ~ fork_length, data = .x)))
# Extract detailed model summaries
    model_summaries <- models %>%
   mutate(tidy_summary = map(model, tidy),
                                     glance_summary = map(model, glance)) %>%
 select(lake, tidy_summary, glance_summary) %>%
 unnest(c(tidy_summary, glance_summary), names_sep = "_")
                            
 # Print model summaries
 print(model_summaries)
 # A tibble: 8 × 18
 # Groups:   lake [4]
# lake    tidy_summary_term tidy_summary_estimate tidy_summary_std.error
 # <fct>   <chr>                             <dbl>                  <dbl>
 #   1 Opeongo (Intercept)                    -3411.                  353.   
# 2 Opeongo fork_length                        9.80                  0.662
# 3 LOTR    (Intercept)                     -401.                  348.   
# 4 LOTR    fork_length                        2.72                  0.862
 # 5 Shirley (Intercept)                    -1293.                  116.   
 # 6 Shirley fork_length                        4.94                  0.302
# 7 Hogan   (Intercept)                    -3448.                  344.   
# 8 Hogan   fork_length                        9.65                  0.615
                            
 # ℹ 14 more variables: tidy_summary_statistic <dbl>, tidy_summary_p.value <dbl>,
#   glance_summary_r.squared <dbl>, glance_summary_adj.r.squared <dbl>,
 #   glance_summary_sigma <dbl>, glance_summary_statistic <dbl>,
#   glance_summary_p.value <dbl>, glance_summary_df <dbl>,
#   glance_summary_logLik <dbl>, glance_summary_AIC <dbl>,
 #   glance_summary_BIC <dbl>, glance_summary_deviance <dbl>,
#   glance_summary_df.residual <int>, glance_summary_nobs <int>
                            
                            
                            
#### look at the SBB data ####
attach(OD_sheet)
 



dim(OD_sheet)  # shows the number of rows x columns
#66 17
class(OD_sheet)  # shows the type or class of object
#"data.frame"
names(OD_sheet)  # shows the column names
# [1] "lake"             "sample_date"      "id"              
# [4] "mass_post_resp_g" "fork_length"      "total_length"    
# [7] "standard_mass"    "sex"              "maturity"        
# [10] "gonad_mass"       "ventricle_mass"   "liver_mass"      
# [13] "Age"              "OD_heart"         "OD_liver"        
# [16] "log_mass"         "log_length"            
str(OD_sheet)  # provides details about the object including 
# Omitting NAs
names(OD_sheet)
OD_sheet %>% na.omit()
                            
#### graphs OD ####                           
                            
SBB_OD1_data <- SBB_OD1 %>%select(Lake, Fork_length, Standard_mass:Liver_mass, OD) %>%
print()
                            
                            
                            
#### OD heart ####
                            
                        
ggplot(data = OD_sheet, mapping = aes(x = lake, y = OD_heart)) + 
geom_boxplot(aes(fill=lake), show.legend = FALSE) +
               theme_classic() #+
 # labs(x = lake, y = OD_heart)
                            
                            
                            
                            
# Create the plot
plot <- ggplot(data = OD_sheet, mapping = aes(x = lake, y = OD_heart)) + 
geom_boxplot(aes(fill = lake), show.legend = FALSE) +
theme_classic() +
labs(x = "Lake", y = "Optical Density of Lipofuscin's Heart") +  # Change axis labels
annotate("text", x = 3, y = 1.3, label = "*", size = 10, color= "red")  + # Add significance mark
annotate("text", x = 4, y = 1.3, label = "*", size = 10, color= "red")  
# Show the plot
print(plot)
# Save the plot as a PNG file
ggsave(filename = "plot_heart_OD_lakes.png", plot = plot, width = 6, height = 4, dpi = 300)
dev.off()
print(plot)
                            
#### anova by lakes ####
# tendency but not significant
# heart
anova_hrt_lakes <- aov(OD_heart ~ lake, data = OD_sheet)
summary(anova_hrt_lakes)
# Df Sum Sq Mean Sq F value Pr(>F)  
# lake         3 0.0674 0.02248   2.207 0.0961 .
# Residuals   62 0.6314 0.01018                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                        
                            
                            
#### linear model OD heart vs lakes ####
# some significance
# linear model
lm_hrt_lakes <- lm(OD_heart ~ lake, data = OD_sheet)
summary(lm_hrt_lakes)
# Call:
#   lm(formula = OD_heart ~ lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211203 -0.065621 -0.005789  0.048009  0.287235 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.06786    0.02447  43.632   <2e-16 ***
#   lakeLOTR    -0.06210    0.03461  -1.794   0.0777 .  
# lakeOpeongo -0.08049    0.03515  -2.290   0.0254 *  
#   lakeShirley -0.07266    0.03515  -2.067   0.0429 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1009 on 62 degrees of freedom
# Multiple R-squared:  0.0965,	Adjusted R-squared:  0.05278 
# F-statistic: 2.207 on 3 and 62 DF,  p-value: 0.09612                   
                            
library(report)
report(lm_hrt_lakes)
                            
# We fitted a linear model (estimated using OLS) to predict OD_heart with
# lake (formula: OD_heart ~ lake). The model explains a statistically not
# significant and weak proportion of variance (R2 = 0.10, F(3, 62) = 2.21,
# p = 0.096, adj. R2 = 0.05). The model's intercept, corresponding to lake
# = hogan, is at 1.07 (95% CI [1.02, 1.12], t(62) = 43.63, p < .001).
# Within this model:
# 
#   - The effect of lake [lotr] is statistically non-significant and
# negative (beta = -0.06, 95% CI [-0.13, 7.09e-03], t(62) = -1.79, p =
# 0.078; Std. beta = -0.60, 95% CI [-1.27, 0.07])
#   - The effect of lake [opeongo] is statistically significant and negative
# (beta = -0.08, 95% CI [-0.15, -0.01], t(62) = -2.29, p = 0.025; Std.
# beta = -0.78, 95% CI [-1.45, -0.10])
#   - The effect of lake [shirley] is statistically significant and negative
# (beta = -0.07, 95% CI [-0.14, -2.40e-03], t(62) = -2.07, p = 0.043; Std.
# beta = -0.70, 95% CI [-1.38, -0.02])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.
# Warning messages:
# 1: Function `data_find()` is deprecated and will be removed in a
#   future release. Please use `extract_column_names()` instead. 
# 2: Function `format_text()` is deprecated and will be removed in a
#   future release. Please use `text_format()` instead. 
                            
# Install and load the broom package
                            
library(broom)
                            
                            
# Generate a tidy summary of the linear model
tidy_summary <- tidy(lm_hrt_lakes)
print(tidy_summary)
# A tibble: 4 × 5
# term        estimate std.error statistic  p.value
# <chr>          <dbl>     <dbl>     <dbl>    <dbl>
#   1 (Intercept)   1.07      0.0245     43.6  2.99e-48
# 2 lakelotr     -0.0621    0.0346     -1.79 7.77e- 2
# 3 lakeopeongo  -0.0805    0.0351     -2.29 2.54e- 2
# 4 lakeshirley  -0.0727    0.0351     -2.07 4.29e- 2
                            
anova(lm_hrt_lakes)
# Analysis of Variance Table
# 
# Response: OD_heart
# Df  Sum Sq  Mean Sq F value  Pr(>F)  
# lake       3 0.06743 0.022478  2.2073 0.09612 .
# Residuals 62 0.63135 0.010183                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
check_model(lm_hrt_lakes)                           
forest_model(lm_hrt_lakes)

# Tidy the model output
tidy_model_lm_hrt_lakes <- tidy(lm_hrt_lakes)

# Create a table
tidy_model_lm_hrt_lakes %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

predictions_lm_hrt_lakes <- OD_sheet %>%
  mutate(predicted_lm_hrt_lakes = predict(lm_hrt_lakes, newdata = .),
         predicted_lm_hrt_lakes = exp(predicted_lm_hrt_lakes)) 

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	1.0678627	0.0244746	43.631502	0.0000000
# lakeLOTR	-0.0620980	0.0346123	-1.794104	0.0776726
# lakeOpeongo	-0.0804877	0.0351489	-2.289905	0.0254423
# lakeShirley	-0.0726596	0.0351489	-2.067192	0.0428972

predictions <- OD_sheet %>%
  mutate(predicted_fork_length = predict(lm_model, newdata = .))
library(ggplot2)

ggplot(OD_sheet, aes(x = Age, y = fork_length, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = Age, y = predicted_fork_length, color = lake)) +
  labs(title = "Fork Length vs Age by Lake",
       x = "Age (years)",
       y = "Fork Length (mm)",
       color = "Lake") +
  theme_minimal()


                            
#### linear model and correlations OD heart vs age ####
                            
#### Age graph heart OD ####
                            
# Create the line graph with regression lines
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart)) + geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
theme_classic() + labs(x = "Age", y = "Optical Density of Lipofuscin's Heart")
                            
# First, install and load the ggpmisc package if you haven't already

library(ggpmisc)
                            
# Create the line graph with regression lines and add the formula
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
 stat_poly_eq(aes(label = ..eq.label..), 
                                           formula = y ~ x, 
                                           parse = TRUE) +  # Adds the formula to the graph
                              theme_classic() +
labs(x = "Age", y = "Optical Density of Lipofuscin's Heart")
                            
                            
# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
 stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), 
                                           formula = y ~ x, 
                                           parse = TRUE) +  # Adds the formula and p-value to the graph
theme_classic() +
labs(x = "Age", y = "Optical Density of Lipofuscin's Heart")
                            
                            
                            
# Create the line graph
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, 
                                group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Heart")
                            
# Create the line graph
 ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, 
                              group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Heart")
 # Create the line graph with regression lines
  ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, 
                                   group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Heart")
                            
                            
                            
# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, group = lake, color = lake)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), 
group = lake, color = lake), formula = y ~ x, 
parse = TRUE) +  # Adds the formula and p-value to the graph
theme_classic() +
labs(x = "Age", y = "Optical Density of Heart")
                            
                            
# Calculate p-values for each lake
p_values <- OD_sheet %>% group_by(lake) %>%
summarise(p_value = cor.test(Age, OD_heart)$p.value)
 # Create the plot
plot <- ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, group = lake, color = lake)) + 
                              geom_point() + 
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Heart")
                            
# Annotate plot with p-values
# Annotate plot with p-values
plot <- plot + geom_text(data = p_values, aes(x = x, y = y, label = paste("p =", round(p_value, 3))), color = "black")
# Show the plot
print(plot)

                            
 # correlation
                            
 # Perform the correlation test
correlation_test <- cor.test(OD_sheet$Age, OD_sheet$OD_heart)
                            
# Print the results
print(correlation_test)
# Pearson's product-moment correlation
 # 
# data:  OD_sheet$Age and OD_sheet$OD_heart
 # t = -0.25399, df = 59, p-value = 0.8004
# alternative hypothesis: true correlation is not equal to 0
 # 95 percent confidence interval:
 #  -0.2825185  0.2206077
# sample estimates:
#       cor 
 # -0.033049 
  # Perform the correlation test
cor.test(OD_sheet$ventricle_mass, OD_sheet$OD_heart)
# 
# Pearson's product-moment correlation
 # 
# data:  OD_sheet$ventricle and OD_sheet$OD_heart
 # t = 0.70701, df = 64, p-value = 0.4821
 # alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1573525  0.3231798
# sample estimates:
#        cor 
# 0.08803299 
                            
 # Remove rows with missing values
OD_sheet_clean <- OD_sheet %>%
               filter(!is.na(Age) & !is.na(OD_heart))
                            
# Calculate correlation and p-values within each lake
results <- OD_sheet %>%
          group_by(lake) %>%
summarise(
      correlation = cor(Age, OD_heart),
                                p_value = cor.test(Age, OD_heart)$p.value
                              )
print(results)
 # A tibble: 4 × 3
 # lake    correlation p_value
  # <chr>         <dbl>   <dbl>
 #   1 Hogan        0.0465   0.859
# 2 LOTR        NA        0.181
 # 3 Opeongo     NA        0.745
# 4 Shirley     -0.0137   0.960
                            
                            
                            
# linear model
lm_hrt_age <- lm(OD_heart ~ Age, data = OD_sheet)
summary(lm_hrt_age)
# Call:
#   lm(formula = OD_heart ~ Age, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.23208 -0.06569 -0.01886  0.08692  0.27114 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.0254750  0.0315117  32.543   <2e-16 ***
#   Age         -0.0007227  0.0024547  -0.294    0.769    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1063 on 61 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.001419,	Adjusted R-squared:  -0.01495 
# F-statistic: 0.08668 on 1 and 61 DF,  p-value: 0.7694                    
                            
report(lm_hrt_age)

#### trying to reinstall report ####
# Install and load the necessary packages
install.packages("report")
install.packages("datawizard")

 library(report)
library(datawizard)

# Fit the linear model
lm_hrt_lake_age <- lm(OD_heart ~ ventricle_mass + lake + Age, data = OD_sheet)

# Generate the report
report(lm_hrt_lake_age)

                            
#### lakes and age for OD heart ####
                            
lm_hrt_lake_age <- lm(OD_heart ~ lake + Age, data = OD_sheet)
summary(lm_hrt_lake_age)
# Call:
#   lm(formula = OD_heart ~ lake + Age, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.210696 -0.069661 -0.005803  0.055137  0.277197 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.095902   0.044523  24.614   <2e-16 ***
#   lakeLOTR    -0.069957   0.040364  -1.733   0.0884 .  
# lakeOpeongo -0.082293   0.036947  -2.227   0.0298 *  
#   lakeShirley -0.074837   0.036149  -2.070   0.0429 *  
#   Age         -0.002028   0.002661  -0.762   0.4489    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1035 on 58 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1014,	Adjusted R-squared:  0.03939 
# F-statistic: 1.636 on 4 and 58 DF,  p-value: 0.1776                 
report(lm_hrt_lake_age)
# We fitted a linear model (estimated using OLS) to predict OD_heart
# with lake and Age (formula: OD_heart ~ lake + Age). The model
# explains a statistically not significant and weak proportion of
# variance (R2 = 0.10, F(4, 58) = 1.64, p = 0.178, adj. R2 = 0.04).
# The model's intercept, corresponding to lake = Hogan and Age = 0,
# is at 1.10 (95% CI [1.01, 1.19], t(58) = 24.61, p < .001). Within
# this model:
# 
#   - The effect of lake [LOTR] is statistically non-significant and
# negative (beta = -0.07, 95% CI [-0.15, 0.01], t(58) = -1.73, p =
# 0.088; Std. beta = -0.66, 95% CI [-1.43, 0.10])
#   - The effect of lake [Opeongo] is statistically significant and
# negative (beta = -0.08, 95% CI [-0.16, -8.34e-03], t(58) = -2.23,
# p = 0.030; Std. beta = -0.78, 95% CI [-1.48, -0.08])
#   - The effect of lake [Shirley] is statistically significant and
# negative (beta = -0.07, 95% CI [-0.15, -2.48e-03], t(58) = -2.07,
# p = 0.043; Std. beta = -0.71, 95% CI [-1.39, -0.02])
#   - The effect of Age is statistically non-significant and negative
# (beta = -2.03e-03, 95% CI [-7.35e-03, 3.30e-03], t(58) = -0.76, p
# = 0.449; Std. beta = -0.11, 95% CI [-0.38, 0.17])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.
check_model(lm_hrt_lake_age)
forest_model(lm_hrt_lake_age)

# Tidy the model output
tidy_model_lm_hrt_lake_age <- tidy(lm_hrt_lake_age)

# Create a table
tidy_model_lm_hrt_lake_age %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	1.0959019	0.0445227	24.6144374	0.0000000
# lakeLOTR	-0.0699568	0.0403641	-1.7331450	0.0883839
# lakeOpeongo	-0.0822930	0.0369468	-2.2273386	0.0298196
# lakeShirley	-0.0748371	0.0361493	-2.0702215	0.0428938
# Age	-0.0020284	0.0026606	-0.7623839	0.4489203

predictions <- OD_sheet %>%
  mutate(predicted_lm = predict(lm_hrt_lake_age, newdata = .))

ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = Age, y = predicted_lm, color = lake)) +
  labs(title = "OD vs Age by Lake",
       x = "Age (years)",
       y = "Optical Density of Heart Tissue",
       color = "Lake") +
  theme_minimal()

                            
lm_hrt_lake_age2 <- lm(OD_heart ~ lake * Age, data = OD_sheet)
summary(lm_hrt_lake_age2)

# Call:
#   lm(formula = OD_heart ~ lake * Age, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211110 -0.069854 -0.005495  0.059582  0.259517 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.0567061  0.0646920  16.334   <2e-16 ***
#   lakeLOTR         0.0574622  0.0902966   0.636    0.527    
# lakeOpeongo     -0.0446913  0.0908729  -0.492    0.625    
# lakeShirley     -0.0567412  0.1045402  -0.543    0.589    
# Age              0.0008071  0.0043128   0.187    0.852    
# lakeLOTR:Age    -0.0146511  0.0087726  -1.670    0.101    
# lakeOpeongo:Age -0.0027033  0.0064547  -0.419    0.677    
# lakeShirley:Age -0.0011806  0.0074807  -0.158    0.875    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1035 on 55 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1464,	Adjusted R-squared:  0.03771 
# F-statistic: 1.347 on 7 and 55 DF,  p-value: 0.2464

predictions2 <- OD_sheet %>%
  mutate(predicted_lm2 = predict(lm_hrt_lake_age2, newdata = .))

ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictions2, aes(x = Age, y = predicted_lm2, color = lake)) +
  labs(title = "OD vs Age by Lake",
       x = "Age (years)", 
       y = "Optical Density of Heart Tissue",
       color = "Lake") +
  theme_minimal()
report(lm_hrt_lake_age2)
# We fitted a linear model (estimated using OLS) to predict OD_heart
# with lake and Age (formula: OD_heart ~ lake * Age). The model
# explains a statistically not significant and moderate proportion
# of variance (R2 = 0.15, F(7, 55) = 1.35, p = 0.246, adj. R2 =
#                0.04). The model's intercept, corresponding to lake = Hogan and
# Age = 0, is at 1.06 (95% CI [0.93, 1.19], t(55) = 16.33, p <
# .001). Within this model:
# 
#   - The effect of lake [LOTR] is statistically non-significant and
# positive (beta = 0.06, 95% CI [-0.12, 0.24], t(55) = 0.64, p =
# 0.527; Std. beta = -1.07, 95% CI [-2.01, -0.13])
#   - The effect of lake [Opeongo] is statistically non-significant
# and negative (beta = -0.04, 95% CI [-0.23, 0.14], t(55) = -0.49, p
# = 0.625; Std. beta = -0.72, 95% CI [-1.44, -3.59e-04])
#   - The effect of lake [Shirley] is statistically non-significant
# and negative (beta = -0.06, 95% CI [-0.27, 0.15], t(55) = -0.54, p
# = 0.589; Std. beta = -0.67, 95% CI [-1.39, 0.05])
#   - The effect of Age is statistically non-significant and positive
# (beta = 8.07e-04, 95% CI [-7.84e-03, 9.45e-03], t(55) = 0.19, p =
# 0.852; Std. beta = 0.04, 95% CI [-0.41, 0.49])
#   - The effect of lake [LOTR] × Age is statistically non-significant
# and negative (beta = -0.01, 95% CI [-0.03, 2.93e-03], t(55) =
# -1.67, p = 0.101; Std. beta = -0.76, 95% CI [-1.68, 0.15])
#   - The effect of lake [Opeongo] × Age is statistically
# non-significant and negative (beta = -2.70e-03, 95% CI [-0.02,
# 0.01], t(55) = -0.42, p = 0.677; Std. beta = -0.14, 95% CI [-0.82,
# 0.53])
#   - The effect of lake [Shirley] × Age is statistically
# non-significant and negative (beta = -1.18e-03, 95% CI [-0.02,
# 0.01], t(55) = -0.16, p = 0.875; Std. beta = -0.06, 95% CI [-0.84,
# 0.72])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.

                            
#### model OD heart vs ventricle mass ####
lm_hrt_ventr <- lm(OD_heart ~ ventricle_mass, data = OD_sheet)
summary(lm_hrt_ventr)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.22246 -0.06808 -0.01495  0.06042  0.28729 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.00106    0.02320  43.142   <2e-16 ***
#   ventricle_mass  0.01252    0.01771   0.707    0.482    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1041 on 64 degrees of freedom
# Multiple R-squared:  0.00775,	Adjusted R-squared:  -0.007754 
# F-statistic: 0.4999 on 1 and 64 DF,  p-value: 0.4821 

lm_hrt_ventr2 <- lm(OD_heart ~ ventricle_mass + lake, data = OD_sheet)
summary(lm_hrt_ventr2)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass + lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211093 -0.065459 -0.005764  0.049134  0.287508 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.065560   0.049881  21.362   <2e-16 ***
#   ventricle_mass  0.001467   0.027617   0.053   0.9578    
# lakeLOTR       -0.060612   0.044723  -1.355   0.1803    
# lakeOpeongo    -0.080738   0.035746  -2.259   0.0275 *  
#   lakeShirley    -0.071099   0.046029  -1.545   0.1276    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
check_model(lm_hrt_ventr2)
forest_model(lm_hrt_ventr2)

# Tidy the model output
tidy_model_lm_hrt_ventr2 <- tidy(lm_hrt_ventr2)

# Create a table
tidy_model_lm_hrt_ventr2 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# cleaning ventricle mass
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(ventricle_mass) & !is.na(ventricle_mass) & is.finite(OD_heart) & !is.na(OD_heart))

lm_hrt_ventr2 <- lm(OD_heart ~ ventricle_mass + lake, data = OD_sheet_clean)

predictionsODV <- OD_sheet_clean %>%
  mutate(predicted_ODV = predict(lm_hrt_ventr2, newdata = .))

ggplot(OD_sheet_clean, aes(x = ventricle_mass, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionsODV, aes(x = ventricle_mass, y = predicted_ODV, color = lake)) +
  labs(title = "Optical Density of Heart vs Ventricle Mass by Lake",
       x = "Ventricle Mass (g)",
       y = "Optical Density of Heart",
       color = "Lake") +
  theme_minimal()

#### model OD heart vs ventricle mass ####
lm_hrt_ventr <- lm(OD_heart ~ ventricle_mass, data = OD_sheet)
summary(lm_hrt_ventr)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.22246 -0.06808 -0.01495  0.06042  0.28729 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.00106    0.02320  43.142   <2e-16 ***
#   ventricle_mass  0.01252    0.01771   0.707    0.482    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1041 on 64 degrees of freedom
# Multiple R-squared:  0.00775,	Adjusted R-squared:  -0.007754 
# F-statistic: 0.4999 on 1 and 64 DF,  p-value: 0.4821 

lm_hrt_ventr2 <- lm(OD_heart ~ ventricle_mass + lake, data = OD_sheet)
summary(lm_hrt_ventr2)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass + lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211093 -0.065459 -0.005764  0.049134  0.287508 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.065560   0.049881  21.362   <2e-16 ***
#   ventricle_mass  0.001467   0.027617   0.053   0.9578    
# lakeLOTR       -0.060612   0.044723  -1.355   0.1803    
# lakeOpeongo    -0.080738   0.035746  -2.259   0.0275 *  
#   lakeShirley    -0.071099   0.046029  -1.545   0.1276    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
check_model(lm_hrt_ventr2)
forest_model(lm_hrt_ventr2)


#### trying to clean the OD heart data ####
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(ventricle_mass) & !is.na(ventricle_mass) & is.finite(OD_heart) & !is.na(OD_heart))

lm_hrt_ventr2 <- lm(OD_heart ~ ventricle_mass + lake, data = OD_sheet_clean)

predictionsODV <- OD_sheet_clean %>%
  mutate(predicted_ODV = predict(lm_hrt_ventr2, newdata = .))

ggplot(OD_sheet_clean, aes(x = ventricle_mass, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionsODV, aes(x = ventricle_mass, y = predicted_ODV, color = lake)) +
  labs(title = "Optical Density of Heart vs Ventricle Mass by Lake",
       x = "Ventricle Mass (g)",
       y = "Optical Density of Heart",
       color = "Lake") +
  theme_minimal()


#### optical density of heart vs ventricle mass and lakes with interaction ####
lm_hrt_ventr3 <- lm(OD_heart ~ ventricle_mass * lake, data = OD_sheet)
summary(lm_hrt_ventr3)
# 
# Call:
#   lm(formula = OD_heart ~ ventricle_mass * lake, data = OD_sheet)\ 
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.210356 -0.062448 -0.005984  0.051472  0.219714 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                 1.08719    0.07001  15.528   <2e-16
# ventricle_mass             -0.01231    0.04184  -0.294   0.7696
# lakeLOTR                    0.14528    0.14601   0.995   0.3239
# lakeOpeongo                -0.11620    0.09903  -1.173   0.2454
# lakeShirley                -0.17941    0.09903  -1.812   0.0752
# ventricle_mass:lakeLOTR    -0.39486    0.22982  -1.718   0.0911
# ventricle_mass:lakeOpeongo  0.02173    0.05625   0.386   0.7007
# ventricle_mass:lakeShirley  0.18508    0.13591   1.362   0.1785
# 
# (Intercept)                ***
#   ventricle_mass                
# lakeLOTR                      
# lakeOpeongo                   
# lakeShirley                .  
# ventricle_mass:lakeLOTR    .  
# ventricle_mass:lakeOpeongo    
# ventricle_mass:lakeShirley    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09996 on 58 degrees of freedom
# Multiple R-squared:  0.1706,	Adjusted R-squared:  0.07049 
# F-statistic: 1.704 on 7 and 58 DF,  p-value: 0.1259
check_model(lm_hrt_ventr3)
forest_model(lm_hrt_ventr3)

# Tidy the model output
tidy_model_lm_hrt_ventr3 <- tidy(lm_hrt_ventr3)

# Create a table
tidy_model_lm_hrt_ventr3 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	1.0871870	0.0700141	15.5281198	0.0000000
# ventricle_mass	-0.0123108	0.0418437	-0.2942085	0.7696489
# lakeLOTR	0.1452751	0.1460114	0.9949573	0.3238900
# lakeOpeongo	-0.1162020	0.0990276	-1.1734311	0.2454193
# lakeShirley	-0.1794070	0.0990323	-1.8116010	0.0752264
# ventricle_mass:lakeLOTR	-0.3948585	0.2298176	-1.7181384	0.0911061
# ventricle_mass:lakeOpeongo	0.0217303	0.0562542	0.3862877	0.7006971
# ventricle_mass:lakeShirley	0.1850837	0.1359070	1.3618401	0.1785147

# cleaning ventricle mass
OD_sheet_clean <- OD_sheet %>%
  filter(is.finite(ventricle_mass) & !is.na(ventricle_mass) & is.finite(OD_heart) & !is.na(OD_heart))

lm_hrt_ventr3 <- lm(OD_heart ~ ventricle_mass * lake, data = OD_sheet_clean)

predictionsODV3 <- OD_sheet_clean %>%
  mutate(predicted_ODV = predict(lm_hrt_ventr3, newdata = .))

ggplot(OD_sheet_clean, aes(x = ventricle_mass, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionsODV3, aes(x = ventricle_mass, y = predicted_ODV, color = lake)) +
  labs(title = "Optical Density of Heart vs Ventricle Mass by Lake",
       x = "Ventricle Mass (g)",
       y = "Optical Density of Heart",
       color = "Lake") +
  theme_minimal()

                           
#### correlation ventr mass and OD ####
                            # Perform correlation test
correlation_OD_ventr <- cor.test(OD_sheet$OD_heart, OD_sheet$ventricle_mass)
                            
# Print the test results
print(correlation_OD_ventr)
# Pearson's product-moment correlation
# 
# data:  OD_sheet$OD_heart and OD_sheet$ventricle_mass
# t = 0.70701, df = 64, p-value = 0.4821
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1573525  0.3231798
# sample estimates:
#        cor 
# 0.08803299                             
                            
library(dplyr)
library(purrr)
                            
                            
# Function to perform correlation test
cor_test <- function(df) {
cor.test(df$OD_heart, df$ventricle_mass)
                            }
                            
# Apply correlation test to each lake
correlation_results <- OD_sheet %>%
            group_by(lake) %>%
            nest() %>%
            mutate(correlation = map(data, cor_test))
                            
# Extract and print results
correlation_results %>%
mutate(correlation_summary = map(correlation, broom::tidy)) %>% select(lake, correlation_summary) %>%
unnest(correlation_summary)
# A tibble: 4 × 9
# Groups:   lake [4]
# lake    estimate statistic p.value parameter conf.low conf.high
# <fct>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl>
#   1 Opeongo   0.110      0.413   0.686        14   -0.408    0.574 
# 2 LOTR     -0.409     -1.74    0.103        15   -0.744    0.0889
# 3 Shirley   0.290      1.13    0.276        14   -0.240    0.687 
# 4 Hogan    -0.0706    -0.274   0.788        15   -0.533    0.424 
# # ℹ 2 more variables: method <chr>, alternative <chr>                         
                            
#### anova OD heart vs lakes ####
anova_hrt_lakes <- aov(OD_heart ~ Lake, data = OD_sheet)
summary(anova_hrt_lakes)
# Df Sum Sq Mean Sq F value Pr(>F)  
# lake         3 0.0674 0.02248   2.207 0.0961 .
# Residuals   62 0.6314 0.01018                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                          
forest_model(lm_hrt_lakes)
                            
# First, install and load the required packages if you haven't already
  install.packages("ggpmisc")
install.packages("ggrepel")
library(ggpmisc)
library(ggrepel)
                            
# Load necessary packages
library(dplyr)
library(broom)
# Filter out non-finite and missing values
OD_sheet_filtered <- OD_sheet %>%
        filter(is.finite(Age), is.finite(OD_heart))
                            
# Fit the linear model for each lake and extract coefficients
equation_labels <- OD_sheet_filtered %>%
    group_by(lake) %>%
    do({
model <- lm(OD_heart ~ Age, data = .)
data.frame(lake = unique(.$lake),intercept = coef(model)[1],
slope = coef(model)[2],p_value = glance(model)$p.value
                                ) }) %>%
  mutate(eq_label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),p_value_label = paste0("p = ", format.pval(p_value, digits = 2)))
                            
# Merge the labels with the original data
labels_data <- OD_sheet_filtered %>%
inner_join(equation_labels, by = "lake") %>%
group_by(lake) %>%summarize(x = max(Age), y = min(OD_heart), eq_label = first(eq_label), p_value_label = first(p_value_label))
# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet_filtered, mapping = aes(x = Age, y = OD_heart, group = lake, color = lake)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
geom_label_repel(data = labels_data, aes(x = x, y = y, label = paste(eq_label, p_value_label, sep = "~~~"), color = lake),
nudge_x = -5, nudge_y = 0.5, size = 3) +  # Positions the text labels
theme_classic() +labs(x = "Age", y = "Optical Density of Heart")

# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, group = lake, color = lake)) + geom_point() +  # Adds points to the line graph
  geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), group = lake, color = lake), formula = y ~ x, parse = TRUE) +  # Adds the formula and p-value to the graph
  theme_classic() +
  labs(x = "Age", y = "Optical Density of Heart")

heart_OD_agebylake <- ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_heart, group = lake, color = lake)) + 
  geom_point() +  # Adds points to the line graph
  geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
  stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), group = lake, color = lake), formula = y ~ x, parse = TRUE) +  # Adds the formula and p-value to the graph
  theme_classic() +
  labs(x = "Age", y = "Optical Density of Heart")
ggsave("heart_OD_agebylake.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(heart_OD_agebylake)                            
                            
#### OD for liver ####
# Load the package
library(ggpmisc)
# Create the line graph with regression lines and add the  ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver)) + 
 geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),formula = y ~ x, parse = TRUE) +  # Adds the formula and p-value to the graph
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Lipofuscin's Liver")
                            
# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), group = lake, color = lake), formula = y ~ x, parse = TRUE) +  # Adds the formula and p-value to the graph
theme_classic() +
labs(x = "Age", y = "Optical Density of Liver")
                            
liver_OD_agebylake <- ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), group = lake, color = lake), formula = y ~ x, parse = TRUE) +  # Adds the formula and p-value to the graph
theme_classic() +
labs(x = "Age", y = "Optical Density of Liver")
ggsave("liver_OD_agebylake.png", plot = final_plot, width = 10, height = 8, dpi = 300)
print(liver_OD_agebylake)


# Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
geom_point() +  # Adds points to the line graph
geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              stat_poly_eq(aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~"), 
                                               group = lake, color = lake), 
                                           formula = y ~ x, 
                                           parse = TRUE) +  # Adds the formula and p-value to the graph
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver") +
                              theme(legend.position = c(0.85, 0.15))  # Positions the legend at the bottom right
                            
                            
  # Calculate the regression equations and p-values for each group
library(dplyr)
library(broom)
                            
 # Fit the linear model for each lake and extract coefficients
 equation_labels <- OD_sheet %>%
              group_by(lake) %>%
               do({
               model <- lm(OD_liver ~ Age, data = .)
                   data.frame(
                                  lake = unique(.$lake),
                                  intercept = coef(model)[1],
                                  slope = coef(model)[2],
                                  p_value = glance(model)$p.value
                                )
                              }) %>%
         mutate(eq_label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),
                                     p_value_label = paste0("p = ", format.pval(p_value, digits = 2)))
                            
 # Merge the labels with the original data
labels_data <- OD_sheet %>%
                  inner_join(equation_labels, by = "lake") %>%
                   group_by(lake) %>%
                  summarize(x = max(Age), y = min(OD_liver), eq_label = first(eq_label), p_value_label = first(p_value_label))
 # Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              geom_text(data = labels_data, aes(x = x, y = y, label = paste(eq_label, p_value_label, sep = "~~~")),
                                        position = position_nudge(x = -0.5, y = 0.5), color = "black", size = 3) +  # Positions the text labels
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
                            
 # Filter out non-finite and missing values
OD_sheet_filtered <- OD_sheet %>%
                              filter(is.finite(Age), is.finite(OD_liver))
                            
 # Fit the linear model for each lake and extract coefficients
 equation_labels <- OD_sheet_filtered %>%
                              group_by(lake) %>%
                              do({
                                model <- lm(OD_liver ~ Age, data = .)
                                data.frame(
                                  lake = unique(.$lake),
                                  intercept = coef(model)[1],
                                  slope = coef(model)[2],
                                  p_value = glance(model)$p.value
                                )
                              }) %>%
                              mutate(eq_label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),
                                     p_value_label = paste0("p = ", format.pval(p_value, digits = 2)))
                            
 # Merge the labels with the original data
 labels_data <- OD_sheet_filtered %>%
                              inner_join(equation_labels, by = "lake") %>%
                              group_by(lake) %>%
                              summarize(x = max(Age), y = min(OD_liver), eq_label = first(eq_label), p_value_label = first(p_value_label))
                            
 # Create the line graph with regression lines and add the formula and p-value
 ggplot(data = OD_sheet_filtered, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              geom_text(data = labels_data, aes(x = x, y = y, label = paste(eq_label, p_value_label, sep = "~~~")),
                                        position = position_nudge(x = -0.5, y = 0.5), color = "black", size = 3) +  # Positions the text labels
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
                            
 # First, install and load the required packages if you haven't already
install.packages("ggpmisc")
install.packages("ggrepel")
library(ggpmisc)
library(ggrepel)
                            
# Load necessary packages
library(dplyr)
library(broom)
                            
# Filter out non-finite and missing values
                             <- OD_sheet %>%
                              filter(is.finite(Age), is.finite(OD_liver))
                            
 # Fit the linear model for each lake and extract coefficients
  equation_labels <- OD_sheet_filtered %>%
                group_by(lake) %>%
                              do({
                                model <- lm(OD_liver ~ Age, data = .)
                                data.frame(
                                  lake = unique(.$lake),
                                  intercept = coef(model)[1],
                                  slope = coef(model)[2],
                                  p_value = glance(model)$p.value
                                )
                              }) %>%
             mutate(eq_label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),
                                     p_value_label = paste0("p = ", format.pval(p_value, digits = 2)))
                            
# Merge the labels with the original data
 labels_data <- OD_sheet_filtered %>%
                              inner_join(equation_labels, by = "lake") %>%
                              group_by(lake) %>%
                              summarize(x = max(Age), y = min(OD_liver), eq_label = first(eq_label), p_value_label = first(p_value_label))
                            
 # Create the line graph with regression lines and add the formula and p-value
 ggplot(data = OD_sheet_filtered, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              geom_label_repel(data = labels_data, aes(x = x, y = y, label = paste(eq_label, p_value_label, sep = "~~~")),
                                               nudge_x = -5, nudge_y = 0.5, color = "black", size = 3) +  # Positions the text labels
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
 # First, install and load the required packages if you haven't already
 install.packages("ggpmisc")
 install.packages("ggrepel")
library(ggpmisc)
library(ggrepel)
# Load necessary packages
library(dplyr)
library(broom)
                            
# Filter out non-finite and missing values
 OD_sheet_filtered <- OD_sheet %>%
             filter(is.finite(Age), is.finite(OD_liver))
                            
 # Fit the linear model for each lake and extract coefficients
equation_labels <- OD_sheet_filtered %>%
                              group_by(lake) %>%
                              do({
                                model <- lm(OD_liver ~ Age, data = .)
                                data.frame(
                                  lake = unique(.$lake),
                                  intercept = coef(model)[1],
                                  slope = coef(model)[2],
                                  p_value = glance(model)$p.value
                                )
                              }) %>%
              mutate(eq_label = paste0("y = ", round(slope, 2), "x + ", round(intercept, 2)),
                                     p_value_label = paste0("p = ", format.pval(p_value, digits = 2)))
                            
 # Merge the labels with the original data
 labels_data <- OD_sheet_filtered %>%
                              inner_join(equation_labels, by = "lake") %>%
                              group_by(lake) %>%
                              summarize(x = max(Age), y = min(OD_liver), eq_label = first(eq_label), p_value_label = first(p_value_label))
                            
 # Create the line graph with regression lines and add the formula and p-value
ggplot(data = OD_sheet_filtered, mapping = aes(x = Age, y = OD_liver, group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              geom_label_repel(data = labels_data, aes(x = x, y = y, label = paste(eq_label, p_value_label, sep = "~~~"), color = lake),
                                               nudge_x = -5, nudge_y = 0.5, size = 3) +  # Positions the text labels
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
                            
 # Create the line graph with regression lines
 ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Lipofuscin's Liver")
                            
# Create the line graph
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, 
                       group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
# Create the line graph
ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, 
                                      group = lake, color = lake)) + 
geom_point() +  # Adds points to the line graph
theme_classic() +
labs(x = "Age", y = "Optical Density of Liver")
                            
# Create the line graph with regression lines
                            ggplot(data = OD_sheet, mapping = aes(x = Age, y = OD_liver, 
                                                                  group = lake, color = lake)) + 
                              geom_point() +  # Adds points to the line graph
                              geom_smooth(method = "lm", se = FALSE) +  # Adds regression lines without confidence intervals
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
                            # Calculate p-values for each lake
                            p_values3 <- OD_sheet %>%
                              group_by(lake) %>%
                              summarise(p_value = cor.test(Age, OD_liver)$p.value)
                            
                            
                            # Create the plot
                            plot3 <- ggplot(data = OD_sheet, 
                                            mapping = aes(x = Age, y = OD_liver, 
                                                          group = lake, color = lake)) + 
                              geom_point() + 
                              theme_classic() +
                              labs(x = "Age", y = "Optical Density of Liver")
                            
                            # Annotate plot with p-values
                            # Annotate plot with p-values
                            plot33 <- plot3 + 
                              geom_text(data = p_values, aes(x = x, y = y, label = paste("p =", round(p_value, 3))), color = "black")
                            
                            # Show the plot
                            print(plot33)
                            # Show the plot
                            print(plot33)
                            
                            # x.lab <- expression("Lakes")
                            ggplot(data = OD_sheet, mapping = aes(x = lake, y = OD_liver)) + 
                              geom_boxplot(aes(fill=lake), show.legend = FALSE) +
                              theme_classic()  +
                              labs(x = "Lake", y = "Optical Density of Lipofuscin's Liver")
                            # labs(x = lake, y = OD_heart)
                            
#### linear models correlation for liver ####
#### lakes and age for OD liver ####
#### linear model OD Liver vs age & lake ####
lm_lvr_age_lake <- lm(OD_liver ~ Age + lake, data = OD_sheet)
summary(lm_lvr_age_lake)
# 
# Call:
#   lm(formula = OD_liver ~ Age + lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48955 -0.05233  0.01598  0.09219  0.24616 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.975060   0.084591  11.527 1.37e-14
# Age         -0.002991   0.004894  -0.611    0.544
# lakeLOTR    -0.003216   0.073958  -0.043    0.966
# lakeOpeongo -0.055273   0.066547  -0.831    0.411
# lakeShirley -0.046624   0.063379  -0.736    0.466
# 
# (Intercept) ***
#   Age            
# lakeLOTR       
# lakeOpeongo    
# lakeShirley    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1562 on 42 degrees of freedom
# (19 observations deleted due to missingness)
# Multiple R-squared:  0.03994,	Adjusted R-squared:  -0.0515 
# F-statistic: 0.4368 on 4 and 42 DF,  p-value: 0.7813
library(report)
report(lm_hrt_age_lake)
# We fitted a linear model (estimated using OLS) to
# predict OD_heart with Age and lake (formula:
#                                       OD_heart ~ Age + lake). The model explains a
# statistically not significant and weak proportion
# of variance (R2 = 0.10, F(4, 58) = 1.64, p =
#                0.178, adj. R2 = 0.04). The model's intercept,
# corresponding to Age = 0 and lake = Hogan, is at
# 1.10 (95% CI [1.01, 1.19], t(58) = 24.61, p <
# .001). Within this model:
# 
#   - The effect of Age is statistically
# non-significant and negative (beta = -2.03e-03,
# 95% CI [-7.35e-03, 3.30e-03], t(58) = -0.76, p =
# 0.449; Std. beta = -0.11, 95% CI [-0.38, 0.17])
#   - The effect of lake [LOTR] is statistically
# non-significant and negative (beta = -0.07, 95% CI
# [-0.15, 0.01], t(58) = -1.73, p = 0.088; Std. beta
# = -0.66, 95% CI [-1.43, 0.10])
#   - The effect of lake [Opeongo] is statistically
# significant and negative (beta = -0.08, 95% CI
# [-0.16, -8.34e-03], t(58) = -2.23, p = 0.030; Std.
# beta = -0.78, 95% CI [-1.48, -0.08])
#   - The effect of lake [Shirley] is statistically
# significant and negative (beta = -0.07, 95% CI
# [-0.15, -2.48e-03], t(58) = -2.07, p = 0.043; Std.
# beta = -0.71, 95% CI [-1.39, -0.02])
# 
# Standardized parameters were obtained by fitting
# the model on a standardized version of the
# dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution
# approximation.
                            
#### linear model OD liver vs age & lake interaction ####
lm_lvr_age_lakei <- lm(OD_liver ~ Age * lake, data = OD_sheet)
summary(lm_lvr_age_lakei)
# 
# Call:
#   lm(formula = OD_liver ~ Age * lake, data = OD_sheet)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48596 -0.03589  0.01626  0.09633  0.22074 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)      0.850803   0.113814   7.475
# Age              0.005379   0.007089   0.759
# lakeLOTR         0.162559   0.169099   0.961
# lakeOpeongo      0.211660   0.171900   1.231
# lakeShirley      0.106015   0.187719   0.565
# Age:lakeLOTR    -0.014078   0.017436  -0.807
# Age:lakeOpeongo -0.019604   0.011748  -1.669
# Age:lakeShirley -0.010829   0.014267  -0.759
# Pr(>|t|)    
# (Intercept)     4.81e-09 ***
#   Age                0.453    
# lakeLOTR           0.342    
# lakeOpeongo        0.226    
# lakeShirley        0.575    
# Age:lakeLOTR       0.424    
# Age:lakeOpeongo    0.103    
# Age:lakeShirley    0.452    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1562 on 39 degrees of freedom
# (19 observations deleted due to missingness)
# Multiple R-squared:  0.1085,	Adjusted R-squared:  -0.0515 
# F-statistic: 0.6781 on 7 and 39 DF,  p-value: 0.6894
forest_model(lm_lvr_age_lakei)
report(lm_lvr_age_lakei)
# Generate predictions
                            predictionslm_lvr_age_lakei <- OD_sheet %>%
                              mutate(predicted_lm_lvr_age_lakei = predict(lm_lvr_age_lakei, newdata = .))
                            
# Create the plot
ggplot(OD_sheet, aes(x = Age, y = OD_liver, color = lake)) +
geom_point() +
geom_line(data = predictionslm_lvr_age_lakei, aes(x = Age, y = predicted_lm_lvr_age_lakei, color = lake)) +
labs(title = "Optical Density Liver vs Age by Lake",
x = "Age (y)",
y = "Optical Density of Lipofuscin's Liver",
color = "Lake") +
theme_minimal()
                            
                            
                            
#### od liver vs age ####
lm_lvr_age <- lm(OD_liver ~ Age, data = OD_sheet)
summary(lm_lvr_age)
# Call:
#   lm(formula = OD_liver ~ Age, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.51049 -0.04096  0.02757  0.08576  0.21507 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.952170   0.053216  17.892   <2e-16 ***
#   Age         -0.003206   0.004128  -0.777    0.441    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.153 on 45 degrees of freedom
# (19 observations deleted due to missingness)
# Multiple R-squared:  0.01323,	Adjusted R-squared:  -0.008701 
# F-statistic: 0.6032 on 1 and 45 DF,  p-value: 0.4414       

library(report)
report(lm_lvr_age)
# We fitted a linear model (estimated using OLS) to predict
# OD_liver with Age (formula: OD_liver ~ Age). The model explains a
# statistically not significant and very weak proportion of
# variance (R2 = 0.01, F(1, 45) = 0.60, p = 0.441, adj. R2 =
#             -8.70e-03). The model's intercept, corresponding to Age = 0, is
# at 0.95 (95% CI [0.84, 1.06], t(45) = 17.89, p < .001). Within
# this model:
# 
#   - The effect of Age is statistically non-significant and negative
# (beta = -3.21e-03, 95% CI [-0.01, 5.11e-03], t(45) = -0.78, p =
# 0.441; Std. beta = -0.12, 95% CI [-0.41, 0.18])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.
                            
check_model(lm_lvr_age)
forest_model(lm_lvr_age)
                            
# Tidy the model output
tidy_model_lm_lvr_age <- tidy(lm_lvr_age)
                            
# Create a table
tidy_model_lm_hrt_age %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)
                            
                            
predictions <- OD_sheet %>% mutate(predicted_lm = predict(lm_hrt_lake_age, newdata = .))
                            
ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
geom_point() + geom_line(data = predictions, aes(x = Age, y = predicted_lm, color = lake)) + labs(title = "OD vs Age by Lake", x = "Age (years)", y = "Optical Density of Heart Tissue", color = "Lake") +
theme_minimal()


#### OD liver vs liver mass only ####

lm_lvr_lmass <- lm(OD_liver ~ liver_mass, data = OD_sheet)
summary(lm_lvr_lmass)
# Call:
#   lm(formula = OD_liver ~ liver_mass, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.52218 -0.04268  0.03023  0.07836  0.22846 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.924713   0.036720  25.183   <2e-16 ***
#   liver_mass  -0.000432   0.002049  -0.211    0.834    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1533 on 47 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.0009447,	Adjusted R-squared:  -0.02031 
# F-statistic: 0.04444 on 1 and 47 DF,  p-value: 0.8339      

library(report)
report(lm_lvr_age)

check_model(lm_lvr_age)
forest_model(lm_lvr_age)

# Tidy the model output
tidy_model_lm_lvr_age <- tidy(lm_lvr_age)

# Create a table
tidy_model_lm_hrt_age %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


predictions <- OD_sheet %>% mutate(predicted_lm = predict(lm_hrt_lake_age, newdata = .))

ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
  geom_point() + geom_line(data = predictions, aes(x = Age, y = predicted_lm, color = lake)) + labs(title = "OD vs Age by Lake", x = "Age (years)", y = "Optical Density of Heart Tissue", color = "Lake") +
  theme_minimal()
#### OD liver vs sex ####
# OD liver vs liver mass only

lm_lvr_sex <- lm(OD_liver ~ sex, data = OD_sheet)
summary(lm_lvr_sex)
# Call:
#   lm(formula = OD_liver ~ sex, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.50546 -0.05275  0.01204  0.08754  0.24254 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.93425    0.03447  27.100   <2e-16 ***
#   sexmale     -0.02879    0.04514  -0.638    0.527    
# sexunknown   0.03425    0.15798   0.217    0.829    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1542 on 46 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.01105,	Adjusted R-squared:  -0.03195 
# F-statistic: 0.257 on 2 and 46 DF,  p-value: 0.7744

#### OD liver vs lake ####
lm_lvr_lake <- lm(OD_liver ~ lake, data = OD_sheet)
summary(lm_lvr_lake)
# 
# Call:
#   lm(formula = OD_liver ~ lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49392 -0.05015  0.01450  0.09908  0.24536 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.93065    0.04304  21.623   <2e-16 ***
#   lakeLOTR     0.01585    0.06212   0.255    0.800    
# lakeOpeongo -0.02802    0.06357  -0.441    0.662    
# lakeShirley -0.03673    0.06087  -0.603    0.549    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1552 on 45 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.01985,	Adjusted R-squared:  -0.04549 
# F-statistic: 0.3038 on 3 and 45 DF,  p-value: 0.8225      

library(report)
report(lm_lvr_lake)
# We fitted a linear model (estimated using OLS) to
# predict OD_liver with liver_mass and lake (formula:
#                                              OD_liver ~ liver_mass * lake). The model explains a
# statistically not significant and weak proportion
# of variance (R2 = 0.04, F(7, 41) = 0.23, p = 0.976,
#              adj. R2 = -0.13). The model's intercept,
# corresponding to liver_mass = 0 and lake = Hogan,
# is at 0.94 (95% CI [0.71, 1.16], t(41) = 8.49, p <
# .001). Within this model:
# 
#   - The effect of liver mass is statistically
# non-significant and negative (beta = -3.02e-04, 95%
# CI [-9.23e-03, 8.62e-03], t(41) = -0.07, p = 0.946;
# Std. beta = -0.02, 95% CI [-0.66, 0.61])
#   - The effect of lake [LOTR] is statistically
# non-significant and negative (beta = -0.08, 95% CI
# [-0.44, 0.29], t(41) = -0.42, p = 0.678; Std. beta
# = 0.75, 95% CI [-1.62, 3.12])
#   - The effect of lake [Opeongo] is statistically
# non-significant and positive (beta = 0.04, 95% CI
# [-0.31, 0.40], t(41) = 0.24, p = 0.809; Std. beta =
# -0.01, 95% CI [-1.20, 1.17])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -0.05, 95% CI
# [-0.37, 0.27], t(41) = -0.30, p = 0.768; Std. beta
# = -0.23, 95% CI [-2.38, 1.93])
#   - The effect of liver mass × lake [LOTR] is
# statistically non-significant and positive (beta =
# 0.01, 95% CI [-0.03, 0.06], t(41) = 0.62, p =
# 0.538; Std. beta = 0.94, 95% CI [-2.12, 4.00])
#   - The effect of liver mass × lake [Opeongo] is
# statistically non-significant and negative (beta =
# -3.13e-03, 95% CI [-0.02, 0.01], t(41) = -0.44, p =
# 0.665; Std. beta = -0.22, 95% CI [-1.25, 0.81])
#   - The effect of liver mass × lake [Shirley] is
# statistically non-significant and positive (beta =
# 8.85e-04, 95% CI [-0.04, 0.04], t(41) = 0.05, p =
# 0.961; Std. beta = 0.06, 95% CI [-2.51, 2.64])
# 
# Standardized parameters were obtained by fitting
# the model on a standardized version of the dataset.
# 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation.

check_model(lm_lvr_lake)
forest_model(lm_lvr_lake)

# Generate predictions using the correct model
predictions_lvr_lake <- OD_sheet %>%
  mutate(predicted_lvr_lake = predict(lm_lvr_lake, newdata = .))
# Create the plot
ggplot(OD_sheet, aes(x = lake, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictions_lvr_lake, aes(x = lake, y = predicted_lvr_lake, color = lake)) +
  labs(title = "OD Liver by Lake",
       x = "Lake",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()



#### OD liver vs liver mass and lake ####
# lake and age

lm_lvr_lvrmass_lake <- lm(OD_liver ~ liver_mass + lake, data = OD_sheet)
summary(lm_lvr_lvrmass_lake)
# Call:
#   lm(formula = OD_liver ~ liver_mass + lake, data = OD_sheet)
# 
# Residuals: 
#   Min       1Q   Median       3Q      Max 
# -0.49414 -0.05824  0.01969  0.10181  0.23434 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.954244   0.086653  11.012 3.14e-14 ***
#   liver_mass  -0.001032   0.003278  -0.315    0.754    
# lakeLOTR    -0.000921   0.082320  -0.011    0.991    
# lakeOpeongo -0.028237   0.064224  -0.440    0.662    
# lakeShirley -0.054055   0.082527  -0.655    0.516    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1568 on 44 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.02205,	Adjusted R-squared:  -0.06685 
# F-statistic: 0.2481 on 4 and 44 DF,  p-value: 0.9093
report(lm_lvr_lvrmass_lake)

check_model(lm_lvr_lvrmass_lake)
forest_model(lm_lvr_lvrmass_lake)

library(knitr)
library(kableExtra)
# Tidy the model output
tidy_model_lm_lvr_lvrmass_lake <- tidy(lm_lvr_lvrmass_lake)

# Create a table
tidy_model_lm_lvr_lvrmass_lake %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

predictions <- OD_sheet %>%
  mutate(predicted_lm = predict(lm_lvr_lvrmass_lake, newdata = .))

ggplot(OD_sheet, aes(x = liver_mass, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = liver_mass, y = predicted_lm, color = lake)) +
  labs(title = "OD vs Liver Mass by Lake",
       x = "Liver Mass (g)",
       y = "Optical Density of Liver Tissue",
       color = "Lake") +
  theme_minimal()
#### OD liver vs liver mass and lake interaction ####
# lake and age

lm_lvr_lvrmass_lakei <- lm(OD_liver ~ liver_mass * lake, data = OD_sheet)
summary(lm_lvr_lvrmass_lakei)
# Call:
#   lm(formula = OD_liver ~ liver_mass * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49380 -0.05975  0.01561  0.09707  0.23493 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)             0.9375490  0.1104852
# liver_mass             -0.0003016  0.0044196
# lakeLOTR               -0.0764219  0.1824825
# lakeOpeongo             0.0427274  0.1759785
# lakeShirley            -0.0471665  0.1591449
# liver_mass:lakeLOTR     0.0132116  0.0212701
# liver_mass:lakeOpeongo -0.0031261  0.0071765
# liver_mass:lakeShirley  0.0008846  0.0179205
# t value Pr(>|t|)    
# (Intercept)              8.486 1.44e-10 ***
#   liver_mass              -0.068    0.946    
# lakeLOTR                -0.419    0.678    
# lakeOpeongo              0.243    0.809    
# lakeShirley             -0.296    0.768    
# liver_mass:lakeLOTR      0.621    0.538    
# liver_mass:lakeOpeongo  -0.436    0.665    
# liver_mass:lakeShirley   0.049    0.961    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1611 on 41 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.03765,	Adjusted R-squared:  -0.1267 
# F-statistic: 0.2292 on 7 and 41 DF,  p-value: 0.9759

report(lm_lvr_lvrmass_lakei)
# We fitted a linear model (estimated using OLS) to
# predict OD_liver with liver_mass and lake (formula:
#                                              OD_liver ~ liver_mass * lake). The model explains a
# statistically not significant and weak proportion
# of variance (R2 = 0.04, F(7, 41) = 0.23, p = 0.976,
#              adj. R2 = -0.13). The model's intercept,
# corresponding to liver_mass = 0 and lake = Hogan,
# is at 0.94 (95% CI [0.71, 1.16], t(41) = 8.49, p <
# .001). Within this model:
# 
#   - The effect of liver mass is statistically
# non-significant and negative (beta = -3.02e-04, 95%
# CI [-9.23e-03, 8.62e-03], t(41) = -0.07, p = 0.946;
# Std. beta = -0.02, 95% CI [-0.66, 0.61])
#   - The effect of lake [LOTR] is statistically
# non-significant and negative (beta = -0.08, 95% CI
# [-0.44, 0.29], t(41) = -0.42, p = 0.678; Std. beta
# = 0.75, 95% CI [-1.62, 3.12])
#   - The effect of lake [Opeongo] is statistically
# non-significant and positive (beta = 0.04, 95% CI
# [-0.31, 0.40], t(41) = 0.24, p = 0.809; Std. beta =
# -0.01, 95% CI [-1.20, 1.17])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -0.05, 95% CI
# [-0.37, 0.27], t(41) = -0.30, p = 0.768; Std. beta
# = -0.23, 95% CI [-2.38, 1.93])
#   - The effect of liver mass × lake [LOTR] is
# statistically non-significant and positive (beta =
# 0.01, 95% CI [-0.03, 0.06], t(41) = 0.62, p =
# 0.538; Std. beta = 0.94, 95% CI [-2.12, 4.00])
#   - The effect of liver mass × lake [Opeongo] is
# statistically non-significant and negative (beta =
# -3.13e-03, 95% CI [-0.02, 0.01], t(41) = -0.44, p =
# 0.665; Std. beta = -0.22, 95% CI [-1.25, 0.81])
#   - The effect of liver mass × lake [Shirley] is
# statistically non-significant and positive (beta =
# 8.85e-04, 95% CI [-0.04, 0.04], t(41) = 0.05, p =
# 0.961; Std. beta = 0.06, 95% CI [-2.51, 2.64])
# 
# Standardized parameters were obtained by fitting
# the model on a standardized version of the dataset.
# 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation.
check_model(lm_lvr_lvrmass_lakei)
forest_model(lm_lvr_lvrmass_lakei)

# Generate predictions
predictionslm_lvr_lvrmass_lakei <- OD_sheet %>%
  mutate(predicted_lm_lvr_lvrmass_lakei = predict(lm_lvr_lvrmass_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = liver_mass, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_lvr_lvrmass_lakei, aes(x = liver_mass, y = predicted_lm_lvr_lvrmass_lakei, color = lake)) +
  labs(title = "OD Liver vs Liver Mass by Lake",
       x = "Liver Mass (g)",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()

### OD liver vs gonad mass and lake interaction ####
# lake and age

lm_lvr_gdmass_lakei <- lm(OD_liver ~ gonad_mass * lake, data = OD_sheet)
summary(lm_lvr_gdmass_lakei)
# Call:
#   lm(formula = OD_liver ~ gonad_mass * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48918 -0.05367  0.01217  0.08924  0.24662 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)             9.264e-01  6.467e-02
# gonad_mass              3.217e-05  3.544e-04
# lakeLOTR                3.587e-02  9.156e-02
# lakeOpeongo            -2.545e-02  1.052e-01
# lakeShirley            -5.184e-02  1.091e-01
# gonad_mass:lakeLOTR    -6.865e-04  1.889e-03
# gonad_mass:lakeOpeongo -1.561e-05  7.617e-04
# gonad_mass:lakeShirley  8.841e-04  3.600e-03
# t value Pr(>|t|)    
# (Intercept)             14.327   <2e-16 ***
#   gonad_mass               0.091    0.928    
# lakeLOTR                 0.392    0.697    
# lakeOpeongo             -0.242    0.810    
# lakeShirley             -0.475    0.637    
# gonad_mass:lakeLOTR     -0.363    0.718    
# gonad_mass:lakeOpeongo  -0.020    0.984    
# gonad_mass:lakeShirley   0.246    0.807    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1622 on 41 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.02458,	Adjusted R-squared:  -0.142 
# F-statistic: 0.1476 on 7 and 41 DF,  p-value: 0.9934

report(lm_lvr_gdmass_lakei)
# We fitted a linear model (estimated using OLS) to
# predict OD_liver with gonad_mass and lake (formula:
#                                              OD_liver ~ gonad_mass * lake). The model explains a
# statistically not significant and weak proportion
# of variance (R2 = 0.02, F(7, 41) = 0.15, p = 0.993,
#              adj. R2 = -0.14). The model's intercept,
# corresponding to gonad_mass = 0 and lake = Hogan,
# is at 0.93 (95% CI [0.80, 1.06], t(41) = 14.33, p <
# .001). Within this model:
# 
#   - The effect of gonad mass is statistically
# non-significant and positive (beta = 3.22e-05, 95%
# CI [-6.84e-04, 7.48e-04], t(41) = 0.09, p = 0.928;
# Std. beta = 0.02, 95% CI [-0.41, 0.44])
#   - The effect of lake [LOTR] is statistically
# non-significant and positive (beta = 0.04, 95% CI
# [-0.15, 0.22], t(41) = 0.39, p = 0.697; Std. beta =
# -0.07, 95% CI [-1.50, 1.35])
#   - The effect of lake [Opeongo] is statistically
# non-significant and negative (beta = -0.03, 95% CI
# [-0.24, 0.19], t(41) = -0.24, p = 0.810; Std. beta
# = -0.17, 95% CI [-1.15, 0.80])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -0.05, 95% CI
# [-0.27, 0.17], t(41) = -0.47, p = 0.637; Std. beta
# = 0.06, 95% CI [-2.38, 2.49])
#   - The effect of gonad mass × lake [LOTR] is
# statistically non-significant and negative (beta =
# -6.86e-04, 95% CI [-4.50e-03, 3.13e-03], t(41) =
# -0.36, p = 0.718; Std. beta = -0.41, 95% CI [-2.68,
# 1.86])
#   - The effect of gonad mass × lake [Opeongo] is
# statistically non-significant and negative (beta =
# -1.56e-05, 95% CI [-1.55e-03, 1.52e-03], t(41) =
# -0.02, p = 0.984; Std. beta = -9.29e-03, 95% CI
# [-0.92, 0.91])
#   - The effect of gonad mass × lake [Shirley] is
# statistically non-significant and positive (beta =
# 8.84e-04, 95% CI [-6.39e-03, 8.15e-03], t(41) =
# 0.25, p = 0.807; Std. beta = 0.53, 95% CI [-3.80,
# 4.85])
# 
# Standardized parameters were obtained by fitting
# the model on a standardized version of the dataset.
# 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation.
check_model(lm_lvr_gdmass_lakei)
forest_model(lm_lvr_gdmass_lakei)

# Generate predictions
predictionslm_lvr_gdmass_lakei <- OD_sheet %>%
  mutate(predicted_lm_lvr_gdmass_lakei = predict(lm_lvr_gdmass_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = gonad_mass, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_lvr_gdmass_lakei, aes(x = gonad_mass, y = predicted_lm_lvr_gdmass_lakei, color = lake)) +
  labs(title = "OD Liver vs Gonad Mass by Lake",
       x = "Gonad Mass (g)",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()



#### Od liver vs liver mass and sex ####

lm_lvr_lvrmass_sex <- lm(OD_liver ~ liver_mass + sex, data = OD_sheet)
summary(lm_lvr_lvrmass_sex)
# Call:
#   lm(formula = OD_liver ~ liver_mass + sex, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.51107 -0.03941  0.01029  0.09147  0.24464 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.959832   0.058748  16.338   <2e-16 ***
#   liver_mass  -0.001262   0.002337  -0.540    0.592    
# sexmale     -0.041367   0.051109  -0.809    0.423    
# sexunknown    0.022195   0.160770   0.138    0.891    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1554 on 45 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.01742,	Adjusted R-squared:  -0.04809 
# F-statistic: 0.2659 on 3 and 45 DF,  p-value: 0.8496
report(lm_lvr_lvrmass_lake)

check_model(lm_lvr_lvrmass_lake)
forest_model(lm_lvr_lvrmass_lake)

library(knitr)
library(kableExtra)
# Tidy the model output
tidy_model_lm_lvr_lvrmass_lake <- tidy(lm_lvr_lvrmass_lake)

# Create a table
tidy_model_lm_lvr_lvrmass_lake %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

predictions <- OD_sheet %>%
  mutate(predicted_lm = predict(lm_lvr_lvrmass_lake, newdata = .))

ggplot(OD_sheet, aes(x = liver_mass, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictions, aes(x = liver_mass, y = predicted_lm, color = lake)) +
  labs(title = "OD vs Liver Mass by Lake",
       x = "Liver Mass (g)",
       y = "Optical Density of Liver Tissue",
       color = "Lake") +
  theme_minimal()
# lake and age
                                                        
lm_lvr_lake_age <- lm(OD_liver ~ lake + Age, data = OD_sheet)
summary(lm_lvr_lake_age)
# Call:
#   lm(formula = OD_liver ~ lake + Age, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48955 -0.05233  0.01598  0.09219  0.24616 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.975060   0.084591  11.527 1.37e-14 ***
#   lakeLOTR    -0.003216   0.073958  -0.043    0.966    
# lakeOpeongo -0.055273   0.066547  -0.831    0.411    
# lakeShirley -0.046624   0.063379  -0.736    0.466    
# Age         -0.002991   0.004894  -0.611    0.544    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1562 on 42 degrees of freedom
# (19 observations deleted due to missingness)
# Multiple R-squared:  0.03994,	Adjusted R-squared:  -0.0515 
# F-statistic: 0.4368 on 4 and 42 DF,  p-value: 0.7813                               
 report(lm_lvr_lake_age)
                   
check_model(lm_lvr_lake_age)
forest_model(lm_lvr_lake_age)

library(knitr)
library(kableExtra)
# Tidy the model output
tidy_model_lm_lvr_lake_age <- tidy(lm_lvr_lake_age)
                            
# Create a table
tidy_model_lm_lvr_lake_age %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)
# 
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	0.9750599	0.0845911	11.5267407	0.0000000
# lakeLOTR	-0.0032157	0.0739584	-0.0434796	0.9655253
# lakeOpeongo	-0.0552732	0.0665470	-0.8305891	0.4109026
# lakeShirley	-0.0466244	0.0633795	-0.7356382	0.4660395
# Age	-0.0029911	0.0048936	-0.6112238	0.5443456                            
             
                            
predictions <- OD_sheet %>%
mutate(predicted_lm = predict(lm_hrt_lake_age, newdata = .))
                            
ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
geom_point() +
geom_line(data = predictions, aes(x = Age, y = predicted_lm, color = lake)) +
labs(title = "OD vs Age by Lake",
x = "Age (years)",
y = "Optical Density of Heart Tissue",
color = "Lake") +
theme_minimal()
                            
                            
lm_hrt_lake_age2 <- lm(OD_heart ~ lake * Age, data = OD_sheet)
summary(lm_hrt_lake_age2)
                            
                           
                            
predictions2 <- OD_sheet %>%
mutate(predicted_lm2 = predict(lm_hrt_lake_age2, newdata = .))
                            
ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
geom_point() +
geom_line(data = predictions2, aes(x = Age, y = predicted_lm2, color = lake)) +
labs(title = "OD vs Age by Lake",
                                   x = "Age (years)", 
                                   y = "Optical Density of Heart Tissue",
color = "Lake") +
theme_minimal()
report(lm_hrt_lake_age2)
        
                            
# Pearson's product-moment correlation
 # data:  OD_sheet$Age and OD_sheet$OD_liver
# t = -1.0041, df = 44, p-value = 0.3208
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4216443  0.1470201
# sample estimates:
#        cor 
# -0.1496657 
# correlation
                            
# Perform the correlation test
cor.test(OD_sheet$Age, OD_sheet$OD_liver)
                            
                            # Remove rows with missing values
                            OD_sheet_clean2 <- OD_sheet %>%
                              filter(!is.na(Age) & !is.na(OD_liver))
                            
                            # Calculate correlation and p-values within each lake
                            results2 <- OD_sheet %>%
                              group_by(lake) %>%
                              summarise(
                                correlation = cor(Age, OD_liver),
                                p_value = cor.test(Age, OD_liver)$p.value
                              )
                            
                            print(results2)
                            # A tibble: 4 × 3
                            # lake    correlation p_value
                            # <chr>         <dbl>   <dbl>
                            #   1 Hogan            NA   0.110
                            # 2 LOTR             NA   0.446
                            # 3 Opeongo          NA   0.194
                            # 4 Shirley          NA   0.714 
                            
                            # linear model
                            lm_lvr_lakes <- lm(OD_liver ~ lake, data = OD_sheet)
                            summary(lm_lvr_lakes)
                            # Call:
                            #   lm(formula = OD_liver ~ lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.49392 -0.05015  0.01450  0.09908  0.24536 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.93065    0.04304  21.623   <2e-16 ***
                            #   lakeLOTR     0.01585    0.06212   0.255    0.800    
                            # lakeOpeongo -0.02802    0.06357  -0.441    0.662    
                            # lakeShirley -0.03673    0.06087  -0.603    0.549    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1552 on 45 degrees of freedom
                            # (17 observations deleted due to missingness)
                            # Multiple R-squared:  0.01985,	Adjusted R-squared:  -0.04549 
                            # F-statistic: 0.3038 on 3 and 45 DF,  p-value: 0.8225
                            
                            # linear model
                            lm_lvr_lakes_age <- lm(OD_liver ~ lake + Age, data = OD_sheet)
                            summary(lm_lvr_lakes_age)
                            
                            # Call:
                            #   lm(formula = OD_liver ~ lake + Age, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.48753 -0.05328  0.01335  0.10597  0.22347 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.995606   0.085957  11.583 1.66e-14 ***
                            #   lakeLOTR     0.002033   0.077815   0.026    0.979    
                            # lakeOpeongo -0.040072   0.067051  -0.598    0.553    
                            # lakeShirley -0.051202   0.064953  -0.788    0.435    
                            # Age         -0.004375   0.004956  -0.883    0.383    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1602 on 41 degrees of freedom
                            # (20 observations deleted due to missingness)
                            # Multiple R-squared:  0.046,	Adjusted R-squared:  -0.04707 
                            # F-statistic: 0.4942 on 4 and 41 DF,  p-value: 0.74
                            
                            # linear model
                            lm_lvr_age <- lm(OD_liver ~ Age, data = OD_sheet)
                            summary(lm_lvr_age)
                            # 
                            # 
                            # 
                            # Call:
                            #   lm(formula = OD_liver ~ Age, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.51512 -0.04396  0.03138  0.08370  0.20428 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.971518   0.056243  17.273   <2e-16 ***
                            #   Age         -0.004338   0.004321  -1.004    0.321    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1566 on 44 degrees of freedom
                            # (20 observations deleted due to missingness)
                            # Multiple R-squared:  0.0224,	Adjusted R-squared:  0.0001816 
                            # F-statistic: 1.008 on 1 and 44 DF,  p-value: 0.3208
                            
                            #### correlation liver mass and OD ####
                            # Perform correlation test
                            correlation_OD_liver <- cor.test(OD_sheet$OD_liver, 
                                                             OD_sheet$liver_mass)
                            
                            # Print the test results
                            print(correlation_OD_liver)
                            # Pearson's product-moment correlation
                            # 
                            # data:  OD_sheet$OD_liver and OD_sheet$liver_mass
                            # t = -0.21081, df = 47, p-value = 0.8339
                            # alternative hypothesis: true correlation is not equal to 0
                            # 95 percent confidence interval:
                            #  -0.3092592  0.2526447
                            # sample estimates:
                            #         cor 
                            # -0.03073546 
                            
                            
                            
                            # Function to perform correlation test
                            cor_test <- function(df) {
                              cor.test(df$OD_liver, df$liver_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test))
                            
                            # Extract and print results
                            correlation_results %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
                            # A tibble: 4 × 9
                            # Groups:   lake [4]
                            # lake  estimate statistic p.value parameter conf.low conf.high method
                            # <fct>    <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr> 
                            #   1 Opeo…   0.110      0.413   0.686        14   -0.408    0.574  Pears…
                            # 2 LOTR   -0.409     -1.74    0.103        15   -0.744    0.0889 Pears…
                            # 3 Shir…   0.290      1.13    0.276        14   -0.240    0.687  Pears…
                            # 4 Hogan  -0.0706    -0.274   0.788        15   -0.533    0.424  Pears…
                            # # ℹ 1 more variable: alternative <chr>
                            
                            
                            #### anova liver mass OD by lake ####
                            anova_livermass_lake <- aov(liver_mass ~ lake, data = OD_sheet)
                            summary(anova_livermass_lake)
                            # Df Sum Sq Mean Sq F value   Pr(>F)    
                            # lake         3   4538  1512.6   32.97 7.29e-13 ***
                            #   Residuals   62   2845    45.9                     
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
                            TukeyHSD(anova_livermass_lake)
                            # Tukey multiple comparisons of means
                            # 95% family-wise confidence level
                            # 
                            # Fit: aov(formula = liver_mass ~ lake, data = OD_sheet)
                            # 
                            # $lake
                            # diff        lwr        upr     p adj
                            # LOTR-Hogan      -16.658647 -22.792626 -10.524668 0.0000000
                            # Opeongo-Hogan    -1.478713  -7.707799   4.750372 0.9230915
                            # Shirley-Hogan   -17.869338 -24.098424 -11.640253 0.0000000
                            # Opeongo-LOTR     15.179934   8.950848  21.409019 0.0000001
                            # Shirley-LOTR     -1.210691  -7.439777   5.018394 0.9556576
                            # Shirley-Opeongo -16.390625 -22.713386 -10.067864 0.0000000
                            
                            
                            #### liver mass linear models ####
                            lm_livermass_lakeage <- lm(liver_mass ~ Age + lake,data = OD_sheet)
                            summary(lm_livermass_lakeage)
                            # Call:
                            #   lm(formula = liver_mass ~ Age + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -15.1771  -2.4621  -0.2071   3.3081  14.4456 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  16.7144     2.7466   6.085 1.10e-07 ***
                            #   Age           0.5493     0.1644   3.341  0.00149 ** 
                            #   lakeLOTR    -14.1690     2.5407  -5.577 7.34e-07 ***
                            #   lakeOpeongo  -0.3453     2.2711  -0.152  0.87971    
                            # lakeShirley -17.2797     2.2220  -7.777 1.80e-10 ***
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 6.359 on 56 degrees of freedom
                            # (5 observations deleted due to missingness)
                            # Multiple R-squared:  0.6838,	Adjusted R-squared:  0.6613 
                            # F-statistic: 30.28 on 4 and 56 DF,  p-value: 2.002e-13
                            
                            lm_OD_livermass_lakeage <- lm(OD_liver ~ liver_mass + Age + lake, data = OD_sheet)
                            summary(lm_OD_livermass_lakeage)
                            # Call:
                            #   lm(formula = OD_liver ~ liver_mass + Age + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.48687 -0.05433  0.01315  0.10489  0.22285 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.988687   0.099776   9.909 2.51e-12 ***
                            #   liver_mass   0.000545   0.003847   0.142    0.888    
                            # Age         -0.004748   0.005666  -0.838    0.407    
                            # lakeLOTR     0.008458   0.090888   0.093    0.926    
                            # lakeOpeongo -0.040984   0.068171  -0.601    0.551    
                            # lakeShirley -0.043285   0.086286  -0.502    0.619    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1622 on 40 degrees of freedom
                            # (20 observations deleted due to missingness)
                            # Multiple R-squared:  0.04648,	Adjusted R-squared:  -0.07271 
                            # F-statistic: 0.3899 on 5 and 40 DF,  p-value: 0.8528
                            
                            lm_OD_livermass_lake <- lm(OD_liver ~ liver_mass + lake, data = OD_sheet)
                            summary(lm_OD_livermass_lake)
                            # Call:
                            #   lm(formula = OD_liver ~ liver_mass + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.49414 -0.05824  0.01969  0.10181  0.23434 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.954244   0.086653  11.012 3.14e-14 ***
                            #   liver_mass  -0.001032   0.003278  -0.315    0.754    
                            # lakeLOTR    -0.000921   0.082320  -0.011    0.991    
                            # lakeOpeongo -0.028237   0.064224  -0.440    0.662    
                            # lakeShirley -0.054055   0.082527  -0.655    0.516    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1568 on 44 degrees of freedom
                            # (17 observations deleted due to missingness)
                            # Multiple R-squared:  0.02205,	Adjusted R-squared:  -0.06685 
                            # F-statistic: 0.2481 on 4 and 44 DF,  p-value: 0.9093
                            
                            lm_OD_livermass <- lm(OD_liver ~ liver_mass, data = OD_sheet)
                            summary(lm_OD_livermass)
                            # Call:
                            #   lm(formula = OD_liver ~ liver_mass, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.52218 -0.04268  0.03023  0.07836  0.22846 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.924713   0.036720  25.183   <2e-16 ***
                            #   liver_mass  -0.000432   0.002049  -0.211    0.834    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1533 on 47 degrees of freedom
                            # (17 observations deleted due to missingness)
                            # Multiple R-squared:  0.0009447,	Adjusted R-squared:  -0.02031 
                            # F-statistic: 0.04444 on 1 and 47 DF,  p-value: 0.8339
                            library(dplyr)
                            library(purrr)
                            library(broom)
                            
                            
                            
                            # Fit linear models for each lake
                            models <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(model = map(data, ~ lm(OD_liver ~ liver_mass, data = .x)))
                            
                            # Extract model summaries
                            model_summaries <- models %>%
                              mutate(summary = map(model, summary),
                                     tidy_summary = map(model, tidy),
                                     glance_summary = map(model, glance)) %>%
                              select(lake, summary, tidy_summary, glance_summary)
                            
                            # Print model summaries
                            print(model_summaries)
                            # did not worked
                            
                            
                            #### Fork Length ####
                            y.lab <- expression("Fork Length")
                            x.lab <- expression("Lakes")
                            ggplot(data = SBB_OD1_data, mapping = aes(x = Llake, y = Fork_length)) + 
                              geom_boxplot(aes(fill=Llake), show.legend = FALSE) +
                              theme_classic()+
                              labs(x = x.lab, y = y.lab)
                            
#### linear model OD heartvs age ####
lm_hrt_age <- aov(OD_heart ~ Age, data = OD_sheet)
summary(lm_hrt_age)
# Df Sum Sq Mean Sq F value Pr(>F)
# Age          1 0.0010 0.00098   0.087  0.769
# Residuals   61 0.6899 0.01131               
# 3 observations deleted due to missingness              
# 5 observations deleted due to missingness
                            
#### linear model OD heart vs age & lake ####
lm_hrt_age_lake <- aov(OD_heart ~ Age + lake, data = OD_sheet)
summary(lm_hrt_age_lake)
# Df Sum Sq Mean Sq F value Pr(>F)
# Age          1 0.0010 0.00098   0.092  0.763
# lake         3 0.0690 0.02302   2.150  0.104
# Residuals   58 0.6208 0.01070               
# 3 observations deleted due to missingness

#### anova model OD heart vs age & lake interaction ####
lm_hrt_age_lakei <- aov(OD_heart ~ Age * lake, data = OD_sheet)
summary(lm_hrt_age_lakei)
# Df Sum Sq Mean Sq F value Pr(>F)
# Age          1 0.0010 0.00098   0.091  0.764
# lake         3 0.0690 0.02302   2.146  0.105
# Age:lake     3 0.0311 0.01036   0.966  0.415
# Residuals   55 0.5897 0.01072               
# 3 observations deleted due to missingness
forest_model(lm_hrt_age_lakei)# did not work
# Generate predictions
predictionslm_hrt_age_lakei <- OD_sheet %>%
  mutate(predicted_lm_hrt_age_lakei = predict(lm_hrt_age_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_age_lakei, aes(x = Age, y = predicted_lm_hrt_age_lakei, color = lake)) +
  labs(title = "OD Heart vs Age by Lake",
       x = "Age (y)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()



#### linear model OD heart vs age & lake interaction ####
lm_hrt_age_lakei <- lm(OD_heart ~ Age * lake, data = OD_sheet)
summary(lm_hrt_age_lakei)
# Call:
#   lm(formula = OD_heart ~ Age * lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211110 -0.069854 -0.005495  0.059582  0.259517 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)      1.0567061  0.0646920  16.334
# Age              0.0008071  0.0043128   0.187
# lakeLOTR         0.0574622  0.0902966   0.636
# lakeOpeongo     -0.0446913  0.0908729  -0.492
# lakeShirley     -0.0567412  0.1045402  -0.543
# Age:lakeLOTR    -0.0146511  0.0087726  -1.670
# Age:lakeOpeongo -0.0027033  0.0064547  -0.419
# Age:lakeShirley -0.0011806  0.0074807  -0.158
# Pr(>|t|)    
# (Intercept)       <2e-16 ***
#   Age                0.852    
# lakeLOTR           0.527    
# lakeOpeongo        0.625    
# lakeShirley        0.589    
# Age:lakeLOTR       0.101    
# Age:lakeOpeongo    0.677    
# Age:lakeShirley    0.875    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1035 on 55 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1464,	Adjusted R-squared:  0.03771 
# F-statistic: 1.347 on 7 and 55 DF,  p-value: 0.2464

forest_model(lm_hrt_age_lakei)# did not work
# Generate predictions
predictionslm_hrt_age_lakei <- OD_sheet %>%
  mutate(predicted_lm_hrt_age_lakei = predict(lm_hrt_age_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = Age, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_age_lakei, aes(x = Age, y = predicted_lm_hrt_age_lakei, color = lake)) +
  labs(title = "Optical Density Heart vs Age by Lake",
       x = "Age (y)",
       y = "Optical Density of Lipofuscin's Heart",
       color = "Lake") +
  theme_minimal()

#### linear model OD heart vs ventricle mass & lake interaction ####
lm_hrt_vm_lakei <- lm(OD_heart ~ ventricle_mass * lake, data = OD_sheet)
summary(lm_hrt_vm_lakei)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass * lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.210356 -0.062448 -0.005984  0.051472  0.219714 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                 1.08719    0.07001
# ventricle_mass             -0.01231    0.04184
# lakeLOTR                    0.14528    0.14601
# lakeOpeongo                -0.11620    0.09903
# lakeShirley                -0.17941    0.09903
# ventricle_mass:lakeLOTR    -0.39486    0.22982
# ventricle_mass:lakeOpeongo  0.02173    0.05625
# ventricle_mass:lakeShirley  0.18508    0.13591
# t value Pr(>|t|)    
# (Intercept)                 15.528   <2e-16 ***
#   ventricle_mass              -0.294   0.7696    
# lakeLOTR                     0.995   0.3239    
# lakeOpeongo                 -1.173   0.2454    
# lakeShirley                 -1.812   0.0752 .  
# ventricle_mass:lakeLOTR     -1.718   0.0911 .  
# ventricle_mass:lakeOpeongo   0.386   0.7007    
# ventricle_mass:lakeShirley   1.362   0.1785    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09996 on 58 degrees of freedom
# Multiple R-squared:  0.1706,	Adjusted R-squared:  0.07049 
# F-statistic: 1.704 on 7 and 58 DF,  p-value: 0.1259

forest_model(lm_hrt_vm_lakei)# did not work
# Generate predictions
predictionslm_hrt_vm_lakei <- OD_sheet %>%
  mutate(predicted_lm_hrt_vm_lakei = predict(lm_hrt_vm_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = ventricle_mass, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_vm_lakei, aes(x = ventricle_mass, y = predicted_lm_hrt_vm_lakei, color = lake)) +
  labs(title = "OD Heart vs Ventricle Mass by Lake",
       x = "Ventricle Mass (g)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()

#### linear model OD heart vs gonad mass & lake interaction ####
lm_hrt_gm_lakei <- lm(OD_heart ~ gonad_mass * lake, data = OD_sheet)
summary(lm_hrt_gm_lakei)
# Call:
#   lm(formula = OD_heart ~ gonad_mass * lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.204974 -0.057707 -0.006784  0.050925  0.275879 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)             1.079e+00  3.838e-02
# gonad_mass             -7.741e-05  2.071e-04
# lakeLOTR               -5.585e-02  5.115e-02
# lakeOpeongo            -7.251e-02  5.898e-02
# lakeShirley            -1.031e-01  6.316e-02
# gonad_mass:lakeLOTR    -4.362e-04  7.126e-04
# gonad_mass:lakeOpeongo -9.416e-05  3.923e-04
# gonad_mass:lakeShirley  9.173e-04  1.855e-03
# t value Pr(>|t|)    
# (Intercept)             28.107   <2e-16 ***
#   gonad_mass              -0.374    0.710    
# lakeLOTR                -1.092    0.279    
# lakeOpeongo             -1.229    0.224    
# lakeShirley             -1.633    0.108    
# gonad_mass:lakeLOTR     -0.612    0.543    
# gonad_mass:lakeOpeongo  -0.240    0.811    
# gonad_mass:lakeShirley   0.494    0.623    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1033 on 58 degrees of freedom
# Multiple R-squared:  0.1145,	Adjusted R-squared:  0.007643 
# F-statistic: 1.072 on 7 and 58 DF,  p-value: 0.3933

forest_model(lm_hrt_gm_lakei)
# Generate predictions
predictionslm_hrt_gm_lakei <- OD_sheet %>%
  mutate(predicted_lm_hrt_gm_lakei = predict(lm_hrt_gm_lakei, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = gonad_mass, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_gm_lakei, aes(x = gonad_mass, y = predicted_lm_hrt_gm_lakei, color = lake)) +
  labs(title = "OD Heart vs Gonad Mass by Lake",
       x = "Gonad Mass (g)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()
                            
#### anova OD heart vs sex ####
                            anova_hrt_sex <- aov(OD ~ Sex, data = SBB_OD1)
                            summary(anova_hrt_sex)
                            # Df Sum Sq  Mean Sq F value Pr(>F)
                            # Sex          3 0.0179 0.005964   0.543  0.655
                            # Residuals   62 0.6809 0.010982 
                            
#### linear model OD heart vs age & lake ####
                            lm_hrt_age_lake_sex <- aov(OD ~ age + Lake + Sex, data = SBB_OD1)
summary(lm_hrt_age_lake_sex)
# Df Sum Sq  Mean Sq F value Pr(>F)
                            # age          1 0.0008 0.000784   0.068  0.796
                            # Lake         3 0.0701 0.023365   2.020  0.122
                            # Sex          3 0.0066 0.002215   0.192  0.902
                            # Residuals   53 0.6129 0.011565               
                            # 5 observations deleted due to missingness
                            
                            #### Standard Mass ####
                            y.lab <- expression("Standard Mass")
                            x.lab <- expression("Lakes")
                            ggplot(data = SBB_OD1_data, mapping = aes(x = Llake, y = Standard_mass)) + 
                              geom_boxplot(aes(fill=Llake), show.legend = FALSE) +
                              theme_classic()+
                              labs(x = x.lab, y = y.lab)
#### Organ mass linear model ####
#### liver mass ####
# Fit the linear model
model_lvrms_mass_lk <- lm(liver_mass ~ 
                      mass_post_resp_g + lake, data = OD_sheet)
                            
# Print the summary of the model
summary(model_lvrms_mass_lk)
# 
# Call:
#   lm(formula   = liver_mass ~ mass_post_resp_g + lake, data = OD_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -9.9913 -1.7912 -0.2295  2.2686 12.3093 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       1.185017   2.443510   0.485    0.629    
# mass_post_resp_g  0.012043   0.001162  10.364 4.45e-15 ***
#   lakeLOTR         -1.896047   2.003985  -0.946    0.348    
# lakeOpeongo       0.062791   1.439235   0.044    0.965    
# lakeShirley      -1.921388   2.101650  -0.914    0.364    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.11 on 61 degrees of freedom
# Multiple R-squared:  0.8604,	Adjusted R-squared:  0.8513 
# F-statistic: 94.02 on 4 and 61 DF,  p-value: < 2.2e-16


report(model_lvrms_mass_lk)
# We fitted a linear model (estimated using OLS) to
# predict liver_mass with mass_post_resp_g and lake
# (formula: liver_mass ~ mass_post_resp_g + lake). The
# model explains a statistically significant and
# substantial proportion of variance (R2 = 0.86, F(4,
#                                                  61) = 94.02, p < .001, adj. R2 = 0.85). The model's
# intercept, corresponding to mass_post_resp_g = 0 and
# lake = Hogan, is at 1.19 (95% CI [-3.70, 6.07], t(61)
# = 0.48, p = 0.629). Within this model:
# 
#   - The effect of mass post resp g is statistically
# significant and positive (beta = 0.01, 95% CI
# [9.72e-03, 0.01], t(61) = 10.36, p < .001; Std. beta =
# 0.85, 95% CI [0.69, 1.02])
#   - The effect of lake [LOTR] is statistically
# non-significant and negative (beta = -1.90, 95% CI
# [-5.90, 2.11], t(61) = -0.95, p = 0.348; Std. beta =
# -0.18, 95% CI [-0.55, 0.20])
#   - The effect of lake [Opeongo] is statistically
# non-significant and positive (beta = 0.06, 95% CI
# [-2.82, 2.94], t(61) = 0.04, p = 0.965; Std. beta =
# 5.89e-03, 95% CI [-0.26, 0.28])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -1.92, 95% CI
# [-6.12, 2.28], t(61) = -0.91, p = 0.364; Std. beta =
# -0.18, 95% CI [-0.57, 0.21])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.

check_model(model_lvrms_mass_lk) # good
forest_model(model_lvrms_mass_lk)
                            
# Tidy the model output
tidy_model_lvrms_mass_lk <- tidy(model_lvrms_mass_lk)
                            
# Create a table
tidy_model_lvrms_mass_lk %>%
kable("html", caption = "Linear Regression Results") %>%
kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	1.1850173	2.443510	0.4849652	0.6294381
# mass_post_resp_g	0.0120430	0.001162	10.3644368	0.0000000
# lakeLOTR	-1.8960473	2.003985	-0.9461385	0.3478109
# lakeOpeongo	0.0627915	1.439235	0.0436284	0.9653431
# lakeShirley	-1.9213879	2.101650	-0.9142283	0.3641966       
# Fit linear models for each lake and extract formulas
formulas10 <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(liver_mass ~ mass_post_resp_g, data = .)) %>%
  mutate(formula = paste0("mass_post_resp_g = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * mass_post_resp_g"))

# Extract the lake and formula into a separate data frame
formulas_df <- formulas10 %>%
  select(lake, formula)

# Create the base plot
lm10 <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Liver Mass vs Mass",
       x = "Mass (g)",
       y = "Liver Mass(g)",
       color = "Lake")

# Add formulas to each facet
lm10 <- lm10 + geom_text(data = formulas_df, aes(x = Inf, y = Inf, label = formula), 
                         hjust = 1.1, vjust = 2, size = 3, color = "black", inherit.aes = FALSE)

# Print the plot
print(lm10)


predictionsPCRhlm <- df2 %>% mutate(predicted_PCRhlm = predict(lmPCRheart4, newdata = .))

ggplot(df2, aes(x = ventricle_mass, y = CT_heart_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionsPCRhlm, aes(x = ventricle_mass, y = predicted_PCRhlm, color = lake)) + labs(title = "CT Heart Ratio vs Ventricle Mass by Lake", x = "Ventricle Mass (g)", y = "CT Heart Tissue", color = "Lake") +
  theme_minimal()
                      

# Load necessary library
library(ggplot2)

# Create the base plot with actual data points
plot_lvrms_mass_lk <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point(size = 3) +  # Dot plot with data points
  
  # Add fitted lines from the linear model
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +  
  
  # Customize the plot
  theme_classic() +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass Post Resp (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme(legend.position = "bottom")

# Print the plot
print(plot_lvrms_mass_lk)


# Generate predictions
predictionslvrmass <- OD_sheet %>%
  mutate(predicted_lvrmass = predict(model_lvrms_mass_lk, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionslvrmass, aes(x = mass_post_resp_g, y = predicted_lvrmass, color = lake)) +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme_minimal()

#### liver mass interaction ####
# Fit the linear model
model_lvrms_mass_lki <- lm(liver_mass ~ 
                            mass_post_resp_g * lake, data = OD_sheet)

# Print the summary of the model
summary(model_lvrms_mass_lki)
# 
# Call:
#   lm(formula = liver_mass ~ mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10.6472  -1.9918  -0.2081   1.6051  13.2304 

# Coefficients:
#   Estimate Std. Error
# (Intercept)                  -0.3871347  3.1250864
# mass_post_resp_g              0.0128618  0.0015423
# lakeLOTR                     -5.9489672  5.7295922
# lakeOpeongo                   5.2080070  4.7925989
# lakeShirley                  -0.4655139  4.7132969
# mass_post_resp_g:lakeLOTR     0.0072844  0.0069403
# mass_post_resp_g:lakeOpeongo -0.0028127  0.0024820
# mass_post_resp_g:lakeShirley -0.0006236  0.0058710
# t value Pr(>|t|)    
# (Intercept)                   -0.124    0.902    
# mass_post_resp_g               8.339 1.66e-11 ***
#   lakeLOTR                      -1.038    0.303    
# lakeOpeongo                    1.087    0.282    
# lakeShirley                   -0.099    0.922    
# mass_post_resp_g:lakeLOTR      1.050    0.298    
# mass_post_resp_g:lakeOpeongo  -1.133    0.262    
# mass_post_resp_g:lakeShirley  -0.106    0.916    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.118 on 58 degrees of freedom
# Multiple R-squared:  0.8668,	Adjusted R-squared:  0.8507 
# F-statistic: 53.92 on 7 and 58 DF,  p-value: < 2.2e-16
library(report)
report(model_lvrms_mass_lki)
# We fitted a linear model (estimated using OLS) to
# predict liver_mass with mass_post_resp_g and lake
# (formula: liver_mass ~ mass_post_resp_g * lake). The
# model explains a statistically significant and
# substantial proportion of variance (R2 = 0.87, F(7,
#                                                  58) = 53.92, p < .001, adj. R2 = 0.85). The model's
# intercept, corresponding to mass_post_resp_g = 0 and
# lake = Hogan, is at -0.39 (95% CI [-6.64, 5.87], t(58)
# = -0.12, p = 0.902). Within this model:
# 
#   - The effect of mass post resp g is statistically
# significant and positive (beta = 0.01, 95% CI
# [9.77e-03, 0.02], t(58) = 8.34, p < .001; Std. beta =
# 0.91, 95% CI [0.69, 1.13])
#   - The effect of lake [LOTR] is statistically
# non-significant and negative (beta = -5.95, 95% CI
# [-17.42, 5.52], t(58) = -1.04, p = 0.303; Std. beta =
# 0.30, 95% CI [-0.48, 1.08])
#   - The effect of lake [Opeongo] is statistically
# non-significant and positive (beta = 5.21, 95% CI
# [-4.39, 14.80], t(58) = 1.09, p = 0.282; Std. beta =
# 0.16, 95% CI [-0.23, 0.54])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -0.47, 95% CI
# [-9.90, 8.97], t(58) = -0.10, p = 0.922; Std. beta =
# -0.12, 95% CI [-0.89, 0.66])
#   - The effect of mass post resp g × lake [LOTR] is
# statistically non-significant and positive (beta =
# 7.28e-03, 95% CI [-6.61e-03, 0.02], t(58) = 1.05, p =
# 0.298; Std. beta = 0.52, 95% CI [-0.47, 1.50])
#   - The effect of mass post resp g × lake [Opeongo] is
# statistically non-significant and negative (beta =
# -2.81e-03, 95% CI [-7.78e-03, 2.16e-03], t(58) =
# -1.13, p = 0.262; Std. beta = -0.20, 95% CI [-0.55,
# 0.15])
#   - The effect of mass post resp g × lake [Shirley] is
# statistically non-significant and negative (beta =
# -6.24e-04, 95% CI [-0.01, 0.01], t(58) = -0.11, p =
# 0.916; Std. beta = -0.04, 95% CI [-0.88, 0.79])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
using a Wald t-distribution approximation.

# Create the base plot with actual data points
plot_lvrms_mass_lki <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point(size = 3) +  # Dot plot with data points
  
  # Add fitted lines from the linear model
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +  
  
  # Customize the plot
  theme_classic() +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass Post Resp (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme(legend.position = "bottom")

# Print the plot
print(plot_lvrms_mass_lki)


# Generate predictions
predictionslvrmassi <- OD_sheet %>%
  mutate(predicted_lvrmassi = predict(model_lvrms_mass_lki, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionslvrmassi, aes(x = mass_post_resp_g, y = predicted_lvrmassi, color = lake)) +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme_minimal()


#### heart mass interaction ####
# Fit the linear model
model_hrtms_mass_lki <- lm(ventricle_mass ~ 
                            mass_post_resp_g * lake, data = OD_sheet)

# Print the summary of the model
summary(model_hrtms_mass_lki) 

# Call:
#   lm(formula = ventricle_mass ~ mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.36560 -0.08004  0.00312  0.06156  0.52042 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                  -3.257e-02  1.245e-01
# mass_post_resp_g              8.345e-04  6.143e-05
# lakeLOTR                      1.971e-01  2.282e-01
# lakeOpeongo                  -3.670e-01  1.909e-01
# lakeShirley                  -6.091e-02  1.877e-01
# mass_post_resp_g:lakeLOTR    -2.695e-04  2.764e-04
# mass_post_resp_g:lakeOpeongo  3.594e-04  9.886e-05
# mass_post_resp_g:lakeShirley  1.718e-04  2.338e-04
# t value Pr(>|t|)    
# (Intercept)                   -0.262  0.79449    
# mass_post_resp_g              13.584  < 2e-16 ***
#   lakeLOTR                       0.864  0.39137    
# lakeOpeongo                   -1.922  0.05949 .  
# lakeShirley                   -0.324  0.74676    
# mass_post_resp_g:lakeLOTR     -0.975  0.33373    
# mass_post_resp_g:lakeOpeongo   3.636  0.00059 ***
#   mass_post_resp_g:lakeShirley   0.734  0.46563    
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.164 on 58 degrees of freedom
# Multiple R-squared:  0.9548,	Adjusted R-squared:  0.9494 
# F-statistic: 175.1 on 7 and 58 DF,  p-value: < 2.2e-16


check_model(model_hrtms_mass_lki) # not working
update.packages("performance")
report(model_hrtms_mass_lki)
# We fitted a linear model (estimated using OLS) to
# predict ventricle_mass with mass_post_resp_g and lake
# (formula: ventricle_mass ~ mass_post_resp_g * lake).
# The model explains a statistically significant and
# substantial proportion of variance (R2 = 0.95, F(7,
#                                                  58) = 175.10, p < .001, adj. R2 = 0.95). The model's
# intercept, corresponding to mass_post_resp_g = 0 and
# lake = Hogan, is at -0.03 (95% CI [-0.28, 0.22], t(58)
# = -0.26, p = 0.794). Within this model:
# 
#   - The effect of mass post resp g is statistically
# significant and positive (beta = 8.35e-04, 95% CI
# [7.12e-04, 9.57e-04], t(58) = 13.58, p < .001; Std.
# beta = 0.86, 95% CI [0.74, 0.99])
#   - The effect of lake [LOTR] is statistically
# non-significant and positive (beta = 0.20, 95% CI
# [-0.26, 0.65], t(58) = 0.86, p = 0.391; Std. beta =
# -0.19, 95% CI [-0.65, 0.26])
#   - The effect of lake [Opeongo] is statistically
# non-significant and negative (beta = -0.37, 95% CI
# [-0.75, 0.02], t(58) = -1.92, p = 0.059; Std. beta =
# 0.11, 95% CI [-0.11, 0.34])
#   - The effect of lake [Shirley] is statistically
# non-significant and negative (beta = -0.06, 95% CI
# [-0.44, 0.31], t(58) = -0.32, p = 0.747; Std. beta =
# 0.21, 95% CI [-0.24, 0.66])
#   - The effect of mass post resp g × lake [LOTR] is
# statistically non-significant and negative (beta =
# -2.69e-04, 95% CI [-8.23e-04, 2.84e-04], t(58) =
# -0.97, p = 0.334; Std. beta = -0.28, 95% CI [-0.85,
# 0.29])
#   - The effect of mass post resp g × lake [Opeongo] is
# statistically significant and positive (beta =
# 3.59e-04, 95% CI [1.62e-04, 5.57e-04], t(58) = 3.64, p
# < .001; Std. beta = 0.37, 95% CI [0.17, 0.58])
#   - The effect of mass post resp g × lake [Shirley] is
# statistically non-significant and positive (beta =
# 1.72e-04, 95% CI [-2.96e-04, 6.40e-04], t(58) = 0.73,
# p = 0.466; Std. beta = 0.18, 95% CI [-0.31, 0.66])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.
forest_model(model_hrtms_mass_lki)

# Tidy the model output
tidy_model_hrtms_mass_lki <- tidy(model_hrtms_mass_lki)

# Create a table
tidy_model_hrtms_mass_lki %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-0.0325738	0.1244771	-0.2616849	0.7944915
# mass_post_resp_g	0.0008345	0.0000614	13.5844145	0.0000000
# lakeLOTR	0.1970872	0.2282187	0.8635889	0.3913710
# lakeOpeongo	-0.3669560	0.1908968	-1.9222742	0.0594875
# lakeShirley	-0.0609139	0.1877381	-0.3244624	0.7467553
# mass_post_resp_g:lakeLOTR	-0.0002695	0.0002764	-0.9747461	0.3337342
# mass_post_resp_g:lakeOpeongo	0.0003594	0.0000989	3.6355713	0.0005902
# mass_post_resp_g:lakeShirley	0.0001718	0.0002339	0.7344560	0.4656299



# Generate predictions
predictionshrtmassi <- OD_sheet %>%
  mutate(predicted_hrtmassi = predict(model_hrtms_mass_lki, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = ventricle_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionshrtmassi, aes(x = mass_post_resp_g, y = predicted_hrtmassi, color = lake)) +
  labs(title = "Ventricle Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Ventricle Mass (g)",
       color = "Lake") +
  theme_minimal()




#### heart mass ####
# Fit the linear model
model_hrtms_mass_lk <- lm(ventricle_mass ~ 
                            mass_post_resp_g + lake, data = OD_sheet)

# Print the summary of the model
summary(model_hrtms_mass_lk)
# Call:
#   lm(formula = ventricle_mass ~ mass_post_resp_g + lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38551 -0.08272 -0.00463  0.06810  0.49534 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -2.785e-01  1.070e-01  -2.602   0.0116 *  
#   mass_post_resp_g  9.626e-04  5.089e-05  18.914  < 2e-16 ***
#   lakeLOTR          1.671e-01  8.777e-02   1.903   0.0617 .  
# lakeOpeongo       2.935e-01  6.304e-02   4.656 1.79e-05 ***
#   lakeShirley       2.110e-01  9.205e-02   2.293   0.0253 *  
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.18 on 61 degrees of freedom
# Multiple R-squared:  0.9428,	Adjusted R-squared:  0.939 
# F-statistic: 251.2 on 4 and 61 DF,  p-value: < 2.2e-16  

report(model_hrtms_mass_lk)
# We fitted a linear model (estimated using OLS) to predict
# ventricle_mass with mass_post_resp_g and lake (formula: ventricle_mass
#                                                ~ mass_post_resp_g + lake). The model explains a statistically
# significant and substantial proportion of variance (R2 = 0.94, F(4, 61)
#                                                     = 251.16, p < .001, adj. R2 = 0.94). The model's intercept,
# corresponding to mass_post_resp_g = 0 and lake = Hogan, is at -0.28
# (95% CI [-0.49, -0.06], t(61) = -2.60, p = 0.012). Within this model:
# 
#   - The effect of mass post resp g is statistically significant and
# positive (beta = 9.63e-04, 95% CI [8.61e-04, 1.06e-03], t(61) = 18.91,
# p < .001; Std. beta = 1.00, 95% CI [0.89, 1.10])
#   - The effect of lake [LOTR] is statistically non-significant and
# positive (beta = 0.17, 95% CI [-8.46e-03, 0.34], t(61) = 1.90, p =
# 0.062; Std. beta = 0.23, 95% CI [-0.01, 0.47])
#   - The effect of lake [Opeongo] is statistically significant and
# positive (beta = 0.29, 95% CI [0.17, 0.42], t(61) = 4.66, p < .001;
# Std. beta = 0.40, 95% CI [0.23, 0.58])
#   - The effect of lake [Shirley] is statistically significant and
# positive (beta = 0.21, 95% CI [0.03, 0.40], t(61) = 2.29, p = 0.025;
# Std. beta = 0.29, 95% CI [0.04, 0.54])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.. The model's intercept,
# corresponding to mass_post_resp_g = 0 and lake = Hogan, is at -0.28
# (95% CI [-0.49, -0.06], t(61) = -2.60, p = 0.012). Within this model:
# 
#   - The effect of mass post resp g is statistically significant and
# positive (beta = 9.63e-04, 95% CI [8.61e-04, 1.06e-03], t(61) = 18.91,
# p < .001; Std. beta = 1.00, 95% CI [0.89, 1.10])
#   - The effect of lake [LOTR] is statistically non-significant and
# positive (beta = 0.17, 95% CI [-8.46e-03, 0.34], t(61) = 1.90, p =
# 0.062; Std. beta = 0.23, 95% CI [-0.01, 0.47])
#   - The effect of lake [Opeongo] is statistically significant and
# positive (beta = 0.29, 95% CI [0.17, 0.42], t(61) = 4.66, p < .001;
# Std. beta = 0.40, 95% CI [0.23, 0.58])
#   - The effect of lake [Shirley] is statistically significant and
# positive (beta = 0.21, 95% CI [0.03, 0.40], t(61) = 2.29, p = 0.025;
# Std. beta = 0.29, 95% CI [0.04, 0.54])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.

check_model(model_hrtms_mass_lk) # not working
update.packages("performance")

forest_model(model_hrtms_mass_lk)

# Tidy the model output
tidy_model_hrtms_mass_lk <- tidy(model_hrtms_mass_lk)

# Create a table
tidy_model_hrtms_mass_lk %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

     
# Fit linear models for each lake and extract formulas
formulas11 <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(ventricle_mass ~ mass_post_resp_g, data = .)) %>%
  mutate(formula = paste0("mass_post_resp_g = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * mass_post_resp_g"))

# Extract the lake and formula into a separate data frame
formulas_11 <- formulas11 %>%
  select(lake, formula)

# Create the base plot
lm11 <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = ventricle_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Ventricle Mass vs Mass",
       x = "Mass (g)",
       y = "Ventricle Mass(g)",
       color = "Lake")

# Add formulas to each facet
lm11 <- lm11 + geom_text(data = formulas_11, aes(x = Inf, y = Inf, label = formula), 
                         hjust = 1.1, vjust = 2, size = 3, color = "black", inherit.aes = FALSE)

# Print the plot
print(lm11)


# Create the base plot with actual data points
plot_hrtms_mass_lk <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = ventricle_mass, color = lake)) +
  geom_point(size = 3) +  # Dot plot with data points
  
  # Add fitted lines from the linear model
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +  
  
  # Customize the plot
  theme_classic() +
  labs(title = "Ventricle Mass vs Mass by Lake",
       x = "Mass Post Resp (g)",
       y = "Ventricle Mass (g)",
       color = "Lake") +
  theme(legend.position = "bottom")

# Print the plot
print(plot_hrtms_mass_lk)

# Check if the variable exists in the dataset
if("mass_post_resp_g" %in% names(OD_sheet)) {
  # Proceed with your mutate function if the variable exists
  OD_sheet <- OD_sheet %>%
    mutate(predicted_hrtmass = predict(model_hrtms_mass_lk, newdata = .))
} else {
  print("Variable 'mass_post_resp_g' not found in the dataset.")
}


predictionshrtmass <- OD_sheet %>% mutate(predicted_hrtmass = predict(model_hrtms_mass_lk, newdata = .))

ggplot(OD_sheet, aes(x = mass_post_resp_g, y = ventricle_mass, color = lake)) +
  geom_point() + geom_line(data = predictionshrtmass, aes(x = ventricle_mass, y = predicted_hrtmass, color = lake)) + labs(title = "Ventricle Mass vs Mass by Lake", x = "Mass (g)", y = "Ventricle Mass (g)", color = "Lake") +
  theme_minimal()



# Generate predictions
your_dataset$predicted_hrtmass <- predict(model_hrtms_mass_lk, newdata = your_dataset)

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = predicted_hrtmass, color=lake)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Predicted Heart Mass vs. Mass Post Resp",
       x = "Mass Post Resp (g)",
       y = "Predicted Heart Mass")

library(dplyr)
library(ggplot2)

# Generate predictions
predictionshrtmass <- OD_sheet %>%
  mutate(predicted_hrtmass = predict(model_hrtms_mass_lk, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = ventricle_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionshrtmass, aes(x = mass_post_resp_g, y = predicted_hrtmass, color = lake)) +
  labs(title = "Ventricle Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Ventricle Mass (g)",
       color = "Lake") +
  theme_minimal()



#### gonad mass linear model ####
# Fit the linear model
model_gndms_mass_lk <- lm(gonad_mass ~ 
                            mass_post_resp_g + lake, data = OD_sheet)

# Print the summary of the model
summary(model_gndms_mass_lk)
# Call:
#   lm(formula = gonad_mass ~ mass_post_resp_g + lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -178.319  -27.161   -5.733   27.994  176.545 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -52.92094   37.60748  -1.407    0.164    
# mass_post_resp_g   0.10068    0.01788   5.630 4.84e-07 ***
#   lakeLOTR          16.35031   30.84286   0.530    0.598    
# lakeOpeongo      -17.68601   22.15093  -0.798    0.428    
# lakeShirley       16.27116   32.34600   0.503    0.617    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 63.25 on 61 degrees of freedom
# Multiple R-squared:  0.5438,	Adjusted R-squared:  0.5138 
# F-statistic: 18.17 on 4 and 61 DF,  p-value: 7.091e-10

check_model(model_gndms_mass_lk) # good
forest_model(model_gndms_mass_lk)

# Tidy the model output
tidy_model_gndms_mass_lk <- tidy(model_gndms_mass_lk)

# Create a table
tidy_model_gndms_mass_lk %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	-52.9209443	37.6074798	-1.4071920	0.1644455
# mass_post_resp_g	0.1006824	0.0178834	5.6299453	0.0000005
# lakeLOTR	16.3503074	30.8428558	0.5301165	0.5979556
# lakeOpeongo	-17.6860150	22.1509313	-0.7984321	0.4277174
# lakeShirley	16.2711633	32.3459954	0.5030349	0.6167513       
# Fit linear models for each lake and extract formulas
formulas12 <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(gonad_mass ~ mass_post_resp_g, data = .)) %>%
  mutate(formula = paste0("mass_post_resp_g = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * mass_post_resp_g"))

# Extract the lake and formula into a separate data frame
formulas_12 <- formulas12 %>%
  select(lake, formula)

# Create the base plot
lm12 <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Gonad Mass vs Mass",
       x = "Mass (g)",
       y = "Gonad Mass(g)",
       color = "Lake")

# Add formulas to each facet
lm12 <- lm12 + geom_text(data = formulas_12, aes(x = Inf, y = Inf, label = formula), 
                         hjust = 1.1, vjust = 2, size = 3, color = "black", inherit.aes = FALSE)

# Print the plot
print(lm12)


# Generate predictions
predictionsgndmass <- OD_sheet %>%
  mutate(predicted_gndmass = predict(model_gndms_mass_lk, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionsgndmass, aes(x = mass_post_resp_g, y = predicted_gndmass, color = lake)) +
  labs(title = "Gonad Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Gonad Mass (g)",
       color = "Lake") +
  theme_minimal()


# Load necessary library
library(ggplot2)

# Create the base plot with actual data points
plot_lvrms_mass_lk <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point(size = 3) +  # Dot plot with data points
  
  # Add fitted lines from the linear model
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +  
  
  # Customize the plot
  theme_classic() +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass Post Resp (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme(legend.position = "bottom")

# Print the plot
print(plot_lvrms_mass_lk)


# Generate predictions
predictionslvrmass <- OD_sheet %>%
  mutate(predicted_lvrmass = predict(model_lvrms_mass_lk, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = liver_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionslvrmass, aes(x = mass_post_resp_g, y = predicted_lvrmass, color = lake)) +
  labs(title = "Liver Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme_minimal()

#### model gonad vs mass w/ interaction #### 
model_gndms_mass_lki <- lm(gonad_mass ~ 
                            mass_post_resp_g * lake, data = OD_sheet)

# Print the summary of the model
summary(model_gndms_mass_lki)
# Call:
#   lm(formula = gonad_mass ~ mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -208.377  -25.266   -2.985   22.909  160.235 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  -102.62870   47.33337  -2.168   0.0343 *  
#   mass_post_resp_g                0.12657    0.02336   5.418 1.21e-06 ***
#   lakeLOTR                       -2.22686   86.78190  -0.026   0.9796    
# lakeOpeongo                   105.81463   72.58995   1.458   0.1503    
# lakeShirley                   106.93615   71.38882   1.498   0.1396    
# mass_post_resp_g:lakeLOTR       0.07248    0.10512   0.689   0.4933    
# mass_post_resp_g:lakeOpeongo   -0.06707    0.03759  -1.784   0.0796 .  
# mass_post_resp_g:lakeShirley   -0.09464    0.08892  -1.064   0.2916    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 62.37 on 58 degrees of freedom
# Multiple R-squared:  0.5783,	Adjusted R-squared:  0.5274 
# F-statistic: 11.36 on 7 and 58 DF,  p-value: 5.699e-09


check_model(model_gndms_mass_lki) # good
forest_model(model_gndms_mass_lki)

# Tidy the model output
tidy_model_gndms_mass_lki <- tidy(model_gndms_mass_lki)

# Create a table
tidy_model_gndms_mass_lki %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

formulas12i <- OD_sheet %>%
  group_by(lake) %>%
  do(model = lm(gonad_mass ~ mass_post_resp_g, data = .)) %>%
  mutate(formula = paste0("mass_post_resp_g = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * mass_post_resp_g"))

# Extract the lake and formula into a separate data frame
formulas_12i <- formulas12i %>%
  select(lake, formula)

# Create the base plot
lm12i <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "Gonad Mass vs Mass",
       x = "Mass (g)",
       y = "Gonad Mass(g)",
       color = "Lake")

# Add formulas to each facet
lm12 <- lm12 + geom_text(data = formulas_12i, aes(x = Inf, y = Inf, label = formula), 
                         hjust = 1.1, vjust = 2, size = 3, color = "black", inherit.aes = FALSE)

# Print the plot
print(lm12i)


# Generate predictions
predictionsgndmassi <- OD_sheet %>%
  mutate(predicted_gndmassi = predict(model_gndms_mass_lki, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionsgndmassi, aes(x = mass_post_resp_g, y = predicted_gndmassi, color = lake)) +
  labs(title = "Gonad Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Gonad Mass (g)",
       color = "Lake") +
  theme_minimal()


# Create the base plot with actual data points
plot_lvrms_mass_lki <- ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) +
  geom_point(size = 3) +  # Dot plot with data points
  
  # Add fitted lines from the linear model
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +  
  
  # Customize the plot
  theme_classic() +
  labs(title = "Gonad Mass vs Mass by Lake",
       x = "Mass Post Resp (g)",
       y = "Liver Mass (g)",
       color = "Lake") +
  theme(legend.position = "bottom")

# Print the plot
print(plot_gndms_mass_lki)


# Generate predictions
predictionsgndmassi <- OD_sheet %>%
  mutate(predicted_gndmassi = predict(model_gndms_mass_lk, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = mass_post_resp_g, y = gonad_mass, color = lake)) +
  geom_point() +
  geom_line(data = predictionsgndmass, aes(x = mass_post_resp_g, y = predicted_gndmass, color = lake)) +
  labs(title = "Gonad Mass vs Mass by Lake",
       x = "Mass (g)",
       y = "Gonad Mass (g)",
       color = "Lake") +
  theme_minimal()

                           
#### Gonad Mass ####
                            y.lab <- expression("Gonad Mass")
                            x.lab <- expression("Lakes")
                            ggplot(data = SBB_OD1_data, mapping = aes(x = Llake, y = Gonad_mass)) + 
                              geom_boxplot(aes(fill=Llake), show.legend = FALSE) +
                              theme_classic()+
                              labs(x = x.lab, y = y.lab)
                            
                            #### Dot Plot with Regression Line ####
                            # OD
                            ggplot(data = SBB_OD1_data, mapping = aes(x= Fork_length, 
                                                                      y=OD, 
                                                                      color=Llake))+
                              geom_point() +
                              geom_smooth(method="lm")
                            
                            # Subset by lakes
                            
                            OPEONGO <- SBB_OD1_data [c(1:12,17:20),]
                            ggplot(data = OPEONGO, mapping = aes(x= Fork_length, y=OD)) +
                              geom_point() +
                              geom_smooth(method="lm") +
                              labs(title= "Opeongo Lake",
                                   y="Optical Density SBB", x = "Fork Length (mm)")
                            
                            #iris %>% group_by(Species) %>% summarise(.)
                            
                            HOGAN <- SBB_OD1_data [c(50:66),]
                            ggplot(data = HOGAN, mapping = aes(x= Fork_length, y=OD)) +
                              geom_point() +
                              geom_smooth(method="lm") +
                              labs(title= "Hogan Lake",
                                   y="Optical Density SBB", x = "Fork Length (mm)")
                            
                            #### regression LOTR ####
                            # graph
                            LOTR <- SBB_OD1_data [c(13:16,21:26,43:49),]
                            ggplot(data = LOTR, mapping = aes(x= Fork_length, y=OD)) +
                              geom_point() +
                              geom_smooth(method="lm") +
                              labs(title= "Lake of Two Rivers Lake",
                                   y="Optical Density SBB", x = "Fork Length (mm)")
                            cor(LOTR$Fork_length, LOTR$OD, method = c("pearson"))
                            #-0.3980014
                            cor.test(LOTR$Fork_length, LOTR$OD, method = c("pearson"))
                            #-0.3980014
                            # Pearson's product-moment correlation
                            # 
                            # data:  LOTR$Fork_length and LOTR$OD
                            # t = -1.6803, df = 15, p-value = 0.1136
                            # alternative hypothesis: true correlation is not equal to 0
                            # 95 percent confidence interval:
                            #  -0.7375541  0.1021926
                            # sample estimates:
                            #        cor 
                            # -0.3980014 
                            
                            #### correlation LORT ####
                            cor.test(LOTR$Fork_length, LOTR$OD, method = c("pearson"))
                            # Pearson's product-moment correlation
                            # 
                            # data:  LOTR$Fork_length and LOTR$OD
                            # t = -1.6803, df = 15, p-value = 0.1136
                            # alternative hypothesis: true correlation is not equal to 0
                            # 95 percent confidence interval:
                            # -0.7375541  0.1021926
                            # sample estimates:
                            # cor 
                            # -0.3980014 
                            # regression LOTR formula
                            summary(lm (OD~ Fork_length, data= LOTR))
                            # Call:
                            #   lm(formula = OD ~ Fork_length, data = LOTR)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.18495 -0.05419 -0.01492  0.04494  0.22555 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  1.5050172  0.2982004   5.047 0.000145 ***
                            #   Fork_length -0.0012396  0.0007377  -1.680 0.113605    
                            # ---
                            #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
                            # 
                            # Residual standard error: 0.1042 on 15 degrees of freedom
                            # Multiple R-squared:  0.1584,	Adjusted R-squared:  0.1023 
                            # F-statistic: 2.823 on 1 and 15 DF,  p-value: 0.1136
                            
                            #### regression Shirley ####
                            SHIRLEY <- SBB_OD1_data [c(27:42),]
                            ggplot(data = SHIRLEY, mapping = aes(x= Fork_length, y=OD)) +
                              geom_point() +
                              geom_smooth(method="lm") +
                              labs(title= "Shirley Lake",
                                   y="Optical Density SBB", x = "Fork Length (mm)")
                            
                            #### correlation Shirley ####
                            cor(SHIRLEY$Fork_length, SHIRLEY$OD, method = c("pearson"))
                            #[1] 0.3372664
                            cor.test(SHIRLEY$Fork_length, SHIRLEY$OD, method = c("pearson"))
                            # Pearson's product-moment correlation
                            # 
                            # data:  SHIRLEY$Fork_length and SHIRLEY$OD
                            # t = 1.3405, df = 14, p-value = 0.2014
                            # alternative hypothesis: true correlation is not equal to 0
                            # 95 percent confidence interval:
                            # -0.1902450  0.7136589
                            # sample estimates:
                            # cor 
                            # 0.3372664 
                            
                            lake_data <- split(SBB_OD1, SBB_OD1$Lake)
                            correlations <- lapply(lake_data, function(df) cor(df$OD, df$age, use="complete.obs"))
                            correlations
                            # $hogan
                            # [1] 0.04653125
                            # 
                            # $lotr
                            # [1] -0.4360697
                            # 
                            # $opeongo
                            # [1] -0.129442
                            # 
                            # $shirley
                            # [1] -0.01374059
                            
                            lake_data <- split(data, data$lake)
                            cor_tests <- lapply(lake_data, function(df) cor.test(df$OD, df$age, use="complete.obs"))
                            results <- lapply(cor_tests, function(test) c(correlation = test$estimate, p.value = test$p.value))
                            results
                            # $hogan
                            # correlation.cor         p.value 
                            # 0.04653125      0.85924539 
                            # 
                            # $lotr
                            # correlation.cor         p.value 
                            # -0.4360697       0.1190598 
                            # 
                            # $opeongo
                            # correlation.cor         p.value 
                            # -0.1294420       0.6591872 
                            # 
                            # $shirley
                            # correlation.cor         p.value 
                            # -0.01374059      0.95971932 
                            
                            #############################
                            ###### Sex Differences ######
                            
                            #### boxplot by lake ####
                            
                            #### OD ####
                            y.labs <- expression("Optical Density SBB")
                            x.labs <- expression("Sex")
                            ggplot(data = SBB_OD1_data, mapping = aes(x = Sex, y = OD)) + 
                              geom_boxplot(aes(fill=Sex)) +
                              theme_classic() +
                              labs(x = x.labs, y = y.labs)
                            
                            #### OD by Lake ####
                            y.labs <- expression("Optical Density SBB")
                            x.labs <- expression("Sex")
                            ggplot(data = SBB_OD1_data, mapping = aes(x = Sex, y = OD)) + 
                              geom_boxplot(aes(fill=Llake)) +
                              theme_classic() +
                              labs(x = x.labs, y = y.labs)
                            
                            # regression Shirley formula
                            summary(lm (OD~ Sex, data= SHIRLEY))
                            # Call:
                            #   lm(formula = OD ~ Sex, data = SHIRLEY)
                            # 
                            # Residuals:
                            #   Min        1Q    Median        3Q       Max 
                            # -0.216062 -0.069234 -0.003188  0.034375  0.254937 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  0.98062    0.06137  15.979  2.2e-10 ***
                            #   Sexmale      0.01944    0.07086   0.274    0.788    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1227 on 14 degrees of freedom
                            # Multiple R-squared:  0.005346,	Adjusted R-squared:  -0.0657 
                            # F-statistic: 0.07524 on 1 and 14 DF,  p-value: 0.7879
                            
                            #### FSA PCA on morfolocical data ####
                            # Install and load the necessary packages
                            install.packages("FSA")
                            library(FSA)
                            
                            # Load your data
                            fish_data <- read.csv("your_data_file.csv")
                            
                            # Remove rows with NAs 
                            fish_data_clean <- na.omit(OD_sheet)
                            
                            # Perform PCA
                            pca_results <- prcomp(fish_data_clean[, c("standard_mass", "fork_length", "Age")], scale. = TRUE)
                            
                            # Plot the PCA results
                            plot(pca_results$x[, 1:2], main = "PCA of Lake Trout Morphology Data", xlab = "PC1", ylab = "PC2")
                            
                            
                            # Install and load the necessary packages
                            install.packages("devtools")
                            devtools::install_github("vqv/ggbiplot")
                            library(ggbiplot)
                            
                            # Perform PCA on the cleaned data (assuming `fish_data_clean` from previous steps)
                            pca_results <- prcomp(fish_data_clean[, c("mass", "length", "age")], scale. = TRUE)
                            
                            # Create a biplot with vectors and variable labels
                            ggbiplot(pca_results, 
                                     obs.scale = 1, 
                                     var.scale = 1, 
                                     groups = fish_data_clean$lake, 
                                     ellipse = TRUE, 
                                     circle = TRUE) + 
                              scale_color_discrete(name = "") +
                              theme_minimal() +
                              labs(title = "PCA Biplot of Fish Morphology Data", x = "PC1", y = "PC2")
                            
                            install.packages("ggfortify")
                            library(ggfortify)
                            library(ggplot2)
                            
                            autoplot(pca_results, data = fish_data_clean, colour = "lake", loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
                              theme_minimal() +
                              labs(title = "PCA Biplot of Fish Morphology Data", x = "PC1", y = "PC2")
                            
                            # Perform PCA on your cleaned data (assuming `fish_data_clean` from previous steps)
                            pca_results <- prcomp(fish_data_clean[, c("mass", "length", "age")], scale. = TRUE)
                            
                            # Eigenvalues (variances of the principal components)
                            eigenvalues <- pca_results$sdev^2
                            
                            # Percentage of variance explained by each principal component
                            variance_explained <- eigenvalues / sum(eigenvalues) * 100
                            
                            # Display the results
                            data.frame(
                              Principal_Component = paste0("PC", 1:length(eigenvalues)),
                              Eigenvalue = round(eigenvalues, 3),
                              Variance_Explained = round(variance_explained, 3)
                            )
                            install.packages("ggfortify") 
                            install.packages("ggrepel")
                            library(ggfortify) 
                            library(ggrepel)
                            
                            
                            # Create a PCA plot with vectors and variable labels, using geom_text_repel
                            autoplot(pca_results, data = fish_data_clean, colour = "lake", loadings = TRUE) +
                              geom_text_repel(aes(label = rownames(pca_results$rotation)), 
                                              size = 3, 
                                              box.padding = unit(0.5, "lines"), 
                                              point.padding = unit(0.5, "lines")) +
                              theme_minimal() +
                              labs(title = "PCA Biplot of Fish Morphology Data", x = "PC1", y = "PC2")
                            
                            # Perform PCA on the cleaned data (assuming `fish_data_clean` from previous steps)
                            pca_results <- prcomp(fish_data_clean[, c("mass", "length", "age")], scale. = TRUE)
                            
                            # Extract the scores for the first two principal components
                            pc_scores <- pca_results$x[, 1:2]
                            
                            # Add the quantitative variable (sites) to the scores data frame
                            pc_scores <- as.data.frame(pc_scores)
                            pc_scores$sites <- fish_data_clean$sites
                            
                            # Calculate the correlation between PC1, PC2, and the quantitative variable (sites)
                            correlation_matrix <- cor(pc_scores, use = "complete.obs")
                            
                            # Display the correlation matrix
                            correlation_matrix
                            
                            # Install and load the necessary package
                            install.packages("Hmisc")
                            library(Hmisc)
                            
                            # Calculate the correlation matrix along with p-values
                            correlation_results <- rcorr(as.matrix(pc_scores[, c("PC1", "PC2", "lake")]), type = "pearson")
                            
                            # Display the correlation coefficients
                            correlation_results$r
                            
                            # Display the p-values
                            correlation_results$P
                            
                            
                            #### anova from age versus site ####
                            
                            
                            # Perform ANOVA for each quantitative variable
                            
                            anova_age <- aov(Age ~ lake, data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(anova_age)
                            # Df Sum Sq Mean Sq F value Pr(>F)  
                            # lake         3  288.3   96.09   3.661 0.0175 *
                            #   Residuals   57 1496.1   26.25                 
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 5 observations deleted due to missingness
                            
                            # Perform Tukey's HSD post-hoc test for age
                            tukey_age <- TukeyHSD(anova_age)
                            
                            # Display the post-hoc test results
                            print(tukey_age)
                            # Tukey multiple comparisons of means
                            # 95% family-wise confidence level
                            # 
                            # Fit: aov(formula = Age ~ lake, data = OD_sheet)
                            # 
                            # $lake
                            # diff         lwr        upr     p adj
                            # LOTR-Hogan      -5.9773756 -10.9728231 -0.9819281 0.0128654
                            # Opeongo-Hogan   -1.7568627  -6.5598873  3.0461618 0.7680477
                            # Shirley-Hogan   -1.0735294  -5.7961432  3.6490843 0.9311458
                            # Opeongo-LOTR     4.2205128  -0.9172212  9.3582469 0.1429005
                            # Shirley-LOTR     4.9038462  -0.1587961  9.9664884 0.0609508
                            # Shirley-Opeongo  0.6833333  -4.1895402  5.5562068 0.9823782
                            # Visualize the Tukey's HSD results
                            plot(tukey_age, las = 1)  # `las = 1` ensures labels are horizontal for better readability
                            
                            # Adjust plot margins
                            par(mar = c(5, 8, 4, 2) + 0.1)  # c(bottom, left, top, right)
                            
                            # Visualize the Tukey's HSD results
                            plot(tukey_age, las = 1)  # `las = 1` ensures labels are horizontal for better readability
                            
                            
                            # Perform ANOVA for each quantitative variable
                            
                            anova_age1 <- aov(Age ~ lake + fork_length + sex , data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(anova_age1)
                            # Df Sum Sq Mean Sq F value  Pr(>F)    
                            # lake         3  288.3    96.1   5.053 0.00375 ** 
                            #   fork_length  1  475.8   475.8  25.018 6.6e-06 ***
                            #   sex          3   12.4     4.1   0.218 0.88340    
                            # Residuals   53 1007.9    19.0                    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 5 observations deleted due to missingness
                            
                            
                            # Perform ANOVA for each quantitative variable
                            
                            anova_age2 <- aov(Age ~ lake + fork_length , data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(anova_age2)
                            # Df Sum Sq Mean Sq F value   Pr(>F)    
                            # lake         3  288.3    96.1   5.274  0.00283 ** 
                            #   fork_length  1  475.8   475.8  26.111 4.04e-06 ***
                            #   Residuals   56 1020.3    18.2                     
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 5 observations deleted due to missingness
                            # Perform Tukey's HSD post-hoc test for age
                            tukey_age2 <- TukeyHSD(anova_age2)
                            
                            # Display the post-hoc test results
                            print(tukey_age1)
                            #Visualize the Tukey's HSD results
                            plot(tukey_age1, las = 1)  # `las = 1` ensures labels are horizontal for better readability
                            
                            # linear model
                            lm_age2 <- lm(Age ~ lake + fork_length , data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(lm_age2)
                            # Call:
                            #   lm(formula = Age ~ lake + fork_length, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min      1Q  Median      3Q     Max 
                            # -8.0655 -2.4053 -0.9133  2.4640  9.2410 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept) -17.56392    6.22910  -2.820  0.00664 ** 
                            #   lakeLOTR      2.66416    2.30939   1.154  0.25356    
                            # lakeOpeongo  -0.52997    1.53105  -0.346  0.73053    
                            # lakeShirley   8.73505    2.42799   3.598  0.00068 ***
                            #   fork_length   0.05645    0.01105   5.110 4.04e-06 ***
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 4.269 on 56 degrees of freedom
                            # (5 observations deleted due to missingness)
                            # Multiple R-squared:  0.4282,	Adjusted R-squared:  0.3873 
                            # F-statistic: 10.48 on 4 and 56 DF,  p-value: 2.075e-06
                            forest_model(lm_age2)
                            
#### fork length vs age ####
                            
                            
anova_fl__age_lake <- aov(log(fork_length) ~ Age + lake , data = OD_sheet)
                            
# Summarize the ANOVA results
                            
summary(anova_fl__age_lake)
# Df Sum Sq Mean Sq F value   Pr(>F)    
                            # Age          1 0.4056  0.4056   50.55 2.28e-09 ***
                            #   lake         3 1.4556  0.4852   60.47  < 2e-16 ***
                            #   Residuals   56 0.4493  0.0080                     
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 5 observations deleted due to missingness
tukey_anova_fl__age_lake <- TukeyHSD(anova_fl__age_lake)
                            
# Display the post-hoc test results
print(anova_fl__age_lake)
                            # Call:
                            #   aov(formula = log(fork_length) ~ Age + lake, data = OD_sheet)
                            # 
                            # Terms:
                            #   Age      lake Residuals
                            # Sum of Squares  0.4056104 1.4555864 0.4493157
                            # Deg. of Freedom         1         3        56
                            # 
                            # Residual standard error: 0.08957396
                            # Estimated effects may be unbalanced
                            # 5 observations deleted due to missingness
                            #Visualize the Tukey's HSD results
                            plot(anova_fl__age_lake, las = 1)  # `las = 1` ensures labels are horizontal for better readability
                            
                            # linear model
                            lm_fl_age_lake <- lm(log(fork_length) ~ Age + lake , data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(lm_fl_age_lake)
                            # Call:
                            #   lm(formula = log(fork_length) ~ Age + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min        1Q    Median        3Q       Max 
                            # -0.201822 -0.049141 -0.007417  0.061800  0.208563 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  6.166231   0.038688 159.383  < 2e-16 ***
                            #   Age          0.010666   0.002316   4.606 2.41e-05 ***
                            #   lakeLOTR    -0.254638   0.035788  -7.115 2.23e-09 ***
                            #   lakeOpeongo -0.018944   0.031991  -0.592    0.556    
                            # lakeShirley -0.360442   0.031299 -11.516  < 2e-16 ***
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.08957 on 56 degrees of freedom
                            # (5 observations deleted due to missingness)
                            # Multiple R-squared:  0.8055,	Adjusted R-squared:  0.7916 
                            # F-statistic: 57.99 on 4 and 56 DF,  p-value: < 2.2e-16
                            check_model(lm_fl_age_lake) #
                            # not too bad, no collinearity
                            forest_model(lm_fl_age_lake)
                            
                            # Tidy the model output
                            tidy_model_lm_fl_age_lake <- tidy(lm_fl_age_lake)
                            library(kableExtra)
                            
                            # Create a table
                            tidy_model_lm_fl_age_lake %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # linear model
                            lm_fl_age_lake_sex <- lm(log(fork_length) ~ Age + lake + sex, data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(lm_fl_age_lake_sex)
                            # Call:
                            #   lm(formula = log(fork_length) ~ Age + lake + sex, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min        1Q    Median        3Q       Max 
                            # -0.198742 -0.050462 -0.007095  0.063048  0.192487 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  6.185501   0.038112 162.296  < 2e-16 ***
                            #   Age          0.009705   0.002220   4.372 5.64e-05 ***
                            #   lakeLOTR    -0.289256   0.036284  -7.972 1.10e-10 ***
                            #   lakeOpeongo -0.046456   0.031985  -1.452   0.1522    
                            # lakeShirley -0.395354   0.033120 -11.937  < 2e-16 ***
                            #   sexmale      0.037187   0.024234   1.535   0.1307    
                            # sexunknown  -0.213407   0.088547  -2.410   0.0194 *  
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.08477 on 54 degrees of freedom
                            # (5 observations deleted due to missingness)
                            # Multiple R-squared:  0.8321,	Adjusted R-squared:  0.8134 
                            # F-statistic: 44.59 on 6 and 54 DF,  p-value: < 2.2e-16
                            check_model(lm_fl_age_lake_sex) #
                            # not too bad, no collinearity
                            forest_model(lm_fl_age_lake_sex)
                            
                            # Tidy the model output
                            tidy_model_lm_fl_age_lake_sex <- tidy(lm_fl_age_lake_sex)
                            library(kableExtra)
                            
                            # Create a table
                            tidy_model_lm_fl_age_lake_sex %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # linear model
                            lm_fl_age_lake <- lm(log(fork_length) ~ Age + lake , data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(lm_fl_age_lake)
                            # Call:
                            #   lm(formula = log(fork_length) ~ Age + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min        1Q    Median        3Q       Max 
                            # -0.201822 -0.049141 -0.007417  0.061800  0.208563 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)  6.166231   0.038688 159.383  < 2e-16 ***
                            #   Age          0.010666   0.002316   4.606 2.41e-05 ***
                            #   lakeLOTR    -0.254638   0.035788  -7.115 2.23e-09 ***
                            #   lakeOpeongo -0.018944   0.031991  -0.592    0.556    
                            # lakeShirley -0.360442   0.031299 -11.516  < 2e-16 ***
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.08957 on 56 degrees of freedom
                            # (5 observations deleted due to missingness)
                            # Multiple R-squared:  0.8055,	Adjusted R-squared:  0.7916 
                            # F-statistic: 57.99 on 4 and 56 DF,  p-value: < 2.2e-16
                            check_model(lm_fl_age_lake) #
                            # not too bad, no collinearity
                            forest_model(lm_fl_age_lake)
                            
                            # Tidy the model output
                            tidy_model_lm_fl_age_lake <- tidy(lm_fl_age_lake)
                            library(kableExtra)
                            
                            # Create a table
                            tidy_model_lm_fl_age_lake %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # linear model
                            lm_fl_age_lake_sex <- lm(log(fork_length) ~ Age + lake + sex, data = OD_sheet)
                            
                            # Summarize the ANOVA results
                            
                            summary(lm_fl_age_lake_sex)
                            
                            check_model(lm_fl_age_lake_sex) #
                            # not too bad, no collinearity
                            forest_model(lm_fl_age_lake_sex)
                            
                            # Tidy the model output
                            tidy_model_lm_fl_age_lake_sex <- tidy(lm_fl_age_lake_sex)
                            library(kableExtra)
                            
                            # Create a table
                            tidy_model_lm_fl_age_lake_sex %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
#### organ mass vs mass by lake ####
# liver vs mass
livermassanc <- lm(liver_mass ~ standard_mass + lake, data=OD_sheet)
                            summary(livermassanc)
                            # 
                            # 
                            # Call:
                            #   lm(formula = liver_mass ~ standard_mass + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -10.8976  -2.0400  -0.3717   2.6634  13.7728 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)    3.238563   2.807475   1.154    0.253    
                            # standard_mass  0.012922   0.001573   8.214 1.89e-11 ***
                            #   lakeLOTR      -3.415065   2.281392  -1.497    0.140    
                            # lakeOpeongo   -0.408986   1.644243  -0.249    0.804    
                            # lakeShirley   -3.605389   2.387921  -1.510    0.136    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 4.706 on 61 degrees of freedom
                            # Multiple R-squared:  0.817,	Adjusted R-squared:  0.805 
                            # F-statistic:  68.1 on 4 and 61 DF,  p-value: < 2.2e-16
                            check_model(livermassanc)
                            forest_model(livermassanc)
                            
                            # Tidy the model output
                            tidy_model_lm_livermassanc <- tidy(livermassanc)
                            library(kableExtra)
                            
                            # Create a table
                            tidy_model_lm_livermassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # ventricle vs mass
                            # ventricle
                            ventriclemassanc <- lm(liver_mass ~ standard_mass + lake, 
                                                   data=OD_sheet)
                            summary(ventriclemassanc)
                            # Call:
                            #   lm(formula = liver_mass ~ standard_mass + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -10.8976  -2.0400  -0.3717   2.6634  13.7728 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)    3.238563   2.807475   1.154    0.253    
                            # standard_mass  0.012922   0.001573   8.214 1.89e-11 ***
                            #   lakeLOTR      -3.415065   2.281392  -1.497    0.140    
                            # lakeOpeongo   -0.408986   1.644243  -0.249    0.804    
                            # lakeShirley   -3.605389   2.387921  -1.510    0.136    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 4.706 on 61 degrees of freedom
                            # Multiple R-squared:  0.817,	Adjusted R-squared:  0.805 
                            # F-statistic:  68.1 on 4 and 61 DF,  p-value: < 2.2e-16
                            check_model(ventriclemassanc)
                            forest_model(ventriclemassanc)
                            
                            # Tidy the model output
                            tidy_model_lm_ventriclemassanc <- tidy(ventriclemassanc)
                            
                            # Create a table
                            tidy_model_lm_ventriclemassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # gonad vs mass
                            # gonad
                            gonadmassanc <- lm(gonad_mass ~ standard_mass + lake, 
                                               data=OD_sheet)
                            summary(gonadmassanc)
                            # Call:
                            #   lm(formula = gonad_mass ~ standard_mass + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min      1Q  Median      3Q     Max 
                            # -167.25  -29.58   -6.26   32.35  204.70 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)   -11.95421   41.33224  -0.289 0.773392    
                            # standard_mass   0.09344    0.02316   4.034 0.000155 ***
                            #   lakeLOTR      -11.30838   33.58713  -0.337 0.737508    
                            # lakeOpeongo   -22.83850   24.20689  -0.943 0.349162    
                            # lakeShirley   -13.91939   35.15547  -0.396 0.693532    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 69.28 on 61 degrees of freedom
                            # Multiple R-squared:  0.4527,	Adjusted R-squared:  0.4168 
                            # F-statistic: 12.61 on 4 and 61 DF,  p-value: 1.536e-07
                            check_model(gonadmassanc)
                            forest_model(gonadmassanc)
                            # Tidy the model output
                            tidy_model_lm_gonadmassanc <- tidy(gonadmassanc)
                            
                            # Create a table
                            tidy_model_lm_gonadmassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # log gonad vs log mass
                            lngonadmassanc <- lm(log(gonad_mass) ~ log(standard_mass) + lake, 
                                                 data=OD_sheet)
                            summary(lngonadmassanc)
                            # Call:
                            #   lm(formula = log(gonad_mass) ~ log(standard_mass) + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min      1Q  Median      3Q     Max 
                            # -3.0856 -0.5220  0.1322  0.7601  1.7663 
                            
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)        -16.1489     3.6746  -4.395 4.51e-05 ***
                            #   log(standard_mass)   2.7491     0.4992   5.507 7.72e-07 ***
                            #   lakeLOTR             1.5357     0.6162   2.492  0.01543 *  
                            #   lakeOpeongo          0.3626     0.4014   0.903  0.36983    
                            # lakeShirley          1.8880     0.6864   2.750  0.00782 ** 
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 1.151 on 61 degrees of freedom
                            # Multiple R-squared:  0.451,	Adjusted R-squared:  0.415 
                            # F-statistic: 12.53 on 4 and 61 DF,  p-value: 1.683e-07
                            check_model(lngonadmassanc)
                            forest_model(lngonadmassanc)
                            
                            # Tidy the model output
                            tidy_model_lm_lngonadmassanc <- tidy(lngonadmassanc)
                            
                            # Create a table
                            tidy_model_lm_lngonadmassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # log ventricle
                            lnventriclemassanc <- lm(log(ventricle_mass) ~ standard_mass + lake, 
                                                     data=OD_sheet)
                            summary(lnventriclemassanc)
                            # Call:
                            #   lm(formula = log(ventricle_mass) ~ standard_mass + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.39225 -0.09775 -0.02565  0.10361  0.45248 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)   -9.194e-01  1.058e-01  -8.694 2.84e-12 ***
                            #   standard_mass  7.994e-04  5.926e-05  13.489  < 2e-16 ***
                            #   lakeLOTR      -1.678e-01  8.594e-02  -1.953  0.05545 .  
                            # lakeOpeongo    1.595e-01  6.194e-02   2.575  0.01246 *  
                            #   lakeShirley   -2.526e-01  8.995e-02  -2.809  0.00667 ** 
                            #   ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.1773 on 61 degrees of freedom
                            # Multiple R-squared:  0.9315,	Adjusted R-squared:  0.9271 
                            # F-statistic: 207.5 on 4 and 61 DF,  p-value: < 2.2e-16
                            check_model(lnventriclemassanc)
                            forest_model(lnventriclemassanc)
                            # Tidy the model output
                            tidy_model_lm_lnventriclemassanc <- tidy(lnventriclemassanc)
                            
                            # Create a table
                            tidy_model_lm_lnventriclemassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            # liver
                            lnlivermassanc <- lm(log(liver_mass) ~ log(standard_mass) + lake, 
                                                 data=OD_sheet)
                            summary(lnlivermassanc)
                            # Call:
                            #   lm(formula = log(liver_mass) ~ log(standard_mass) + lake, data = OD_sheet)
                            # 
                            # Residuals:
                            #   Min       1Q   Median       3Q      Max 
                            # -0.56055 -0.22650 -0.00544  0.21123  0.52065 
                            # 
                            # Coefficients:
                            #   Estimate Std. Error t value Pr(>|t|)    
                            # (Intercept)        -4.87926    0.87319  -5.588 5.69e-07 ***
                            #   log(standard_mass)  1.08691    0.11862   9.163 4.53e-13 ***
                            #   lakeLOTR           -0.11615    0.14643  -0.793    0.431    
                            # lakeOpeongo         0.01871    0.09538   0.196    0.845    
                            # lakeShirley        -0.11413    0.16312  -0.700    0.487    
                            # ---
                            #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                            # 
                            # Residual standard error: 0.2735 on 61 degrees of freedom
                            # Multiple R-squared:  0.8724,	Adjusted R-squared:  0.864 
                            # F-statistic: 104.2 on 4 and 61 DF,  p-value: < 2.2e-16
                            check_model(lnlivermassanc)
                            forest_model(lnlivermassanc)
                            # Tidy the model output
                            tidy_model_lm_lnlivermassanc <- tidy(lnlivermassanc)
                            
                            # Create a table
                            tidy_model_lm_lnlivermassanc %>%
                              kable("html", caption = "Linear Regression Results") %>%
                              kable_styling(full_width = FALSE)
                            
                            #### correlation organ vs mass per lake separate ####
                            library(dplyr)
                            library(purrr)
                            library(broom)
                            
                            # liver
                            # Function to perform correlation test
                            cor_test1 <- function(OD_sheet) {
                              cor.test(OD_sheet$liver_mass, OD_sheet$standard_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results1 <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test1))
                            
                            # Extract and print results
                            correlation_results1 %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
                            # A tibble: 4 × 9
                            # Groups:   lake [4]
                            # lake    estimate statistic   p.value parameter conf.low conf.high method alternative
                            # <chr>      <dbl>     <dbl>     <dbl>     <int>    <dbl>     <dbl> <chr>  <chr>      
                            #   1 Opeongo    0.578      2.65 0.0191           14    0.115     0.834 Pears… two.sided  
                            # 2 LOTR       0.809      5.33 0.0000843        15    0.537     0.929 Pears… two.sided  
                            # 3 Shirley    0.814      5.24 0.000125         14    0.533     0.933 Pears… two.sided  
                            # 4 Hogan      0.822      5.60 0.0000511        15    0.565     0.934 Pears… two.sided  
                            
                            
                            # ventricle
                            # Function to perform correlation test
                            cor_test2 <- function(OD_sheet) {
                              cor.test(OD_sheet$ventricle_mass, OD_sheet$standard_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results2 <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test2))
                            
                            # Extract and print results
                            correlation_results2 %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
                            # A tibble: 4 × 9
                            # Groups:   lake [4]
                            # lake    estimate statistic   p.value parameter conf.low conf.high method alternative
                            # <chr>      <dbl>     <dbl>     <dbl>     <int>    <dbl>     <dbl> <chr>  <chr>      
                            #   1 Opeongo    0.954     11.8    1.11e-8        14    0.868     0.984 Pears… two.sided  
                            # 2 LOTR       0.799      5.15   1.18e-4        15    0.518     0.925 Pears… two.sided  
                            # 3 Shirley    0.960     12.9    3.69e-9        14    0.887     0.986 Pears… two.sided  
                            # 4 Hogan      0.955     12.4    2.64e-9        15    0.876     0.984 Pears… two.sided  
                            
                            
                            
                            # gonad
                            # Function to perform correlation test
                            cor_test3 <- function(OD_sheet) {
                              cor.test(OD_sheet$gonad_mass, OD_sheet$standard_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results3 <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test3))
                            
                            # Extract and print results
                            correlation_results3 %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
                            A tibble: 4 × 9
                            # Groups:   lake [4]
                            # lake    estimate statistic p.value parameter conf.low conf.high method   alternative
                            # <chr>      <dbl>     <dbl>   <dbl>     <int>    <dbl>     <dbl> <chr>    <chr>      
                            #   1 Opeongo    0.276      1.07 0.301          14   -0.255     0.679 Pearson… two.sided  
                            # 2 LOTR       0.668      3.48 0.00337        15    0.276     0.870 Pearson… two.sided  
                            # 3 Shirley    0.341      1.36 0.196          14   -0.186     0.716 Pearson… two.sided  
                            # 4 Hogan      0.564      2.64 0.0184         15    0.114     0.822 Pearson… two.sided  
                            
                            
                            #### differently weighted linear model ####
                            # Define weights (for illustration, you might have your own logic for weighting)
                            OD_sheet$weights <- c(1, 1.5, 2, 2.5, 3, 2, 1.5, 1)
                            
                            # Fit weighted linear models for each lake
                            formulas_weighted <- OD_sheet %>%
                              group_by(lake) %>%
                              do(model = lm(fork_length ~ log(Age), data = ., weights = weights)) %>%
                              mutate(formula = paste0("fork_length = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * log(Age)"))
                            
                            # Create the base plot with weights
                            lalm_weighted <- ggplot(OD_sheet, aes(x = log(Age), y = fork_length, color = lake)) + 
                              geom_point(aes(size = weights), alpha = 0.7) +  # Adjust point size based on weights
                              geom_smooth(method = "lm", aes(weight = weights), se = FALSE) +  # Weighted regression line
                              facet_wrap(~ lake) +  # Create separate panels for each lake
                              theme_classic() +
                              labs(title = "Fork Length vs Age (Weighted)",
                                   x = "ln Age (y)",
                                   y = "Fork Length (mm)",
                                   color = "Lake")
                            
                            # Add weighted formulas to each facet
                            lalm_weighted + geom_text(data = formulas_weighted, aes(x = Inf, y = Inf, label = formula), 
                                                      hjust = 1.1, vjust = 2, size = 3, color = "black")
                            
#### liver vs mass differently weighted linear model ####
                            # Define weights (for illustration, you might have your own logic for weighting)
                            OD_sheet$weights <- c(1, 1.5, 2, 2.5, 3, 2, 1.5, 1)
                            
                            # Fit weighted linear models for each lake
                            formulas_weighted <- OD_sheet %>%
                              group_by(lake) %>%
                              do(model = lm(fork_length ~ log(Age), data = ., weights = weights)) %>%
                              mutate(formula = paste0("fork_length = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * log(Age)"))
                            
                            # Create the base plot with weights
                            lalm_weighted <- ggplot(OD_sheet, aes(x = log(Age), y = fork_length, color = lake)) + 
                              geom_point(aes(size = weights), alpha = 0.7) +  # Adjust point size based on weights
                              geom_smooth(method = "lm", aes(weight = weights), se = FALSE) +  # Weighted regression line
                              facet_wrap(~ lake) +  # Create separate panels for each lake
                              theme_classic() +
                              labs(title = "Fork Length vs Age (Weighted)",
                                   x = "ln Age (y)",
                                   y = "Fork Length (mm)",
                                   color = "Lake")
                            
                            # Add weighted formulas to each facet
                            lalm_weighted + geom_text(data = formulas_weighted, aes(x = Inf, y = Inf, label = formula), 
                                                      hjust = 1.1, vjust = 2, size = 3, color = "black")
                            
                            
                            #### writing the report ####
                            #install.packages("quarto")
                            library(quarto)
                            
                            
#### Analysis without the outlier ####
view(OD_sheet)
# Print the row to confirm the outlier
print(OD_sheet[25, ])
                            # lake sample_date     id mass_post_resp_g fork_length total_length
                            # 25 LOTR   9/13/2018 pol_25              588         485          525
                            # standard_mass  sex maturity gonad_mass ventricle_mass liver_mass Age
                            # 25           532 male   mature     16.529          0.475      4.759   6
                            # OD_heart OD_liver
                            # 25   1.0065       NA
                            
                            # Remove the outlier row (e.g., row 10)
OD_sheet_no_outlier <- OD_sheet[-25, ]
                            
#### FL vs mass correlation  ####

# Perform correlation test
cor__noout_fl_stmass <- cor.test(OD_sheet_no_outlier$fork_length, 
                                                          OD_sheet_no_outlier$standard_mass)
                            
# Print the test results
print(cor__noout_fl_stmass)
                            # Pearson's product-moment correlation
                            # 
                            # data:  OD_sheet_no_outlier$fork_length and OD_sheet_no_outlier$standard_mass
                            # t = 36.694, df = 63, p-value < 2.2e-16
                            # alternative hypothesis: true correlation is not equal to 0
                            # 95 percent confidence interval:
                            #  0.9630807 0.9861985
                            # sample estimates:
                            #       cor 
                            # 0.9773951 
                            
                            # Function to perform correlation test
                            cor_test <- function(OD_sheet_no_outlier) {
                              cor.test(OD_sheet_no_outlier$fork_length, 
                                       OD_sheet_no_outlier$standard_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test))
                            
                            # Extract and print results
                            correlation_results %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
                            
                            # A tibble: 4 × 9
                            # Groups:   lake [4]
                            # lake    estimate statistic  p.value parameter conf.low conf.high method       
                            # <chr>      <dbl>     <dbl>    <dbl>     <int>    <dbl>     <dbl> <chr>        
                            #   1 Opeongo    0.969     14.6  7.35e-10        14    0.910     0.989 Pearson's pr…
                            # 2 LOTR       0.647      3.29 4.95e- 3        15    0.242     0.860 Pearson's pr…
                            # 3 Shirley    0.978     17.4  6.88e-11        14    0.935     0.992 Pearson's pr…
                            # 4 Hogan      0.965     14.3  3.93e-10        15    0.904     0.988 Pearson's pr…
                            # # ℹ 1 more variable: alternative <chr>
                            
                            #### linear models per lake ####
                            library(dplyr)
                            library(purrr)
                            library(broom)
                            
                            #### standard mass vs fl per separate lake ####
                            # Fit linear models for each lake
                            models_no <- OD_sheet_no_outlier %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(model = map(data, ~ lm(standard_mass ~ fork_length, data = .x)))
                            
                            # Extract model summaries
                            model_summaries_no <- models_no %>%
                              mutate(summary = map(model, glance)) %>%
                              unnest(summary)
                            
                            # Print model summaries
                            print(model_summaries_no)
                            # A tibble: 4 × 15
                            # Groups:   lake [4]
                            # lake    data     model  r.squared adj.r.squared sigma statistic  p.value    df
                            # <chr>   <list>   <list>     <dbl>         <dbl> <dbl>     <dbl>    <dbl> <dbl>
                            #   1 Opeongo <tibble> <lm>       0.938         0.934 124.       213. 7.35e-10     1
                            # 2 LOTR    <tibble> <lm>       0.912         0.906  33.8      145. 8.79e- 9     1
                            # 3 Shirley <tibble> <lm>       0.956         0.953  35.1      304. 6.88e-11     1
                            # 4 Hogan   <tibble> <lm>       0.931         0.927 149.       204. 3.93e-10     1
                            # ℹ 6 more variables: logLik <dbl>, AIC <dbl>, BIC <dbl>, deviance <dbl>,
                            #   df.residual <int>, nobs <int>
                            
                            # Extract detailed model summaries 
                            model_summaries_no <- models_no %>% mutate(tidy_summary = map(model, tidy), 
                                                                       glance_summary = map(model, glance)) %>% 
                              select(lake, tidy_summary, glance_summary) %>% unnest(c(tidy_summary, glance_summary)) 
                            # Print model summaries print(model_summaries)
                            # did not work
                            # one more time
                            
                            # Fit linear models for each lake
                            models <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(model = map(data, ~ lm(mass_post_resp_g ~ fork_length, data = .x)))
                            
                            # Extract detailed model summaries
                            model_summaries_no <- models_no %>%
                              mutate(tidy_summary = map(model, tidy),
                                     glance_summary = map(model, glance)) %>%
                              select(lake, tidy_summary, glance_summary) %>%
                              unnest(c(tidy_summary, glance_summary), names_sep = "_")
                            
                            # Print model summaries
                            print(model_summaries_no)
                            # A tibble: 8 × 18
                            # Groups:   lake [4]
                            # lake    tidy_summary_term tidy_summary_estimate tidy_summary_std.error
                            # <chr>   <chr>                             <dbl>                  <dbl>
                            #   1 Opeongo (Intercept)                    -3046.                  316.   
                            # 2 Opeongo fork_length                        8.65                  0.593
                            # 3 LOTR    (Intercept)                     -825.                  119.   
                            # 4 LOTR    fork_length                        3.61                  0.299
                            # 5 Shirley (Intercept)                    -1103.                   93.9  
                            # 6 Shirley fork_length                        4.26                  0.245
                            # 7 Hogan   (Intercept)                    -2778.                  311.   
                            # 8 Hogan   fork_length                        7.93                  0.556
                            
                            # ℹ 14 more variables: tidy_summary_statistic <dbl>,
                            #   tidy_summary_p.value <dbl>, glance_summary_r.squared <dbl>,
                            #   glance_summary_adj.r.squared <dbl>, glance_summary_sigma <dbl>,
                            #   glance_summary_statistic <dbl>, glance_summary_p.value <dbl>,
                            #   glance_summary_df <dbl>, glance_summary_logLik <dbl>,
                            #   glance_summary_AIC <dbl>, glance_summary_BIC <dbl>,
                            #   glance_summary_deviance <dbl>, glance_summary_df.residual <int>, …
                            
                            #### comparison to 
                            # Function to perform correlation test
                            cor_test <- function(df) {
                              cor.test(df$fork_length, df$standard_mass)
                            }
                            
                            # Apply correlation test to each lake
                            correlation_results <- OD_sheet %>%
                              group_by(lake) %>%
                              nest() %>%
                              mutate(correlation = map(data, cor_test))
                            
                            # Extract and print results
                            correlation_results %>%
                              mutate(correlation_summary = map(correlation, broom::tidy)) %>%
                              select(lake, correlation_summary) %>%
                              unnest(correlation_summary)
#### as relative mass organ stuff ####
                            
#### Liver relative mass ####
# Create a new column with the ratio of liver_mass to body_mass
OD_sheet$liver_body_ratio <- OD_sheet$liver_mass / OD_sheet$mass_post_resp_g
# Preview the new column
head(OD_sheet$liver_body_ratio)
                            
# View the dataset structure
str(OD_sheet)
# Avoid division by zero or missing values
QC_sheet$liver_body_ratio <- ifelse(QC_sheet$body_mass == 0 | is.na(QC_sheet$body_mass), NA, QC_sheet$liver_mass / QC_sheet$body_mass)

OD_sheet <- read.csv("C:/Users/user/Desktop/CH3/OD_sheet.csv")
# > names(OD_sheet)
# [1] "lake"             "sample_date"      "id"              
# [4] "mass_post_resp_g" "fork_length"      "total_length"    
# [7] "standard_mass"    "sex"              "maturity"        
# [10] "gonad_mass"       "ventricle_mass"   "liver_mass"      
# [13] "Age"              "OD_heart"         "OD_liver"                            
# Define the columns for the numerators
numerator_cols <- c("liver_mass", "gonad_mass", "ventricle_mass")

# Create ratios using apply() and store in new columns
OD_sheet[paste0(numerator_cols, "_mass_post_resp_g")] <- 
  sapply(numerator_cols, function(col) OD_sheet[[col]] / OD_sheet$mass_post_resp_g)

head(OD_sheet)

# names(OD_sheet)
# [1] "lake"                           
# [2] "sample_date"                    
# [3] "id"                              
# [4] "mass_post_resp_g"               
# [5] "fork_length"                    
# [6] "total_length"                   
# [7] "standard_mass"                  
# [8] "sex"                            
# [9] "maturity"                       
# [10] "gonad_mass"                     
# [11] "ventricle_mass"                 
# [12] "liver_mass"                     
# [13] "Age"                            
# [14] "OD_heart"                       
# [15] "OD_liver"                       
# [16] "liver_mass_mass_post_resp_g"    
# [17] "gonad_mass_mass_post_resp_g"    
# [18] "ventricle_mass_mass_post_resp_g"


#### linear model OD heart vs relative gonad mass & lake interaction ####
lm_hrt_gm_lakeii <- lm(OD_heart ~ gonad_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_hrt_gm_lakeii)
# Call:
#   lm(formula = OD_heart ~ gonad_mass_mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.211702 -0.057198 -0.004226  0.050393  0.278638 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)                              1.08267    0.04014  26.971
# gonad_mass_mass_post_resp_g             -0.23096    0.48873  -0.473
# lakeLOTR                                -0.05756    0.05543  -1.038
# lakeOpeongo                             -0.07320    0.06153  -1.190
# lakeShirley                             -0.08078    0.06528  -1.237
# gonad_mass_mass_post_resp_g:lakeLOTR    -0.22361    0.83564  -0.268
# gonad_mass_mass_post_resp_g:lakeOpeongo -0.13728    0.81069  -0.169
# gonad_mass_mass_post_resp_g:lakeShirley  0.06134    1.22988   0.050
# Pr(>|t|)    
# (Intercept)                               <2e-16 ***
#   gonad_mass_mass_post_resp_g                0.638    
# lakeLOTR                                   0.303    
# lakeOpeongo                                0.239    
# lakeShirley                                0.221    
# gonad_mass_mass_post_resp_g:lakeLOTR       0.790    
# gonad_mass_mass_post_resp_g:lakeOpeongo    0.866    
# gonad_mass_mass_post_resp_g:lakeShirley    0.960    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1034 on 58 degrees of freedom
# Multiple R-squared:  0.1121,	Adjusted R-squared:  0.004952 
# F-statistic: 1.046 on 7 and 58 DF,  p-value: 0.4096


forest_model(lm_hrt_gm_lakeii)
# Generate predictions
predictionslm_hrt_gm_lakeii <- OD_sheet %>%
  mutate(predicted_lm_hrt_gm_lakeii = predict(lm_hrt_gm_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = gonad_mass_mass_post_resp_g, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_gm_lakeii, aes(x = gonad_mass_mass_post_resp_g, y = predicted_lm_hrt_gm_lakeii, color = lake)) +
  labs(title = "OD Heart vs Relative Gonad Mass by Lake",
       x = "Relative Gonad Mass (g)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()

#### linear model OD heart vs relative ventricle mass & lake interaction ####
lm_hrt_vn_lakeii <- lm(OD_heart ~ ventricle_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_hrt_vn_lakeii)
# Call:
#   lm(formula = OD_heart ~ ventricle_mass_mass_post_resp_g * lake, 
#      data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.19503 -0.05981 -0.01108  0.04643  0.28824 
# 
# Coefficients:
#   Estimate
# (Intercept)                                  1.055e+00
# ventricle_mass_mass_post_resp_g              1.567e+01
# lakeLOTR                                    -7.773e-02
# lakeOpeongo                                 -9.663e-03
# lakeShirley                                  2.913e-02
# ventricle_mass_mass_post_resp_g:lakeLOTR     1.944e+01
# ventricle_mass_mass_post_resp_g:lakeOpeongo -7.678e+01
# ventricle_mass_mass_post_resp_g:lakeShirley -1.220e+02
# Std. Error t value
# (Intercept)                                  2.130e-01   4.953
# ventricle_mass_mass_post_resp_g              2.580e+02   0.061
# lakeLOTR                                     2.838e-01  -0.274
# lakeOpeongo                                  2.872e-01  -0.034
# lakeShirley                                  2.883e-01   0.101
# ventricle_mass_mass_post_resp_g:lakeLOTR     3.450e+02   0.056
# ventricle_mass_mass_post_resp_g:lakeOpeongo  3.271e+02  -0.235
# ventricle_mass_mass_post_resp_g:lakeShirley  3.457e+02  -0.353
# Pr(>|t|)    
# (Intercept)                                 6.67e-06 ***
#   ventricle_mass_mass_post_resp_g                0.952    
# lakeLOTR                                       0.785    
# lakeOpeongo                                    0.973    
# lakeShirley                                    0.920    
# ventricle_mass_mass_post_resp_g:lakeLOTR       0.955    
# ventricle_mass_mass_post_resp_g:lakeOpeongo    0.815    
# ventricle_mass_mass_post_resp_g:lakeShirley    0.725    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.104 on 58 degrees of freedom
# Multiple R-squared:  0.1017,	Adjusted R-squared:  -0.006763 
# F-statistic: 0.9376 on 7 and 58 DF,  p-value: 0.4849

forest_model(lm_hrt_vn_lakeii)
# Generate predictions
predictionslm_hrt_vn_lakeii <- OD_sheet %>%
  mutate(predicted_lm_hrt_vn_lakeii = predict(lm_hrt_vn_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = ventricle_mass_mass_post_resp_g, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_vn_lakeii, aes(x = ventricle_mass_mass_post_resp_g, y = predicted_lm_hrt_vn_lakeii, color = lake)) +
  labs(title = "OD Heart vs Relative Ventricle Mass by Lake",
       x = "Relative Ventricle Mass (g)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()

#### linear model OD heart vs relative liver mass & lake interaction ####
lm_hrt_lv_lakeii <- lm(OD_heart ~ liver_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_hrt_lv_lakeii)
# Call:
#   lm(formula = OD_heart ~ liver_mass_mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18929 -0.05262 -0.00825  0.04459  0.28801 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                               1.2713     0.1326
# liver_mass_mass_post_resp_g             -16.2295    10.3935
# lakeLOTR                                 -0.2721     0.1666
# lakeOpeongo                              -0.2725     0.1731
# lakeShirley                              -0.1936     0.1703
# liver_mass_mass_post_resp_g:lakeLOTR     16.8429    13.8404
# liver_mass_mass_post_resp_g:lakeOpeongo  15.3482    13.3274
# liver_mass_mass_post_resp_g:lakeShirley   8.4616    14.2601
# t value Pr(>|t|)    
# (Intercept)                               9.586 1.45e-13 ***
#   liver_mass_mass_post_resp_g              -1.562    0.124    
# lakeLOTR                                 -1.634    0.108    
# lakeOpeongo                              -1.575    0.121    
# lakeShirley                              -1.137    0.260    
# liver_mass_mass_post_resp_g:lakeLOTR      1.217    0.229    
# liver_mass_mass_post_resp_g:lakeOpeongo   1.152    0.254    
# liver_mass_mass_post_resp_g:lakeShirley   0.593    0.555    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1017 on 58 degrees of freedom
# Multiple R-squared:  0.1422,	Adjusted R-squared:  0.03862 
# F-statistic: 1.373 on 7 and 58 DF,  p-value: 0.234

forest_model(lm_hrt_lv_lakeii)
# Generate predictions
predictionslm_hrt_lv_lakeii <- OD_sheet %>%
  mutate(predicted_lm_hrt_lv_lakeii = predict(lm_hrt_lv_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = liver_mass_mass_post_resp_g, y = OD_heart, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_hrt_lv_lakeii, aes(x = liver_mass_mass_post_resp_g, y = predicted_lm_hrt_lv_lakeii, color = lake)) +
  labs(title = "OD Heart vs Relative Liver Mass by Lake",
       x = "Relative Liver Mass (g)",
       y = "OD Heart",
       color = "Lake") +
  theme_minimal()

#### linear model OD liver vs relative gonad mass & lake interaction ####

lm_lvr_gm_lakeii <- lm(OD_liver ~ gonad_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_lvr_gm_lakeii)
# Call:
#   lm(formula = OD_liver ~ gonad_mass_mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.49360 -0.06062  0.01695  0.09584  0.26976 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                              0.924241   0.068273
# gonad_mass_mass_post_resp_g              0.112329   0.904060
# lakeLOTR                                 0.048284   0.097311
# lakeOpeongo                             -0.065241   0.114963
# lakeShirley                              0.002976   0.118712
# gonad_mass_mass_post_resp_g:lakeLOTR    -0.852170   1.718732
# gonad_mass_mass_post_resp_g:lakeOpeongo  0.724032   1.758651
# gonad_mass_mass_post_resp_g:lakeShirley -1.033197   2.550261
# t value Pr(>|t|)    
# (Intercept)                              13.537   <2e-16 ***
#   gonad_mass_mass_post_resp_g               0.124    0.902    
# lakeLOTR                                  0.496    0.622    
# lakeOpeongo                              -0.567    0.573    
# lakeShirley                               0.025    0.980    
# gonad_mass_mass_post_resp_g:lakeLOTR     -0.496    0.623    
# gonad_mass_mass_post_resp_g:lakeOpeongo   0.412    0.683    
# gonad_mass_mass_post_resp_g:lakeShirley  -0.405    0.687    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1612 on 41 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.03695,	Adjusted R-squared:  -0.1275 
# F-statistic: 0.2248 on 7 and 41 DF,  p-value: 0.9772

forest_model(lm_lvr_gm_lakeii)
# Generate predictions
predictionslm_lvr_gm_lakeii <- OD_sheet %>%
  mutate(predicted_lm_lvr_gm_lakeii = predict(lm_lvr_gm_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = gonad_mass_mass_post_resp_g, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_lvr_gm_lakeii, aes(x = gonad_mass_mass_post_resp_g, y = predicted_lm_lvr_gm_lakeii, color = lake)) +
  labs(title = "OD Liver vs Relative Gonad Mass by Lake",
       x = "Relative Gonad Mass (g)",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()

#### linear model OD heart vs relative liver mass & lake interaction ####
lm_lv_lv_lakeii <- lm(OD_liver ~ liver_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_lv_lv_lakeii)
# Call:
#   lm(formula = OD_liver ~ liver_mass_mass_post_resp_g * lake, data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.40339 -0.06257  0.01023  0.09323  0.25980 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                               0.92355    0.21330
# liver_mass_mass_post_resp_g               0.59399   17.49009
# lakeLOTR                                 -0.05391    0.28169
# lakeOpeongo                              -0.28301    0.28917
# lakeShirley                               0.27305    0.27422
# liver_mass_mass_post_resp_g:lakeLOTR      6.99125   24.83384
# liver_mass_mass_post_resp_g:lakeOpeongo  20.71442   23.32049
# liver_mass_mass_post_resp_g:lakeShirley -29.83537   23.79696
# t value Pr(>|t|)    
# (Intercept)                               4.330 9.38e-05 ***
#   liver_mass_mass_post_resp_g               0.034    0.973    
# lakeLOTR                                 -0.191    0.849    
# lakeOpeongo                              -0.979    0.333    
# lakeShirley                               0.996    0.325    
# liver_mass_mass_post_resp_g:lakeLOTR      0.282    0.780    
# liver_mass_mass_post_resp_g:lakeOpeongo   0.888    0.380    
# liver_mass_mass_post_resp_g:lakeShirley  -1.254    0.217    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1529 on 41 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.1335,	Adjusted R-squared:  -0.01442 
# F-statistic: 0.9025 on 7 and 41 DF,  p-value: 0.5139

forest_model(lm_lv_lv_lakeii)
# Generate predictions
predictionslm_lv_lv_lakeii <- OD_sheet %>%
  mutate(predicted_lm_lv_lv_lakeii = predict(lm_lv_lv_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = liver_mass_mass_post_resp_g, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_lv_lv_lakeii, aes(x = liver_mass_mass_post_resp_g, y = predicted_lm_lv_lv_lakeii, color = lake)) +
  labs(title = "OD Liver vs Relative Liver Mass by Lake",
       x = "Relative Liver Mass (g)",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()

#### linear model OD liver vs relative ventricle mass & lake interaction ####
lm_lv_vn_lakeii <- lm(OD_liver ~ ventricle_mass_mass_post_resp_g * lake, data = OD_sheet)
summary(lm_lv_vn_lakeii)
# Call:
#   lm(formula = OD_liver ~ ventricle_mass_mass_post_resp_g * lake, 
#      data = OD_sheet)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.44210 -0.05371  0.01519  0.09046  0.24720 
# 
# Coefficients:
#   Estimate
# (Intercept)                                    0.6265
# ventricle_mass_mass_post_resp_g              364.5697
# lakeLOTR                                       0.5368
# lakeOpeongo                                    0.3973
# lakeShirley                                    0.5894
# ventricle_mass_mass_post_resp_g:lakeLOTR    -629.7290
# ventricle_mass_mass_post_resp_g:lakeOpeongo -489.8427
# ventricle_mass_mass_post_resp_g:lakeShirley -745.0604
# Std. Error t value
# (Intercept)                                     0.3887   1.612
# ventricle_mass_mass_post_resp_g               462.9226   0.788
# lakeLOTR                                        0.5034   1.066
# lakeOpeongo                                     0.5409   0.735
# lakeShirley                                     0.5625   1.048
# ventricle_mass_mass_post_resp_g:lakeLOTR      603.4368  -1.044
# ventricle_mass_mass_post_resp_g:lakeOpeongo   602.5112  -0.813
# ventricle_mass_mass_post_resp_g:lakeShirley   665.1737  -1.120
# Pr(>|t|)
# (Intercept)                                    0.115
# ventricle_mass_mass_post_resp_g                0.435
# lakeLOTR                                       0.292
# lakeOpeongo                                    0.467
# lakeShirley                                    0.301
# ventricle_mass_mass_post_resp_g:lakeLOTR       0.303
# ventricle_mass_mass_post_resp_g:lakeOpeongo    0.421
# ventricle_mass_mass_post_resp_g:lakeShirley    0.269
# 
# Residual standard error: 0.1591 on 41 degrees of freedom
# (17 observations deleted due to missingness)
# Multiple R-squared:  0.06172,	Adjusted R-squared:  -0.09847 
# F-statistic: 0.3853 on 7 and 41 DF,  p-value: 0.9056

forest_model(lm_lv_vn_lakeii)
# Generate predictions
predictionslm_lv_vn_lakeii <- OD_sheet %>%
  mutate(predicted_lm_lv_vn_lakeii = predict(lm_lv_vn_lakeii, newdata = .))

# Create the plot
ggplot(OD_sheet, aes(x = ventricle_mass_mass_post_resp_g, y = OD_liver, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_lv_vn_lakeii, aes(x = ventricle_mass_mass_post_resp_g, y = predicted_lm_lv_vn_lakeii, color = lake)) +
  labs(title = "OD Liver vs Relative Ventricle Mass by Lake",
       x = "Relative Ventricle Mass (g)",
       y = "OD Liver",
       color = "Lake") +
  theme_minimal()


#### correlation across tissue OD ####

# Pearson correlation test
result <- cor.test(OD_sheet$OD_heart, OD_sheet$OD_liver, method = "pearson")

# Display the result
print(result)
# Pearson's product-moment correlation
# 
# data:  OD_sheet$OD_heart and OD_sheet$OD_liver
# t = 1.444, df = 47, p-value = 0.1554
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.07970628  0.46061067
# sample estimates:
#       cor 
# 0.2061099 
