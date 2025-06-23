
#### Packages used ####
library(tidyverse)
library(FSA)
library(nlstools)
library(plotrix)

# An R script to analyze a panel of metabolites using linear models
library(tidyverse)
library(sjstats)
library(performance)
library(forestmodel)
library(emmeans)
library(pwr)
library(forestmodel)
library(emmeans)
library(ggpubr)
library(rstatix)
library(FSAdata) # for data
library(FSA)     # for vbFuns(), vbStarts(), confint.bootCase()
library(car)     # for Boot()
library(dplyr)   # for filter(), mutate()
library(ggplot2) 
setwd("E:/PCR")
setwd("G:/PCR")
PCR_lvr_RBC <- read.csv("D:/PCR/PCR_sheet_lvr_RBC.csv", quote="")
PCR_lvr_RBC <- read.csv("G:/PCR/PCR_sheet_lvr_RBC.csv", quote="")
 View(PCR_lvr_RBC)
 names(PCR_lvr_RBC)
 PCR_liver_heart_RBCs <- read.csv("E:/PCR/PCR_liver_heart_RBCs.csv", quote="")
 PCR_liver_heart_RBCs <- read.csv("G:/PCR/PCR_liver_heart_RBCs.csv", quote="")
    View(PCR_liver_heart_RBCs)
    
PCR_liver_heart_RBCs <- read.csv("C:/Users/user/Desktop/CH3/PCR_liver_heart_RBCs.csv", quote="")
View(PCR_liver_heart_RBCs)
names(PCR_liver_heart_RBCs)    
# [1] "lake"           "run_date"       "sample_date"   
# [4] "pol_id"         "mass_pre_resp"  "mass_post_resp"
# [7] "fork_length"    "total_length"   "standard_mass" 
# [10] "sex"            "maturity"       "age"           
# [13] "CT_OX_liver"    "FSH_CT_liver"   "TELB_CT_liver" 
# [16] "T.S_FSH_liver"  "T.S_OX_liver"   "OX_RBCs"       
# [19] "FSH_RBCs"       "TELB_CT_RBCs"   "T.S_FSH_RBCs"  
# [22] "T.S_OX_RBCs"    "OX_Ct_Heart"    "FSH_Ct_Heart"  
# [25] "TELB_Ct_Heart"  "T.S_FSH_Heart"  "T.S_OX_Heart"  
# [28] "gonad_mass"     "ventricle_mass" "liver_mass"     
 names(PCR_liver_heart_RBCs)
 # [1] "lake"           "run_date"       "sample_date"    "pol_id"        
 # [5] "mass_pre_resp"  "mass_post_resp" "fork_length"    "total_length"  
 # [9] "standard_mass"  "sex"            "maturity"       "age"           
 # [13] "CT_OX_liver"    "FSH_CT_liver"   "TELB_CT_liver"  "T.S_FSH_liver" 
 # [17] "T.S_OX_liver"   "OX_RBCs"        "FSH_RBCs"       "TELB_CT_RBCs"  
 # [21] "T.S_FSH_RBCs"   "T.S_OX_RBCs"    "OX_Ct_Heart"    "FSH_Ct_Heart"  
 #[25] "TELB_Ct_Heart"  "T.S_FSH_Heart"  "T.S_OX_Heart"  
 summary(PCR_liver_heart_RBCs)
 # lake             run_date         sample_date       
 # Length:31          Length:31          Length:31         
 # Class :character   Class :character   Class :character  
 # Mode  :character   Mode  :character   Mode  :character  
 # 
 # 
 # 
 # 
 # pol_id          mass_pre_resp    mass_post_resp   fork_length   
 # Length:31          Min.   :0.3800   Min.   :0.240   Min.   :332.0  
 # Class :character   1st Qu.:0.7025   1st Qu.:0.600   1st Qu.:404.2  
 # Mode  :character   Median :1.0850   Median :0.815   Median :486.0  
 # Mean   :1.3017   Mean   :1.201   Mean   :474.2  
 # 3rd Qu.:1.8625   3rd Qu.:1.683   3rd Qu.:532.2  
 # Max.   :2.9900   Max.   :2.850   Max.   :624.0  
 # NA's   :1        NA's   :1       NA's   :1      
 #  total_length   standard_mass        sex              maturity        
 # Min.   :369.0   Min.   : 314.0   Length:31          Length:31         
 # 1st Qu.:445.2   1st Qu.: 579.2   Class :character   Class :character  
 # Median :530.5   Median : 889.5   Mode  :character   Mode  :character  
 # Mean   :520.0   Mean   :1097.3                                        
 # 3rd Qu.:583.5   3rd Qu.:1507.2                                        
 # Max.   :679.0   Max.   :2537.0                                        
 # NA's   :1       NA's   :1                                             
 #      age         CT_OX_liver     FSH_CT_liver   TELB_CT_liver  
 # Min.   : 5.00   Min.   :23.14   Min.   :22.30   Min.   : 9.08  
 # 1st Qu.: 7.00   1st Qu.:24.07   1st Qu.:23.14   1st Qu.: 9.80  
 # Median :10.00   Median :24.57   Median :23.45   Median :10.67  
 # Mean   :11.17   Mean   :24.49   Mean   :23.47   Mean   :10.86  
 # 3rd Qu.:15.00   3rd Qu.:24.86   3rd Qu.:23.66   3rd Qu.:11.48  
 # Max.   :22.00   Max.   :25.83   Max.   :24.88   Max.   :13.38  
 # NA's   :1       NA's   :1       NA's   :1       NA's   :1      
 # T.S_FSH_liver    T.S_OX_liver      OX_RBCs         FSH_RBCs     
 # Min.   : 1319   Min.   : 2885   Min.   :18.83   Min.   : 5.873  
 # 1st Qu.: 3974   1st Qu.: 8237   1st Qu.:23.39   1st Qu.: 6.607  
 # Median : 7650   Median :13163   Median :23.72   Median : 6.962  
 # Mean   : 7507   Mean   :15892   Mean   :23.58   Mean   : 7.239  
 # 3rd Qu.:10255   3rd Qu.:20592   3rd Qu.:24.09   3rd Qu.: 7.680  
 # Max.   :16203   Max.   :57444   Max.   :25.83   Max.   :10.256  
 # NA's   :1       NA's   :1       NA's   :3       NA's   :3       
 #  TELB_CT_RBCs     T.S_FSH_RBCs     T.S_OX_RBCs      OX_Ct_Heart   
 # Min.   : 5.873   Min.   :0.1630   Min.   :  2073   Min.   :23.74  
 # 1st Qu.: 6.586   1st Qu.:0.4928   1st Qu.: 79084   1st Qu.:24.39  
 # Median : 6.962   Median :1.0200   Median :106387   Median :24.84  
 # Mean   : 7.212   Mean   :1.4465   Mean   :115028   Mean   :24.87  
 # 3rd Qu.: 7.642   3rd Qu.:1.6292   3rd Qu.:151418   3rd Qu.:25.41  
 # Max.   :10.256   Max.   :5.3855   Max.   :263595   Max.   :25.97  
 # NA's   :3        NA's   :3        NA's   :3        NA's   :1      
 #  FSH_Ct_Heart   TELB_Ct_Heart   T.S_FSH_Heart     T.S_OX_Heart   
 # Min.   :23.06   Min.   :5.910   Min.   : 35041   Min.   : 59843  
 # 1st Qu.:23.97   1st Qu.:7.246   1st Qu.: 57928   1st Qu.: 82100  
 # Median :24.19   Median :8.121   Median : 84041   Median :110320  
 # Mean   :24.32   Mean   :8.020   Mean   : 90160   Mean   :137117  
 # 3rd Qu.:24.92   3rd Qu.:8.872   3rd Qu.:107511   3rd Qu.:158677  
 # Max.   :25.34   Max.   :9.401   Max.   :266095   Max.   :463262  
 # NA's   :1       NA's   :1       NA's   :1        NA's   :1   
summary(df2)
df2<- PCR_liver_heart_RBCs[1:30,]
# lake             run_date         sample_date       
# Length:30          Length:30          Length:30         
# Class :character   Class :character   Class :character  
# Mode  :character   Mode  :character   Mode  :character  
# 
# 
# 
# 
# pol_id          mass_pre_resp    mass_post_resp   fork_length   
# Length:30          Min.   :0.3800   Min.   :0.240   Min.   :332.0  
# Class :character   1st Qu.:0.7025   1st Qu.:0.600   1st Qu.:404.2  
# Mode  :character   Median :1.0850   Median :0.815   Median :486.0  
# Mean   :1.3017   Mean   :1.201   Mean   :474.2  
# 3rd Qu.:1.8625   3rd Qu.:1.683   3rd Qu.:532.2  
# Max.   :2.9900   Max.   :2.850   Max.   :624.0  
# 
# total_length   standard_mass        sex              maturity        
# Min.   :369.0   Min.   : 314.0   Length:30          Length:30         
# 1st Qu.:445.2   1st Qu.: 579.2   Class :character   Class :character  
# Median :530.5   Median : 889.5   Mode  :character   Mode  :character  
# Mean   :520.0   Mean   :1097.3                                        
# 3rd Qu.:583.5   3rd Qu.:1507.2                                        
# Max.   :679.0   Max.   :2537.0                                        
# 
# age         CT_OX_liver     FSH_CT_liver   TELB_CT_liver  
# Min.   : 5.00   Min.   :23.14   Min.   :22.30   Min.   : 9.08  
# 1st Qu.: 7.00   1st Qu.:24.07   1st Qu.:23.14   1st Qu.: 9.80  
# Median :10.00   Median :24.57   Median :23.45   Median :10.67  
# Mean   :11.17   Mean   :24.49   Mean   :23.47   Mean   :10.86  
# 3rd Qu.:15.00   3rd Qu.:24.86   3rd Qu.:23.66   3rd Qu.:11.48  
# Max.   :22.00   Max.   :25.83   Max.   :24.88   Max.   :13.38  
# 
# T.S_FSH_liver    T.S_OX_liver      OX_RBCs         FSH_RBCs     
# Min.   : 1319   Min.   : 2885   Min.   :18.83   Min.   : 5.873  
# 1st Qu.: 3974   1st Qu.: 8237   1st Qu.:23.39   1st Qu.: 6.607  
# Median : 7650   Median :13163   Median :23.72   Median : 6.962  
# Mean   : 7507   Mean   :15892   Mean   :23.58   Mean   : 7.239  
# 3rd Qu.:10255   3rd Qu.:20592   3rd Qu.:24.09   3rd Qu.: 7.680  
# Max.   :16203   Max.   :57444   Max.   :25.83   Max.   :10.256  
# NA's   :2       NA's   :2       
# TELB_CT_RBCs     T.S_FSH_RBCs     T.S_OX_RBCs      OX_Ct_Heart   
# Min.   : 5.873   Min.   :0.1630   Min.   :  2073   Min.   :23.74  
# 1st Qu.: 6.586   1st Qu.:0.4928   1st Qu.: 79084   1st Qu.:24.39  
# Median : 6.962   Median :1.0200   Median :106387   Median :24.84  
# Mean   : 7.212   Mean   :1.4465   Mean   :115028   Mean   :24.87  
# 3rd Qu.: 7.642   3rd Qu.:1.6292   3rd Qu.:151418   3rd Qu.:25.41  
# Max.   :10.256   Max.   :5.3855   Max.   :263595   Max.   :25.97  
# NA's   :2        NA's   :2        NA's   :2                       
#   FSH_Ct_Heart   TELB_Ct_Heart   T.S_FSH_Heart     T.S_OX_Heart   
#  Min.   :23.06   Min.   :5.910   Min.   : 35041   Min.   : 59843  
#  1st Qu.:23.97   1st Qu.:7.246   1st Qu.: 57928   1st Qu.: 82100  
#  Median :24.19   Median :8.121   Median : 84041   Median :110320  
#  Mean   :24.32   Mean   :8.020   Mean   : 90160   Mean   :137117  
#  3rd Qu.:24.92   3rd Qu.:8.872   3rd Qu.:107511   3rd Qu.:158677  
#  Max.   :25.34   Max.   :9.401   Max.   :266095   Max.   :463262  
names(df2)
# [1] "lake"           "run_date"       "sample_date"    "pol_id"         "mass_pre_resp" 
# [6] "mass_post_resp" "fork_length"    "total_length"   "standard_mass"  "sex"           
# [11] "maturity"       "age"            "CT_OX_liver"    "FSH_CT_liver"   "TELB_CT_liver" 
# [16] "T.S_FSH_liver"  "T.S_OX_liver"   "OX_RBCs"        "FSH_RBCs"       "TELB_CT_RBCs"  
# [21] "T.S_FSH_RBCs"   "T.S_OX_RBCs"    "OX_Ct_Heart"    "FSH_Ct_Heart"   "TELB_Ct_Heart" 
# [26] "T.S_FSH_Heart"  "T.S_OX_Heart"  

# averaging the data for uniting both genes

df2 <- df2 %>%
  mutate(
    CT_liver_ratio = (T.S_FSH_liver + T.S_OX_liver) / 2,
    CT_RBC_ratio = (T.S_FSH_RBCs + T.S_OX_RBCs) / 2,
    CT_heart_ratio = (T.S_FSH_Heart + T.S_OX_Heart) / 2
  )
names(df2)
# [1] "lake"           "run_date"       "sample_date"    "pol_id"        
# [5] "mass_pre_resp"  "mass_post_resp" "fork_length"    "total_length"  
# [9] "standard_mass"  "sex"            "maturity"       "age"           
# [13] "CT_OX_liver"    "FSH_CT_liver"   "TELB_CT_liver"  "T.S_FSH_liver" 
# [17] "T.S_OX_liver"   "OX_RBCs"        "FSH_RBCs"       "TELB_CT_RBCs"  
# [21] "T.S_FSH_RBCs"   "T.S_OX_RBCs"    "OX_Ct_Heart"    "FSH_Ct_Heart"  
# [25] "TELB_Ct_Heart"  "T.S_FSH_Heart"  "T.S_OX_Heart"   "CT_liver_ratio"
# [29] "CT_RBC_ratio"   "CT_heart_ratio"

#### Liver FSH Ratios Graphs ####
# correlation graph between individuals
livercor <- cor.test(df2$T.S_FSH_liver, df2$T.S_OX_liver, method = "pearson")
livercor
# Pearson's product-moment correlation
# 
# data:  df2$T.S_FSH_liver and df2$T.S_OX_liver
# t = 6.3207, df = 28, p-value = 7.755e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5616604 0.8830930
# sample estimates:
#       cor 
# 0.7667738 
#[1] 0.7667738
livercor$p.value
#[1] 7.7551e-07

# without lakes as groups
g1<- ggplot(data = df2, aes(age, T.S_FSH_liver)) + geom_point() + geom_smooth(method="lm")   
g1 

FSH_liver_agevsCT_dotplot <- g1+ ggtitle("Age Vs T/S Ratio FSH/Telb", subtitle="Liver Tissue") + xlab("Age") + ylab("Ct Value Ratio")
FSH_liver_agevsCT_dotplot
png("FSH_liver_agevsCT_dotplot.png")

dev.off()
# adding formula to line in graph
modelg1 <- lm(T.S_FSH_liver ~ age, data = df2)
formula <- paste("y = ", round(coef(modelg1)[1], 2), " + ", round(coef(modelg1)[2], 2), "x", sep = "")


FSH_liver_agevsCT_dotplot_formula <- g1 + 
  ggtitle("Age Vs T/S Ratio FSH/Telb", subtitle="Liver Tissue") + 
  xlab("Age") + 
  ylab("Ct Value Ratio") + 
  annotate("text", x = Inf, y = Inf, label = formula, hjust = 1.1, vjust = 2, size = 5, color = "blue")
FSH_liver_agevsCT_dotplot_formula


# with p value
modelg1 <- lm(T.S_FSH_liver ~ age, data = df2)
summary_modelg1 <- summary(modelg1)
p_valueg1 <- summary_modelg1$coefficients[2, 4]
formulag1 <- paste("y = ", round(coef(modelg1)[1], 2), " + ", round(coef(modelg1)[2], 2), "x", sep = "")
p_value_textg1 <- paste("p-value = ", round(p_valueg1, 4), sep = "")

FSH_liver_agevsCT_dotplotg1 <- g1 + 
  ggtitle("Age Vs T/S Ratio FSH/Telb", subtitle="Liver Tissue") + 
  xlab("Age") + 
  ylab("Ct Value Ratio") + 
  annotate("text", x = Inf, y = Inf, label = formulag1, hjust = 1.1, vjust = 2, size = 5, color = "blue") + 
  annotate("text", x = Inf, y = Inf, label = p_value_textg1, hjust = 1.1, vjust = 3, size = 5, color = "red")
FSH_liver_agevsCT_dotplotg1

png("FSH_liver_agevsCT_dotplotg1.png")

dev.off()

# with lakes as groups
gg1 <- ggplot(df2, aes(x=age, y=T.S_FSH_liver)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio FSH/Telb", subtitle="Liver Tissue", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =15000)
plot(gg1)

#### Liver OX Ratios Graphs ####


# without lakes as groups
g2<- ggplot(data = df2, aes(age, T.S_OX_liver)) + geom_point() + geom_smooth(method="lm")   
g2 
g2+ ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="Liver Tissue") + xlab("Age") + ylab("Ct Value Ratio")
# with lakes as groups
gg2 <- ggplot(df2, aes(x=age, y=T.S_OX_liver)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio OX/Telb", subtitle="Liver Tissue", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =50000)
plot(gg2)

# with p value
modelg2 <- lm(T.S_OX_liver ~ age, data = df2)
summary_modelg2 <- summary(modelg2)
p_valueg2 <- summary_modelg1$coefficients[2, 4]
formulag2 <- paste("y = ", round(coef(modelg2)[1], 2), " + ", round(coef(modelg2)[2], 2), "x", sep = "")
p_value_textg2 <- paste("p-value = ", round(p_valueg2, 4), sep = "")

OX_liver_agevsCT_dotplotg2 <- g2 + 
  ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="Liver Tissue") + 
  xlab("Age") + 
  ylab("Ct Value Ratio") + 
  annotate("text", x = Inf, y = Inf, label = formulag2, hjust = 1.1, vjust = 2, size = 5, color = "blue") + 
  annotate("text", x = Inf, y = Inf, label = p_value_textg2, hjust = 1.1, vjust = 3, size = 5, color = "red")
OX_liver_agevsCT_dotplotg2

png("OX_liver_agevsCT_dotplotg2.png")

dev.off()

#### RBC OX Ratios Graphs ####
# correlation graph between individuals
RBCcor <- cor.test(df2$T.S_FSH_RBCs, df2$T.S_OX_RBCs, method = "pearson")
RBCcor
# Pearson's product-moment correlation
# 
# data:  df2$T.S_FSH_RBCs and df2$T.S_OX_RBCs
# t = 5.5448, df = 26, p-value = 8.028e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5004265 0.8701897
# sample estimates:
#       cor 
# 0.7360792
# without lakes as groups
g3<- ggplot(data = PCR_lvr_RBC, aes(age, T.S_Ratio_OX_RBCs)) + geom_point() + geom_smooth(method="lm")   
g3 
g3+ ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="RBCs") + xlab("Age") + ylab("Ct Value Ratio")
# with lakes as groups
gg3 <- ggplot(PCR_lvr_RBC, aes(x=age, y=T.S_Ratio_OX_RBCs)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio OX/Telb", subtitle="RBCs", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =250000)
plot(gg3)

#### RBC FSH Ratios Graphs ####
# without lakes as groups
g4<- ggplot(data = df2, aes(age, T.S_FSH_RBCs)) + geom_point() + geom_smooth(method="lm")   
g4 + ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="RBCs") + xlab("Age") + ylab("Ct Value Ratio")
# with lakes as groups
gg4 <- ggplot(df2, aes(x=age, y=T.S_FSH_RBCs)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio FSH/Telb", subtitle="RBCs", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =5)
plot(gg4)

#### heart OX Ratios Graphs ####
# correlation graph between individuals
heartcor <- cor.test(df2$T.S_FSH_Heart, df2$T.S_OX_Heart, method = "pearson")
heartcor
# Pearson's product-moment correlation
# 
# data:  df2$T.S_FSH_Heart and df2$T.S_OX_Heart
# t = 4.5436, df = 28, p-value = 9.648e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3804867 0.8194100
# sample estimates:
#       cor 
# 0.6514559 

# without lakes as groups
g5<- ggplot(data = df2, aes(age, T.S_OX_Heart)) + geom_point() + geom_smooth(method="lm")   
g5 
g5+ ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="heart") + xlab("Age") + ylab("Ct Value Ratio")
# with lakes as groups
gg5 <- ggplot(df2, aes(x=age, y=T.S_OX_Heart)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio OX/Telb", subtitle="RBCs", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =250000)
plot(gg5)

#### RBC FSH Ratios Graphs ####
# without lakes as groups
g4<- ggplot(data = df2, aes(age, T.S_FSH_RBCs)) + geom_point() + geom_smooth(method="lm")   
g4 + ggtitle("Age Vs T/S Ratio OX/Telb", subtitle="RBCs") + xlab("Age") + ylab("Ct Value Ratio")
# with lakes as groups
gg4 <- ggplot(df2, aes(x=age, y=T.S_FSH_RBCs)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio FSH/Telb", subtitle="RBCs", y="Ct Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =5)
plot(gg4)




#### Linear Models ####
#### Linear Models average CT ####
# "CT_liver_ratio"
# [29] "CT_RBC_ratio"   "CT_heart_ratio"


#### ct lver vs age ####
lmAgePCRliver1 = lm(CT_liver_ratio~age, data = df2) #Create the linear regression
summary(lmAgePCRliver1)
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -10066  -4033  -1006   2842  22315 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  13790.7     3593.0   3.838 0.000647 ***
#   age           -187.3      297.5  -0.630 0.534075    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7503 on 28 degrees of freedom
# Multiple R-squared:  0.01396,	Adjusted R-squared:  -0.02126 
# F-statistic: 0.3964 on 1 and 28 DF,  p-value: 0.5341


# not significant 

lmAgePCRliver3 = lm(CT_liver_ratio~age+lake+sex, data = df2) #Create the linear regression
summary(lmAgePCRliver3)
# Call:
#   lm(formula = CT_liver_ratio ~ age + lake + sex, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10454.2  -3818.2    402.8   3134.7  17351.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   8780.2     3604.1   2.436   0.0230 *
#   age           -153.2      254.8  -0.601   0.5536  
# lakelotr      6866.1     3648.5   1.882   0.0726 .
# lakeopeongo  -1981.3     3607.7  -0.549   0.5882  
# lakeshirley   5374.8     3820.0   1.407   0.1728  
# sexmale       4189.6     2872.7   1.458   0.1582  
# sexunknown   -4146.1     6582.7  -0.630   0.5350  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6149 on 23 degrees of freedom
# Multiple R-squared:  0.4559,	Adjusted R-squared:  0.3139 
# F-statistic: 3.212 on 6 and 23 DF,  p-value: 0.01934

# with interaction 

lmAgePCRliver3i = lm(CT_liver_ratio~age*lake*sex, data = df2) #Create the linear regression
summary(lmAgePCRliver3i)
# Call:
#   lm(formula = CT_liver_ratio ~ age * lake * sex, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -7434  -2930      0   1605  15695 
# 
# Coefficients: (10 not defined because of singularities)
# Estimate Std. Error t value
# (Intercept)                  8054.14    7109.23   1.133
# age                           -52.24     593.91  -0.088
# lakelotr                     7437.47   12857.27   0.578
# lakeopeongo                -12566.78   12620.88  -0.996
# lakeshirley                -33034.24   24014.71  -1.376
# sexmale                     16384.58   14757.14   1.110
# sexunknown                  -4227.55    6980.48  -0.606
# age:lakelotr                  196.70    1067.15   0.184
# age:lakeopeongo              1522.25    1473.42   1.033
# age:lakeshirley              1658.05    1427.36   1.162
# age:sexmale                 -1406.46    1195.12  -1.177
# age:sexunknown                    NA         NA      NA
# lakelotr:sexmale            -3752.29    8726.97  -0.430
# lakeopeongo:sexmale         -3981.47   10748.36  -0.370
# lakeshirley:sexmale         24365.92   12227.55   1.993
# lakelotr:sexunknown               NA         NA      NA
# lakeopeongo:sexunknown            NA         NA      NA
# lakeshirley:sexunknown            NA         NA      NA
# age:lakelotr:sexmale              NA         NA      NA
# age:lakeopeongo:sexmale           NA         NA      NA
# age:lakeshirley:sexmale           NA         NA      NA
# age:lakelotr:sexunknown           NA         NA      NA
# age:lakeopeongo:sexunknown        NA         NA      NA
# age:lakeshirley:sexunknown        NA         NA      NA
# Pr(>|t|)  
# (Intercept)                  0.2739  
# age                          0.9310  
# lakelotr                     0.5710  
# lakeopeongo                  0.3342  
# lakeshirley                  0.1879  
# sexmale                      0.2833  
# sexunknown                   0.5533  
# age:lakelotr                 0.8561  
# age:lakeopeongo              0.3169  
# age:lakeshirley              0.2624  
# age:sexmale                  0.2565  
# age:sexunknown                   NA  
# lakelotr:sexmale             0.6730  
# lakeopeongo:sexmale          0.7159  
# lakeshirley:sexmale          0.0636 .
# lakelotr:sexunknown              NA  
# lakeopeongo:sexunknown           NA  
# lakeshirley:sexunknown           NA  
# age:lakelotr:sexmale             NA  
# age:lakeopeongo:sexmale          NA  
# age:lakeshirley:sexmale          NA  
# age:lakelotr:sexunknown          NA  
# age:lakeopeongo:sexunknown       NA  
# age:lakeshirley:sexunknown       NA  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6269 on 16 degrees of freedom
# Multiple R-squared:  0.6066,	Adjusted R-squared:  0.2869 
# F-statistic: 1.897 on 13 and 16 DF,  p-value: 0.1124

predictions_liver <- df2 %>%
  mutate(predicted_liver_ratio = predict(lmAgePCRliver3i, newdata = .))

forest_model(lmAgePCRliver3i)


lmAgePCRliver2 = lm(CT_liver_ratio~age+lake, data = df2) #Create the linear regression
summary(lmAgePCRliver2)
# Call:
#   lm(formula = CT_liver_ratio ~ age + lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -12993.1  -3275.7   -139.4   2520.4  17966.1 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   8380.6     3472.2   2.414  0.02345 * 
#   age           -116.0      255.7  -0.454  0.65383   
# lakelotr      9913.3     3163.9   3.133  0.00437 **
#   lakeopeongo    959.8     3153.6   0.304  0.76338   
# lakeshirley   8903.8     3161.8   2.816  0.00935 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6225 on 25 degrees of freedom
# Multiple R-squared:  0.3939,	Adjusted R-squared:  0.2969 
# F-statistic: 4.061 on 4 and 25 DF,  p-value: 0.01134

report(lmAgePCRliver2)
# We fitted a linear model (estimated using OLS) to predict CT_liver_ratio
# with age and lake (formula: CT_liver_ratio ~ age + lake). The model explains
# a statistically significant and substantial proportion of variance (R2 =
#                                                                       We fitted a linear model (estimated using OLS) to predict CT_liver_ratio
#                                                                     with age and lake (formula: CT_liver_ratio ~ age + lake). The model explains
# #                                                                     a statistically significant and substantial proportion of variance (R2 =
# #                                                                                                                                           0.39, F(4, 25) = 4.06, p = 0.011, adj. R2 = 0.30). The model's intercept,
# corresponding to age = 0 and lake = hogan, is at 8380.63 (95% CI [1229.56,
# 15531.71], t(25) = 2.41, p = 0.023). Within this model:
# 
#   - The effect of age is statistically non-significant and negative (beta =
# -116.04, 95% CI [-642.59, 410.51], t(25) = -0.45, p = 0.654; Std. beta =
# -0.07, 95% CI [-0.41, 0.26])
#   - The effect of lake [lotr] is statistically significant and positive (beta
# = 9913.32, 95% CI [3397.19, 16429.46], t(25) = 3.13, p = 0.004; Std. beta =
# 1.34, 95% CI [0.46, 2.21])
#   - The effect of lake [opeongo] is statistically non-significant and positive
# (beta = 959.80, 95% CI [-5535.07, 7454.67], t(25) = 0.30, p = 0.763; Std.
# beta = 0.13, 95% CI [-0.75, 1.00])
#   - The effect of lake [shirley] is statistically significant and positive
# (beta = 8903.83, 95% CI [2391.94, 15415.72], t(25) = 2.82, p = 0.009; Std.
# beta = 1.20, 95% CI [0.32, 2.08])
# 
# Standardized parameters were obtained by fitting the model on a standardized
# version of the dataset. 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation. CI [3397.19, 16429.46], t(25) = 3.13, p = 0.004; Std. beta =
# 1.34, 95% CI [0.46, 2.21])
#   - The effect of lake [opeongo] is statistically non-significant and positive
# (beta = 959.80, 95% CI [-5535.07, 7454.67], t(25) = 0.30, p = 0.763; Std.
# beta = 0.13, 95% CI [-0.75, 1.00])
#   - The effect of lake [shirley] is statistically significant and positive
# (beta = 8903.83, 95% CI [2391.94, 15415.72], t(25) = 2.82, p = 0.009; Std.
# beta = 1.20, 95% CI [0.32, 2.08])
# 
# Standardized parameters were obtained by fitting the model on a standardized
# version of the dataset. 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation.

anova(lmAgePCRliver2)
# Analysis of Variance Table
# 
# Response: CT_liver_ratio
# Df    Sum Sq   Mean Sq F value   Pr(>F)   
# age        1  22310983  22310983  0.5757 0.455091   
# lake       3 607243458 202414486  5.2231 0.006147 **
#   Residuals 25 968836279  38753451                    
# ---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
check_model(lmAgePCRliver2)
forest_model(lmAgePCRliver2)
confint(lmAgePCRliver2)
# 2.5 %     97.5 %
#   (Intercept)  1229.5597 15531.7070
# age          -642.5877   410.5095
# lakelotr     3397.1876 16429.4569
# lakeopeongo -5535.0747  7454.6729
# lakeshirley  2391.9427 15415.7188

levels(df2$lake)
table(df2$lake)
df2$lake <- droplevels(df2$lake)
df2$lake <- as.factor(df2$lake)
# formulaslma <- df2 %>%
#   group_by(lake) %>%
#   do(model = lm(CT_liver_ratio ~ age, data = .)) %>%
#   mutate(formula = paste0("lake = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * age"))
# # Create the base plot
# mlma <- ggplot(df2, aes(x = age, y = CT_liver_ratio, color = lake)) + 
#   geom_point(size = 3) + 
#   geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
#   facet_wrap(~ lake) +  # Create separate panels for each lake
#   theme_classic() +
#   labs(title = "CT Liver Ratio vs Age",
#        x = "Age (y)",
#        y = "Ct Liver Ratio",
#        color = "Lake")
# mlma
# # Add formulas to each facet
# mlma + geom_text(data = formulaslm, aes(x = Inf, y = Inf, label = formulaslm), 
#                 hjust = 1.1, vjust = 2, size = 3, color = "black")

formulaslm <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ age, data = .)) %>%
  mutate(formula = paste0("CT_liver_ratio = ", round(coef(model)[1], 2), 
                          " + ", round(coef(model)[2], 2), " * age"))

# Create the base plot
mlma <- ggplot(df2, aes(x = age, y = CT_liver_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "CT Liver Ratio vs Age",
       x = "Age (y)",
       y = "Ct Liver Ratio",
       color = "Lake")

# Add formulas to each facet
mlma <- mlma + geom_text(data = formulaslm, aes(x = Inf, y = Inf, label = formula), 
                         hjust = 1.1, vjust = 2, size = 3, color = "black", inherit.aes = FALSE)

# Print the plot
print(mlma)

# Tidy the model output


# Tidy the model output
tidy_model_PCRlm2 <- tidy(lmAgePCRliver2)
library(kableExtra)
# Create a table
tidy_model_PCRlm2 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	8380.6334	3472.1728	2.4136568	0.0234479
# age	-116.0391	255.6634	-0.4538746	0.6538348
# lakelotr	9913.3222	3163.8809	3.1332792	0.0043739
# lakeopeongo	959.7991	3153.5578	0.3043544	0.7633770
# lakeshirley	8903.8307	3161.8190	2.8160470	0.0093506

predictionsPCRlm2 <- df2 %>% mutate(predicted_PCRlm2 = predict(lmAgePCRliver2, newdata = .))

ggplot(df2, aes(x = age, y = CT_liver_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionsPCRlm2, aes(x = age, y = predicted_PCRlm2, color = lake)) + labs(title = "CT Liver Ratio vs Age by Lake", x = "Age (y)", y = "CT Liver Tissue", color = "Lake") +
  theme_minimal()

#### CT liver vs age by lake separated ####
# Group by lakes and run linear models
model_liver_age_bylake <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ age, data = .))

# To view the models
model_liver_age_bylake

# View summaries of each model
model_summaries_liver_age_bylake <- model_liver_age_bylake %>%
  rowwise() %>%
  mutate(summary = list(summary(model)))

model_summaries_liver_age_bylake$summary

# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -3811  -2243   1257   1389   4855 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  6576.72    3396.49   1.936    0.094 .
# age            49.63     295.63   0.168    0.871  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3250 on 7 degrees of freedom
# Multiple R-squared:  0.00401,	Adjusted R-squared:  -0.1383 
# F-statistic: 0.02818 on 1 and 7 DF,  p-value: 0.8714
# 
# 
# [[2]]
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 6522 -2802 -5258  4559  8557 -6752 -4827 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  22710.5     6455.5   3.518    0.017 *
#   age           -591.7      635.2  -0.931    0.394  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6942 on 5 degrees of freedom
# Multiple R-squared:  0.1479,	Adjusted R-squared:  -0.02256 
# F-statistic: 0.8676 on 1 and 5 DF,  p-value: 0.3944
# 
# 
# [[3]]
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 1688 -4417  1961  3856 -4047 -1351  2310 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   5332.6     3569.7   1.494    0.195
# age            214.0      271.9   0.787    0.467
# 
# Residual standard error: 3595 on 5 degrees of freedom
# Multiple R-squared:  0.1103,	Adjusted R-squared:  -0.06766 
# F-statistic: 0.6198 on 1 and 5 DF,  p-value: 0.4668
# 
# 
# [[4]]
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1      2      3      4      5      6      7 
# -12455  -5214  -4052  -1295   1302  17916   3798 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  18748.7    10430.4   1.798    0.132
# age           -233.8      777.7  -0.301    0.776
# 
# Residual standard error: 10370 on 5 degrees of freedom
# Multiple R-squared:  0.01776,	Adjusted R-squared:  -0.1787 
# F-statistic: 0.09041 on 1 and 5 DF,  p-value: 0.7758
# View summaries of each model and include lake name
model_summaries_liver_age_bylake <- model_liver_age_bylake %>%
  rowwise() %>%
  mutate(lake = lake, summary = list(summary(model)))

# Display the results
model_summaries_liver_age_bylake %>% select(lake, summary)
# Print the summaries with lake names
for (i in 1:nrow(model_summaries_liver_age_bylake)) {
  cat("Lake:", model_summaries_liver_age_bylake$lake[i], "\n")
  print(model_summaries_liver_age_bylake$summary[[i]])
  cat("\n")
}

# Lake: 1 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -3811  -2243   1257   1389   4855 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  6576.72    3396.49   1.936    0.094 .
# age            49.63     295.63   0.168    0.871  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3250 on 7 degrees of freedom
# Multiple R-squared:  0.00401,	Adjusted R-squared:  -0.1383 
# F-statistic: 0.02818 on 1 and 7 DF,  p-value: 0.8714
# 
# 
# Lake: 2 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 6522 -2802 -5258  4559  8557 -6752 -4827 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  22710.5     6455.5   3.518    0.017 *
#   age           -591.7      635.2  -0.931    0.394  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6942 on 5 degrees of freedom
# Multiple R-squared:  0.1479,	Adjusted R-squared:  -0.02256 
# F-statistic: 0.8676 on 1 and 5 DF,  p-value: 0.3944
# 
# 
# Lake: 3 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 1688 -4417  1961  3856 -4047 -1351  2310 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   5332.6     3569.7   1.494    0.195
# age            214.0      271.9   0.787    0.467
# 
# Residual standard error: 3595 on 5 degrees of freedom
# Multiple R-squared:  0.1103,	Adjusted R-squared:  -0.06766 
# F-statistic: 0.6198 on 1 and 5 DF,  p-value: 0.4668
# 
# 
# Lake: 4 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1      2      3      4      5      6      7 
# -12455  -5214  -4052  -1295   1302  17916   3798 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  18748.7    10430.4   1.798    0.132
# age           -233.8      777.7  -0.301    0.776
# 
# Residual standard error: 10370 on 5 degrees of freedom
# Multiple R-squared:  0.01776,	Adjusted R-squared:  -0.1787 
# F-statistic: 0.09041 on 1 and 5 DF,  p-value: 0.7758

# Sort data frame by lake
df2 <- df2 %>% arrange(lake)

# Group by lakes and run linear models
model_liver_age_bylake <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ age, data = .))

# View summaries of each model and include lake name
model_summaries_liver_age_bylake <- model_liver_age_bylake %>%
  rowwise() %>%
  mutate(lake = lake, summary = list(summary(model)))

# Print the summaries with lake names
for (i in 1:nrow(model_summaries_liver_age_bylake)) {
  cat("Lake:", model_summaries_liver_age_bylake$lake[i], "\n")
  print(model_summaries_liver_age_bylake$summary[[i]])
  cat("\n")
}
# Lake: 1 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -3811  -2243   1257   1389   4855 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  6576.72    3396.49   1.936    0.094 .
# age            49.63     295.63   0.168    0.871  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3250 on 7 degrees of freedom
# Multiple R-squared:  0.00401,	Adjusted R-squared:  -0.1383 
# F-statistic: 0.02818 on 1 and 7 DF,  p-value: 0.8714
# 
# 
# Lake: 2 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 6522 -2802 -5258  4559  8557 -6752 -4827 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  22710.5     6455.5   3.518    0.017 *
#   age           -591.7      635.2  -0.931    0.394  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6942 on 5 degrees of freedom
# Multiple R-squared:  0.1479,	Adjusted R-squared:  -0.02256 
# F-statistic: 0.8676 on 1 and 5 DF,  p-value: 0.3944
# 
# 
# Lake: 3 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 1688 -4417  1961  3856 -4047 -1351  2310 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   5332.6     3569.7   1.494    0.195
# age            214.0      271.9   0.787    0.467
# 
# Residual standard error: 3595 on 5 degrees of freedom
# Multiple R-squared:  0.1103,	Adjusted R-squared:  -0.06766 
# F-statistic: 0.6198 on 1 and 5 DF,  p-value: 0.4668
# 
# 
# Lake: 4 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ age, data = .)
# 
# Residuals:
#   1      2      3      4      5      6      7 
# -12455  -5214  -4052  -1295   1302  17916   3798 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  18748.7    10430.4   1.798    0.132
# age           -233.8      777.7  -0.301    0.776
# 
# Residual standard error: 10370 on 5 degrees of freedom
# Multiple R-squared:  0.01776,	Adjusted R-squared:  -0.1787 
# F-statistic: 0.09041 on 1 and 5 DF,  p-value: 0.7758


#### CT liver vs age and lake with interaction ####
# not significant for any of the variables
# including the intercept
lmAgePCRliver2i = lm(CT_liver_ratio~age*lake, data = df2) #Create the linear regression
summary(lmAgePCRliver2i)
# Call:
#   lm(formula = CT_liver_ratio ~ age * lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -12454.5  -3988.5    -18.9   2780.6  17915.6 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)      6576.72    6746.27   0.975   0.3402  
# age                49.63     587.19   0.085   0.9334  
# lakelotr        16133.82    9031.02   1.786   0.0878 .
# lakeopeongo     -1244.14    9306.69  -0.134   0.8949  
# lakeshirley     12171.94    9364.61   1.300   0.2071  
# age:lakelotr     -641.30     832.94  -0.770   0.4495  
# age:lakeopeongo   164.39     763.65   0.215   0.8315  
# age:lakeshirley  -283.47     761.14  -0.372   0.7131  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6456 on 22 degrees of freedom
# Multiple R-squared:  0.4263,	Adjusted R-squared:  0.2438 
# F-statistic: 2.335 on 7 and 22 DF,  p-value: 0.06068
report(lmAgePCRliver2i)
# We fitted a linear model (estimated using OLS) to
# predict CT_liver_ratio with age and lake (formula:
#                                             CT_liver_ratio ~ age * lake). The model explains a
# statistically not significant and substantial
# proportion of variance (R2 = 0.43, F(7, 22) = 2.34, p =
#                           0.061, adj. R2 = 0.24). The model's intercept,
# corresponding to age = 0 and lake = hogan, is at
# 6576.72 (95% CI [-7414.19, 20567.62], t(22) = 0.97, p =
# 0.340). Within this model:
# 
#   - The effect of age is statistically non-significant
# and positive (beta = 49.63, 95% CI [-1168.13, 1267.38],
# t(22) = 0.08, p = 0.933; Std. beta = 0.03, 95% CI
# [-0.74, 0.80])
#   - The effect of lake [lotr] is statistically
# non-significant and positive (beta = 16133.82, 95% CI
# [-2595.38, 34863.01], t(22) = 1.79, p = 0.088; Std.
# beta = 1.21, 95% CI [0.25, 2.17])
#   - The effect of lake [opeongo] is statistically
# non-significant and negative (beta = -1244.14, 95% CI
# [-20545.03, 18056.74], t(22) = -0.13, p = 0.895; Std.
# beta = 0.08, 95% CI [-0.84, 1.00])
#   - The effect of lake [shirley] is statistically
# non-significant and positive (beta = 12171.94, 95% CI
# [-7249.07, 31592.95], t(22) = 1.30, p = 0.207; Std.
# beta = 1.21, 95% CI [0.29, 2.14])
#   - The effect of age × lake [lotr] is statistically
# non-significant and negative (beta = -641.30, 95% CI
# [-2368.72, 1086.12], t(22) = -0.77, p = 0.450; Std.
# beta = -0.40, 95% CI [-1.49, 0.69])
#   - The effect of age × lake [opeongo] is statistically
# non-significant and positive (beta = 164.39, 95% CI
# [-1419.32, 1748.11], t(22) = 0.22, p = 0.832; Std. beta
# = 0.10, 95% CI [-0.90, 1.10])
#   - The effect of age × lake [shirley] is statistically
# non-significant and negative (beta = -283.47, 95% CI
# [-1861.98, 1295.03], t(22) = -0.37, p = 0.713; Std.
# beta = -0.18, 95% CI [-1.17, 0.82])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.


# Generate predictions
predictionslmAgePCRliver2i <- df2 %>%
  mutate(predicted_lmAgePCRliver2i = predict(lmAgePCRliver2i, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmAgePCRliver2i, aes(x = age, y = predicted_lmAgePCRliver2i, color = lake)) +
  labs(title = "Telomere Ratio Liver vs Age by Lake",
       x = "Age (y)",
       y = "Telomere Liver Ratio",
       color = "Lake") +
  theme_minimal()

#### CT liver vs age and lake and sex with interaction ####
# not significant for any of the variables
# including the intercept
lmAgePCRliver6i = lm(CT_liver_ratio~age*lake*sex, data = df2) #Create the linear regression
summary(lmAgePCRliver6i)
# Call:
#   lm(formula = CT_liver_ratio ~ age * lake * sex, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -7434  -2930      0   1605  15695 
# 
# Coefficients: (10 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)                  8054.14    7109.23   1.133   0.2739
# age                           -52.24     593.91  -0.088   0.9310
# lakelotr                     7437.47   12857.27   0.578   0.5710
# lakeopeongo                -12566.78   12620.88  -0.996   0.3342
# lakeshirley                -33034.24   24014.71  -1.376   0.1879
# sexmale                     16384.58   14757.14   1.110   0.2833
# sexunknown                  -4227.55    6980.48  -0.606   0.5533
# age:lakelotr                  196.70    1067.15   0.184   0.8561
# age:lakeopeongo              1522.25    1473.42   1.033   0.3169
# age:lakeshirley              1658.05    1427.36   1.162   0.2624
# age:sexmale                 -1406.46    1195.12  -1.177   0.2565
# age:sexunknown                    NA         NA      NA       NA
# lakelotr:sexmale            -3752.29    8726.97  -0.430   0.6730
# lakeopeongo:sexmale         -3981.47   10748.36  -0.370   0.7159
# lakeshirley:sexmale         24365.92   12227.55   1.993   0.0636
# lakelotr:sexunknown               NA         NA      NA       NA
# lakeopeongo:sexunknown            NA         NA      NA       NA
# lakeshirley:sexunknown            NA         NA      NA       NA
# age:lakelotr:sexmale              NA         NA      NA       NA
# age:lakeopeongo:sexmale           NA         NA      NA       NA
# age:lakeshirley:sexmale           NA         NA      NA       NA
# age:lakelotr:sexunknown           NA         NA      NA       NA
# age:lakeopeongo:sexunknown        NA         NA      NA       NA
# age:lakeshirley:sexunknown        NA         NA      NA       NA
# 
# (Intercept)                 
# age                         
# lakelotr                    
# lakeopeongo                 
# lakeshirley                 
# sexmale                     
# sexunknown                  
# age:lakelotr                
# age:lakeopeongo             
# age:lakeshirley             
# age:sexmale                 
# age:sexunknown              
# lakelotr:sexmale            
# lakeopeongo:sexmale         
# lakeshirley:sexmale        .
# lakelotr:sexunknown         
# lakeopeongo:sexunknown      
# lakeshirley:sexunknown      
# age:lakelotr:sexmale        
# age:lakeopeongo:sexmale     
# age:lakeshirley:sexmale     
# age:lakelotr:sexunknown     
# age:lakeopeongo:sexunknown  
# age:lakeshirley:sexunknown  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6269 on 16 degrees of freedom
# Multiple R-squared:  0.6066,	Adjusted R-squared:  0.2869 
# F-statistic: 1.897 on 13 and 16 DF,  p-value: 0.1124
forest_model(lmAgePCRliver6i)
report(lmAgePCRliver6i)
# We fitted a linear model (estimated using OLS) to predict
# CT_liver_ratio with age, lake and sex (formula: CT_liver_ratio ~
#                                          age * lake * sex). The model explains a statistically not
# significant and substantial proportion of variance (R2 = 0.61,
#                                                     F(13, 16) = 1.90, p = 0.112, adj. R2 = 0.29). The model's
# intercept, corresponding to age = 0, lake = hogan and sex =
# female, is at 8054.14 (95% CI [-7016.75, 23125.04], t(16) = 1.13,
# p = 0.274). Within this model:
# 
#   - The effect of age is statistically non-significant and negative
# (beta = -52.24, 95% CI [-1311.27, 1206.80], t(16) = -0.09, p =
# 0.931; Std. beta = -0.03, 95% CI [-0.83, 0.76])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 7437.47, 95% CI [-19818.73, 34693.66], t(16) =
# 0.58, p = 0.571; Std. beta = 1.30, 95% CI [-0.14, 2.73])
#   - The effect of lake [opeongo] is statistically non-significant
# and negative (beta = -12566.78, 95% CI [-39321.85, 14188.29],
# t(16) = -1.00, p = 0.334; Std. beta = 0.60, 95% CI [-1.56, 2.75])
#   - The effect of lake [shirley] is statistically non-significant
# and negative (beta = -33034.24, 95% CI [-83943.16, 17874.67],
# t(16) = -1.38, p = 0.188; Std. beta = -1.96, 95% CI [-4.84,
# 0.93])
#   - The effect of sex [male] is statistically non-significant and
# positive (beta = 16384.58, 95% CI [-14899.16, 47668.33], t(16) =
# 1.11, p = 0.283; Std. beta = 0.09, 95% CI [-1.82, 2.01])
#   - The effect of sex [unknown] is statistically non-significant
# and negative (beta = -4227.55, 95% CI [-19025.50, 10570.41],
# t(16) = -0.61, p = 0.553; Std. beta = -0.57, 95% CI [-2.56,
# 1.42])
#   - The effect of age × lake [lotr] is statistically
# non-significant and positive (beta = 196.70, 95% CI [-2065.56,
# 2458.96], t(16) = 0.18, p = 0.856; Std. beta = 0.12, 95% CI
# [-1.30, 1.55])
#   - The effect of age × lake [opeongo] is statistically
# non-significant and positive (beta = 1522.25, 95% CI [-1601.25,
# 4645.76], t(16) = 1.03, p = 0.317; Std. beta = 0.96, 95% CI
# [-1.01, 2.93])
#   - The effect of age × lake [shirley] is statistically
# non-significant and positive (beta = 1658.05, 95% CI [-1367.83,
# 4683.93], t(16) = 1.16, p = 0.262; Std. beta = 1.05, 95% CI
# [-0.86, 2.95])
#   - The effect of age × sex [male] is statistically non-significant
# and negative (beta = -1406.46, 95% CI [-3940.00, 1127.08], t(16)
# = -1.18, p = 0.256; Std. beta = -0.89, 95% CI [-2.49, 0.71])
#   - The effect of age × sex [unknown] is statistically
# non-significant and negative (beta = -3752.29, 95% CI [-22252.65,
# 14748.06], t(16) = -0.43, p = 0.673; Std. beta = -0.51, 95% CI
# [-3.00, 1.99])
#   - The effect of lake [lotr] × sex [male] is statistically
# non-significant and negative (beta = -3981.47, 95% CI [-26766.99,
# 18804.04], t(16) = -0.37, p = 0.716; Std. beta = -0.54, 95% CI
# [-3.61, 2.53])
#   - The effect of lake [opeongo] × sex [male] is statistically
# non-significant and positive (beta = 24365.92, 95% CI [-1555.34,
# 50287.17], t(16) = 1.99, p = 0.064; Std. beta = 3.28, 95% CI
# [-0.21, 6.77])
#   - The effect of lake [lotr] × sex [unknown] is statistically
# non-significant and negative (beta = -52.24, 95% CI [-1311.27,
# 1206.80], t(16) = -0.09, p = 0.931; Std. beta = -0.03, 95% CI
# [-0.83, 0.76])
#   - The effect of lake [opeongo] × sex [unknown] is statistically
# non-significant and positive (beta = 7437.47, 95% CI [-19818.73,
# 34693.66], t(16) = 0.58, p = 0.571; Std. beta = 1.30, 95% CI
# [-0.14, 2.73])
#   - The effect of lake [shirley] × sex [unknown] is statistically
# non-significant and negative (beta = -12566.78, 95% CI
# [-39321.85, 14188.29], t(16) = -1.00, p = 0.334; Std. beta =
# 0.60, 95% CI [-1.56, 2.75])
#   - The effect of (age × lake [lotr]) × sex [male] is statistically
# non-significant and negative (beta = -33034.24, 95% CI
# [-83943.16, 17874.67], t(16) = -1.38, p = 0.188; Std. beta =
# -1.96, 95% CI [-4.84, 0.93])
#   - The effect of (age × lake [opeongo]) × sex [male] is
# statistically non-significant and positive (beta = 16384.58, 95%
# CI [-14899.16, 47668.33], t(16) = 1.11, p = 0.283; Std. beta =
# 0.09, 95% CI [-1.82, 2.01])
#   - The effect of (age × lake [shirley]) × sex [male] is
# statistically non-significant and negative (beta = -4227.55, 95%
# CI [-19025.50, 10570.41], t(16) = -0.61, p = 0.553; Std. beta =
# -0.57, 95% CI [-2.56, 1.42])
#   - The effect of (age × lake [lotr]) × sex [unknown] is
# statistically non-significant and positive (beta = 196.70, 95% CI
# [-2065.56, 2458.96], t(16) = 0.18, p = 0.856; Std. beta = 0.12,
# 95% CI [-1.30, 1.55])
#   - The effect of (age × lake [opeongo]) × sex [unknown] is
# statistically non-significant and positive (beta = 1522.25, 95%
# CI [-1601.25, 4645.76], t(16) = 1.03, p = 0.317; Std. beta =
# 0.96, 95% CI [-1.01, 2.93])
#   - The effect of (age × lake [shirley]) × sex [unknown] is
# statistically non-significant and positive (beta = 1658.05, 95%
# CI [-1367.83, 4683.93], t(16) = 1.16, p = 0.262; Std. beta =
# 1.05, 95% CI [-0.86, 2.95])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.

# Generate predictions
predictions_liver <- df2 %>%
  mutate(predicted_liver_ratio = predict(lmAgePCRliver6i, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_liver_ratio, color = lake, shape = sex)) +
  geom_point() +
  geom_line(data = predictions_liver, aes(x = age, y = predicted_liver_ratio, color = lake, linetype = sex)) +
  labs(title = "CT Liver Ratio vs Age by Lake and Sex",
       x = "Age",
       y = "CT Liver Ratio",
       color = "Lake",
       shape = "Sex",
       linetype = "Sex") +
  theme_minimal()

#### liver ct vs liver mass and lake ####

lmPCRliver4 = lm(CT_liver_ratio~ liver_mass+lake, data = df2) #Create the linear regression
summary(lmPCRliver4)
# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass + lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13683.6  -2477.1   -161.6   2391.8  17780.7 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   4379.1     5293.7   0.827   0.4159  
# liver_mass     110.5      196.7   0.562   0.5791  
# lakelotr     12097.4     4737.1   2.554   0.0171 *
#   lakeopeongo   1203.5     3206.1   0.375   0.7105  
# lakeshirley  10784.2     4819.1   2.238   0.0344 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6212 on 25 degrees of freedom
# Multiple R-squared:  0.3965,	Adjusted R-squared:  0.2999 
# F-statistic: 4.106 on 4 and 25 DF,  p-value: 0.0108
check_model(lmPCRliver4)
forest_model(lmPCRliver4)
confint(lmPCRliver4)
# 2.5 %     97.5 %
#   (Intercept) -6523.5401 15281.7460
# liver_mass   -294.5671   515.6513
# lakelotr     2341.1515 21853.7301
# lakeopeongo -5399.5909  7806.5906
# lakeshirley   859.1333 20709.3203

# Fit linear models for each lake and extract formulas
formulaslm <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ liver_mass, data = .)) %>%
  mutate(formula = paste0("CT_liver_ratio = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * liver_mass"))

# Create a separate data frame with lake and formula
formulaslm_df <- formulaslm %>%
  select(lake, formula)
# Create the base plot
mlm <- ggplot(df2, aes(x = liver_mass, y = CT_liver_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "CT Liver Ratio vs Liver Mass",
       x = "Liver Mass (g)",
       y = "CT Liver Ratio",
       color = "Lake")
mlm
# Add formulas to each facet
mlm + geom_text(data = formulaslm_df, aes(x = Inf, y = Inf, label = formula), 
                hjust = 1.1, vjust = 2, size = 3, color = "black")


# Tidy the model output
tidy_model_PCRlm <- tidy(lmPCRliver4)
library(kableExtra)
# Create a table
tidy_model_PCRlm %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# 
# term	estimate	std.error	statistic	p.value
# (Intercept)	4379.1030	5293.731	0.8272243	0.4159364
# liver_mass	110.5421	196.699	0.5619862	0.5791297
# lakelotr	12097.4408	4737.124	2.5537523	0.0171351
# lakeopeongo	1203.4998	3206.102	0.3753779	0.7105456
# lakeshirley	10784.2268	4819.086	2.2378158	0.0343825

predictionsPCRlm <- df2 %>% mutate(predicted_PCRlm = predict(lmPCRliver4, newdata = .))

ggplot(df2, aes(x = liver_mass, y = CT_liver_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionsPCRlm, aes(x = liver_mass, y = predicted_PCRlm, color = lake)) + labs(title = "CT Liver Ratio vs Liver Mass by Lake", x = "Liver Mass (g)", y = "CT Liver Tissue", color = "Lake") +
  theme_minimal()
report(lmPCRliver4)
# We fitted a linear model (estimated using OLS) to predict
# CT_liver_ratio with liver_mass and lake (formula: CT_liver_ratio ~
#                                            liver_mass + lake). The model explains a statistically significant
# and substantial proportion of variance (R2 = 0.40, F(4, 25) =
#                                           4.11, p = 0.011, adj. R2 = 0.30). The model's intercept,
# corresponding to liver_mass = 0 and lake = hogan, is at 4379.10
# (95% CI [-6523.54, 15281.75], t(25) = 0.83, p = 0.416). Within
# this model:
# 
#   - The effect of liver mass is statistically non-significant and
# positive (beta = 110.54, 95% CI [-294.57, 515.65], t(25) = 0.56, p
# = 0.579; Std. beta = 0.16, 95% CI [-0.41, 0.72])
#   - The effect of lake [lotr] is statistically significant and
# positive (beta = 12097.44, 95% CI [2341.15, 21853.73], t(25) =
# 2.55, p = 0.017; Std. beta = 1.63, 95% CI [0.32, 2.94])
#   - The effect of lake [opeongo] is statistically non-significant
# and positive (beta = 1203.50, 95% CI [-5399.59, 7806.59], t(25) =
# 0.38, p = 0.711; Std. beta = 0.16, 95% CI [-0.73, 1.05])
#   - The effect of lake [shirley] is statistically significant and
# positive (beta = 10784.23, 95% CI [859.13, 20709.32], t(25) =
# 2.24, p = 0.034; Std. beta = 1.45, 95% CI [0.12, 2.79])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
approximation.


# Sort data frame by lake
df2 <- df2 %>% arrange(lake)

# Group by lakes and run linear models
model_lmPCRliver4 <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ liver_mass, data = .))

# View summaries of each model and include lake name
model_summaries_lmPCRliver4 <- model_lmPCRliver4 %>%
  rowwise() %>%
  mutate(lake = lake, summary = list(summary(model)))

# Print the summaries with lake names
for (i in 1:nrow(model_summaries_lmPCRliver4)) {
  cat("Lake:", model_summaries_lmPCRliver4$lake[i], "\n")
  print(model_summaries_lmPCRliver4$summary[[i]])
  cat("\n")
}
# Lake: 1 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass, data = .)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4130.4 -1960.7    61.7  2023.4  4264.5 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   4035.7     2888.7   1.397    0.205
# liver_mass     124.4      109.5   1.137    0.293
# 
# Residual standard error: 2992 on 7 degrees of freedom
# Multiple R-squared:  0.1558,	Adjusted R-squared:  0.0352 
# F-statistic: 1.292 on 1 and 7 DF,  p-value: 0.2931
# 
# 
# Lake: 2 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass, data = .)
# 
# Residuals:
#   1       2       3       4       5       6       7 
# 6806.7 -1338.5 -2479.8   101.5 10713.0 -9881.7 -3921.2 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  16484.1     7308.3   2.256   0.0738 .
# liver_mass     109.4     1006.1   0.109   0.9176  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7511 on 5 degrees of freedom
# Multiple R-squared:  0.002359,	Adjusted R-squared:  -0.1972 
# F-statistic: 0.01182 on 1 and 5 DF,  p-value: 0.9176
# 
# 
# Lake: 3 
# 
# Call:
#   lm( formula = CT_liver_ratio ~ liver_mass, data = .)
# 
# Residuals:
#   1     2     3     4     5     6     7 
# 2313 -5458  1676  4698 -2098 -2367  1236 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  6219.65    6432.90   0.967    0.378
# liver_mass     80.56     295.18   0.273    0.796
# 
# Residual standard error: 3783 on 5 degrees of freedom
# Multiple R-squared:  0.01468,	Adjusted R-squared:  -0.1824 
# F-statistic: 0.07448 on 1 and 5 DF,  p-value: 0.7958
# 
# 
# Lake: 4 
# 
# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass, data = .)
# 
# Residuals:
#   1      2      3      4      5      6      7 
# -13423  -4934  -2594   -540   1836  18164   1490 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 16270.31   12410.09   1.311    0.247
# liver_mass    -69.69    1915.30  -0.036    0.972
# 
# Residual standard error: 10460 on 5 degrees of freedom
# Multiple R-squared:  0.0002647,	Adjusted R-squared:  -0.1997 
# F-statistic: 0.001324 on 1 and 5 DF,  p-value: 0.9724

#### CT liver vs lake ####

lmPCRliver = lm(CT_liver_ratio~lake, data = df2) #Create the linear regression
summary(lmPCRliver)
# Call:
#   lm(formula = CT_liver_ratio ~ lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13523.6  -2679.1     46.5   2045.5  18015.8 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   7117.1     2043.1   3.483  0.00177 **
#   lakelotr     10099.4     3088.9   3.270  0.00303 **
#   lakeopeongo    814.3     3088.9   0.264  0.79416   
# lakeshirley   8725.2     3088.9   2.825  0.00897 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6129 on 26 degrees of freedom
# Multiple R-squared:  0.3889,	Adjusted R-squared:  0.3184 
# F-statistic: 5.515 on 3 and 26 DF,  p-value: 0.004559

check_model(lmPCRliver4)
forest_model(lmPCRliver4)
confint(lmPCRliver4)


# Fit linear models for each lake and extract formulas
formulaslm <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_liver_ratio ~ liver_mass, data = .)) %>%
  mutate(formula = paste0("liver_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * liver_mass"))

# Create the base plot
mlm <- ggplot(df2, aes(x = liver_mass, y = CT_liver_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "CT Liver Ratio vs Liver Mass",
       x = "Liver Mass (g)",
       y = "Ct Liver Ratio",
       color = "Lake")
mlm
# Add formulas to each facet
mlm + geom_text(data = formulaslm, aes(x = Inf, y = Inf, label = formulaslm), 
                hjust = 1.1, vjust = 2, size = 3, color = "black")
# Tidy the model output
tidy_model_PCRlm <- tidy(lmPCRliver4)
library(kableExtra)
# Create a table
tidy_model_PCRlm %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)


predictionsPCRlm <- df2 %>% mutate(predicted_PCRlm = predict(lmPCRliver4, newdata = .))

ggplot(df2, aes(x = liver_mass, y = CT_liver_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionsPCRlm, aes(x = liver_mass, y = predicted_PCRlm, color = lake)) + labs(title = "CT Liver Ratio vs Liver Mass by Lake", x = "Liver Mass (g)", y = "CT Liver Tissue", color = "Lake") +
  theme_minimal()

#### CT liver vs lake & sex ####

lmPCRlivers = lm(CT_liver_ratio~lake + sex, data = df2) #Create the linear regression
summary(lmPCRlivers)
# Call:
#   lm(formula = CT_liver_ratio ~ lake + sex, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11310.1  -3281.4    346.9   2841.9  17427.0 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     7065       2174   3.250   0.0034 **
#   lakelotr        7207       3556   2.027   0.0539 . 
# lakeopeongo    -2078       3556  -0.584   0.5644   
# lakeshirley     5244       3763   1.394   0.1762   
# sexmale         4122       2832   1.455   0.1585   
# sexunknown     -3657       6445  -0.567   0.5757   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6067 on 24 degrees of freedom
# Multiple R-squared:  0.4473,	Adjusted R-squared:  0.3322 
# F-statistic: 3.885 on 5 and 24 DF,  p-value: 0.01012


#### ct liver vs gonad mass and lake with interaction ####
lmPCRliver4ig = lm(CT_liver_ratio~ gonad_mass*lake, data = df2) #Create the linear regression
summary(lmPCRliver4ig)
# Call:
#   lm(formula = CT_liver_ratio ~ gonad_mass * lake, data = df2) 
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10146.0  -3848.5   -525.9   1825.6  18808.7 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)             6722.958   3240.622   2.075
# gonad_mass                 2.619     16.224   0.161
# lakelotr               10536.731   5248.712   2.007
# lakeopeongo             3786.035   7576.612   0.500
# lakeshirley            13886.767   5440.246   2.553
# gonad_mass:lakelotr       -4.190    122.675  -0.034
# gonad_mass:lakeopeongo   -34.042     79.789  -0.427
# gonad_mass:lakeshirley  -202.159    153.252  -1.319
# Pr(>|t|)  
# (Intercept)              0.0499 *
#   gonad_mass               0.8732  
# lakelotr                 0.0571 .
# lakeopeongo              0.6222  
# lakeshirley              0.0181 *
#   gonad_mass:lakelotr      0.9731  
# gonad_mass:lakeopeongo   0.6738  
# gonad_mass:lakeshirley   0.2007  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6393 on 22 degrees of freedom
# Multiple R-squared:  0.4375,	Adjusted R-squared:  0.2585 
# F-statistic: 2.445 on 7 and 22 DF,  p-value: 0.05146
report(lmPCRliver4ig)
# We fitted a linear model (estimated using OLS) to
# predict CT_liver_ratio with gonad_mass and lake
# (formula: CT_liver_ratio ~ gonad_mass * lake). The
# model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.44, F(7, 22)
#                                     = 2.44, p = 0.051, adj. R2 = 0.26). The model's
# intercept, corresponding to gonad_mass = 0 and lake =
# hogan, is at 6722.96 (95% CI [2.32, 13443.60], t(22) =
# 2.07, p = 0.050). Within this model:
# 
#   - The effect of gonad mass is statistically
# non-significant and positive (beta = 2.62, 95% CI
# [-31.03, 36.27], t(22) = 0.16, p = 0.873; Std. beta =
# 0.03, 95% CI [-0.39, 0.46])
#   - The effect of lake [lotr] is statistically
# non-significant and positive (beta = 10536.73, 95% CI
# [-348.43, 21421.89], t(22) = 2.01, p = 0.057; Std. beta
# = 1.38, 95% CI [-0.54, 3.29])
#   - The effect of lake [opeongo] is statistically
# non-significant and positive (beta = 3786.03, 95% CI
# [-11926.90, 19498.97], t(22) = 0.50, p = 0.622; Std.
# beta = 0.16, 95% CI [-0.81, 1.13])
#   - The effect of lake [shirley] is statistically
# significant and positive (beta = 13886.77, 95% CI
# [2604.39, 25169.15], t(22) = 2.55, p = 0.018; Std. beta
# = -0.21, 95% CI [-2.64, 2.22])
#   - The effect of gonad mass × lake [lotr] is
# statistically non-significant and negative (beta =
# -4.19, 95% CI [-258.60, 250.22], t(22) = -0.03, p =
# 0.973; Std. beta = -0.05, 95% CI [-3.25, 3.14])
#   - The effect of gonad mass × lake [opeongo] is
# statistically non-significant and negative (beta =
# -34.04, 95% CI [-199.51, 131.43], t(22) = -0.43, p =
# 0.674; Std. beta = -0.43, 95% CI [-2.51, 1.65])
#   - The effect of gonad mass × lake [shirley] is
# statistically non-significant and negative (beta =
# -202.16, 95% CI [-519.98, 115.67], t(22) = -1.32, p =
# 0.201; Std. beta = -2.54, 95% CI [-6.53, 1.45])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.
forest_model(lmPCRliver4ig)

# Generate predictions
predictionslmPCRliver4ig <- df2 %>%
  mutate(predicted_lmPCRliver4ig = predict(lmPCRliver4ig, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmPCRliver4ig, aes(x = gonad_mass, y = predicted_lmPCRliver4ig, color = lake)) +
  labs(title = "CT Liver Ratio vs Gonad Mass by Lake",
       x = "Gonad Mass (g)",
       y = "CT Liver Ratio",
       color = "Lake") +
  theme_minimal()

#### ct liver vs liver mass and lake with interaction ####
lmPCRliver4il = lm(CT_liver_ratio~ liver_mass*lake, data = df2) #Create the linear regression
summary(lmPCRliver4il)
# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass * lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13422.6  -2565.1   -239.2   1976.6  18164.0 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)             4035.70    6387.95   0.632
# liver_mass               124.41     242.04   0.514
# lakelotr               12448.45    9069.97   1.372
# lakeopeongo             2183.95   12939.40   0.169
# lakeshirley            12234.61   10121.52   1.209
# liver_mass:lakelotr      -15.00     918.87  -0.016
# liver_mass:lakeopeongo   -43.85     570.26  -0.077
# liver_mass:lakeshirley  -194.10    1235.63  -0.157
# Pr(>|t|)
# (Intercept)               0.534
# liver_mass                0.612
# lakelotr                  0.184
# lakeopeongo               0.868
# lakeshirley               0.240
# liver_mass:lakelotr       0.987
# liver_mass:lakeopeongo    0.939
# liver_mass:lakeshirley    0.877
# 
# Residual standard error: 6617 on 22 degrees of freedom
# Multiple R-squared:  0.3973,	Adjusted R-squared:  0.2055 
# F-statistic: 2.072 on 7 and 22 DF,  p-value: 0.09074
report(lmPCRliver4il)
# We fitted a linear model (estimated using OLS) to predict CT_liver_ratio
# with liver_mass and lake (formula: CT_liver_ratio ~ liver_mass * lake).
# The model explains a statistically not significant and substantial
# proportion of variance (R2 = 0.40, F(7, 22) = 2.07, p = 0.091, adj. R2 =
#                           0.21). The model's intercept, corresponding to liver_mass = 0 and lake =
# hogan, is at 4035.70 (95% CI [-9212.09, 17283.49], t(22) = 0.63, p =
# 0.534). Within this model:
# 
#   - The effect of liver mass is statistically non-significant and positive
# (beta = 124.41, 95% CI [-377.55, 626.37], t(22) = 0.51, p = 0.612; Std.
# beta = 0.17, 95% CI [-0.53, 0.88])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 12448.45, 95% CI [-6361.51, 31258.40], t(22) = 1.37, p
# = 0.184; Std. beta = 1.65, 95% CI [-0.78, 4.07])
#   - The effect of lake [opeongo] is statistically non-significant and
# positive (beta = 2183.95, 95% CI [-24650.71, 29018.62], t(22) = 0.17, p
# = 0.868; Std. beta = 0.20, 95% CI [-1.21, 1.61])
#   - The effect of lake [shirley] is statistically non-significant and
# positive (beta = 12234.61, 95% CI [-8756.14, 33225.37], t(22) = 1.21, p
# = 0.240; Std. beta = 1.25, 95% CI [-2.08, 4.57])
#   - The effect of liver mass × lake [lotr] is statistically
# non-significant and negative (beta = -15.00, 95% CI [-1920.62, 1890.62],
# t(22) = -0.02, p = 0.987; Std. beta = -0.02, 95% CI [-2.70, 2.66])
#   - The effect of liver mass × lake [opeongo] is statistically
# non-significant and negative (beta = -43.85, 95% CI [-1226.49, 1138.79],
# t(22) = -0.08, p = 0.939; Std. beta = -0.06, 95% CI [-1.72, 1.60])
#   - The effect of liver mass × lake [shirley] is statistically
# non-significant and negative (beta = -194.10, 95% CI [-2756.63,
# 2368.43], t(22) = -0.16, p = 0.877; Std. beta = -0.27, 95% CI [-3.87,
# 3.33])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.
forest_model(lmPCRliver4il)

# Generate predictions using the correct model
predictionslmPCRliver4il <- df2 %>%
  mutate(predicted_lmPCRliver4il = predict(lmPCRliver4il, newdata = .))

# Create the plot with the corrected variable
ggplot(df2, aes(x = liver_mass, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmPCRliver4il, aes(x = liver_mass, y = predicted_lmPCRliver4il, color = lake)) +
  labs(title = "CT Liver Ratio vs Liver Mass by Lake",
       x = "Liver Mass (g)",
       y = "CT Liver Ratio",
       color = "Lake") +
  theme_minimal()




#### ct RBC vs gonad mass and lake with interaction ####
lmPCRRBC4ig = lm(CT_RBC_ratio~ gonad_mass*lake, data = df2) #Create the linear regression
summary(lmPCRRBC4ig)
# Call:
#   lm(formula = CT_RBC_ratio ~ gonad_mass * lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -58782 -22803  -2194  18726  57274 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)            41374.61   17416.98   2.376
# gonad_mass               -37.59      88.47  -0.425
# lakelotr               47684.97   28694.56   1.662
# lakeopeongo            31653.90   40586.24   0.780
# lakeshirley            40101.77   29162.60   1.375
# gonad_mass:lakelotr     -969.96     656.99  -1.476
# gonad_mass:lakeopeongo    43.63     427.43   0.102
# gonad_mass:lakeshirley  -739.61     820.50  -0.901
# Pr(>|t|)  
# (Intercept)              0.0276 *
#   gonad_mass               0.6755  
# lakelotr                 0.1121  
# lakeopeongo              0.4446  
# lakeshirley              0.1843  
# gonad_mass:lakelotr      0.1554  
# gonad_mass:lakeopeongo   0.9197  
# gonad_mass:lakeshirley   0.3781  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 34220 on 20 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2959,	Adjusted R-squared:  0.04949 
# F-statistic: 1.201 on 7 and 20 DF,  p-value: 0.347
report(lmPCRRBC4ig)
# We fitted a linear model (estimated using OLS) to
# predict CT_RBC_ratio with gonad_mass and lake
# (formula: CT_RBC_ratio ~ gonad_mass * lake). The
# model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.30, F(7,
#                                                  20) = 1.20, p = 0.347, adj. R2 = 0.05). The model's
# intercept, corresponding to gonad_mass = 0 and lake
# = hogan, is at 41374.61 (95% CI [5043.43,
# 77705.80], t(20) = 2.38, p = 0.028). Within this
# model:
# 
#   - The effect of gonad mass is statistically
# non-significant and negative (beta = -37.59, 95% CI
# [-222.14, 146.96], t(20) = -0.42, p = 0.675; Std.
# beta = -0.10, 95% CI [-0.58, 0.38])
#   - The effect of lake [lotr] is statistically
# non-significant and positive (beta = 47684.97, 95%
# CI [-12170.83, 1.08e+05], t(20) = 1.66, p = 0.112;
# Std. beta = -0.66, 95% CI [-2.75, 1.44])
#   - The effect of lake [opeongo] is statistically
# non-significant and positive (beta = 31653.90, 95%
# CI [-53007.51, 1.16e+05], t(20) = 0.78, p = 0.445;
# Std. beta = 0.99, 95% CI [-0.14, 2.13])
#   - The effect of lake [shirley] is statistically
# non-significant and positive (beta = 40101.77, 95%
# CI [-20730.34, 1.01e+05], t(20) = 1.38, p = 0.184;
# Std. beta = -0.39, 95% CI [-3.02, 2.23])
#   - The effect of gonad mass × lake [lotr] is
# statistically non-significant and negative (beta =
# -969.96, 95% CI [-2340.41, 400.49], t(20) = -1.48,
# p = 0.155; Std. beta = -2.54, 95% CI [-6.13, 1.05])
#   - The effect of gonad mass × lake [opeongo] is
# statistically non-significant and positive (beta =
# 43.63, 95% CI [-847.96, 935.23], t(20) = 0.10, p =
# 0.920; Std. beta = 0.11, 95% CI [-2.22, 2.45])
#   - The effect of gonad mass × lake [shirley] is
# statistically non-significant and negative (beta =
# -739.61, 95% CI [-2451.13, 971.92], t(20) = -0.90,
# p = 0.378; Std. beta = -1.94, 95% CI [-6.42, 2.55])
# 
# Standardized parameters were obtained by fitting
# the model on a standardized version of the dataset.
# 95% Confidence Intervals (CIs) and p-values were
# computed using a Wald t-distribution approximation.
forest_model(lmPCRRBC4ig)

# Generate predictions
predictionslmPCRRBC4ig <- df2 %>%
  mutate(predicted_lmPCRRBC4ig = predict(lmPCRRBC4ig, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmPCRRBC4ig, aes(x = gonad_mass, y = predicted_lmPCRRBC4ig, color = lake)) +
  labs(title = "CT RBC Ratio vs Gonad Mass by Lake",
       x = "Gonad Mass (g)",
       y = "CT RBC Ratio",
       color = "Lake") +
  theme_minimal()


#### heart age ####

#### CT heart vs age by lake separated ####
# Group by lakes and run linear models
model_heart_age_bylake <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_heart_ratio ~ age, data = .))

# To view the models
model_heart_age_bylake

# View summaries of each model
model_summaries_heart_age_bylake <- model_heart_age_bylake %>%
  rowwise() %>%
  mutate(summary = list(summary(model)))

model_summaries_heart_age_bylake$summary
# [1]]
# 
# Call:
#   lm(formula = CT_heart_ratio ~ age, data = .)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -33957  -8175  -1612  13369  31348 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   110374      20724   5.326  0.00109 **
#   age            -2987       1804  -1.656  0.14172   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 19830 on 7 degrees of freedom
# Multiple R-squared:  0.2815,	Adjusted R-squared:  0.1788 
# F-statistic: 2.742 on 1 and 7 DF,  p-value: 0.1417
# 
# 
# [[2]]
# 
# Call:
#   lm(formula = CT_heart_ratio ~ age, data = .)
# 
# Residuals:
#   1       2       3       4       5       6       7 
# 38112 -102027  -39162  -18671  183480   22382  -84114 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   248785      97523   2.551   0.0512 .
# age           -11264       9596  -1.174   0.2933  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 104900 on 5 degrees of freedom
# Multiple R-squared:  0.2161,	Adjusted R-squared:  0.05927 
# F-statistic: 1.378 on 1 and 5 DF,  p-value: 0.2933
# 
# 
# [[3]]
# 
# Call:
#   lm(formula = CT_heart_ratio ~ age, data = .)
# 
# Residuals:
#   1      2      3      4      5      6      7 
# -32718 -55963  15307 113710 -59702   3246  16119 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept) 130051.4    64663.7   2.011    0.100
# age            142.6     4924.4   0.029    0.978
# 
# Residual standard error: 65120 on 5 degrees of freedom
# Multiple R-squared:  0.0001676,	Adjusted R-squared:  -0.1998 
# F-statistic: 0.0008383 on 1 and 5 DF,  p-value: 0.978
# 
# 
# [[4]]
# 
# Call:
#   lm(formula = CT_heart_ratio ~ age, data = .)
# 
# Residuals:
#   1         2         3         4         5         6 
# -5312.04  51129.76  20539.14 -35794.38 -13873.50 -16768.11 
# 7 
# 79.14 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    63602      31232   2.036   0.0973 .
# age             3810       2329   1.636   0.1627  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 31050 on 5 degrees of freedom
# Multiple R-squared:  0.3487,	Adjusted R-squared:  0.2185 
# F-statistic: 2.677 on 1 and 5 DF,  p-value: 0.1627

#### CT heart vs age and lake with interaction ####
# not significant for any of the variables
# including the intercept
lmAgePCRheart2i = lm(CT_heart_ratio~age*lake, data = df2) #Create the linear regression
summary(lmAgePCRheart2i)

# Call:
#   lm(formula = CT_heart_ratio ~ age * lake, data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -102027  -29206   -1673   15916  183480 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)       110374      64476   1.712    0.101
# age                -2987       5612  -0.532    0.600
# lakelotr          138411      86312   1.604    0.123
# lakeopeongo        19677      88947   0.221    0.827
# lakeshirley       -46772      89501  -0.523    0.606
# age:lakelotr       -8278       7961  -1.040    0.310
# age:lakeopeongo     3130       7298   0.429    0.672
# age:lakeshirley     6797       7274   0.934    0.360
# 
# Residual standard error: 61700 on 22 degrees of freedom
# Multiple R-squared:  0.319,	Adjusted R-squared:  0.1023 
# F-statistic: 1.472 on 7 and 22 DF,  p-value: 0.2284

report(lmAgePCRheart2i)
# We fitted a linear model (estimated using OLS) to predict
# CT_heart_ratio with age and lake (formula: CT_heart_ratio ~ age *
#                                     lake). The model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.32, F(7, 22) = 1.47, p
#                                     = 0.228, adj. R2 = 0.10). The model's intercept, corresponding to
# age = 0 and lake = hogan, is at 1.10e+05 (95% CI [-23341.50,
# 2.44e+05], t(22) = 1.71, p = 0.101). Within this model:
# 
#   - The effect of age is statistically non-significant and negative
# (beta = -2986.91, 95% CI [-14625.34, 8651.53], t(22) = -0.53, p =
# 0.600; Std. beta = -0.21, 95% CI [-1.05, 0.62])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 1.38e+05, 95% CI [-40589.85, 3.17e+05], t(22) =
# 1.60, p = 0.123; Std. beta = 0.71, 95% CI [-0.34, 1.75])
#   - The effect of lake [opeongo] is statistically non-significant
# and positive (beta = 19677.42, 95% CI [-1.65e+05, 2.04e+05],
# t(22) = 0.22, p = 0.827; Std. beta = 0.84, 95% CI [-0.16, 1.84])
#   - The effect of lake [shirley] is statistically non-significant
# and negative (beta = -46772.38, 95% CI [-2.32e+05, 1.39e+05],
# t(22) = -0.52, p = 0.606; Std. beta = 0.45, 95% CI [-0.56, 1.46])
#   - The effect of age × lake [lotr] is statistically
# non-significant and negative (beta = -8277.60, 95% CI [-24787.07,
# 8231.88], t(22) = -1.04, p = 0.310; Std. beta = -0.60, 95% CI
# [-1.78, 0.59])
#   - The effect of age × lake [opeongo] is statistically
# non-significant and positive (beta = 3129.48, 95% CI [-12006.57,
# 18265.53], t(22) = 0.43, p = 0.672; Std. beta = 0.23, 95% CI
# [-0.86, 1.31])
#   - The effect of age × lake [shirley] is statistically
# non-significant and positive (beta = 6797.32, 95% CI [-8288.91,
# 21883.55], t(22) = 0.93, p = 0.360; Std. beta = 0.49, 95% CI
# [-0.60, 1.57])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.
forest_model(lmAgePCRheart2i)

# Generate predictions
predictionslmAgePCRheart2i <- df2 %>%
  mutate(predicted_lmAgePCRheart2i = predict(lmAgePCRheart2i, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmAgePCRheart2i, aes(x = age, y = predicted_lmAgePCRheart2i, color = lake)) +
  labs(title = "Telomere Ratio Heart vs Age by Lake",
       x = "Age (y)",
       y = "Telomere Heart Ratio",
       color = "Lake") +
  theme_minimal()


#### ct heart vs ventricle mass and lake ####
lmPCRheart4 = lm(CT_heart_ratio~ ventricle_mass+lake, data = df2) #Create the linear regression
summary(lmPCRheart4)
# Call:
#   lm(formula = CT_heart_ratio ~ ventricle_mass + lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -95897 -31262  -3923  17073 221315 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)       61651      51631   1.194   0.2437  
# ventricle_mass    11556      33552   0.344   0.7334  
# lakelotr          76224      43143   1.767   0.0895 .
# lakeopeongo       51386      33042   1.555   0.1325  
# lakeshirley       43342      43816   0.989   0.3320  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 63900 on 25 degrees of freedom
# Multiple R-squared:  0.1699,	Adjusted R-squared:  0.03711 
# F-statistic: 1.279 on 4 and 25 DF,  p-value: 0.3046
check_model(lmPCRheart4)
forest_model(lmPCRheart4)
confint(lmPCRheart4)
# 2.5 %    97.5 %
#   (Intercept)    -44684.68 167987.31
# ventricle_mass -57545.36  80656.86
# lakelotr       -12630.00 165077.18
# lakeopeongo    -16665.38 119438.01
# lakeshirley    -46897.39 133581.84
report(lmPCRheart4)
# We fitted a linear model (estimated using OLS) to predict
# CT_heart_ratio with ventricle_mass and lake (formula:
#                                                CT_heart_ratio ~ ventricle_mass + lake). The model explains a
# statistically not significant and moderate proportion of
# variance (R2 = 0.17, F(4, 25) = 1.28, p = 0.305, adj. R2 =
#             0.04). The model's intercept, corresponding to ventricle_mass
# = 0 and lake = hogan, is at 61651.32 (95% CI [-44684.68,
# 1.68e+05], t(25) = 1.19, p = 0.244). Within this model:
# 
#   - Th e effect of ventricle mass is statistically
# non-significant and positive (beta = 11555.75, 95% CI
# [-57545.36, 80656.86], t(25) = 0.34, p = 0.733; Std. beta =
# 0.11, 95% CI [-0.54, 0.76])
#   - The effect of lake [lotr] is statistically non-significant
# and positive (beta = 76223.59, 95% CI [-12630.00, 1.65e+05],
# t(25) = 1.77, p = 0.089; Std. beta = 1.17, 95% CI [-0.19,
# 2.53])
#   - The effect of lake [opeongo] is statistically
# non-significant and positive (beta = 51386.31, 95% CI
# [-16665.38, 1.19e+05], t(25) = 1.56, p = 0.132; Std. beta =
# 0.79, 95% CI [-0.26, 1.83])
#   - The effect of lake [shirley] is statistically
# non-significant and positive (beta = 43342.23, 95% CI
# [-46897.39, 1.34e+05], t(25) = 0.99, p = 0.332; Std. beta =
# 0.67, 95% CI [-0.72, 2.05])
# 
# Standardized parameters were obtained by fitting the model on
# a standardized version of the dataset. 95% Confidence
# Intervals (CIs) and p-values were computed using a Wald
# t-distribution approximation.
# Fit linear models for each lake and extract formulas
formulashlm <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_heart_ratio ~ ventricle_mass, data = .)) %>%
  mutate(formula = paste0("ventricle_mass = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * ventricle_mass"))


formulashlm
# Create the base plot
hlm <- ggplot(df2, aes(x = ventricle_mass, y = CT_heart_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "CT Heart Ratio vs Liver Mass",
       x = "Ventricle Mass (g)",
       y = "Ct Heart Ratio",
       color = "Lake")
hlm
# Add formulas to each facet
hlm + geom_text(data = formulashlm, aes(x = Inf, y = Inf, label = formulashlm), 
                hjust = 1.1, vjust = 2, size = 3, color = "black")
# Tidy the model output
tidy_model_PCRhlm <- tidy(lmPCRheart4)
library(kableExtra)
# Create a table
tidy_model_PCRhlm %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)
# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	61651.32	51630.98	1.1940761	0.2436555
# ventricle_mass	11555.75	33551.74	0.3444157	0.7334141
# lakelotr	76223.59	43142.47	1.7667876	0.0894722
# lakeopeongo	51386.31	33042.20	1.5551722	0.1324747
# lakeshirley	43342.23	43815.45	0.9891996	0.3320436

predictionsPCRhlm <- df2 %>% mutate(predicted_PCRhlm = predict(lmPCRheart4, newdata = .))

ggplot(df2, aes(x = ventricle_mass, y = CT_heart_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionsPCRhlm, aes(x = ventricle_mass, y = predicted_PCRhlm, color = lake)) + labs(title = "CT Heart Ratio vs Ventricle Mass by Lake", x = "Ventricle Mass (g)", y = "CT Heart Tissue", color = "Lake") +
  theme_minimal()

#### ct heart vs ventricle mass and lake with interaction ####
lmPCRheart4i = lm(CT_heart_ratio~ ventricle_mass*lake, data = df2) #Create the linear regression
summary(lmPCRheart4i)
# Call:
#   lm(formula = CT_heart_ratio ~ ventricle_mass * lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -71178 -33794 -10059  15812 193928 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)                  127998      78111   1.639
# ventricle_mass               -35775      53745  -0.666
# lakelotr                     220112     149084   1.476
# lakeopeongo                  -87387     107860  -0.810
# lakeshirley                  -24820     107480  -0.231
# ventricle_mass:lakelotr     -337615     234766  -1.438
# ventricle_mass:lakeopeongo    91979      69160   1.330
# ventricle_mass:lakeshirley    50846     145895   0.349
# Pr(>|t|)
# (Intercept)                   0.116
# ventricle_mass                0.513
# lakelotr                      0.154
# lakeopeongo                   0.427
# lakeshirley                   0.820
# ventricle_mass:lakelotr       0.164
# ventricle_mass:lakeopeongo    0.197
# ventricle_mass:lakeshirley    0.731
# 
# Residual standard error: 61880 on 22 degrees of freedom
# Multiple R-squared:  0.3152,	Adjusted R-squared:  0.09726 
# F-statistic: 1.446 on 7 and 22 DF,  p-value: 0.2375

forest_model(lmPCRheart4i)
confint(lmPCRheart4i)
# 2.5 %    97.5 %
#   (Intercept)                 -33992.85 289989.65
# ventricle_mass             -147234.80  75684.94
# lakelotr                    -89068.67 529292.72
# lakeopeongo                -311073.46 136300.21
# lakeshirley                -247718.81 198079.23
# ventricle_mass:lakelotr    -824490.83 149260.12
# ventricle_mass:lakeopeongo  -51450.79 235408.65
# ventricle_mass:lakeshirley -251721.04 353413.10
report(lmPCRheart4i)
# We fitted a linear model (estimated using OLS) to
# predict CT_heart_ratio with ventricle_mass and lake
# (formula: CT_heart_ratio ~ ventricle_mass * lake). The
# model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.32, F(7, 22)
#                                     = 1.45, p = 0.238, adj. R2 = 0.10). The model's
# intercept, corresponding to ventricle_mass = 0 and lake
# = hogan, is at 1.28e+05 (95% CI [-33992.85, 2.90e+05],
# t(22) = 1.64, p = 0.116). Within this model:
# 
#   - The effect of ventricle mass is statistically
# non-significant and negative (beta = -35774.93, 95% CI
# [-1.47e+05, 75684.94], t(22) = -0.67, p = 0.513; Std.
# beta = -0.34, 95% CI [-1.38, 0.71])
#   - The effect of lake [lotr] is statistically
# non-significant and positive (beta = 2.20e+05, 95% CI
# [-89068.67, 5.29e+05], t(22) = 1.48, p = 0.154; Std.
# beta = -2.05, 95% CI [-5.87, 1.78])
#   - The effect of lake [opeongo] is statistically
# non-significant and negative (beta = -87386.63, 95% CI
# [-3.11e+05, 1.36e+05], t(22) = -0.81, p = 0.427; Std.
# beta = 0.14, 95% CI [-1.27, 1.55])
#   - The effect of lake [shirley] is statistically
# non-significant and negative (beta = -24819.79, 95% CI
# [-2.48e+05, 1.98e+05], t(22) = -0.23, p = 0.820; Std.
# beta = 0.44, 95% CI [-2.13, 3.01])
#   - The effect of ventricle mass × lake [lotr] is
# statistically non-significant and negative (beta =
# -3.38e+05, 95% CI [-8.24e+05, 1.49e+05], t(22) = -1.44,
# p = 0.164; Std. beta = -3.16, 95% CI [-7.73, 1.40])
#   - The effect of ventricle mass × lake [opeongo] is
# statistically non-significant and positive (beta =
# 91978.93, 95% CI [-51450.79, 2.35e+05], t(22) = 1.33, p
# = 0.197; Std. beta = 0.86, 95% CI [-0.48, 2.21])
#   - The effect of ventricle mass × lake [shirley] is
# statistically non-significant and positive (beta =
# 50846.03, 95% CI [-2.52e+05, 3.53e+05], t(22) = 0.35, p
# = 0.731; Std. beta = 0.48, 95% CI [-2.36, 3.31])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.

# Generate predictions
predictionslmPCRheart4i <- df2 %>%
  mutate(predicted_lmPCRheart4i = predict(lmPCRheart4i, newdata = .))

# Create the plot
ggplot(df2, aes(x = ventricle_mass, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmPCRheart4i, aes(x = ventricle_mass, y = predicted_lmPCRheart4i, color = lake)) +
  labs(title = "CT Heart Ratio vs Ventricle Mass by Lake",
       x = "Ventricle Mass (g)",
       y = "CT Heart Ratio",
       color = "Lake") +
  theme_minimal()

#### ct heart vs gonad mass and lake with interaction ####
lmPCRheart4ig = lm(CT_heart_ratio~ gonad_mass*lake, data = df2) #Create the linear regression
summary(lmPCRheart4ig)
# Call:
#   lm(formula = CT_heart_ratio ~ gonad_mass * lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -74292 -30515  -6569  22012 193331 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)             84591.1    31065.0   2.723
# gonad_mass                -44.8      155.5  -0.288
# lakelotr               127552.6    50314.8   2.535
# lakeopeongo             96259.1    72630.3   1.325
# lakeshirley             14800.7    52150.8   0.284
# gonad_mass:lakelotr     -2423.4     1176.0  -2.061
# gonad_mass:lakeopeongo   -553.4      764.9  -0.723
# gonad_mass:lakeshirley    529.0     1469.1   0.360
# Pr(>|t|)  
# (Intercept)              0.0124 *
#   gonad_mass               0.7760  
# lakelotr                 0.0189 *
#   lakeopeongo              0.1987  
# lakeshirley              0.7792  
# gonad_mass:lakelotr      0.0513 .
# gonad_mass:lakeopeongo   0.4770  
# gonad_mass:lakeshirley   0.7222  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 61280 on 22 degrees of freedom
# Multiple R-squared:  0.3282,	Adjusted R-squared:  0.1145 
# F-statistic: 1.536 on 7 and 22 DF,  p-value: 0.2072
report(lmPCRheart4ig)
# We fitted a linear model (estimated using OLS) to
# predict CT_heart_ratio with gonad_mass and lake
# (formula: CT_heart_ratio ~ gonad_mass * lake). The
# model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.33, F(7, 22)
#                                     = 1.54, p = 0.207, adj. R2 = 0.11). The model's
# intercept, corresponding to gonad_mass = 0 and lake =
# hogan, is at 84591.13 (95% CI [20166.28, 1.49e+05],
# t(22) = 2.72, p = 0.012). Within this model:
# 
#   - The effect of gonad mass is statistically
# non-significant and negative (beta = -44.80, 95% CI
# [-367.34, 277.75], t(22) = -0.29, p = 0.776; Std. beta
# = -0.06, 95% CI [-0.53, 0.40])
#   - The effect of lake [lotr] is statistically
# significant and positive (beta = 1.28e+05, 95% CI
# [23206.08, 2.32e+05], t(22) = 2.54, p = 0.019; Std.
# beta = -0.88, 95% CI [-2.97, 1.21])
#   - The effect of lake [opeongo] is statistically
# non-significant and positive (beta = 96259.10, 95% CI
# [-54366.94, 2.47e+05], t(22) = 1.33, p = 0.199; Std.
# beta = 0.83, 95% CI [-0.23, 1.89])
#   - The effect of lake [shirley] is statistically
# non-significant and positive (beta = 14800.67, 95% CI
# [-93353.57, 1.23e+05], t(22) = 0.28, p = 0.779; Std.
# beta = 0.85, 95% CI [-1.81, 3.50])
#   - The effect of gonad mass × lake [lotr] is
# statistically non-significant and negative (beta =
# -2423.38, 95% CI [-4862.21, 15.45], t(22) = -2.06, p =
# 0.051; Std. beta = -3.47, 95% CI [-6.96, 0.02])
#   - The effect of gonad mass × lake [opeongo] is
# statistically non-significant and negative (beta =
# -553.38, 95% CI [-2139.61, 1032.85], t(22) = -0.72, p =
# 0.477; Std. beta = -0.79, 95% CI [-3.06, 1.48])
#   - The effect of gonad mass × lake [shirley] is
# statistically non-significant and positive (beta =
# 528.96, 95% CI [-2517.74, 3575.66], t(22) = 0.36, p =
# 0.722; Std. beta = 0.76, 95% CI [-3.60, 5.12])
# 
# Standardized parameters were obtained by fitting the
# model on a standardized version of the dataset. 95%
# Confidence Intervals (CIs) and p-values were computed
# using a Wald t-distribution approximation.
forest_model(lmPCRheart4ig)

# Generate predictions
predictionslmPCRheart4ig <- df2 %>%
  mutate(predicted_lmPCRheart4ig = predict(lmPCRheart4ig, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslmPCRheart4ig, aes(x = gonad_mass, y = predicted_lmPCRheart4ig, color = lake)) +
  labs(title = "CT Heart Ratio vs Gonad Mass by Lake",
       x = "Gonad Mass (g)",
       y = "CT Heart Ratio",
       color = "Lake") +
  theme_minimal()

#### ct heart vs age and lake and sex ####
lmPCRheart4 = lm(CT_heart_ratio~ age*lake*sex, data = df2) #Create the linear regression
summary(lmPCRheart4)
# Call:
#   lm(formula = CT_heart_ratio ~ age * lake * sex, data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -109499  -23339       0   13986  153941 
# 
# Coefficients: (10 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)                  97710.8    73143.6   1.336    0.200
# age                          -2096.1     6110.5  -0.343    0.736
# lakelotr                      -967.3   132282.4  -0.007    0.994
# lakeopeongo                 -40547.7   129850.3  -0.312    0.759
# lakeshirley                -219576.8   247076.2  -0.889    0.387
# sexmale                     104322.7   151829.4   0.687    0.502
# sexunknown                   36884.4    71818.9   0.514    0.615
# age:lakelotr                  -832.7    10979.4  -0.076    0.940
# age:lakeopeongo               9805.3    15159.3   0.647    0.527
# age:lakeshirley              16503.9    14685.5   1.124    0.278
# age:sexmale                 -10412.6    12296.0  -0.847    0.410
# age:sexunknown                    NA         NA      NA       NA
# lakelotr:sexmale             89720.0    89787.7   0.999    0.333
# lakeopeongo:sexmale          16947.0   110584.9   0.153    0.880
# lakeshirley:sexmale          79874.5   125803.6   0.635    0.534
# lakelotr:sexunknown               NA         NA      NA       NA
# lakeopeongo:sexunknown            NA         NA      NA       NA
# lakeshirley:sexunknown            NA         NA      NA       NA
# age:lakelotr:sexmale              NA         NA      NA       NA
# age:lakeopeongo:sexmale           NA         NA      NA       NA
# age:lakeshirley:sexmale           NA         NA      NA       NA
# age:lakelotr:sexunknown           NA         NA      NA       NA
# age:lakeopeongo:sexunknown        NA         NA      NA       NA
# age:lakeshirley:sexunknown        NA         NA      NA       NA
# 
# Residual standard error: 64500 on 16 degrees of freedom
# Multiple R-squared:  0.4588,	Adjusted R-squared:  0.01899 
# F-statistic: 1.043 on 13 and 16 DF,  p-value: 0.4613
check_model(lmPCRheart4)
forest_model(lmPCRheart4)
confint(lmPCRheart4)

# Generate predictions
predictions_heart <- df2 %>%
  mutate(predicted_heart_ratio = predict(lmPCRheart4, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_heart_ratio, color = lake, shape = sex)) +
  geom_point() +
  geom_line(data = predictions_heart, aes(x = age, y = predicted_heart_ratio, color = lake, linetype = sex)) +
  labs(title = "CT Liver Heart vs Age by Lake and Sex",
       x = "Age",
       y = "CT Heart Ratio",
       color = "Lake",
       shape = "Sex",
       linetype = "Sex") +
  theme_minimal()


#### ct heart vs age ####
lmAgePCRheart1 = lm(CT_heart_ratio~age, data = df2) #Create the linear regression
summary(lmAgePCRheart1)
# Call:
#   lm(formula = CT_heart_ratio ~ age, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -66098 -44147 -18592  24317 241159 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   134995      31438   4.294  0.00019 ***
#   age            -1913       2603  -0.735  0.46855    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 65650 on 28 degrees of freedom
# Multiple R-squared:  0.01892,	Adjusted R-squared:  -0.01612 
# F-statistic:  0.54 on 1 and 28 DF,  p-value: 0.4685

report(lmAgePCRheart1)
# We fitted a linear model (estimated using OLS) to predict CT_heart_ratio
# with age (formula: CT_heart_ratio ~ age). The model explains a
# statistically not significant and very weak proportion of variance (R2 =
#                                                                       0.02, F(1, 28) = 0.54, p = 0.469, adj. R2 = -0.02). The model's
# intercept, corresponding to age = 0, is at 1.35e+05 (95% CI [70597.97,
# 1.99e+05], t(28) = 4.29, p < .001). Within this model:
# 
#   - The effect of age is statistically non-significant and negative (beta
# = -1912.56, 95% CI [-7243.93, 3418.82], t(28) = -0.73, p = 0.469; Std.
# beta = -0.14, 95% CI [-0.52, 0.25])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.
# p-values were computed using a Wald t-distribution approximation.

lmAgePCRheart3 = lm(CT_heart_ratio~age+lake+sex, data = df2) #Create the linear regression
summary(lmAgePCRheart3)
# Call:
#   lm(formula = CT_heart_ratio ~ age + lake + sex, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -83157 -40950  -5625  27157 203690 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    86283      37097   2.326   0.0292 *
#   age            -1638       2622  -0.625   0.5384  
# lakelotr       44561      37555   1.187   0.2475  
# lakeopeongo    36838      37135   0.992   0.3315  
# lakeshirley    10772      39321   0.274   0.7866  
# sexmale        39972      29569   1.352   0.1896  
# sexunknown     44647      67757   0.659   0.5165  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 63300 on 23 degrees of freedom
# Multiple R-squared:  0.2508,	Adjusted R-squared:  0.05533 
# F-statistic: 1.283 on 6 and 23 DF,  p-value: 0.3037
lmAgePCRheart2 = lm(CT_heart_ratio~age+lake, data = df2) #Create the linear regression
summary(lmAgePCRheart2)
# Call:
#   lm(formula = CT_heart_ratio ~ age + lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -82930 -32016  -6612  19810 214927 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    96295      35430   2.718   0.0118 *
#   age            -1694       2609  -0.649   0.5221  
# lakelotr       63621      32284   1.971   0.0599 .
# lakeopeongo    56057      32179   1.742   0.0938 .
# lakeshirley    35718      32263   1.107   0.2788  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 63520 on 25 degrees of freedom
# Multiple R-squared:  0.1798,	Adjusted R-squared:  0.04858 
# F-statistic:  1.37 on 4 and 25 DF,  p-value: 0.2726

#### RBC ####

#### rbc vs age ####
lmAgePCRRBC1 = lm(CT_RBC_ratio~age, data = df2) #Create the linear regression
summary(lmAgePCRRBC1)
# Call:
#   lm(formula = CT_RBC_ratio ~ age, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -55228 -24920  -7049  23260  76393 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    80965      17167   4.716 7.11e-05 ***
#   age            -2058       1395  -1.476    0.152    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 34360 on 26 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.07728,	Adjusted R-squared:  0.04179 
# F-statistic: 2.177 on 1 and 26 DF,  p-value: 0.1521
report(lmAgePCRRBC1)

#### ct RBCs vs age lake and sex ####
lmAgePCRRBC3 = lm(CT_RBC_ratio~age+lake+sex, data = df2) #Create the linear regression
summary(lmAgePCRRBC3)
# Call:
#   lm(formula = CT_RBC_ratio ~ age + lake + sex, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -62282 -21251    608  17849  60420 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  63430.2    20813.2   3.048  0.00612 **
#   age          -2431.9     1417.7  -1.715  0.10100   
# lakelotr     25324.0    22253.3   1.138  0.26794   
# lakeopeongo  42775.7    20543.2   2.082  0.04974 * 
#   lakeshirley  33484.0    21820.2   1.535  0.13983   
# sexmale      -4412.5    16961.0  -0.260  0.79728   
# sexunknown    -615.1    36561.4  -0.017  0.98674   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 33820 on 21 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2781,	Adjusted R-squared:  0.07179 
# F-statistic: 1.348 on 6 and 21 DF,  p-value: 0.2806

check_model(lmAgePCRRBC3)
forest_model(lmAgePCRRBC3)
confint(lmAgePCRRBC3)

# 2.5 %      97.5 %
#   (Intercept)  20146.91856 106713.5747
# age          -5380.13130    516.3806
# lakelotr    -20954.38542  71602.3141
# lakeopeongo     53.69187  85497.7167
# lakeshirley -11893.52295  78861.5493
# sexmale     -39684.74494  30859.7529
# sexunknown  -76648.73222  75418.6272

report(lmAgePCRRBC3)
# We fitted a linear model (estimated using OLS) to predict CT_RBC_ratio
# with age, lake and sex (formula: CT_RBC_ratio ~ age + lake + sex). The
# model explains a statistically not significant and substantial
# proportion of variance (R2 = 0.28, F(6, 21) = 1.35, p = 0.281, adj. R2 =
#                           0.07). The model's intercept, corresponding to age = 0, lake = hogan and
# sex = female, is at 63430.25 (95% CI [20146.92, 1.07e+05], t(21) = 3.05,
# p = 0.006). Within this model:
# 
#   - The effect of age is statistically non-significant and negative (beta
# = -2431.88, 95% CI [-5380.13, 516.38], t(21) = -1.72, p = 0.101; Std.
# beta = -0.33, 95% CI [-0.73, 0.07])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 25323.96, 95% CI [-20954.39, 71602.31], t(21) = 1.14, p
# = 0.268; Std. beta = 0.72, 95% CI [-0.60, 2.04])
#   - The effect of lake [opeongo] is statistically significant and positive
# (beta = 42775.70, 95% CI [53.69, 85497.72], t(21) = 2.08, p = 0.050;
# Std. beta = 1.22, 95% CI [1.53e-03, 2.44])
#   - The effect of lake [shirley] is statistically non-significant and
# positive (beta = 33484.01, 95% CI [-11893.52, 78861.55], t(21) = 1.53, p
# = 0.140; Std. beta = 0.95, 95% CI [-0.34, 2.25])
#   - The effect of sex [male] is statistically non-significant and negative
# (beta = -4412.50, 95% CI [-39684.74, 30859.75], t(21) = -0.26, p =
# 0.797; Std. beta = -0.13, 95% CI [-1.13, 0.88])
#   - The effect of sex [unknown] is statistically non-significant and
# negative (beta = -615.05, 95% CI [-76648.73, 75418.63], t(21) = -0.02, p
# = 0.987; Std. beta = -0.02, 95% CI [-2.18, 2.15])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.

forest_model(lmAgePCRRBC3)
confint(lmAgePCRRBC3)

# Generate predictions
predictions_RBC <- df2 %>%
  mutate(predicted_RBC_ratio = predict(lmAgePCRRBC3, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake, shape = sex)) +
  geom_point() +
  geom_line(data = predictions_RBC, aes(x = age, y = predicted_RBC_ratio, color = lake, linetype = sex)) +
  labs(title = "CT Ratio RBC vs Age by Lake and Sex",
       x = "Age",
       y = "CT RBC Ratio",
       color = "Lake",
       shape = "Sex",
       linetype = "Sex") +
  theme_minimal()

#### ct RBCs vs age lake and sex interaction ####
lmAgePCRRBC3i = lm(CT_RBC_ratio~age*lake*sex, data = df2) #Create the linear regression
summary(lmAgePCRRBC3i)
# Call:
#   lm(formula = CT_RBC_ratio ~ age * lake * sex, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -62286 -17812      0  19197  63748 
# 
# Coefficients: (11 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)                   62159      44067   1.411    0.179
# age                           -1960       3592  -0.546    0.593
# lakelotr                     123517     181502   0.681    0.507
# lakeopeongo                   64981      86376   0.752    0.464
# lakeshirley                  112651     141836   0.794    0.439
# sexmale                      -89713     140284  -0.640    0.532
# sexunknown                    -3123      42502  -0.073    0.942
# age:lakelotr                  -8333      11523  -0.723    0.481
# age:lakeopeongo               -3272      11123  -0.294    0.773
# age:lakeshirley               -4983      10126  -0.492    0.630
# age:sexmale                    5128       9947   0.516    0.614
# age:sexunknown                   NA         NA      NA       NA
# lakelotr:sexmale              17655      65820   0.268    0.792
# lakeopeongo:sexmale           30772      84572   0.364    0.721
# lakeshirley:sexmale              NA         NA      NA       NA
# lakelotr:sexunknown              NA         NA      NA       NA
# lakeopeongo:sexunknown           NA         NA      NA       NA
# lakeshirley:sexunknown           NA         NA      NA       NA
# age:lakelotr:sexmale             NA         NA      NA       NA
# age:lakeopeongo:sexmale          NA         NA      NA       NA
# age:lakeshirley:sexmale          NA         NA      NA       NA
# age:lakelotr:sexunknown          NA         NA      NA       NA
# age:lakeopeongo:sexunknown       NA         NA      NA       NA
# age:lakeshirley:sexunknown       NA         NA      NA       NA
# 
# Residual standard error: 37590 on 15 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3628,	Adjusted R-squared:  -0.1469 
# F-statistic: 0.7118 on 12 and 15 DF,  p-value: 0.7203

forest_model(lmAgePCRRBC3i)
confint(lmAgePCRRBC3i)

# Generate predictions
predictions_RBC <- df2 %>%
  mutate(predicted_RBC_ratio = predict(lmAgePCRRBC3i, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake, shape = sex)) +
  geom_point() +
  geom_line(data = predictions_RBC, aes(x = age, y = predicted_RBC_ratio, color = lake, linetype = sex)) +
  labs(title = "CT Ratio RBC vs Age by Lake and Sex",
       x = "Age",
       y = "CT RBC Ratio",
       color = "Lake",
       shape = "Sex",
       linetype = "Sex") +
  theme_minimal()


#### rbc vs age and lake ####
lmAgePCRRBC2 = lm(CT_RBC_ratio~age+lake, data = df2) #Create the linear regression
summary(lmAgePCRRBC2)
# Call:
#   lm(formula = CT_RBC_ratio ~ age + lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -62910 -20862     85  17220  59154 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    62729      18699   3.355  0.00274 **
#   age            -2425       1344  -1.804  0.08437 . 
# lakelotr       22283      17549   1.270  0.21685   
# lakeopeongo    40245      16820   2.393  0.02528 * 
#   lakeshirley    30321      16859   1.798  0.08524 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 32360 on 23 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2757,	Adjusted R-squared:  0.1498 
# F-statistic: 2.189 on 4 and 23 DF,  p-value: 0.1021
check_model(lmAgePCRRBC2)
forest_model(lmAgePCRRBC2)
confint(lmAgePCRRBC2)
# 2.5 %      97.5 %
#   (Intercept)  24046.907 101410.8455
# age          -5206.570    356.0785
# lakelotr    -14018.723  58584.8479
# lakeopeongo   5450.009  75039.5673
# lakeshirley  -4555.508  65197.2004
report(lmAgePCRRBC2)
# We fitted a linear model (estimated using OLS) to predict CT_RBC_ratio
# with age and lake (formula: CT_RBC_ratio ~ age + lake). The model
# explains a statistically not significant and substantial proportion of
# variance (R2 = 0.28, F(4, 23) = 2.19, p = 0.102, adj. R2 = 0.15). The
# model's intercept, corresponding to age = 0 and lake = hogan, is at
# 62728.88 (95% CI [24046.91, 1.01e+05], t(23) = 3.35, p = 0.003). Within
# this model:
# 
#   - The effect of age is statistically non-significant and negative (beta
# = -2425.25, 95% CI [-5206.57, 356.08], t(23) = -1.80, p = 0.084; Std.
# beta = -0.33, 95% CI [-0.70, 0.05])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 22283.06, 95% CI [-14018.72, 58584.85], t(23) = 1.27, p
# = 0.217; Std. beta = 0.63, 95% CI [-0.40, 1.67])
#   - The effect of lake [opeongo] is statistically significant and positive
# (beta = 40244.79, 95% CI [5450.01, 75039.57], t(23) = 2.39, p = 0.025;
# Std. beta = 1.15, 95% CI [0.16, 2.14])
#   - The effect of lake [shirley] is statistically non-significant and
# positive (beta = 30320.85, 95% CI [-4555.51, 65197.20], t(23) = 1.80, p
# = 0.085; Std. beta = 0.86, 95% CI [-0.13, 1.86])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals (CIs) and
# p-values were computed using a Wald t-distribution approximation.

# Fit linear models for each lake and extract formulas
formula_lmAgePCRRBC2 <- df2 %>%
  group_by(lake) %>%
  do(model = lm(CT_RBC_ratio ~ age, data = .)) %>%
  mutate(formula = paste0("CT_RBC_ratio = ", round(coef(model)[1], 2), " + ", round(coef(model)[2], 2), " * age"))

formula_lmAgePCRRBC2
# A tibble: 4 × 3
# Rowwise: 
# lake    model  formula                                  
# <fct>   <list> <chr>                                    
#   1 hogan   <lm>   CT_RBC_ratio = 56746.63 + -1881.4 * age  
# 2 lotr    <lm>   CT_RBC_ratio = 117622.36 + -5741.56 * age
# 3 opeongo <lm>   CT_RBC_ratio = 91244.54 + -1459.32 * age 
# 4 shirley <lm>   CT_RBC_ratio = 84648.16 + -1749.26 * age 
# Create the base plot
lmalRBC <- ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  facet_wrap(~ lake) +  # Create separate panels for each lake
  theme_classic() +
  labs(title = "CT RBC Ratio vs Age",
       x = "Age (y)",
       y = "CT RBC Ratio",
       color = "Lake")

lmalRBC

# Add formulas to each facet
lmalRBC + geom_text(data = formula_lmAgePCRRBC2, aes(x = Inf, y = Inf, label = formula), 
                    hjust = 1.1, vjust = 2, size = 3, color = "black")



# Tidy the model output
tidy_model_lmAgePCRRBC2 <- tidy(lmAgePCRRBC2)

# Create a table
tidy_model_lmAgePCRRBC2 %>%
  kable("html", caption = "Linear Regression Results") %>%
  kable_styling(full_width = FALSE)

# Linear Regression Results
# term	estimate	std.error	statistic	p.value
# (Intercept)	62728.876	18699.068	3.354653	0.0027437
# age	-2425.246	1344.507	-1.803818	0.0843748
# lakelotr	22283.063	17548.475	1.269801	0.2168531
# lakeopeongo	40244.788	16819.980	2.392678	0.0252834
# lakeshirley	30320.846	16859.414	1.798452	0.0852443

# plot with different slopes
# Create the base plot
plot_lmAgePCRRBC2 <- ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line
  theme_classic() +
  labs(title = "CT RBC Ratio vs Age by Lake",
       x = "Age (years)",
       y = "CT RBC Ratio",
       color = "Lake")

plot_lmAgePCRRBC2

# Create the base plot with separate regression lines for each lake
plot_lmAgePCRRBC2 <- ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake, group = lake)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = FALSE) +  # Add a linear regression line for each lake
  theme_classic() +
  labs(title = "CT RBC Ratio vs Age by Lake",
       x = "Age (years)",
       y = "CT RBC Ratio",
       color = "Lake")

plot_lmAgePCRRBC2

# graph with same slope
df2$predicted_CT_RBC_ratio <- predict(lmAgePCRRBC2, newdata = df2)
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) + 
  geom_point(size = 3) + 
  geom_line(aes(y = predicted_CT_RBC_ratio), linetype = "dashed", size = 1) +  # Add predicted values
  theme_classic() +
  labs(title = "CT RBC Ratio vs Age by Lake (with Predictions)",
       x = "Age (years)",
       y = "CT RBC Ratio",
       color = "Lake")


predictionslmAgePCRRBC2 <- df2 %>% mutate(predicted_lmAgePCRRBC2 = predict(lmAgePCRRBC2, newdata = .))

ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() + geom_line(data = predictionslmAgePCRRBC2, aes(x = ventricle_mass, y = predicted_lmAgePCRRBC2, color = lake)) + labs(title = "CT RBC Ratio vs Age by Lake", x = "Age (y)", y = "CT RBC", color = "Lake") +
  theme_minimal()

df2 <- df2 %>% mutate(predicted_lmAgePCRRBC2 = predict(lmAgePCRRBC2, newdata = df2))

ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(aes(y = predicted_lmAgePCRRBC2), linetype = "dashed", size = 1) +
  labs(title = "CT RBC Ratio vs Age by Lake", x = "Age (years)", y = "CT RBC Ratio", color = "Lake") +
  theme_minimal()

ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(aes(y = predicted_lmAgePCRRBC2), size = 1) +
  labs(title = "CT RBC Ratio vs Age by Lake", x = "Age (years)", y = "CT RBC Ratio", color = "Lake") +
  theme_minimal()


#### rbc vs age and lake interaction ###

lmAgePCRRBC5 = lm(CT_RBC_ratio~age*lake, data = df2) #Create the linear regression
summary(lmAgePCRRBC5)
# Call:
#   lm(formula = CT_RBC_ratio ~ age * lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -62620 -17802   1771  19107  58326 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)      56746.6    35814.2   1.584    0.129
# age              -1881.4     3071.0  -0.613    0.547
# lakelotr         60875.7    49956.3   1.219    0.237
# lakeopeongo      34497.9    48976.5   0.704    0.489
# lakeshirley      27901.5    49275.4   0.566    0.578
# age:lakelotr     -3860.2     4474.9  -0.863    0.399
# age:lakeopeongo    422.1     3988.0   0.106    0.917
# age:lakeshirley    132.1     3974.9   0.033    0.974
# 
# Residual standard error: 33640 on 20 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3195,	Adjusted R-squared:  0.08128 
# F-statistic: 1.341 on 7 and 20 DF,  p-value: 0.282

report(lmAgePCRRBC5)
# We fitted a linear model (estimated using OLS) to predict
# CT_RBC_ratio with age and lake (formula: CT_RBC_ratio ~ age *
#                                   lake). The model explains a statistically not significant and
# substantial proportion of variance (R2 = 0.32, F(7, 20) = 1.34, p
#                                     = 0.283, adj. R2 = 0.08). The model's intercept, corresponding to
# age = 0 and lake = hogan, is at 56746.63 (95% CI [-17960.54,
# 1.31e+05], t(20) = 1.58, p = 0.129). Within this model:
# 
#   - The effect of age is statistically non-significant and negative
# (beta = -1881.40, 95% CI [-8287.49, 4524.68], t(20) = -0.61, p =
# 0.547; Std. beta = -0.25, 95% CI [-1.12, 0.61])
#   - The effect of lake [lotr] is statistically non-significant and
# positive (beta = 60875.73, 95% CI [-43331.35, 1.65e+05], t(20) =
# 1.22, p = 0.237; Std. beta = 0.48, 95% CI [-0.64, 1.60])
#   - The effect of lake [opeongo] is statistically non-significant
# and positive (beta = 34497.92, 95% CI [-67665.28, 1.37e+05],
# t(20) = 0.70, p = 0.489; Std. beta = 1.12, 95% CI [0.08, 2.16])
#   - The effect of lake [shirley] is statistically non-significant
# and positive (beta = 27901.53, 95% CI [-74885.12, 1.31e+05],
# t(20) = 0.57, p = 0.578; Std. beta = 0.84, 95% CI [-0.21, 1.89])
#   - The effect of age × lake [lotr] is statistically
# non-significant and negative (beta = -3860.15, 95% CI [-13194.69,
# 5474.38], t(20) = -0.86, p = 0.399; Std. beta = -0.52, 95% CI
# [-1.78, 0.74])
#   - The effect of age × lake [opeongo] is statistically
# non-significant and positive (beta = 422.09, 95% CI [-7896.64,
# 8740.82], t(20) = 0.11, p = 0.917; Std. beta = 0.06, 95% CI
# [-1.07, 1.18])
#   - The effect of age × lake [shirley] is statistically
# non-significant and positive (beta = 132.15, 95% CI [-8159.32,
# 8423.62], t(20) = 0.03, p = 0.974; Std. beta = 0.02, 95% CI
# [-1.10, 1.14])
# 
# Standardized parameters were obtained by fitting the model on a
# standardized version of the dataset. 95% Confidence Intervals
# (CIs) and p-values were computed using a Wald t-distribution
# approximation.
forest_model(lmAgePCRRBC5)

# Generate predictions
predictionslmAgePCRRBC5 <- df2 %>%
  mutate(predicted_lmAgePCRRBC5 = predict(lmAgePCRRBC5, newdata = .))

# Create the plot
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(data = lmAgePCRRBC5, aes(x = age, y = predicted_lmAgePCRRBC5, color = lake)) +
  labs(title = "CT RBC Ratio vs Age by Lake",
       x = "Age (y)",
       y = "CT Ratio RBCs",
       color = "Lake") +
  theme_minimal()

# Generate predictions and add them to df2
df2 <- df2 %>%
  mutate(predicted_lmAgePCRRBC5 = predict(lmAgePCRRBC5, newdata = df2))

# Create the plot
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(aes(x = age, y = predicted_lmAgePCRRBC5, color = lake)) +
  labs(title = "Telomere RBC Ratio vs Age by Lake",
       x = "Age (y)",
       y = "Telomere Ratio RBCs",
       color = "Lake") +
  theme_minimal()



# Create the linear regression model
lmAgePCRRBC5 <- lm(CT_RBC_ratio ~ age * lake, data = df2)

# Add predicted values to the dataframe
df2$predicted_lmAgePCRRBC5 <- predict(lmAgePCRRBC5, newdata = df2)

# Create the plot
ggplot(df2, aes(x = age, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(aes(y = predicted_lmAgePCRRBC5)) +
  labs(title = "CT RBC Ratio vs Age by Lake",
       x = "Age (y)",
       y = "CT Ratio RBCs",
       color = "Lake") +
  theme_minimal()


#### rbc vs lake ####
lmAgePCRRBC4 = lm(CT_RBC_ratio~lake, data = df2) #Create the linear regression
summary(lmAgePCRRBC4)
# Call:
#   lm(formula = CT_RBC_ratio ~ lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -61871 -19971  -5255  18395  70637 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    36051      11967   3.012  0.00602 **
#   lakelotr       25113      18281   1.374  0.18222   
# lakeopeongo    37473      17519   2.139  0.04281 * 
#   lakeshirley    26856      17519   1.533  0.13835   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 33850 on 24 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1733,	Adjusted R-squared:  0.06993 
# F-statistic: 1.677 on 3 and 24 DF,  p-value: 0.1986
check_model(lmAgePCRRBC4)
forest_model(lmAgePCRRBC4)
confint(lmAgePCRRBC4)
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


#### liver ? ####

lmAgePCRliver4 = lm(CT_liver_ratio~lake, data = df2) #Create the linear regression
summary(lmAgePCRliver4)
# Call:
#   lm(formula = CT_liver_ratio ~ lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -13523.6  -2679.1     46.5   2045.5  18015.8 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   7117.1     2043.1   3.483  0.00177 **
#   lakelotr     10099.4     3088.9   3.270  0.00303 **
#   lakeopeongo    814.3     3088.9   0.264  0.79416   
# lakeshirley   8725.2     3088.9   2.825  0.00897 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6129 on 26 degrees of freedom
# Multiple R-squared:  0.3889,	Adjusted R-squared:  0.3184 
# F-statistic: 5.515 on 3 and 26 DF,  p-value: 0.004559
check_model(lmAgePCRliver4)
forest_model(lmAgePCRliver4)
confint(lmAgePCRliver4)


lmAgePCRheart4 = lm(CT_heart_ratio~lake, data = df2) #Create the linear regression
summary(lmAgePCRheart4)
# Call:
#   lm(formula = CT_heart_ratio ~ lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -94304 -31181  -4738  15633 220492 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    77850      20937   3.718  0.00097 ***
#   lakelotr       66336      31654   2.096  0.04600 *  
#   lakeopeongo    53933      31654   1.704  0.10034    
# lakeshirley    33110      31654   1.046  0.30519    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 62810 on 26 degrees of freedom
# Multiple R-squared:  0.166,	Adjusted R-squared:  0.06975 
# F-statistic: 1.725 on 3 and 26 DF,  p-value: 0.1865
check_model(lmAgePCRheart4)
forest_model(lmAgePCRheart4)
confint(lmAgePCRheart4)
# 2.5 %    97.5 %
#   (Intercept)  34813.227 120886.59
# lakelotr      1270.746 131401.43
# lakeopeongo -11132.545 118998.14
# lakeshirley -31955.641  98175.05

#### graphs on TS ratio average ####

# RBC

ggplot(data = df2, mapping = aes(x = lake, y = CT_RBC_ratio)) + 
  geom_boxplot(aes(fill=lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "T/S Ratio for RBCs on Lake Trout")

# RBC
ggplot(data = df2, mapping = aes(x = lake, y = CT_RBC_ratio)) + 
  geom_boxplot(aes(fill = lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "Telomere Ratio for RBCs") +
  annotate("text", x = 3, y = 130000, label = "*", color = "red", size = 15)# +  # Adding text annotation
  #annotate("rect", xmin = 1.5, xmax = 2.5, ymin = 1.2, ymax = 1.8, alpha = 0.2, fill = "blue")  # Adding a rectangle annotation
ggplot(data = df2, mapping = aes(x = lake, y = CT_RBC_ratio)) + 
  geom_boxplot(aes(fill = lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "Telomere Ratio for RBCs") +
  annotate("text", x = 3, y = 130000, label = "*", color = "red", size = 15) +  # Adding text annotation
  scale_y_continuous(limits = c(0, 150000))  # Adjust y-axis limits


# with lakes as groups
ggRBC <- ggplot(df2, aes(x=age, y=CT_RBC_ratio)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio RBCs", subtitle="RBCs", y="T/S Value Ratio", x="Age") +
  stat_cor(method = "pearson", label.x = 15, label.y =140000)
plot(ggRBC)

# liver

ggplot(data = df2, mapping = aes(x = lake, y = CT_liver_ratio)) + 
  geom_boxplot(aes(fill=lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "T/S Ratio for Liver on Lake Trout") 

ggplot(data = df2, mapping = aes(x = lake, y = CT_liver_ratio)) + 
  geom_boxplot(aes(fill=lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "Telomere Ratio for Liver") +
  annotate("text", x = 2, y = 30000, label = "*", color = "red", size = 15) + # +  # Adding text annotation
  annotate("text", x = 4, y = 30000, label = "*", color = "red", size = 15)# +  # Adding text annotation


# with lakes as groups
ggliver <- ggplot(df2, aes(x=age, y=CT_liver_ratio)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio Liver", subtitle="Liver Tissue", y="T/S Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =30000)
plot(ggliver)

# with lakes as groups
ggliver <- ggplot(df2, aes(x = age, y = CT_liver_ratio)) + 
  geom_point(aes(col = lake), size = 2) +  # Set color to vary based on lake
  geom_smooth(method = "lm", size = 1) + 
  labs(title = "Age Vs T/S Ratio Liver", subtitle = "Liver Tissue", y = "T/S Value Ratio", x = "Age") +
  stat_cor(method = "pearson", label.x = 15, label.y = 30000) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = "black"),
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
plot(ggliver)

# with lakes as groups
ggliver <- ggplot(df2, aes(x = age, y = CT_liver_ratio)) + 
  geom_point(aes(col = lake), size = 2) +  # Set color to vary based on lake
  labs(title = "Age Vs T/S Ratio Liver", subtitle = "Liver Tissue", y = "T/S Value Ratio", x = "Age") +
  stat_cor(method = "pearson", label.x = 15, label.y = 30000) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white", color = "black"),
    legend.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
plot(ggliver)



# heart

ggplot(data = df2, mapping = aes(x = lake, y = CT_heart_ratio)) + 
  geom_boxplot(aes(fill=lake), show.legend = FALSE) +
  theme_classic() +
  labs(x = "Lake", y = "Telomere Ratio for Heart")

# with lakes as groups
ggheart <- ggplot(df2, aes(x=age, y=CT_heart_ratio)) + 
  geom_point(aes(col=lake), size=2) +  # Set color to vary based on lake.
  geom_smooth(method="lm", size=1) + 
  labs(title="Age Vs T/S Ratio Heart", subtitle="Heart Tissue", y="T/S Value Ratio", x="Age")+
  stat_cor(method = "pearson", label.x = 15, label.y =300000)
plot(ggheart)


#Liver DNA by lake FSH
lmAgePCR = lm(T.S_Ratio_FSH_liver~age, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCR) #Review the results
# 
# Call:
#   lm(formula = T.S_Ratio_FSH_liver ~ age, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6630.5 -2648.1  -341.1  2895.5  8185.0 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   8692.2     1927.6   4.509 0.000106 ***
#   age           -106.1      159.6  -0.665 0.511405    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4025 on 28 degrees of freedom
# Multiple R-squared:  0.01555,	Adjusted R-squared:  -0.0196 
# F-statistic: 0.4424 on 1 and 28 DF,  p-value: 0.5114
check_model(lmAgePCR)
forest_model(lmAgePCR)
confint(lmAgePCR)
anova(lmAgePCR)
lmAgePCRlake = lm(T.S_Ratio_FSH_liver~age+lake, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlake) #Review the results
# Call:
#   lm(formula = T.S_Ratio_FSH_liver ~ age + lake, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6601.0 -2075.0   424.4  1763.7  4913.2 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5212.13    1865.40   2.794  0.00984 ** 
#   age           -49.23     137.35  -0.358  0.72302    
# lakelotr     6373.22    1699.77   3.749  0.00094 ***
#   lakeopeongo  1839.77    1694.23   1.086  0.28788    
# lakeshirley  3977.98    1698.66   2.342  0.02746 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3344 on 25 degrees of freedom
# Multiple R-squared:  0.3931,	Adjusted R-squared:  0.296 
# F-statistic: 4.049 on 4 and 25 DF,  p-value: 0.0115
check_model(lmAgePCRlake)
forest_model(lmAgePCRlakeSex)
confint(lmAgePCRlakeSex)
# 2.5 %    97.5 %
#   (Intercept)  1555.7684 9397.9663
# age          -350.0412  204.3111
# lakelotr      491.2635 8430.1188
# lakeopeongo -3930.3482 3919.8537
# lakeshirley -2389.7921 5922.3441
# sexmale      -511.2687 5739.5043
# sexunknown  -9842.6191 4480.8221
anova(lmAgePCRlakeSex)
# Analysis of Variance Table
# 
# Response: T.S_Ratio_FSH_liver
# Df    Sum Sq  Mean Sq F value   Pr(>F)   
# age        1   7167284  7167284  0.6853 0.416290   
# lake       3 173983422 57994474  5.5448 0.005157 **
#   sex        2  39071244 19535622  1.8678 0.177152   
# Residuals 23 240563949 10459302                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
emmeans (lmAgePCRlake,  ~ lake, adjust = "tukey")
# lake    emmean   SE df lower.CL upper.CL
# hogan     4662 1115 25     1669     7655
# lotr     11036 1290 25     7574    14498
# opeongo   6502 1271 25     3091     9913
# shirley   8640 1276 25     5217    12064
# 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates

lmAgePCRlakeSex = lm(T.S_Ratio_FSH_liver~age+lake+sex, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlakeSex) #Review the results

#Call:
#   lm(formula = T.S_Ratio_FSH_liver ~ age + lake + sex, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -6819.3 -1915.1   377.1  2610.9  4088.7 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  5476.867   1895.480   2.889  0.00827 **
#   age           -72.865    133.988  -0.544  0.59180   
# lakelotr     4460.691   1918.842   2.325  0.02927 * 
#   lakeopeongo    -5.247   1897.415  -0.003  0.99782   
# lakeshirley  1766.276   2009.065   0.879  0.38841   
# sexmale      2614.118   1510.828   1.730  0.09698 . 
# sexunknown  -2680.899   3462.014  -0.774  0.44660   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3234 on 23 degrees of freedom
# Multiple R-squared:  0.4779,	Adjusted R-squared:  0.3417 
# F-statistic: 3.509 on 6 and 23 DF,  p-value: 0.013
check_model(lmAgePCRlakeSex)
emmeans (lmAgePCRlakeSex,  ~ sex, adjust = "tukey")
# Note: adjust = "tukey" was changed to "sidak"
# because "tukey" is only appropriate for one set of pairwise comparisons
# sex     emmean   SE df lower.CL upper.CL
# female    6219 1091 23     3410     9027
# male      8833  850 23     6644    11022
# unknown   3538 3492 23    -5453    12528
# 
# Results are averaged over the levels of: lake 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 3 estimates


#Liver DNA by lake OX
lmAgePCR2 = lm(T.S_Ratio_OX_liver~age, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCR2) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_liver ~ age, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -13501  -6179  -2103   3578  41776 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  18889.2     5605.3   3.370  0.00221 **
#   age           -268.4      464.1  -0.578  0.56763   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 11700 on 28 degrees of freedom
# Multiple R-squared:  0.01181,	Adjusted R-squared:  -0.02349 
# F-statistic: 0.3345 on 1 and 28 DF,  p-value: 0.5676
check_model(lmAgePCR2)
lmAgePCRlake2 = lm(T.S_Ratio_OX_liver~age+lake, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlake2) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_liver ~ age + lake, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -19385  -5185     98   2596  34259 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) 11549.13    5530.22   2.088   0.0471 *
#   age          -182.84     407.20  -0.449   0.6573  
# lakelotr    13453.42    5039.19   2.670   0.0131 *
#   lakeopeongo    79.83    5022.75   0.016   0.9874  
# lakeshirley 13829.68    5035.91   2.746   0.0110 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 9915 on 25 degrees of freedom
# Multiple R-squared:  0.3669,	Adjusted R-squared:  0.2656 
# F-statistic: 3.621 on 4 and 25 DF,  p-value: 0.01844
check_model(lmAgePCRlake2)
emmeans (lmAgePCRlake2,  ~ lake, adjust = "tukey")
# Note: adjust = "tukey" was changed to "sidak"
# because "tukey" is only appropriate for one set of pairwise comparisons
# lake    emmean   SE df lower.CL upper.CL
# hogan     9507 3307 25      634    18381
# lotr     22961 3825 25    12697    33224
# opeongo   9587 3769 25     -525    19699
# shirley  23337 3783 25    13187    33487
# 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 

lmAgePCRlakeSex2 = lm(T.S_Ratio_OX_liver~age+lake+sex, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlakeSex2) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_liver ~ age + lake + sex, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -14212  -4834   -409   3856  33414 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  12083.5     5824.0   2.075   0.0494 *
#   age           -233.5      411.7  -0.567   0.5761  
# lakelotr      9271.4     5895.8   1.573   0.1295  
# lakeopeongo  -3957.4     5829.9  -0.679   0.5040  
# lakeshirley   8983.3     6173.0   1.455   0.1591  
# sexmale       5765.0     4642.1   1.242   0.2268  
# sexunknown   -5611.3    10637.2  -0.528   0.6029  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 9937 on 23 degrees of freedom
# Multiple R-squared:  0.4149,	Adjusted R-squared:  0.2623 
# F-statistic: 2.719 on 6 and 23 DF,  p-value: 0.03816
check_model(lmAgePCRlakeSex2)
emmeans (lmAgePCRlakeSex2,  ~ lake, adjust = "tukey")
# Note: adjust = "tukey" was changed to "sidak"
# because "tukey" is only appropriate for one set of pairwise comparisons
# lake    emmean   SE df lower.CL upper.CL
# hogan     9527 4235 23    -1913    20968
# lotr     18799 5529 23     3863    33734
# opeongo   5570 5379 23    -8961    20101
# shirley  18511 5608 23     3362    33660
# 
# Results are averaged over the levels of: sex 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates
emmeans (lmAgePCRlakeSex2,  ~ sex, adjust = "tukey")
# sex     emmean    SE df lower.CL upper.CL
# female   13051  3352 23     4421    21680
# male     18816  2613 23    12089    25542
# unknown   7439 10730 23   -20185    35064
# 
# Results are averaged over the levels of: lake 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 3 estimates
forest_model(lmAgePCRlakeSex2)
confint(lmAgePCRlakeSex2)
# 2.5 %     97.5 %
#   (Intercept)     35.73338 24131.3522
# age          -1085.13118   618.1491
# lakelotr     -2924.86041 21467.7441
# lakeopeongo -16017.54943  8102.6622
# lakeshirley  -3786.46577 21753.0665
# sexmale      -3837.92597 15367.9453
# sexunknown  -27616.10161 16393.5224
anova(lmAgePCRlakeSex2)
# Analysis of Variance Table
# 
# Response: T.S_Ratio_OX_liver
# Df     Sum Sq   Mean Sq F value  Pr(>F)  
# age        1   45829138  45829138  0.4641 0.50250  
# lake       3 1378243933 459414644  4.6527 0.01103 *
#   sex        2  186643980  93321990  0.9451 0.40322  
# Residuals 23 2271073130  98742310                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#RBC DNA by lake OX
lmAgePCR3 = lm(T.S_Ratio_OX_RBCs~age, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCR3) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -110456  -49839  -14097   46519  152784 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   161928      34333   4.716 7.11e-05 ***
#   age            -4117       2790  -1.476    0.152    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 68710 on 26 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.07728,	Adjusted R-squared:  0.04179 
# F-statistic: 2.178 on 1 and 26 DF,  p-value: 0.1521
check_model(lmAgePCR3)
lmAgePCRlake3 = lm(T.S_Ratio_OX_RBCs~age+lake, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlake3) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age + lake, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -125819  -41724     169   34440  118307 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   125456      37398   3.355  0.00274 **
#   age            -4850       2689  -1.804  0.08437 . 
# lakelotr       44566      35096   1.270  0.21685   
# lakeopeongo    80488      33640   2.393  0.02528 * 
#   lakeshirley    60641      33718   1.798  0.08524 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 64730 on 23 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2757,	Adjusted R-squared:  0.1498 
# F-statistic: 2.189 on 4 and 23 DF,  p-value: 0.1021
check_model(lmAgePCRlake3)
emmeans (lmAgePCRlake3,  ~ lake, adjust = "tukey")
# lake    emmean    SE df lower.CL upper.CL
# hogan    70196 22908 23     8313   132078
# lotr    114762 26755 23    42489   187034
# opeongo 150684 24547 23    84375   216993
# shirley 130837 24622 23    64326   197349
# 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates
lmAgePCRlakeSex3 = lm(T.S_Ratio_OX_RBCs~age+lake+sex, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlakeSex3) #Review the results
# 
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age + lake + sex, data = PCR_lvr_RBC)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -124564  -42502    1215   35699  120840 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   126859      41626   3.048  0.00612 **
#   age            -4864       2835  -1.715  0.10099   
# lakelotr       50647      44506   1.138  0.26794   
# lakeopeongo    85550      41086   2.082  0.04974 * 
#   lakeshirley    66967      43640   1.535  0.13982   
# sexmale        -8825      33921  -0.260  0.79728   
# sexunknown     -1231      73122  -0.017  0.98673   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 67630 on 21 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2781,	Adjusted R-squared:  0.07179 
# F-statistic: 1.348 on 6 and 21 DF,  p-value: 0.2806
check_model(lmAgePCRlakeSex3)
emmeans (lmAgePCRlakeSex3,  ~ lake, adjust = "tukey")
# lake    emmean    SE df lower.CL upper.CL
# hogan    68096 29574 21   -12434   148626
# lotr    118743 40409 21     8708   228779
# opeongo 153646 36916 21    53122   254169
# shirley 135063 38665 21    29776   240350
# 
# Results are averaged over the levels of: sex 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
emmeans (lmAgePCRlakeSex3,  ~ sex, adjust = "tukey")
forest_model(lmAgePCRlakeSex3)
# sex     emmean    SE df lower.CL upper.CL
# female  122239 25241 21    56774   187703
# male    113414 17939 21    66888   159940
# unknown 121008 73732 21   -70222   312238
# 
# Results are averaged over the levels of: lake 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 3 estimates
confint(lmAgePCRlakeSex3)
# 2.5 %     97.5 %
#   (Intercept)   40293.7824 213424.598
# age          -10760.1456   1032.708
# lakelotr     -41908.1551 143202.575
# lakeopeongo     107.1363 170992.722
# lakeshirley  -23786.5803 157720.947
# sexmale      -79368.0722  61718.889
# sexunknown  -153295.7299 150834.604
anova(lmAgePCRlakeSex3)
# Analysis of Variance Table
# 
# Response: T.S_Ratio_OX_RBCs
# Df     Sum Sq    Mean Sq F value Pr(>F)
# age        1 1.0281e+10 1.0281e+10  2.2479 0.1487
# lake       3 2.6402e+10 8.8007e+09  1.9242 0.1566
# sex        2 3.0954e+08 1.5477e+08  0.0338 0.9668
# Residuals 21 9.6047e+10 4.5737e+09    

#RBC DNA by lake FSH
lmAgePCR4 = lm(T.S_Ratio_OX_RBCs~age, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCR4) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -110456  -49839  -14097   46519  152784 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   161928      34333   4.716 7.11e-05 ***
#   age            -4117       2790  -1.476    0.152    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 68710 on 26 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.07728,	Adjusted R-squared:  0.04179 
# F-statistic: 2.178 on 1 and 26 DF,  p-value: 0.1521
check_model(lmAgePCR4)
lmAgePCRlake4 = lm(T.S_Ratio_OX_RBCs~age+lake, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlake4) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age + lake, data = PCR_lvr_RBC)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -125819  -41724     169   34440  118307 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   125456      37398   3.355  0.00274 **
#   age            -4850       2689  -1.804  0.08437 . 
# lakelotr       44566      35096   1.270  0.21685   
# lakeopeongo    80488      33640   2.393  0.02528 * 
#   lakeshirley    60641      33718   1.798  0.08524 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 64730 on 23 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2757,	Adjusted R-squared:  0.1498 
# F-statistic: 2.189 on 4 and 23 DF,  p-value: 0.1021
check_model(lmAgePCRlake4)
emmeans (lmAgePCRlake4,  ~ lake, adjust = "tukey")
# lake    emmean    SE df lower.CL upper.CL
# hogan    70196 22908 23     8313   132078
# lotr    114762 26755 23    42489   187034
# opeongo 150684 24547 23    84375   216993
# shirley 130837 24622 23    64326   197349
# 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
lmAgePCRlakeSex4 = lm(T.S_Ratio_OX_RBCs~age+lake+sex, data = PCR_lvr_RBC) #Create the linear regression
summary(lmAgePCRlakeSex4) #Review the results
# Call:
#   lm(formula = T.S_Ratio_OX_RBCs ~ age + lake + sex, data = PCR_lvr_RBC)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -124564  -42502    1215   35699  120840 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   126859      41626   3.048  0.00612 **
#   age            -4864       2835  -1.715  0.10099   
# lakelotr       50647      44506   1.138  0.26794   
# lakeopeongo    85550      41086   2.082  0.04974 * 
#   lakeshirley    66967      43640   1.535  0.13982   
# sexmale        -8825      33921  -0.260   0.79728   
# sexunknown     -1231      73122  -0.017  0.98673   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 67630 on 21 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2781,	Adjusted R-squared:  0.07179 
# F-statistic: 1.348 on 6 and 21 DF,  p-value: 0.2806
check_model(lmAgePCRlakeSex4)
emmeans (lmAgePCRlakeSex4,  ~ lake, adjust = "tukey")
# lake    emmean    SE df lower.CL upper.CL
# hogan    68096 29574 21   -12434   148626
# lotr    118743 40409 21     8708   228779
# opeongo 153646 36916 21    53122   254169
# shirley 135063 38665 21    29776   240350
# 
# Results are averaged over the levels of: sex 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
emmeans (lmAgePCRlakeSex2,  ~ sex, adjust = "tukey")
# sex     emmean    SE df lower.CL upper.CL
# female   13051  3352 23     4421    21680
# male     18816  2613 23    12089    25542
# unknown   7439 10730 23   -20185    35064
# 
# Results are averaged over the levels of: lake 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 3 estimates 
forest_model(lmAgePCRlakeSex2)
confint(lmAgePCRlakeSex2)
# 2.5 %     97.5 %
#   (Intercept)     35.73338 24131.3522
# age          -1085.13118   618.1491
# lakelotr     -2924.86041 21467.7441
# lakeopeongo -16017.54943  8102.6622
# lakeshirley  -3786.46577 21753.0665
# sexmale      -3837.92597 15367.9453
# sexunknown  -27616.10161 16393.5224
anova(lmAgePCRlakeSex2)
# Analysis of Variance Table
# 
# Response: T.S_Ratio_OX_liver
# Df     Sum Sq   Mean Sq F value  Pr(>F)  
# age        1   45829138  45829138  0.4641 0.50250  
# lake       3 1378243933 459414644  4.6527 0.01103 *
#   sex        2  186643980  93321990  0.9451 0.40322  
# Residuals 23 2271073130  98742310                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Bar Plots ####


# Standard deviation of the mean as error bar
# Barplot
# Barplot
ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_FSH_liver)) + 
  geom_bar(stat = "identity")

ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_FSH_liver)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=T.S_Ratio_FSH_live-sd, ymax=T.S_Ratio_FSH_live+sd), width=.2,
                position=position_dodge(.9))
# Change box plot colors by groups
bp<-ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_FSH_liver, fill=lake)) +
  geom_boxplot()
bp

bp1 <- ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_FSH_liver)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=T.S_Ratio_FSH_liver-sd, ymax=T.S_Ratio_FSH_liver+sd), width=.2,
                position=position_dodge(.9))

bp1 + scale_fill_brewer(palette="Paired") + theme_minimal()
bp1

bp1<-ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_OX_liver, fill=lake)) +
  geom_boxplot()
bp1

bp2<-ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_OX_RBCs, fill=lake)) +
  geom_boxplot()
bp2

bp3<-ggplot(PCR_lvr_RBC, aes(x=lake, y=T.S_Ratio_FSH_RBCs, fill=lake)) +
  geom_boxplot()
bp3

#### QC ####
QC_sheet <- read.csv("C:/Users/user/Desktop/CH3/QC_sheet.csv")
View(QC_sheet)
names(QC_sheet)

QC_sheet$OD_heart <- as.numeric(QC_sheet$OD_heart)
QC_sheet$heart_TS <- as.numeric(QC_sheet$heart_TS)


# Perform correlation test
cor.test(QC_sheet$OD_heart, QC_sheet$Heart_TS, method = "pearson")

QC_sheet$OD_heart <- as.numeric(as.character(QC_sheet$OD_heart))
QC_sheet$heart_TS <- as.numeric(as.character(QC_sheet$heart_TS))


# Inspect unique values in heart_TS
unique(QC_sheet$heart_TS)

# Check if there are any non-numeric entries
QC_sheet[!is.na(as.numeric(as.character(QC_sheet$heart_TS))), ]
# Replace non-numeric entries with NA
QC_sheet$heart_TS <- suppressWarnings(as.numeric(as.character(QC_sheet$heart_TS)))

# Check the class and unique values
class(QC_sheet$heart_TS)
unique(QC_sheet$heart_TS)

names(QC_sheet)

# linear model
#### heart QC ####
# Perform correlation test with 0 intercept
model <- lm(Heart_TS ~ 0 + OD_heart, data = QC_sheet)
summary(model)
# Call:
#   lm(formula = Heart_TS ~ 0 + OD_heart, data = QC_sheet)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -65004 -40858 -18137  21741 253430 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# OD_heart   110530      12034   9.185 6.08e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 66390 on 28 degrees of freedom
# (37 observations deleted due to missingness)
# Multiple R-squared:  0.7508,	Adjusted R-squared:  0.7419 
# F-statistic: 84.36 on 1 and 28 DF,  p-value: 6.075e-10

# Perform correlation test with 0 intercept
model9 <- lm(Heart_TS ~ OD_heart, data = QC_sheet)
summary(model9)

# Call:
#   lm(formula = Heart_TS ~ OD_heart, data = QC_sheet)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -61229 -30095  -5928  24455 137384 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    19073      86090   0.222    0.827
# OD_heart       93169      84628   1.101    0.285
# 
# Residual standard error: 46950 on 19 degrees of freedom
# (28 observations deleted due to missingness)
# Multiple R-squared:  0.05997,	Adjusted R-squared:  0.01049 
# F-statistic: 1.212 on 1 and 19 DF,  p-value: 0.2847

# Perform correlation test with 0 intercept
model22 <- lm(liver_TS ~ 0 + OD_liver, data = QC_sheet)
summary(model22)

# Inspect unique values
unique(QC_sheet$liver_TS)
unique(QC_sheet$OD_liver)

# Check for NA, NaN, or Inf
sum(is.na(QC_sheet$liver_TS))    # Count of NA values in liver_TS
sum(is.na(QC_sheet$OD_liver))    # Count of NA values in OD_liver
any(is.nan(QC_sheet$liver_TS))   # Check for NaN values
any(is.nan(QC_sheet$OD_liver))   # Check for NaN values
any(is.infinite(QC_sheet$liver_TS))   # Check for Inf values
any(is.infinite(QC_sheet$OD_liver))   # Check for Inf values


# Remove rows with NA, NaN, or Inf in either column
QC_sheet <- QC_sheet[!is.na(QC_sheet$liver_TS) & 
                       !is.na(QC_sheet$OD_liver) & 
                       !is.nan(QC_sheet$liver_TS) & 
                       !is.nan(QC_sheet$OD_liver) & 
                       !is.infinite(QC_sheet$liver_TS) & 
                       !is.infinite(QC_sheet$OD_liver), ]


# Check data types
class(QC_sheet$liver_TS)
class(QC_sheet$OD_liver)

# Convert to numeric if necessary
QC_sheet$liver_TS <- as.numeric(QC_sheet$liver_TS)
QC_sheet$OD_liver <- as.numeric(QC_sheet$OD_liver)

#### liver QC ####
model22 <- lm(liver_TS ~ 0 + OD_liver, data = QC_sheet)
summary(model22)
# Call:
#   lm(formula = liver_TS ~ 0 + OD_liver, data = QC_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -9531.4 -4586.8  -490.9  3356.1 20503.8 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# OD_liver    13361       1646   8.117 9.31e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6834 on 20 degrees of freedom
# (28 observations deleted due to missingness)
# Multiple R-squared:  0.7671,	Adjusted R-squared:  0.7555 
# F-statistic: 65.89 on 1 and 20 DF,  p-value: 9.308e-08

#### liver QC ####
model222 <- lm(liver_TS ~ OD_liver, data = QC_sheet)
# Call:
#   lm(formula = liver_TS ~ OD_liver, data = QC_sheet)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -9492.2 -4621.4  -410.9  3058.8 20567.1 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)    712.6     9476.5   0.075    0.941
# OD_liver     12584.7    10460.3   1.203    0.244
# 
# Residual standard error: 7010 on 19 degrees of freedom
# (28 observations deleted due to missingness)
# Multiple R-squared:  0.07079,	Adjusted R-squared:  0.02188 
# F-statistic: 1.447 on 1 and 19 DF,  p-value: 0.2437
summary(model222)



# Basic scatterplot
plot(QC_sheet$OD_liver, QC_sheet$liver_TS, 
     main = "Regression of liver_TS on OD_liver", 
     xlab = "OD_liver", ylab = "liver_TS", 
     pch = 4, col = "blue")
# Add regression line with intercept set to zero
abline(model22, col = "red", lwd = 2)
# Add enhancements to the plot
legend("topleft", legend = c("Data Points", "Regression Line"), 
       col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1), lwd = 2)
grid(col = "gray", lty = "dotted")

# Scatterplot with axes limits that include (0,0)
plot(QC_sheet$OD_liver, QC_sheet$liver_TS, 
     main = "Regression of liver_TS on OD_liver", 
     xlab = "OD_liver", ylab = "liver_TS", 
     pch = 16, col = "blue", 
     xlim = c(0, max(QC_sheet$OD_liver, na.rm = TRUE)),
     ylim = c(0, max(QC_sheet$liver_TS, na.rm = TRUE)))

# Add regression line with intercept set to zero
abline(model22, col = "red", lwd = 2)

plot(QC_sheet$OD_heart, QC_sheet$Heart_TS, 
  main = "Regression of heart_TS on OD_Heart", 
  xlab = "OD_heart", ylab = "Heart_TS", 
  pch = 4, col = "blue", 
  xlim = c(0, max(QC_sheet$OD_heart, na.rm = TRUE)),
  ylim = c(0, max(QC_sheet$Heart_TS, na.rm = TRUE)))

# Add regression line with intercept set to zero
abline(model, col = "red", lwd = 2)

# Calculate the mean and difference
mean_values <- (QC_sheet$Heart_TS + QC_sheet$OD_heart) / 2
differences <- QC_sheet$Heart_TS - QC_sheet$OD_heart

# Create the Bland-Altman plot
plot(mean_values, differences,
     main = "Bland-Altman Plot",
     xlab = "Mean of Heart_TS and OD_heart",
     ylab = "Difference (Heart_TS - OD_heart)",
     pch = 16, col = "blue")

# Add the mean difference line (bias)
abline(h = mean(differences), col = "red", lwd = 2)

# Add limits of agreement (±1.96 SD of differences)
abline(h = mean(differences) + 1.96 * sd(differences), col = "green", lty = 2)
abline(h = mean(differences) - 1.96 * sd(differences), col = "green", lty = 2)
legend("topleft", legend = c("Mean Difference", "Limits of Agreement"),
       col = c("red", "green"), lty = c(1, 2), lwd = 2)
grid(col = "gray", lty = "dotted")

# Check the differences
head(differences)
summary(differences)
differences <- differences[!is.na(differences)]
mean_value <- mean(differences)
mean_value  # Verify the value
abline(h = mean(differences), col = "red", lwd = 2)

plot(mean_values, differences,
     main = "Bland-Altman Plot",
     xlab = "Mean of Heart_TS and OD_heart",
     ylab = "Difference (Heart_TS - OD_heart)",
     pch = 16, col = "blue",
     ylim = c(min(differences) - 1, max(differences) + 1))  # Expand y-axis range
abline(h = mean(differences), col = "red", lwd = 2)

#### relative organ mass linear models ####

# Define the columns for the numerators
numerator_cols <- c("liver_mass", "gonad_mass", "ventricle_mass")

# Create ratios using apply() and store in new columns
df2[paste0(numerator_cols, "_mass_post_resp_g")] <- 
  sapply(numerator_cols, function(col) OD_sheet[[col]] / df2$mass_post_resp)

head(df2)
names(df2)

#### linear model heart df2 vs relative gonad mass & lake interaction ####
lm_PCRh_gm_lakeii <- lm(CT_heart_ratio ~ gonad_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRh_gm_lakeii)
# 
# Call:
#   lm(formula = CT_heart_ratio ~ gonad_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -103691  -22935   -5213   15789  205593 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                              81136.8    28286.4
# gonad_mass_mass_post_resp_g               -116.2      644.6
# lakelotr                                110076.5    51218.4
# lakeopeongo                              68520.9    72255.3
# lakeshirley                              -8021.2    68438.9
# gonad_mass_mass_post_resp_g:lakelotr     -3570.8     2815.2
# gonad_mass_mass_post_resp_g:lakeopeongo  -6186.9    21801.5
# gonad_mass_mass_post_resp_g:lakeshirley   3023.1     4447.5
# t value Pr(>|t|)   
# (Intercept)                               2.868  0.00893 **
#   gonad_mass_mass_post_resp_g              -0.180  0.85856   
# lakelotr                                  2.149  0.04288 * 
#   lakeopeongo                               0.948  0.35327   
# lakeshirley                              -0.117  0.90776   
# gonad_mass_mass_post_resp_g:lakelotr     -1.268  0.21791   
# gonad_mass_mass_post_resp_g:lakeopeongo  -0.284  0.77923   
# gonad_mass_mass_post_resp_g:lakeshirley   0.680  0.50378   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 64890 on 22 degrees of freedom
# Multiple R-squared:  0.2469,	Adjusted R-squared:  0.007232 
# F-statistic:  1.03 on 7 and 22 DF,  p-value: 0.4386

forest_model(lm_PCRh_gm_lakeii)
# Generate predictions
predictionslm_PCRh_gm_lakeii <- df2 %>%
  mutate(predicted_lm_PCRh_gm_lakeii = predict(lm_PCRh_gm_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass_mass_post_resp_g, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRh_gm_lakeii, aes(x = gonad_mass_mass_post_resp_g, y = predicted_lm_PCRh_gm_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Gonad Mass by Lake",
       x = "Relative Gonad Mass (g)",
       y = " Relative Telomere Length Heart",
       color = "Lake") +
  theme_minimal()


#### linear model heart df2 vs relative ventricle mass & lake interaction ####
lm_PCRh_vn_lakeii <- lm(CT_heart_ratio ~ ventricle_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRh_vn_lakeii)
# Call:
#   lm(formula = CT_heart_ratio ~ ventricle_mass_mass_post_resp_g * 
#        lake, data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -100207  -30299    -846   27709  114980 
# 
# Coefficients:
#   Estimate
# (Intercept)                                  72017.98
# ventricle_mass_mass_post_resp_g                 50.01
# lakelotr                                    -28797.29
# lakeopeongo                                  63546.47
# lakeshirley                                  45005.34
# ventricle_mass_mass_post_resp_g:lakelotr      1375.45
# ventricle_mass_mass_post_resp_g:lakeopeongo   -401.28
# ventricle_mass_mass_post_resp_g:lakeshirley    -93.67
# Std. Error t value
# (Intercept)                                   27827.19   2.588
# ventricle_mass_mass_post_resp_g                 174.51   0.287
# lakelotr                                      48082.56  -0.599
# lakeopeongo                                   53750.72   1.182
# lakeshirley                                   42589.38   1.057
# ventricle_mass_mass_post_resp_g:lakelotr        494.59   2.781
# ventricle_mass_mass_post_resp_g:lakeopeongo    3779.01  -0.106
# ventricle_mass_mass_post_resp_g:lakeshirley     245.64  -0.381
# Pr(>|t|)  
# (Intercept)                                   0.0168 *
#   ventricle_mass_mass_post_resp_g               0.7771  
# lakelotr                                      0.5553  
# lakeopeongo                                   0.2497  
# lakeshirley                                   0.3021  
# ventricle_mass_mass_post_resp_g:lakelotr      0.0109 *
#   ventricle_mass_mass_post_resp_g:lakeopeongo   0.9164  
# ventricle_mass_mass_post_resp_g:lakeshirley   0.7066  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 56940 on 22 degrees of freedom
# Multiple R-squared:  0.4201,	Adjusted R-squared:  0.2356 
# F-statistic: 2.277 on 7 and 22 DF,  p-value: 0.0663

forest_model(lm_PCRh_vn_lakeii)
# Generate predictions
predictionslm_PCRh_vn_lakeii <- df2 %>%
  mutate(predicted_lm_PCRh_vn_lakeii = predict(lm_PCRh_vn_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = ventricle_mass_mass_post_resp_g, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRh_vn_lakeii, aes(x = ventricle_mass_mass_post_resp_g, y = predicted_lm_PCRh_vn_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Ventricle Mass by Lake",
       x = "Relative Ventricle Mass",
       y = " Relative Telomere Length Heart",
       color = "Lake") +
  theme_minimal()




#### linear model heart df2 vs relative liver mass & lake interaction ####
lm_PCRh_lv_lakeii <- lm(CT_heart_ratio ~ liver_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRh_lv_lakeii)
# Call:
#   lm(formula = CT_heart_ratio ~ liver_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -89866 -23415  -8963  22455 210970 
# 
# Coefficients:
#   Estimate Std. Error
# (Intercept)                              78876.7    29209.5
# liver_mass_mass_post_resp_g               -101.2     1936.3
# lakelotr                                 14582.8    58441.1
# lakeopeongo                              23427.3    65766.9
# lakeshirley                              66416.9    54483.2
# liver_mass_mass_post_resp_g:lakelotr      1852.0     2467.0
# liver_mass_mass_post_resp_g:lakeopeongo   2238.5     4340.9
# liver_mass_mass_post_resp_g:lakeshirley  -1165.4     2410.5
# t value Pr(>|t|)  
# (Intercept)                               2.700   0.0131 *
#   liver_mass_mass_post_resp_g              -0.052   0.9588  
# lakelotr                                  0.250   0.8053  
# lakeopeongo                               0.356   0.7251  
# lakeshirley                               1.219   0.2357  
# liver_mass_mass_post_resp_g:lakelotr      0.751   0.4608  
# liver_mass_mass_post_resp_g:lakeopeongo   0.516   0.6112  
# liver_mass_mass_post_resp_g:lakeshirley  -0.483   0.6335  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 64840 on 22 degrees of freedom
# Multiple R-squared:  0.2479,	Adjusted R-squared:  0.008568 
# F-statistic: 1.036 on 7 and 22 DF,  p-value: 0.4351

forest_model(lm_PCRh_lv_lakeii)
# Generate predictions
predictionslm_PCRh_lv_lakeii <- df2 %>%
  mutate(predicted_lm_PCRh_lv_lakeii = predict(lm_PCRh_lv_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = liver_mass_mass_post_resp_g, y = CT_heart_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRh_lv_lakeii, aes(x = liver_mass_mass_post_resp_g, y = predicted_lm_PCRh_lv_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Liver Mass by Lake",
       x = "Relative Liver Mass",
       y = " Relative Telomere Length Heart",
       color = "Lake") +
  theme_minimal()

#### linear model liver df2 vs relative liver mass & lake interaction ####
lm_PCRl_lv_lakeii <- lm(CT_liver_ratio ~ liver_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRl_lv_lakeii)

# Call:
#   lm(formula = CT_liver_ratio ~ liver_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8578.1 -2871.1   518.2  1718.6 12913.3 
# 
# Coefficients:
#   Estimate Std. Error t value
# (Intercept)                              8384.03    2547.16   3.292
# liver_mass_mass_post_resp_g              -124.86     168.85  -0.739
# lakelotr                                 5113.91    5096.24   1.003
# lakeopeongo                              1873.11    5735.07   0.327
# lakeshirley                             -1410.13    4751.10  -0.297
# liver_mass_mass_post_resp_g:lakelotr      253.21     215.13   1.177
# liver_mass_mass_post_resp_g:lakeopeongo   -43.76     378.54  -0.116
# liver_mass_mass_post_resp_g:lakeshirley   452.02     210.20   2.150
# Pr(>|t|)   
# (Intercept)                              0.00333 **
#   liver_mass_mass_post_resp_g              0.46744   
# lakelotr                                 0.32655   
# lakeopeongo                              0.74705   
# lakeshirley                              0.76940   
# liver_mass_mass_post_resp_g:lakelotr     0.25177   
# liver_mass_mass_post_resp_g:lakeopeongo  0.90902   
# liver_mass_mass_post_resp_g:lakeshirley  0.04277 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5655 on 22 degrees of freedom
# Multiple R-squared:  0.5599,	Adjusted R-squared:  0.4199 
# F-statistic: 3.999 on 7 and 22 DF,  p-value: 0.005774
forest_model(lm_PCRl_lv_lakeii)
# Generate predictions
predictionslm_PCRl_lv_lakeii <- df2 %>%
  mutate(predicted_lm_PCRl_lv_lakeii = predict(lm_PCRl_lv_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = liver_mass_mass_post_resp_g, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRl_lv_lakeii, aes(x = liver_mass_mass_post_resp_g, y = predicted_lm_PCRl_lv_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Liver Mass by Lake",
       x = "Relative Liver Mass",
       y = " Relative Telomere Length Liver",
       color = "Lake") +
  theme_minimal()

#### linear model liver df2 vs relative gonad mass & lake interaction ####
lm_PCRl_gn_lakeii <- lm(CT_liver_ratio ~ gonad_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRl_gn_lakeii)
# Call:
#   lm(formula = CT_liver_ratio ~ gonad_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11933.3  -2358.0    215.7   1765.1  18617.8 
# 
# Coefficients:
#   Estimate
# (Intercept)                              7980.27
# gonad_mass_mass_post_resp_g               -30.52
# lakelotr                                12737.16
# lakeopeongo                             -2850.03
# lakeshirley                              3930.11
# gonad_mass_mass_post_resp_g:lakelotr     -243.96
# gonad_mass_mass_post_resp_g:lakeopeongo  1018.27
# gonad_mass_mass_post_resp_g:lakeshirley   332.53
# Std. Error t value
# (Intercept)                                2783.04   2.867
# gonad_mass_mass_post_resp_g                  63.42  -0.481
# lakelotr                                   5039.27   2.528
# lakeopeongo                                7109.05  -0.401
# lakeshirley                                6733.56   0.584
# gonad_mass_mass_post_resp_g:lakelotr        276.98  -0.881
# gonad_mass_mass_post_resp_g:lakeopeongo    2145.00   0.475
# gonad_mass_mass_post_resp_g:lakeshirley     437.58   0.760
# Pr(>|t|)   
# (Intercept)                              0.00895 **
#   gonad_mass_mass_post_resp_g              0.63508   
# lakelotr                                 0.01917 * 
#   lakeopeongo                              0.69236   
# lakeshirley                              0.56539   
# gonad_mass_mass_post_resp_g:lakelotr     0.38795   
# gonad_mass_mass_post_resp_g:lakeopeongo  0.63967   
# gonad_mass_mass_post_resp_g:lakeshirley  0.45536   
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6384 on 22 degrees of freedom
# Multiple R-squared:  0.439,	Adjusted R-squared:  0.2605 
# F-statistic:  2.46 on 7 and 22 DF,  p-value: 0.05031
forest_model(lm_PCRl_gn_lakeii)
# Generate predictions
predictionslm_PCRl_gn_lakeii <- df2 %>%
  mutate(predicted_lm_PCRl_gn_lakeii = predict(lm_PCRl_gn_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass_mass_post_resp_g, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRl_gn_lakeii, aes(x = gonad_mass_mass_post_resp_g, y = predicted_lm_PCRl_gn_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Gonad Mass by Lake",
       x = "Relative Gonad Mass",
       y = " Relative Telomere Length Liver",
       color = "Lake") +
  theme_minimal()

#### linear model liver df2 vs relative ventricle mass & lake interaction ####
lm_PCRl_vn_lakeii <- lm(CT_liver_ratio ~ ventricle_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRl_vn_lakeii)
# Call:
#   lm(formula = CT_liver_ratio ~ ventricle_mass_mass_post_resp_g * 
#        lake, data = df2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -11230.3  -3048.9   -624.7   2357.8  14158.6 
# 
# Coefficients:
#   Estimate
# (Intercept)                                  6658.216
# ventricle_mass_mass_post_resp_g                 3.935
# lakelotr                                     8117.939
# lakeopeongo                                  2551.360
# lakeshirley                                 14070.857
# ventricle_mass_mass_post_resp_g:lakelotr       30.518
# ventricle_mass_mass_post_resp_g:lakeopeongo  -122.661
# ventricle_mass_mass_post_resp_g:lakeshirley   -39.123
# Std. Error
# (Intercept)                                   2983.259
# ventricle_mass_mass_post_resp_g                 18.709
# lakelotr                                      5154.770
# lakeopeongo                                   5762.433
# lakeshirley                                   4565.864
# ventricle_mass_mass_post_resp_g:lakelotr        53.024
# ventricle_mass_mass_post_resp_g:lakeopeongo    405.135
# ventricle_mass_mass_post_resp_g:lakeshirley     26.334
# t value
# (Intercept)                                   2.232
# ventricle_mass_mass_post_resp_g               0.210
# lakelotr                                      1.575
# lakeopeongo                                   0.443
# lakeshirley                                   3.082
# ventricle_mass_mass_post_resp_g:lakelotr      0.576
# ventricle_mass_mass_post_resp_g:lakeopeongo  -0.303
# ventricle_mass_mass_post_resp_g:lakeshirley  -1.486
# Pr(>|t|)   
# (Intercept)                                  0.03613 * 
#   ventricle_mass_mass_post_resp_g              0.83535   
# lakelotr                                     0.12957   
# lakeopeongo                                  0.66226   
# lakeshirley                                  0.00545 **
#   ventricle_mass_mass_post_resp_g:lakelotr     0.57076   
# ventricle_mass_mass_post_resp_g:lakeopeongo  0.76491   
# ventricle_mass_mass_post_resp_g:lakeshirley  0.15156   
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6104 on 22 degrees of freedom
# Multiple R-squared:  0.4872,	Adjusted R-squared:  0.324 
# F-statistic: 2.986 on 7 and 22 DF,  p-value: 0.0232
forest_model(lm_PCRl_vn_lakeii)
# Generate predictions
predictionslm_PCRl_vn_lakeii <- df2 %>%
  mutate(predicted_lm_PCRl_vn_lakeii = predict(lm_PCRl_vn_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = ventricle_mass_mass_post_resp_g, y = CT_liver_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRl_vn_lakeii, aes(x = ventricle_mass_mass_post_resp_g, y = predicted_lm_PCRl_vn_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Ventricle Mass by Lake",
       x = "Relative Ventricle Mass",
       y = " Relative Telomere Length Liver",
       color = "Lake") +
  theme_minimal()

#### linear model RBC df2 vs relative ventricle mass & lake interaction ####
lm_PCRr_vn_lakeii <- lm(CT_RBC_ratio ~ ventricle_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRr_vn_lakeii) # Opeongo significant
# some lake effects of lakes on the telomere length for the RBCs

# Call:
#   lm(formula = CT_RBC_ratio ~ ventricle_mass_mass_post_resp_g * 
#        lake, data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -49363 -18066   -757  15295  68932 
# 
# Coefficients:
#   Estimate
# (Intercept)                                 22015.936
# ventricle_mass_mass_post_resp_g               118.969
# lakelotr                                    47946.080
# lakeopeongo                                 71914.921
# lakeshirley                                 25045.067
# ventricle_mass_mass_post_resp_g:lakelotr     -250.848
# ventricle_mass_mass_post_resp_g:lakeopeongo -2014.464
# ventricle_mass_mass_post_resp_g:lakeshirley    -4.864
# Std. Error
# (Intercept)                                  17345.139
# ventricle_mass_mass_post_resp_g                105.156
# lakelotr                                     29317.820
# lakeopeongo                                  32676.972
# lakeshirley                                  26035.188
# ventricle_mass_mass_post_resp_g:lakelotr       304.230
# ventricle_mass_mass_post_resp_g:lakeopeongo   2275.744
# ventricle_mass_mass_post_resp_g:lakeshirley    147.973
# t value
# (Intercept)                                   1.269
# ventricle_mass_mass_post_resp_g               1.131
# lakelotr                                      1.635
# lakeopeongo                                   2.201
# lakeshirley                                   0.962
# ventricle_mass_mass_post_resp_g:lakelotr     -0.825
# ventricle_mass_mass_post_resp_g:lakeopeongo  -0.885
# ventricle_mass_mass_post_resp_g:lakeshirley  -0.033
# Pr(>|t|)  
# (Intercept)                                   0.2189  
# ventricle_mass_mass_post_resp_g               0.2713  
# lakelotr                                      0.1176  
# lakeopeongo                                   0.0397 *
#   lakeshirley                                   0.3476  
# ventricle_mass_mass_post_resp_g:lakelotr      0.4194  
# ventricle_mass_mass_post_resp_g:lakeopeongo   0.3866  
# ventricle_mass_mass_post_resp_g:lakeshirley   0.9741  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 34290 on 20 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.2931,	Adjusted R-squared:  0.04567 
# F-statistic: 1.185 on 7 and 20 DF,  p-value: 0.3552
forest_model(lm_PCRr_vn_lakeii)
# Generate predictions
predictionslm_PCRr_vn_lakeii <- df2 %>%
  mutate(predicted_lm_PCRr_vn_lakeii = predict(lm_PCRr_vn_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = ventricle_mass_mass_post_resp_g, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRr_vn_lakeii, aes(x = ventricle_mass_mass_post_resp_g, y = predicted_lm_PCRr_vn_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Ventricle Mass by Lake",
       x = "Relative Ventricle Mass",
       y = " Relative Telomere Length RBCs",
       color = "Lake") +
  theme_minimal()

#### linear model RBC df2 vs relative liver mass & lake interaction ####
lm_PCRr_lv_lakeii <- lm(CT_RBC_ratio ~ liver_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRr_lv_lakeii) # no significance o effects
# Call:
#   lm(formula = CT_RBC_ratio ~ liver_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -56317 -22071  -1864  20591  57069 
# 
# Coefficients:
#   Estimate
# (Intercept)                              29601.8
# liver_mass_mass_post_resp_g                600.8
# lakelotr                                 -9202.6
# lakeopeongo                              43941.7
# lakeshirley                              42957.9
# liver_mass_mass_post_resp_g:lakelotr       931.3
# liver_mass_mass_post_resp_g:lakeopeongo   -602.2
# liver_mass_mass_post_resp_g:lakeshirley   -956.9
# Std. Error t value
# (Intercept)                                16326.7   1.813
# liver_mass_mass_post_resp_g                 1027.9   0.585
# lakelotr                                   31368.1  -0.293
# lakeopeongo                                34975.2   1.256
# lakeshirley                                29144.3   1.474
# liver_mass_mass_post_resp_g:lakelotr        1340.6   0.695
# liver_mass_mass_post_resp_g:lakeopeongo     2283.7  -0.264
# liver_mass_mass_post_resp_g:lakeshirley     1274.5  -0.751
# Pr(>|t|)  
# (Intercept)                               0.0849 .
# liver_mass_mass_post_resp_g               0.5654  
# lakelotr                                  0.7723  
# lakeopeongo                               0.2235  
# lakeshirley                               0.1561  
# liver_mass_mass_post_resp_g:lakelotr      0.4953  
# liver_mass_mass_post_resp_g:lakeopeongo   0.7947  
# liver_mass_mass_post_resp_g:lakeshirley   0.4615  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 34040 on 20 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3033,	Adjusted R-squared:  0.05952 
# F-statistic: 1.244 on 7 and 20 DF,  p-value: 0.326
forest_model(lm_PCRr_lv_lakeii)
# Generate predictions
predictionslm_PCRr_lv_lakeii <- df2 %>%
  mutate(predicted_lm_PCRr_lv_lakeii = predict(lm_PCRr_lv_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = liver_mass_mass_post_resp_g, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRr_lv_lakeii, aes(x = liver_mass_mass_post_resp_g, y = predicted_lm_PCRr_lv_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Liver Mass by Lake",
       x = "Relative Liver Mass",
       y = " Relative Telomere Length RBCs",
       color = "Lake") +
  theme_minimal()

#### linear model RBC df2 vs relative gonad mass & lake interaction ####
lm_PCRr_gn_lakeii <- lm(CT_RBC_ratio ~ gonad_mass_mass_post_resp_g * lake, data = df2)
summary(lm_PCRr_gn_lakeii) # no lake effects or mass effects
# Call:
#   lm(formula = CT_RBC_ratio ~ gonad_mass_mass_post_resp_g * lake, 
#      data = df2)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -60680 -19134  -1709  10300  68280 
# 
# Coefficients:
#   Estimate
# (Intercept)                              32312.7
# gonad_mass_mass_post_resp_g                120.5
# lakelotr                                 35940.6
# lakeopeongo                             -20695.2
# lakeshirley                              22816.5
# gonad_mass_mass_post_resp_g:lakelotr      -640.1
# gonad_mass_mass_post_resp_g:lakeopeongo  21709.1
# gonad_mass_mass_post_resp_g:lakeshirley    476.9
# Std. Error t value
# (Intercept)                                16067.3   2.011
# gonad_mass_mass_post_resp_g                  345.6   0.349
# lakelotr                                   29213.5   1.230
# lakeopeongo                                38224.0  -0.541
# lakeshirley                                36262.5   0.629
# gonad_mass_mass_post_resp_g:lakelotr        1513.8  -0.423
# gonad_mass_mass_post_resp_g:lakeopeongo    11372.8   1.909
# gonad_mass_mass_post_resp_g:lakeshirley     2321.4   0.205
# Pr(>|t|)  
# (Intercept)                               0.0580 .
# gonad_mass_mass_post_resp_g               0.7310  
# lakelotr                                  0.2329  
# lakeopeongo                               0.5942  
# lakeshirley                               0.5363  
# gonad_mass_mass_post_resp_g:lakelotr      0.6769  
# gonad_mass_mass_post_resp_g:lakeopeongo   0.0707 .
# gonad_mass_mass_post_resp_g:lakeshirley   0.8393  
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 33850 on 20 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3111,	Adjusted R-squared:  0.06998 
# F-statistic:  1.29 on 7 and 20 DF,  p-value: 0.3048
forest_model(lm_PCRr_gn_lakeii)
# Generate predictions
predictionslm_PCRr_gn_lakeii <- df2 %>%
  mutate(predicted_lm_PCRr_gn_lakeii = predict(lm_PCRr_gn_lakeii, newdata = .))

# Create the plot
ggplot(df2, aes(x = gonad_mass_mass_post_resp_g, y = CT_RBC_ratio, color = lake)) +
  geom_point() +
  geom_line(data = predictionslm_PCRr_gn_lakeii, aes(x = gonad_mass_mass_post_resp_g, y = predicted_lm_PCRr_gn_lakeii, color = lake)) +
  labs(title = "Telomere Length vs Relative Gonad Mass by Lake",
       x = "Relative Gonad Mass",
       y = " Relative Telomere Length RBCs",
       color = "Lake") +
  theme_minimal()

#### comparison CT values across tissues ####


# Subset the data for three specific variables
subset_data <- df2[, c("CT_liver_ratio", "CT_RBC_ratio", "CT_heart_ratio")]
cor_matrix <- cor(subset_data, use = "complete.obs")

# Correlation matrix
cor_matrix <- cor(subset_data)

# Display as a table
library(knitr)
kable(cor_matrix, caption = "Correlation Matrix for Tissue Relative Telomere Length")

# Table: Correlation Matrix for Tissue Relative Telomere Length
# 
# |               | CT_liver_ratio| CT_RBC_ratio| CT_heart_ratio|
#   |:--------------|--------------:|------------:|--------------:|
#   |CT_liver_ratio |      1.0000000|   -0.1228103|      0.3921196|
#   |CT_RBC_ratio   |     -0.1228103|    1.0000000|      0.1331776|
#   |CT_heart_ratio |      0.3921196|    0.1331776|      1.0000000|
# Using Hmisc for correlation with significance
library(Hmisc)
rcorr_matrix <- rcorr(as.matrix(data))

# Correlation coefficients
print(cor_matrix$r)

# P-values
print(rcorr_matrix$P)

# Install and load Hmisc package
install.packages("Hmisc")  # Run this if you haven't already installed Hmisc
library(Hmisc)



# Calculate correlations and p-values
results <- rcorr(as.matrix(subset_data))

# Correlation coefficients
cor_matrix <- results$r

# P-values
p_values <- results$P

# View results
print(cor_matrix)  # Correlation coefficients
# CT_liver_ratio CT_RBC_ratio CT_heart_ratio
# CT_liver_ratio      1.0000000   -0.1228103      0.3806350
# CT_RBC_ratio       -0.1228103    1.0000000      0.1331776
# CT_heart_ratio      0.3806350    0.1331776      1.0000000
print(p_values)    # P-values
# CT_liver_ratio CT_RBC_ratio CT_heart_ratio
# CT_liver_ratio             NA    0.5335560     0.03797932
# CT_RBC_ratio       0.53355600           NA     0.49929480
# CT_heart_ratio     0.03797932    0.4992948             NA

