#setwd("C:/Users/yinliao/OneDrive - Texas Tech University/Research/Prognosis Analysis - Cancer/Codes/Ying/Matlab/Bayesian clustering - R/CaseStudy/HPV/2Groups/Oropharynx")

############# Library #################

library(readxl)
library(writexl)
library(tictoc)

#1: Kruskal-Wallis test (one ANOVA on ranks)
#2: Pairwise comparisons using Wilcoxon rank sum test with continuity correction
tic()
#========================== HDP-WMM ==========================#

#-------------- Load data --------------#
Dat <- read_excel("data_censor_Original_seed_3_Chi.xlsx", col_names = TRUE)
ClusterEst <- read_excel("LabelEst_Overall_HDP.xlsx", col_names = TRUE)

Dat$ClusterEst <- as.factor(ClusterEst$HDP.cluster.est) 

Num_C = max(ClusterEst$HDP.cluster.est) # number of clusters
Num_C
Num_cluster = 1:Num_C  # number of obs for each cluster
for (i in 1:Num_C){
  Num_cluster[i] = sum(ClusterEst$HDP.cluster.est==i)
}

#-------------- age (continuous) ------------#

kruskal.test(Age ~ ClusterEst, data = Dat)
pairwise.wilcox.test(Dat$Age, Dat$ClusterEst,
                     p.adjust.method = "BH")

#-------------- age (categorical) ------------#
#var = Dat$`Age Standard for Survival (15-44,45-54,55-64,65-74,75+)`
#var = as.factor(var)
#levels(var)
#Dat$Age.discrete = var
#kruskal.test(Age.discrete ~ ClusterEst, data = Dat)
#pairwise.wilcox.test(Dat$Age.discrete, Dat$ClusterEst,
#                     p.adjust.method = "BH")

#-------------- sex ------------#
var = Dat$Sex
var = as.factor(var)
levels(var)
Dat$Sex = var
kruskal.test(Sex ~ ClusterEst, data = Dat)

#-------------- marital status ------------#
var = Dat$`Marital status at diagnosis`
var = as.factor(var)
levels(var)
Dat$Marital = var
kruskal.test(Marital ~ ClusterEst, data = Dat)

#-------------- race ------------#
var = Dat$`Race recode (White, Black, Other)`
var = as.factor(var)
levels(var)
Dat$Race = var
kruskal.test(Race ~ ClusterEst, data = Dat)

#-------------- Summary Stage ------------#
var = Dat$`Combined Summary Stage (2004+)`
var = as.factor(var)
levels(var)
Dat$Stage = var
kruskal.test(Stage ~ ClusterEst, data = Dat)

#-------------- Treatment ------------#
var = Dat$Treatment
var = as.factor(var)
levels(var)
Dat$Treatment = var
kruskal.test(Treatment ~ ClusterEst, data = Dat)







#========================== DP-WMM-All ==========================#
#-------------- Load data --------------#
Dat <- read_excel("data_censor_Original_seed_3_Chi.xlsx", col_names = TRUE)
ClusterEst <- read_excel("LabelEst_Overall_DP_All.xlsx", col_names = TRUE)

Dat$ClusterEst <- as.factor(ClusterEst$DP1.cluster.est) 

Num_C = max(ClusterEst$DP1.cluster.est) # number of clusters
Num_C
Num_cluster = 1:Num_C  # number of obs for each cluster
for (i in 1:Num_C){
  Num_cluster[i] = sum(ClusterEst$DP1.cluster.est==i)
}

#-------------- age (continuous) ------------#

kruskal.test(Age ~ ClusterEst, data = Dat)

#-------------- age (categorical) ------------#
#var = Dat$`Age Standard for Survival (15-44,45-54,55-64,65-74,75+)`
#var = as.factor(var)
#levels(var)
#Dat$Age.discrete = var
#kruskal.test(Age.discrete ~ ClusterEst, data = Dat)

#-------------- sex ------------#
var = Dat$Sex
var = as.factor(var)
levels(var)
Dat$Sex = var
kruskal.test(Sex ~ ClusterEst, data = Dat)

#-------------- marital status ------------#
var = Dat$`Marital status at diagnosis`
var = as.factor(var)
levels(var)
Dat$Marital = var
kruskal.test(Marital ~ ClusterEst, data = Dat)

#-------------- race ------------#
var = Dat$`Race recode (White, Black, Other)`
var = as.factor(var)
levels(var)
Dat$Race = var
kruskal.test(Race ~ ClusterEst, data = Dat)

#-------------- Summary Stage ------------#
var = Dat$`Combined Summary Stage (2004+)`
var = as.factor(var)
levels(var)
Dat$Stage = var
kruskal.test(Stage ~ ClusterEst, data = Dat)

#-------------- Treatment ------------#
var = Dat$Treatment
var = as.factor(var)
levels(var)
Dat$Treatment = var
kruskal.test(Treatment ~ ClusterEst, data = Dat)






#========================== DP-WMM-Each ==========================#
#-------------- Load data --------------#
Dat <- read_excel("data_censor_Original_seed_3_Chi.xlsx", col_names = TRUE)
ClusterEst <- read_excel("LabelEst_Overall_DP_Each.xlsx", col_names = TRUE)

Dat$ClusterEst <- as.factor(ClusterEst$DP2.cluster.est) 

Num_C = max(ClusterEst$DP2.cluster.est) # number of clusters
Num_C
Num_cluster = 1:Num_C  # number of obs for each cluster
for (i in 1:Num_C){
  Num_cluster[i] = sum(ClusterEst$DP2.cluster.est==i)
}

#-------------- age (continuous) ------------#

kruskal.test(Age ~ ClusterEst, data = Dat)

#-------------- age (categorical) ------------#
#var = Dat$`Age Standard for Survival (15-44,45-54,55-64,65-74,75+)`
#var = as.factor(var)
#levels(var)
#Dat$Age.discrete = var
#kruskal.test(Age.discrete ~ ClusterEst, data = Dat)

#-------------- sex ------------#
var = Dat$Sex
var = as.factor(var)
levels(var)
Dat$Sex = var
kruskal.test(Sex ~ ClusterEst, data = Dat)

#-------------- marital status ------------#
var = Dat$`Marital status at diagnosis`
var = as.factor(var)
levels(var)
Dat$Marital = var
kruskal.test(Marital ~ ClusterEst, data = Dat)

#-------------- race ------------#
var = Dat$`Race recode (White, Black, Other)`
var = as.factor(var)
levels(var)
Dat$Race = var
kruskal.test(Race ~ ClusterEst, data = Dat)

#-------------- Summary Stage ------------#
var = Dat$`Combined Summary Stage (2004+)`
var = as.factor(var)
levels(var)
Dat$Stage = var
kruskal.test(Stage ~ ClusterEst, data = Dat)

#-------------- Treatment ------------#
var = Dat$Treatment
var = as.factor(var)
levels(var)
Dat$Treatment = var
kruskal.test(Treatment ~ ClusterEst, data = Dat)

toc()
