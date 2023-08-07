
#setwd("C:/Users/yinliao/OneDrive - Texas Tech University/Research/Prognosis Analysis - Cancer/Codes/Ying/Matlab/IISE-Revision-1/results/SimulationStudy/seed_3/censor_5")

############################# Library #################################

# library(mcclust) # (2009) Process an MCMC Sample of Clusterings
library(mcclust.ext) #(2015) Point estimation and credible balls for Bayesian cluster analysis
library(readxl)
library(fields) # image.plot
library(latex2exp)
library(tictoc)

tic()
########################### Parameters ################################
nj = 200 # number of samples in each group
J = 2 # number of groups

n = nj * J
########################### Computation ##############################
#------------------------ True: 0.0 -----------------------#
idx.true <- read_excel("idx_true_0.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.0 = Label.true
psm_hc.true.0 = psm_hc.true

#------------------------ HDP-WMM: 0.0 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_0.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.0 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.0 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_0.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.0 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.0 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_0.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.0 = DP2.psm_hc




#------------------------ True: 0.1 -----------------------#
idx.true <- read_excel("idx_true_1.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.1 = Label.true
psm_hc.true.1 = psm_hc.true

#------------------------ HDP-WMM: 0.1 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_1.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.1 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.1 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_1.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.1 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.1 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_1.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.1 = DP2.psm_hc




#------------------------ True: 0.2 -----------------------#
idx.true <- read_excel("idx_true_2.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.2 = Label.true
psm_hc.true.2 = psm_hc.true

#------------------------ HDP-WMM: 0.2 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_2.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.2 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.2 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_2.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.2 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.2 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_2.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.2 = DP2.psm_hc



#------------------------ True: 0.3 -----------------------#
idx.true <- read_excel("idx_true_3.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.3 = Label.true
psm_hc.true.3 = psm_hc.true

#------------------------ HDP-WMM: 0.3 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_3.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.3 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.3 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_3.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.3 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.3 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_3.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.3 = DP2.psm_hc



#------------------------ True: 0.4 -----------------------#
idx.true <- read_excel("idx_true_4.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.4 = Label.true
psm_hc.true.4 = psm_hc.true

#------------------------ HDP-WMM: 0.4 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_4.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.4 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.4 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_4.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.4 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.4 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_4.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.4 = DP2.psm_hc



#------------------------ True: 0.5 -----------------------#
idx.true <- read_excel("idx_true_5.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]

Label.true.5 = Label.true
psm_hc.true.5 = psm_hc.true

#------------------------ HDP-WMM: 0.5 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM_5.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm
HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

HDP.psm_hc.5 = HDP.psm_hc

#------------------------ DP_WMM_All: 0.5 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All_5.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm
DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

DP1.psm_hc.5 = DP1.psm_hc

#------------------------ DP_WMM_Each: 0.5 -----------------------#
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each_5.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm
DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

DP2.psm_hc.5 = DP2.psm_hc





########################### Final figure ################################

#---------------------- Line Parameters -------------------------------#
table(Label.true.0[1,1:n])
image(1:n,1:n,psm_hc.true.0, ylab="", 
      xlab="(a) True clustering") # ,cex.axis=2
line1 = 136
line2 = line1 + 186
abline(v=line1,lty=2) 
abline(v=line2,lty=2) 

table(Label.true.1[1,1:n])
image(1:n,1:n,psm_hc.true.1, ylab="", 
      xlab="(a) True clustering") # ,cex.axis=2
line1 = 134
line2 = line1 + 182
abline(v=line1,lty=2) 
abline(v=line2,lty=2)

table(Label.true.2[1,1:n])
image(1:n,1:n,psm_hc.true.2, ylab="", 
      xlab="(a) True clustering") # ,cex.axis=2
line1 = 134
line2 = line1 + 182
abline(v=line1,lty=2) 
abline(v=line2,lty=2)

dev.off()

#---------------------- Figure Parameters -------------------------------#
#tiff("clustering_seed_3_censor_all.tif", width = 8.4, height = 15.6, units = 'in', res = 300)
png("clustering_seed_3_censor_all.png", width = 8, height = 7.8, units = 'in', res = 600)
par(mfrow=c(6,4),oma=c( 1,1,1,3),mar=c(0.5,4,1.5,1)) # mar=c(2,4.5,2,2)
# oma (outer): margin of 4 spaces width at right hand side; 
# mar (inner): bottom, left, top, and right

cex = 2
cex.label = 1 # 1.5
colorTable<- designer.colors(51, c( "white","yellow", "red") )
windowsFonts(A = windowsFont("Times New Roman"))

#---------------------- True clustering: 0.0 -------------------------------#
image(1:n,1:n,psm_hc.true.0, col=colorTable, ylab=TeX('$r_c=$0.0'), family="A",xlab="",
      main = "True",xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",, cex.lab=cex.label

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.0 -------------------------------#
image(1:n,1:n,HDP.psm_hc.0, col=colorTable, ylab="", family="A",xlab="",
      main = "HDP-WMM",xaxt='n',yaxt='n') # xlab="(b) HDP-WMM", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.0 -------------------------------#
image(1:n,1:n,DP1.psm_hc.0, col=colorTable, ylab="", family="A",xlab="",
      main = "DP-WMM-All",xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.0 -------------------------------#
image(1:n,1:n,DP2.psm_hc.0, col=colorTable, ylab="", family="A",xlab="",
      main = "DP-WMM-Each",xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)  




 
#---------------------- True clustering: 0.1 -------------------------------#
image(1:n,1:n,psm_hc.true.1, col=colorTable, ylab=TeX('$r_c=0.1$'), family="A",xlab="",
      xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.1 -------------------------------#
image(1:n,1:n,HDP.psm_hc.1, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(b) HDP-WMM",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.1 -------------------------------#
image(1:n,1:n,DP1.psm_hc.1, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.1 -------------------------------#
image(1:n,1:n,DP2.psm_hc.1, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)




#---------------------- True clustering: 0.2 -------------------------------#
image(1:n,1:n,psm_hc.true.2, col=colorTable, ylab=TeX('$r_c=0.2$'), family="A",xlab="",
      xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.2 -------------------------------#
image(1:n,1:n,HDP.psm_hc.2, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(b) HDP-WMM",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.2 -------------------------------#
image(1:n,1:n,DP1.psm_hc.2, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.2 -------------------------------#
image(1:n,1:n,DP2.psm_hc.2, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)





#---------------------- True clustering: 0.3 -------------------------------#
image(1:n,1:n,psm_hc.true.3, col=colorTable, ylab=TeX('$r_c=0.3$'), family="A",xlab="",
      xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.3 -------------------------------#
image(1:n,1:n,HDP.psm_hc.3, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(b) HDP-WMM",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.3 -------------------------------#
image(1:n,1:n,DP1.psm_hc.3, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.3 -------------------------------#
image(1:n,1:n,DP2.psm_hc.3, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)





#---------------------- True clustering: 0.4 -------------------------------#
image(1:n,1:n,psm_hc.true.4, col=colorTable, ylab=TeX('$r_c=0.4$'), family="A",xlab="",
      xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.4 -------------------------------#
image(1:n,1:n,HDP.psm_hc.4, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(b) HDP-WMM",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.4 -------------------------------#
image(1:n,1:n,DP1.psm_hc.4, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.4 -------------------------------#
image(1:n,1:n,DP2.psm_hc.4, col=colorTable, ylab="", family="A",xlab="",
      xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)






#---------------------- True clustering: 0.5 -------------------------------#
image(1:n,1:n,psm_hc.true.5, col=colorTable, ylab=TeX('$r_c=0.5$'), family="A",xlab="",
      xaxt='n',yaxt='n') # ,cex.axis=2, xlab="(a) True clustering",cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#---------------------- HDP-WMM: 0.5 -------------------------------#
image(1:n,1:n,HDP.psm_hc.5, col=colorTable, ylab="", family="A",xlab="",
     xaxt='n',yaxt='n') # xlab="(b) HDP-WMM", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


#---------------------- DP_WMM_All: 0.5 -------------------------------#
image(1:n,1:n,DP1.psm_hc.5, col=colorTable, ylab="", family="A",xlab="",
     xaxt='n',yaxt='n') # xlab="(c) DP-WMM-All", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 



#---------------------- DP_WMM_Each: 0.5 -------------------------------#
image(1:n,1:n,DP2.psm_hc.5, col=colorTable, ylab="", family="A",xlab="",
     xaxt='n',yaxt='n') # xlab="(d) DP-WMM-Each", cex.lab=cex.label,

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)

# oma (outer): margin of 4 spaces width at right hand side; 
# mar (inner): bottom, left, top, and right
par(oma=c( 1,0,1,1))# reset margin to be much smaller.
image.plot(1:n,1:n,psm_hc.true.0, legend.only=TRUE, col=colorTable,family="A") 

dev.off()

toc()
