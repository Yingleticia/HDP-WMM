
#setwd("C:/Users/yinliao/OneDrive - Texas Tech University/Research/Prognosis Analysis - Cancer/Codes/Ying/Matlab/IISE-Revision-1/results/SimulationStudy/seed_1")

############################# Library #################################

# library(mcclust) # (2009) Process an MCMC Sample of Clusterings
library(mcclust.ext) #(2015) Point estimation and credible balls for Bayesian cluster analysis
library(readxl)
library(fields) # image.plot
library(tictoc)

tic()
########################### Parameters ################################
nj = 200 # number of samples in each group
J = 2 # number of groups

n = nj * J
########################### True Labels ################################

idx.true <- read_excel("idx_true.xlsx", col_names = FALSE)
Label.true <- matrix(unlist(idx.true), ncol = n)
all.true.psm = comp.psm(Label.true) 
plotpsm(all.true.psm)
hc.true = hclust(as.dist(1-all.true.psm), method = "complete", members = NULL)
psm_hc.true = all.true.psm
idx_order = hc.true$order
psm_hc.true[1:n,] = psm_hc.true[idx_order,]
psm_hc.true[,1:n] = psm_hc.true[,idx_order]


########################### HDP-WMM ################################
# B rows and N columns: B posterior samples of the clustering of the N data points
HDP.ClusterDraws <- read_excel("ClusterDraws_HDP_WMM.xlsx", col_names = FALSE)
HDP.Matrix_draw <- matrix(unlist(HDP.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
HDP.all.psm = comp.psm(HDP.Matrix_draw)
plotpsm(HDP.all.psm)
HDP.all.VI = minVI(HDP.all.psm,HDP.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(HDP.all.VI)

HDP.hc = hclust(as.dist(1-HDP.all.psm), method = "complete", members = NULL)
HDP.psm_hc = HDP.all.psm

#HDP.psm_hc[1:n,]=HDP.psm_hc[HDP.hc$order,]
#HDP.psm_hc[,1:n]=HDP.psm_hc[,HDP.hc$order]

HDP.psm_hc[1:n,]=HDP.psm_hc[idx_order,]
HDP.psm_hc[,1:n]=HDP.psm_hc[,idx_order]

#---------------------- HDP-WMM (accuracy) -------------------------------#
#HDP.cluster.est = HDP.all.VI$cl[1,]
# overall accuracy
#sum(HDP.cluster.est==Label.true[1,1:n])/n
# group-specific accuracy
#acc_HDP = rep(0,1,J)
#for (i in 1:J){
#  idx = ((i-1)*nj+1):(i*nj)
#  acc_HDP[i] = sum(HDP.cluster.est[idx]==Label.true[1,idx])/nj
#}
#acc_HDP

#---------------------- HDP-WMM (adjusted Rand index) -------------------------------#
HDP.cluster.est = HDP.all.VI$cl[1,]
# overall arandi
arandi(HDP.cluster.est, Label.true[1,1:n])
# group-specific arandi
arandi_HDP = rep(0,1,J)
for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  arandi_HDP[i] = arandi(HDP.cluster.est[idx],Label.true[1,idx])
}
arandi_HDP


#---------------------- HDP-WMM (clustering) -------------------------------#
# overall clustering
# true 
table(Label.true[1,1:n])
# estimate
table(HDP.cluster.est)

for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  paste(cat('Group ',as.character(i))) 
  print(table(Label.true[1,idx]))
  print(table(HDP.cluster.est[idx]))
}



########################### DP_WMM_All ################################
# B rows and N columns: B posterior samples of the clustering of the N data points
DP1.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_All.xlsx", col_names = FALSE)
DP1.Matrix_draw <- matrix(unlist(DP1.ClusterDraws), ncol = n)

# Bayesian Clustering Analysis #
DP1.all.psm = comp.psm(DP1.Matrix_draw)
plotpsm(DP1.all.psm)
DP1.all.VI = minVI(DP1.all.psm,DP1.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP1.all.VI)

DP1.hc = hclust(as.dist(1-DP1.all.psm), method = "complete", members = NULL)
DP1.psm_hc = DP1.all.psm

#DP1.psm_hc[1:n,]=DP1.psm_hc[DP1.hc$order,]
#DP1.psm_hc[,1:n]=DP1.psm_hc[,DP1.hc$order]

DP1.psm_hc[1:n,]=DP1.psm_hc[idx_order,]
DP1.psm_hc[,1:n]=DP1.psm_hc[,idx_order]

#---------------------- DP_WMM_All (adjusted Rand index) -------------------------------#
DP1.cluster.est = DP1.all.VI$cl[1,]
# overall arandi
arandi(DP1.cluster.est, Label.true[1,1:n])

# group-specific arandi
arandi_DP1 = rep(0,1,J)
for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  arandi_DP1[i] = arandi(DP1.cluster.est[idx],Label.true[1,idx])
}
arandi_DP1


#---------------------- DP_WMM_All (clustering) -------------------------------#
# overall clustering
# true 
table(Label.true[1,1:n])
# estimate
table(DP1.cluster.est)

# group1 clustering
for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  paste(cat('Group ',as.character(i))) 
  print(table(Label.true[1,idx]))
  print(table(DP1.cluster.est[idx]))
}



########################### DP_WMM_Each ################################
# B rows and N columns: B posterior samples of the clustering of the N data points
DP2.ClusterDraws <- read_excel("ClusterDraws_DP_WMM_Each.xlsx", col_names = FALSE)
DP2.Matrix_draw <- matrix(unlist(DP2.ClusterDraws), ncol = n)


# Bayesian Clustering Analysis #
DP2.all.psm = comp.psm(DP2.Matrix_draw)
plotpsm(DP2.all.psm)
DP2.all.VI = minVI(DP2.all.psm,DP2.Matrix_draw,method=("all"),include.greedy=TRUE) # Time costly at this step
summary(DP2.all.VI)

DP2.hc = hclust(as.dist(1-DP2.all.psm), method = "complete", members = NULL)
DP2.psm_hc = DP2.all.psm

#DP2.psm_hc[1:n,]=DP2.psm_hc[DP2.hc$order,]
#DP2.psm_hc[,1:n]=DP2.psm_hc[,DP2.hc$order]

DP2.psm_hc[1:n,]=DP2.psm_hc[idx_order,]
DP2.psm_hc[,1:n]=DP2.psm_hc[,idx_order]

#---------------------- DP_WMM_Each (adjusted Rand index) -------------------------------#
DP2.cluster.est = DP2.all.VI$cl[1,]

# overall arandi
arandi(DP2.cluster.est, Label.true[1,1:n])

# group-specific arandi
arandi_DP2 = rep(0,1,J)
for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  arandi_DP2[i] = arandi(DP2.cluster.est[idx],Label.true[1,idx])
}
arandi_DP2

#---------------------- DP_WMM_Each (clustering) -------------------------------#
# overall clustering
# true 
table(Label.true[1,1:n])
# estimate
table(DP2.cluster.est)


for (i in 1:J){
  idx = ((i-1)*nj+1):(i*nj)
  paste(cat('Group ',as.character(i))) 
  print(table(Label.true[1,idx]))
  print(table(DP2.cluster.est[idx]))
}




########################### Final figure ################################

#---------------------- Line Parameters -------------------------------#
table(Label.true[1,1:n])
image(1:n,1:n,psm_hc.true, ylab="", 
      xlab="(a) True clustering") # ,cex.axis=2
line1 = 134
line2 = line1 + 182
abline(v=line1,lty=2)
abline(v=line2,lty=2)

dev.off()

#---------------------- Figure Parameters -------------------------------#
#tiff("clustering_seed_3_censor_0.tif", width = 8.4, height = 2.6, units = 'in', res = 300)
png("clustering_seed_3_censor_0.png", width = 8.4, height = 2.6, units = 'in', res = 300)
par(mfcol=c(1,4),oma=c( 1,1,1,3),mar=c(4,2,2,2)) 
# oma (outer): margin of 4 spaces width at right hand side; 
# mar (inner): bottom, left, top, and right

cex = 2
cex.label = 1.3
colorTable<- designer.colors(51, c( "white","yellow", "red") )
windowsFonts(A = windowsFont("Times New Roman"))

#---------------------- True clustering (figure) -------------------------------#
image(1:n,1:n,psm_hc.true, col=colorTable, ylab="", family="A",
           xlab="(a) True clustering",cex.lab=cex.label) # ,cex.axis=2

#image.plot(1:n,1:n,psm_hc.true, col=colorTable, ylab="", 
#           xlab="(a) True clustering",cex.lab=cex.label) # ,cex.axis=2


abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 

#points( 151/2, 151/2, pch="3", cex=cex)
#points( 151+169/2, 151+169/2, pch="1", cex=cex)
#points( 151+169+80/2, 151+169+80/2, pch="2", cex=cex)


#---------------------- HDP-WMM (figure) -------------------------------#
image(1:n,1:n,HDP.psm_hc, col=colorTable, ylab="", family="A",
           xlab="(b) HDP-WMM",cex.lab=cex.label)

#image.plot(1:n,1:n,HDP.psm_hc, col=colorTable, ylab="", 
#           xlab="(b) HDP-WMM",cex.lab=cex.label)

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


# estimate
#table(HDP.cluster.est)

#abline(v=0) 
#abline(h=400) 
### !!! edits needed !!! ###
#abline(v=155,lty=2) 
#abline(h=155,lty=2) 
#abline(v=155+85,lty=2) # 122+64=186
#abline(h=155+85,lty=2)

#points( 155/2, 155/2, pch="3", cex=cex)
#points( 155+85/2, 155+85/2, pch="2", cex=cex)
#points( 155+85+160/2, 155+85+160/2, pch="1", cex=cex)


#---------------------- DP_WMM_All (figure) -------------------------------#
image(1:n,1:n,DP1.psm_hc, col=colorTable, ylab="", family="A",
           xlab="(c) DP-WMM-All",cex.lab=cex.label)

#image.plot(1:n,1:n,DP1.psm_hc, col=colorTable, ylab="", 
#           xlab="(c) DP-WMM-All",cex.lab=cex.label)

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2) 


# estimate
#DP1.cluster.est = DP1.all.VI$cl[1,]
#table(DP1.cluster.est)

#abline(v=0) 
#abline(h=400) 
### !!! edits needed !!! ###
#abline(v=155,lty=2) 
#abline(h=155,lty=2) 

#points( 155/2, 155/2, pch="2", cex=cex)
#points( 155+245/2, 155+245/2, pch="1", cex=cex)


#---------------------- DP_WMM_Each (figure) -------------------------------#
image(1:n,1:n,DP2.psm_hc, col=colorTable, ylab="", family="A",
           xlab="(d) DP-WMM-Each",cex.lab=cex.label)

#image.plot(1:n,1:n,DP2.psm_hc, col=colorTable, ylab="", 
#           xlab="(d) DP-WMM-Each",cex.lab=cex.label)

abline(v=0) 
abline(h=n) 
abline(v=line1,lty=2) 
abline(h=line1,lty=2) 
abline(v=line2,lty=2) 
abline(h=line2,lty=2)  

# estimate
#table(DP2.cluster.est)

#abline(v=0) 
#abline(h=400) 
### !!! edits needed !!! ###
#abline(v=47,lty=2) 
#abline(h=47,lty=2) 
#abline(v=47 + 144,lty=2) 
#abline(h=47 + 144,lty=2)
#abline(v=47 + 144 + 9,lty=2) 
#abline(h=47 + 144 + 9,lty=2)
#abline(v=47 + 144 + 9 + 110,lty=2) 
#abline(h=47 + 144 + 9 + 110,lty=2)

#points( 47/2, 47/2, pch="5", cex=cex)
#points( 47+144/2, 47+144/2, pch="3", cex=cex)
#points( 47+144+9/2, 47+144+9/2, pch="4", cex=cex)
#points( 47+144+9+110/2, 47+144+9+110/2, pch="2", cex=cex)
#points( 47+144+9+110+90/2, 47+144+9+110+90/2, pch="1", cex=cex)

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot(1:n,1:n,psm_hc.true, legend.only=TRUE, col=colorTable,family="A") 

dev.off()
toc()
