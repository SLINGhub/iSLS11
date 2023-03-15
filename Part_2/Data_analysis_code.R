#################################
#### Synopsis
#################################

dir.create("Output")

### Read concentration data, plaque volume data, and sample injection sequence table
### (forcing the software to avoid converting special characters)
qdata = read.delim("data/SPERFECT_concentration_data.txt", 
                   header=TRUE, as.is=TRUE, check.names=FALSE)
pdata = read.delim("data/SPERFECT_plaque_data.txt", 
                   header=T, as.is=T, check.names=F)
sdata = read.delim("data/SPERFECT_MS_injection_sequence.txt", 
                   header=T, as.is=T, check.names=F) 

### Note that plaque volumes measured by computed tomography coronary angiography (CTCA)
### and not all individuals were included. 

############################################################
### Task 1. Quick practice to write functions 
### Goal: filter and match subjects across the data sets
############################################################

## Concentration data distribution check first
tmp = t(qdata[,-1])
par(mfrow=c(1,1))
boxplot(log2(tmp)) ## global distribution check
apply(log2(tmp), 2, function(x) mean(x, na.rm=TRUE))
## 376th column
rr = qdata$SampleID[376] 
ss = sdata$`Subject ID`[sdata$`Sample ID` == rr] 
# this individual has only three data points and we have to remove one.
# will later be disqualified from the analysis. 
qdata = qdata[-376, ] ## remove the row
sdata = sdata[sdata$`Sample ID` != rr, ]

## Parse sample ID into batch and run order variables
sdata$Batch = NA
sdata$OrderInBatch = NA
for(i in 1:nrow(sdata)) {
  tmp = strsplit(sdata$`Sample ID`[i], "_")[[1]]
  sdata$Batch[i] = gsub("batch", "", tmp[2])
  sdata$OrderInBatch[i] = tmp[3]
}
sdata$Batch = as.numeric(sdata$Batch)
sdata$OrderInBatch = as.numeric(sdata$OrderInBatch)

usub = unique(sdata$`Subject ID`)
sdata$SIDnum = match(sdata$`Subject ID`, usub) ## numerical ID of subjects

## Use the MS injection sequence data (sdata) as the anchoring point
## between concentration data and plaque volume data
counts = table(sdata$`Subject ID`)
patients = names(counts)
## or patients = unique(sdata$`Subject ID`)

## At this point we decide to remove patients with only two visits
remove.id = patients[counts < 3] ## two individuals will be lost: C103, C131
remove.index = sdata$`Subject ID` %in% remove.id
sdata = sdata[-remove.index, ]

## We now remove subjects without plaque volume data (outcomes)
## Plaque data (smaller, one time point)
all(pdata$`PATIENT CODE` %in% sdata$`Subject ID`)  ## evaluates to TRUE
mid = sdata$`Subject ID` %in% pdata$`PATIENT CODE`
sdata = sdata[mid, ]

## finally match the subjects with the concentration table using the sample ID as key
all(sdata$`Sample ID` %in% qdata$SampleID) ## should evaluate to TRUE
mid = match(sdata$`Sample ID`, qdata$SampleID)
qdata = qdata[mid, ]
rownames(qdata) = qdata$SampleID
qdata = qdata[,-1]

## now all data points are aligned

## we can export tables too.
lipids = colnames(qdata)[-1]
ltab = data.frame(Lipids=lipids, Class=NA, SubClass=NA,
                  stringsAsFactors=F, check.names=F)
write.table(ltab, "data/Lipid_table.txt", sep="\t", quote=F, row.names=F, na="")
## change file name, edit in Excel
ltab = read.delim("data/Lipid_table_editted.txt", header=T, as.is=T, check.names=F)

#############################################################################
### Task 2. Checking drift and batch effects
###   - Make a copy of the data and re-arrange in order of analysis sequence
###   - Plot the concentrations of individual species, with batch indicators
#############################################################################
ord = order(sdata$`Analysis Sequence`)
sdata2 = sdata[ord, ]
qdata2 = qdata[ord, ]

ticks = diff(as.numeric(sdata2$Batch))
ticks = which(ticks > 0)
nticks = length(ticks)

pdf("Output/analysis_sequence_plots.pdf", height=8, width=10, useDingbats=F)

par(mfrow=c(3,1))
for(k in 1:ncol(qdata2)) {
  plot(qdata2[,k], pch=19, cex=.5, col=sdata2$SIDbatch, 
       xlab="Analysis Sequence", ylab="Concentrations", main=colnames(qdata2)[k])
  for(j in 1:nticks) abline(v=ticks[j]+0.5, lty=2, col=2)
}

dev.off()

pdf("Output/analysis_sequence_plots_by_subject.pdf", height=8, width=10, useDingbats=F)

par(mfrow=c(3,1))
for(k in 1:ncol(qdata2)) {
  plot(qdata2[,k], pch=19, cex=.5, col=sdata2$SIDnum, 
       xlab="Analysis Sequence", ylab="Concentrations", main=colnames(qdata2)[k])
  for(j in 1:nticks) abline(v=ticks[j]+0.5, lty=2, col=2)
}

dev.off()

####### Lastly, we should check missing data
####### Track by rows and columns
any(qdata <= 0)
sum(is.na(qdata))
col.ct = apply(qdata, 2, function(x) sum(is.na(x)))
col.ct[col.ct > 0]
row.ct = apply(qdata, 1, function(x) sum(is.na(x)))
row.ct[row.ct > 0]
## Turns out the batch5_65 sample had S1P values missing!

## Decide to work with KNN imputation
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("impute")
library(impute)

qdata.tmp = log2(qdata)
qdata.new = impute.knn(as.matrix(qdata.tmp))$data
### Important note: In the original implementation,
### they find K nearest "genes" positioned across the rows of the matrix
### In lipidomics, we don't have that many "genes" (lipids), but better to
### find nearest K "samples" with similar patterns in other lipid patterns
### So we put samples in the rows, lipids in the columns
### (contrary to the instruction in impute.knn function)

### Check the imputed values for S1P species
gg = grep("S1P", colnames(qdata.new))
sid = ifelse(rownames(qdata.new) == "LT_batch5_65", 2, 1)  ## red dot
cid = ifelse(rownames(qdata.new) == "LT_batch5_65", 2, 0.2)
par(mfrow=c(3,2))
for(k in gg) {
  plot(qdata.new[,k], pch=19, cex=cid, col=sid, 
       xlab="Analysis Sequence", ylab="Concentrations", main=colnames(qdata2)[k])
}


################################################################
### Task 3. Visualization of entire data in
###  (i) Heatmaps with or without per individual normalization
###  (ii) Dimension reduction: PCA / PLS-DA
################################################################

## for downstream analysis, let's integrate plaque data (pdata) into the sample-level data (sdata)
mm = match(sdata$`Subject ID`, pdata$`PATIENT CODE`)
sdata = data.frame(sdata, pdata[mm, ],  ## you can use merge, dplyr::join family functions
                    stringsAsFactors=F, check.names=F)
q.med = apply(qdata.new, 2, median)
qdata.norm = sweep(qdata.new, 2, q.med) ## median normalized data, different from qdata.new (log2 transformed, imputed)

### Before moving onto downstream analysis, let's re-order the individuals
ord = order(sdata$`total lipid plaq vol index`, sdata$SIDnum)
sdata = sdata[ord, ]
qdata = qdata[ord, ]
qdata.new = qdata.new[ord, ]
qdata.norm = qdata.norm[ord, ]

## Complex heatmap with Plaque volume as row-wise information
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#install.packages("circlize")
#or devtools::install_github("jokergoo/circlize")
library(circlize)
library(ComplexHeatmap) 

set.seed(12345679)
row_ha = HeatmapAnnotation(Total=anno_barplot(sdata$`Total Plaque vol index`),
                           Calcified=anno_barplot(sdata$`total calc plaq vol index`),
                           LipidRich=anno_barplot(sdata$`total lipid plaq vol index`),
                           Fibrotic=anno_barplot(sdata$`total fibrot plaq vol index`),
                           annotation_name_rot = 270,
                           which="row", border=TRUE, show_legend=TRUE)

col_ha = HeatmapAnnotation(LipidClass=ltab$Class,
                           which="col", border=TRUE, show_legend=TRUE)

pdf("Output/heatmap_lipids.pdf", height=30, width=30)
Heatmap(as.matrix(qdata.norm), name = "Normalized Levels", 
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        column_title = "",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        left_annotation = row_ha,
        top_annotation = col_ha, 
        clustering_distance_rows = "pearson",
        clustering_method_rows = "average",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "average"
)
dev.off()

## Principal Component Analysis (PCA)
## Here we will plot the overall PCA plot first
## and track each individual's trajectory along the time course

tmp = qdata.norm
tmp.pca = prcomp(tmp, scale.=TRUE)
vv = tmp.pca$sdev^2
vv = round(vv / sum(vv) * 100, 2)
print(vv)

#install.packages("scales")
#install.packages("devtools")
#devtools::install_github("r-lib/scales")
library(scales)
### This is going to be a multi-page deck of plots

pdf("Output/PCAplot.pdf", height=12, width=11, useDingbats = FALSE)

  par(mfrow=c(3,3))
  THRES = max(abs(tmp.pca$x[,1]))
  dot.col = alpha(sdata$SIDnum, 0.5)
  dot.size = sdata$`Total Plaque vol index`
  dot.size = dot.size / mean(dot.size)  ## "normal" dot size is 1

  ##### OVERALL PLOT FIRST
  plot(tmp.pca$x[,1], tmp.pca$x[,2], 
       col=dot.col, pch=19, 
       cex = dot.size,  ### dot size proportional to total plaque volume
       xlab="PC1 (23.4%)", ylab="PC2 (21.3%)", main="SPERFECT", 
       xlim=c(-THRES,THRES), ylim=c(-THRES,THRES))
  abline(v=0, lty=2)
  abline(h=0, lty=2)

  #### Here we draw individual trajectory
  #### Per individual --> for loop
  subjects = unique(sdata$`Subject ID`)
  nsubjects = length(subjects)
  mm = match(subjects, sdata$`Subject ID`) ## mapping back to expanded table
  subject.col = sdata$SIDnum[mm]
  proj = tmp.pca$x[,1:2] ### PC1 and PC2 coordinates only

  for(i in 1:nsubjects) {
  
    ## background plot
    plot(tmp.pca$x[,1], tmp.pca$x[,2], 
         col=dot.col, cex=.5,
         xlab="PC1 (23.4%)", ylab="PC2 (21.3%)", main=subjects[i], 
         xlim=c(-THRES,THRES), ylim=c(-THRES,THRES))
    abline(v=0, lty=2)
    abline(h=0, lty=2)
  
    ## get data point indices for the corresponding individual
    coord = which(sdata$`Subject ID` == subjects[i])
    nn = length(coord)
  
    ## arrows
    for(k in 2:nn) {
      x0 = proj[coord[k-1],1]
      y0 = proj[coord[k-1],2]
      x1 = proj[coord[k],1]
      y1 = proj[coord[k],2]
      arrows(x0,y0,x1,y1, col=subject.col[i], length=0.07, lwd=1.5)    
    }

    ## text labels
    for(k in 1:nn) {
      x1 = proj[coord[k],1]
      y1 = proj[coord[k],2]
      text(x1, y1, subjects[i], cex=.4, col=subject.col[i])    
    }  
  
  }

dev.off()


################################################################
## Task 4. Clustering of subjects based on outcomes
################################################################

library(circlize)
library(ComplexHeatmap)

set.seed(12345679)
pdata2 = pdata[,-1]
rownames(pdata2) = pdata[,1]
for(k in 1:ncol(pdata2)) {
  pdata2[,k] = pdata2[,k] / mean(pdata2[,k])
}

pdf("Output/heatmap_subjects_outcome.pdf", height=10, width=5)
Heatmap(as.matrix(pdata2), name = "Normalized Volumes", 
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        column_title = "Plaque volume",
        col = colorRamp2(c(0, 3), c("white", "red")),
        clustering_distance_rows = "pearson",
        clustering_method_rows = "average",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "average"
)
dev.off()

## Let's use hierarchical clustering as it seems to make sense 
## (this may differ slightly from the data published in the paper)

dd = as.dist(1-cor(t(pdata2)))
hc=hclust(dd, method="average")
par(mfrow=c(1,1))
plot(hc, main="Pearson / Average Linkage", xlab="", sub="", cex=.9)
clus = cutree(hc, 3)
C1 = names(clus)[clus == 1] ## Lipid-rich plaque
C2 = names(clus)[clus == 2] ## Low plaque
C3 = names(clus)[clus == 3] ## Calcified plaque
### Try to interpret this

sdata$Cluster = NA
sdata$Cluster[sdata$`Subject ID` %in% C1] = "Lipid"
sdata$Cluster[sdata$`Subject ID` %in% C2] = "Low"
sdata$Cluster[sdata$`Subject ID` %in% C3] = "Calcified"
sdata$Cluster = factor(sdata$Cluster, levels=c("Low","Lipid","Calcified"))

########################################################################
### Task 5. Within-individual and between-group variability
### Exploratory analysis without formal LME-based decomposition
### (beyond the scope) via CV_g and CV_w approximation
########################################################################

### First, let's write a function to compute
### the coefficient of variation (CoV%) from log-transformed data
### Show Canchola et al paper
cov.logged = function(x, logbase=2) {
  var.x = var(x, na.rm=TRUE)
  cv = logbase^(log(logbase) * var.x) - 1  #log function in R is natural log
  sqrt(cv) * 100
}

### Intra-individual variability (CoV_w)
subjects = unique(sdata$`Subject ID`)
nsubjects = length(subjects)
lipids = colnames(qdata.new)
nlipids = length(lipids)
COVtab = matrix(NA, nlipids, nsubjects) ### placeholder to keep CV% values
rownames(COVtab) = lipids
colnames(COVtab) = subjects

for(i in 1:nlipids) {
  for(j in 1:nsubjects) {
    rid = which(sdata$`Subject ID` == subjects[j])
    COVtab[i,j] = cov.logged(qdata.new[rid,i])
  }
}

groups = sdata$Cluster[match(subjects, sdata$`Subject ID`)]
ord = order(groups)
COVtab = COVtab[,ord]
groups = groups[ord]

set.seed(12345679)
row_ha = HeatmapAnnotation(Class=ltab$Class,
                           annotation_name_rot = 270,
                           which="row", border=TRUE, show_legend=TRUE)
col_ha = HeatmapAnnotation(Group=groups,
                           which="col", border=TRUE, show_legend=TRUE)

pdf("Output/COV_w.pdf", height=25, width=15)
Heatmap(as.matrix(COVtab), name = "CoV", 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        column_title = "Plaque Group",
        col = colorRamp2(c(20,70), c("white", "red")),
        left_annotation = row_ha,
        top_annotation = col_ha,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "average",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "average"
)
dev.off()
## Try to match the lipids on top with the CV_w in Figure 3 of Tan et al. 

### How to estimate inter-individual? 
### We can approximate this by calculating CoV_g at the averaged log2 data
mm = match(subjects, sdata$`Subject ID`)
qdata.avg = qdata.new[mm, ]

for(i in 1:nsubjects) {
  rid = which(sdata$`Subject ID` == subjects[i])
  qdata.avg[i, ] = apply(qdata.new[rid, ], 2, mean)  ## averaging over replicates
}

COVg_vec = rep(NA, nlipids)
for(i in 1:nlipids) {
  COVg_vec[i] = cov.logged(qdata.avg[,i])
}

COVw_vec = apply(COVtab, 1, mean)

## colors for lipid classes
cl = as.numeric(factor(ltab$Class)) * 3
cl = alpha(cl, 0.6) ### opacity

plot(COVw_vec, COVg_vec, xlab="CoV% (intra)", ylab="CoV% (inter)", 
     xlim=c(0,100), ylim=c(0,100), 
     pch=19, col=cl)
abline(0,1.6,lty=2)  ## shows CoV_g > CoV_w
abline(0,1,lty=2)
## create a legend
uclass = unique(ltab$Class)
ucolor = cl[match(uclass,ltab$Class)]
legend("bottomright", legend=uclass, col=ucolor, pch=19, ncol=2)

## Note: this approximation-based analysis has an issue: 
## - It can over-estimate the CoV_g and under-estimate CoV_w (because they were averaged) - either or both.
## - The two variabilities have to be simultaneously deconvoluted, not separately. 
## - Purely analytical variability (CoV_a) was also ignored, but this can only be estimated
##   using repeated injections of a QC sample
##
## Here is the formal way to obtain the proper variance parameters 
## from the data (usually done from a population-scale data set)

library(nlme)
sid = sdata$`Subject ID`
sigma_g = rep(NA, nlipids)
sigma_w = rep(NA, nlipids)
names(sigma_g) = names(sigma_w) = ltab$Lipids
for(j in 1:nlipids) {
  if(j %% 20 == 0) print(j) ## print progress every 20th lipid
  y = qdata.new[,j]
  lme.fit = lme(y ~ 1, random = ~1|sid)
  vv = as.numeric(VarCorr(lme.fit)[,2])
  sigma_g[j] = vv[1]
  sigma_w[j] = vv[2]
}

## another function that directly translates standard deviation into CoV
cov.logged.sd = function(sd.x, logbase=2) {
  cv = logbase^(log(logbase) * sd.x^2) - 1  #log function in R is natural log
  sqrt(cv) * 100
}

COVw_vec = cov.logged.sd(sigma_w) ## intra-individual CoV (+ analytical)
COVg_vec = cov.logged.sd(sigma_g) ## inter-individual CoV

plot(COVw_vec, COVg_vec, xlab="CoV% (intra)", ylab="CoV% (inter)", 
     xlim=c(0,100), ylim=c(0,100), 
     pch=19, col=cl)
abline(0,1.2,lty=2)
abline(0,1,lty=1)
abline(0,1/1.2,lty=2)
legend("topleft", legend=uclass, col=ucolor, pch=19, ncol=6)
text(COVw_vec, COVg_vec, labels = ltab$Lipids, cex=.4, col=cl)
## This recovers the main Figure 3 in Tan et al. (except for a few species)
## Discussion point: which analytes would be the best clinical markers 
## if (hypothetically) all of them are candidates?

########################################################################
### Task 6. Association analysis 
### Individual variability (STDEV) ~ Plaque group (overall and each plaque type)
########################################################################

## Create a new data set where 
mm = match(subjects, sdata$`Subject ID`)
qdata.sd = qdata.new[mm, ]

for(i in 1:nsubjects) {
  rid = which(sdata$`Subject ID` == subjects[i])
  qdata.sd[i, ] = apply(qdata.new[rid, ], 2, sd)  ## averaging over replicates
}
rownames(qdata.sd) = subjects

Cluster = sdata$Cluster[mm]
TotalPL = sdata$`Total Plaque vol index`[mm]
CalcPL = sdata$`total calc plaq vol index`[mm]
LipidPL = sdata$`total lipid plaq vol index`[mm]
FibroPL = sdata$`total fibrot plaq vol index`[mm]

## Association test between in-person SD and plaques
tertiles = function(x) {
  qpt = quantile(x, c(0,1/3,2/3,1), na.rm=TRUE)
  nx = length(x)
  y = rep(1, nx)
  for(i in 1:nx) {
    for(k in 2:3) {
      if(x[i] > qpt[k] & x[i] <= qpt[k+1]) y[i] = k
    }
  }
  #y = factor(y, levels=c(1:3))
  y
}

quartiles = function(x) {
  qpt = quantile(x, c(0,0.25,0.5,0.75,1), na.rm=TRUE)
  nx = length(x)
  y = rep(1, nx)
  for(i in 1:nx) {
    for(k in 2:4) {
      if(x[i] > qpt[k] & x[i] <= qpt[k+1]) y[i] = k
    }
  }
  #y = factor(y, levels=c(1:4))
  y
}
quartiles(LipidPL) ## test

## Kruskal-Wallis (Non-parametric ANOVA)
pval.T = pval.F = pval.L = pval.C = rep(NA, nlipids)
names(pval.T) = names(pval.L) = names(pval.C) = colnames(qdata.sd)

for(j in 1:nlipids) {
  
  ## Total
  tmp.test = lm(qdata.sd[,j] ~ quartiles(TotalPL))
  tmp.test.0 = lm(qdata.sd[,j] ~ 1)
  pval.T[j] = anova(tmp.test, tmp.test.0)$`Pr(>F)`[2]
  ## Calcified
  tmp.test = lm(qdata.sd[,j] ~ tertiles(CalcPL))
  tmp.test.0 = lm(qdata.sd[,j] ~ 1)
  pval.C[j] = anova(tmp.test, tmp.test.0)$`Pr(>F)`[2]
  ## Lipid-rich
  tmp.test = lm(qdata.sd[,j] ~ quartiles(LipidPL))
  tmp.test.0 = lm(qdata.sd[,j] ~ 1)
  pval.L[j] = anova(tmp.test, tmp.test.0)$`Pr(>F)`[2]
  ## Fibrotic
  tmp.test = lm(qdata.sd[,j] ~ quartiles(FibroPL))
  tmp.test.0 = lm(qdata.sd[,j] ~ 1)
  pval.F[j] = anova(tmp.test, tmp.test.0)$`Pr(>F)`[2]
  
}


## With heatmap -- significant lipids only
qdata.sd.norm = qdata.sd
for(k in 1:ncol(qdata.sd)) qdata.sd.norm[,k] = qdata.sd.norm[,k] / mean(qdata.sd.norm[,k])

### Lipid-rich plque
ord = order(LipidPL)
col_ha = HeatmapAnnotation(Cluster=Cluster,
                           Total=anno_barplot(TotalPL),
                           Calcified=anno_barplot(CalcPL),
                           LipidRich=anno_barplot(LipidPL),
                           Fibrotic=anno_barplot(FibroPL),
                           which="col", border=TRUE, show_legend=TRUE)

sel = pval.L <= 0.1
row_ha = HeatmapAnnotation(LipidClass=ltab$Class[sel],
                           which="row", border=TRUE, show_legend=TRUE)
set.seed(12345679)

pdf("Output/SDassociation_heatmap_LipidRich.pdf", height=6, width=9)
X = as.matrix(qdata.sd.norm)[,sel]
Heatmap(t(X), name = "CoV", 
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 8),
        column_title = "Plaque Group",
        col = colorRamp2(c(0.5,2), c("white", "red")),
        left_annotation = row_ha,
        top_annotation = col_ha,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "average",
        clustering_distance_columns = "pearson",
        clustering_method_columns = "average"
)
dev.off()

### Results are slightly different from the published results
### beause we cannot reveal gender and age in this data set 
### (due to the risk of re-identifiability). 

