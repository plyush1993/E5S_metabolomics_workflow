##################################################################################
# Primary step: initial settings
##################################################################################

# set working directory
setwd("C:/...")

# required packages
library(caret)
library(tuple)
library(factoextra)
library(FactoMineR)
library(rafalib)
library(data.table)
library(parallel)
library(doParallel)
library(RSEIS)
library(ggsci)
library(cowplot)
library(ggplot2)

# stop parallel processing
stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()

# start parallel processing
fc <- as.numeric(detectCores(logical = F))
cl <- makePSOCKcluster(fc+1)
registerDoParallel(cl)

##################################################################################
# 1st step: Normalization
##################################################################################


# EigenMS normalization
library(ProteoMM)
library(data.table)
setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T)) # load data table with ordered sample name (order as in experimental batches) with Label
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

m_logInts = as.data.frame(t(data[,-1]))
m_logInts = convert_log2(m_logInts)
grps = as.factor(c(data[,1]))
m_prot.info = cbind(colnames(data)[-1], colnames(data)[-1])
m_ints_eig1 = eig_norm1(m=m_logInts,treatment=grps,prot.info=m_prot.info)
m_ints_norm1 = eig_norm2(rv=m_ints_eig1)

data_em <- data
data_em[,-1] <- as.data.frame(t(m_ints_norm1$norm_m))
data_em <- cbind(rownames(data), data_em)
fwrite(data_em, "dsr.csv")


# Combat normalization
library(sva)
library(data.table)
setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

data_combat <- as.matrix(t(data[,-1]))

# for batch in ComBat
pheno_data <-as.data.frame(fread(input = ".csv", header=T)) # load pheno_data with columns for Batch and Label, Order (order as in data file)
batch <- pheno_data$Batch

combat_adjdata1 = ComBat(dat=data_combat, batch=batch, mod=NULL, par.prior=T, prior.plots=F, mean.only = F)

data_c <- data
data_c[,-1] <- as.data.frame(t(combat_adjdata1))
data_c <- cbind(rownames(data), data_c)
fwrite(data_c, ".csv")


# VSN normalization
library(vsn)
library(data.table)
setwd("C:/...")
data <-as.data.frame(fread(input = ".csv", header=T))
row.names(data) <- data[,1]
data <- data[,-1]
colnames(data)[1] <- c("Label")

data_vsn <- as.matrix(t(data[,-1]))

VSN <- function(data) {
  vsn.model <- vsn2(data)
  vsn.data <- predict(vsn.model, data)
  return(vsn.data)
}

vsn_data <- as.data.frame(VSN(data_vsn))

data_v <- data
data_v[,-1] <- as.data.frame(t(vsn_data))
data_v <- cbind(rownames(data), data_v)
fwrite(data_v, ".csv")

####################################################################################
# 2nd step: Half Minimum Missing Value Imputation + Univariate Filtering (HM + UVF)
####################################################################################

# Load data
#by data.table package
dsr <-as.data.frame(fread(input = ".csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
colnames(dsr)[1] <- c("Label")

# Half minimum missing value imputation
# by apply
HM_apply <- function(data) {
  result <- data
  result <- apply(result, 2, function(x) {
    x[is.na(x)] <- min(x, na.rm = T)/2
    x
  })
  return(result)
}

dsr <- as.data.frame(cbind(Label = dsr[,1], HM_apply(dsr[,-1]))) # exclude Label column

# Univariate filtering
uvf <- function(x){
  
  norm_homog_tests <- function(x) {
    xx <- x[,-1]
    # normality test
    norm.test <- apply(xx, 2, function(t) shapiro.test(t)$p.value)
    # homogeneity test
    homog.test <- apply(xx, 2, function(t) bartlett.test(t,g = x[,1])$p.value)
    return(as.data.frame(cbind(norm.test, homog.test)))}
  
  res_tests <- norm_homog_tests(x)
  
  
  wilcox_test <- function(x,y) {
    xx <- x[,-1]
    wx.t <- as.vector(which(y[,1] < 0.05))
    wilcox_test <- list()
    ifelse(identical(wx.t, integer(0)), return (wilcox_test <- 1), wx.t)
    wilcox_test <- lapply(as.data.frame(xx[,wx.t]),  function(t) as.vector(pairwise.wilcox.test(x = t, g =  x[,1], p.adjust.method ="BH", paired=F)$p.value))
    names(wilcox_test) <- (colnames(x)[-1])[wx.t]
    return(as.list(wilcox_test))}
  
  wx.t.res <- wilcox_test(x, res_tests)
  
  
  welch_test <- function(x,y) {
    xx <- x[,-1]
    wl.t <- as.vector(which(y[,1] > 0.05 & y[,2] < 0.05))
    welch_test <- list()
    ifelse(identical(wl.t, integer(0)), return (welch_test <- 1), wl.t)
    welch_test <- lapply(as.data.frame(xx[,wl.t]), function(t) as.vector(pairwise.t.test(x = t, g = x[,1], p.adjust.method = "BH", pool.sd = F)$p.value))
    names(welch_test) <- (colnames(x)[-1])[wl.t]
    return(as.list(welch_test))}
  
  wl.t.res <- welch_test(x, res_tests)
  
  
  student_test <- function(x,y) {
    xx <- x[,-1]
    st.t <- as.vector(which(y[,1] > 0.05 & y[,2] > 0.05))
    student_test <- list()
    ifelse(identical(st.t, integer(0)), return (student_test <- 1), st.t)
    student_test <- lapply(as.data.frame(xx[,st.t]), function(t) as.vector(pairwise.t.test(x = t, g = x[,1], p.adjust.method = "BH", pool.sd = T)$p.value))
    names(student_test) <- (colnames(x)[-1])[st.t]
    return(as.list(student_test))}
  
  st.t.res <- student_test(x, res_tests)
  
  filt_p_val <- function(x, y, z, w){
    
    #x = ds
    #y = wx.t.res
    #z = wl.t.res
    #w = st.t.res
    
    wx.t.n <- names(y)
    wx.t.res2 <-as.data.frame(y)
    colnames(wx.t.res2) <- wx.t.n
    wx.t.res2[is.na(wx.t.res2)] <- max(wx.t.res2, na.rm = T)
    
    wl.t.n <- names(z)
    wl.t.res2 <- as.data.frame(z)
    colnames(wl.t.res2) <- wl.t.n
    wl.t.res2[is.na(wl.t.res2)] <- max(wl.t.res2, na.rm = T)
    
    st.t.n <- names(w)
    st.t.res2 <- as.data.frame(w)
    colnames(st.t.res2) <- st.t.n
    st.t.res2[is.na(st.t.res2)] <- max(st.t.res2, na.rm = T)
    
    wxx <- apply(wx.t.res2, 2, min)
    wll <- apply(wl.t.res2, 2, min)
    stt <- apply(st.t.res2, 2, min)
    wxx2 <- as.data.frame(wxx)
    wll2 <- as.data.frame(wll)
    stt2 <- as.data.frame(stt)
    
    wxx3 <- rownames(wxx2)[which(wxx2 <= 0.05)]
    wll3 <- rownames(wll2)[which(wll2 <= 0.05)]
    stt3 <- rownames(stt2)[which(stt2 <= 0.05)]
    aff <- c(wxx3, wll3, stt3)
    
    ds_fil <- cbind(x[,1], x[, aff])
    return(ds_fil)
  }
  return(filt_p_val(x, wx.t.res, wl.t.res, st.t.res))
}

ds <- uvf(dsr)
colnames(ds)[1] <- "Label"

####################################################################################
# 3rd step: Machine Learning + Stable Feature Extraction (ML+SFE)
####################################################################################

#repeated cross validation

trainControl <- trainControl(method="repeatedcv", number=5, repeats=3)
metric <- "Accuracy"

# KNN
set.seed(1234)
fit.knn <- train(Label~ ., data=ds, method="knn", metric=metric, trControl=trainControl, tuneLength = 5)

# SVM
set.seed(1234)
fit.svm <- train(Label~ ., data=ds, method="svmRadial", metric=metric, trControl=trainControl, tuneLength = 5)

# PLS
set.seed(1234)
fit.pls <- train(Label~ ., data=ds, method="pls", metric=metric, trControl=trainControl, tuneLength = 5)

# RF
set.seed(1234)
fit.rf <- train(Label~ ., data=ds, method="rf", metric=metric, trControl=trainControl, tuneLength = 5)

# only accuracy for all models
results <- resamples(list(knn=fit.knn, svm=fit.svm, rf=fit.rf, pls=fit.pls), trControl = trainControl, metric=metric)
results_df <- as.data.frame(results)
results_df_fin <- apply(results_df[,-5], 2, mean)
results_df_fin

# sort by model specific value
set.seed(1234)
Imp.rf <- varImp(fit.rf, scale = FALSE)
Imp.rf <- Imp.rf$importance
rownames(Imp.rf) <- gsub("`", '', rownames(Imp.rf))

set.seed(1234)
Imp.pls <- varImp(fit.pls, scale = FALSE)
Imp.pls <- Imp.pls$importance
rownames(Imp.pls) <- gsub("`", '', rownames(Imp.pls))

set.seed(1234)
Imp.svm <- varImp(fit.svm, scale = FALSE)
Imp.svm <- Imp.svm$importance
rownames(Imp.svm) <- gsub("`", '', rownames(Imp.svm))

set.seed(1234)
Imp.knn <- varImp(fit.knn, scale = FALSE)
Imp.knn <- Imp.knn$importance
rownames(Imp.knn) <- gsub("`", '', rownames(Imp.knn))

# creating of list with top = n important features by model
n_model = 50
set.seed(1234)

Imp.rf[,c("sum")] <- apply(X = Imp.rf, 1, FUN = sum)
Imp.rf <- Imp.rf[order(Imp.rf$sum, decreasing = T),]
Imp.rf <- rownames(Imp.rf)[c(1:n_model)]

Imp.pls[,c("sum")] <- apply(X = Imp.pls, 1, FUN = sum)
Imp.pls <- Imp.pls[order(Imp.pls$sum, decreasing = T),]
Imp.pls <- rownames(Imp.pls)[c(1:n_model)]

Imp.svm[,c("sum")] <- apply(X = Imp.svm, 1, FUN = sum)
Imp.svm <- Imp.svm[order(Imp.svm$sum, decreasing = T),]
Imp.svm <- rownames(Imp.svm)[c(1:n_model)]

Imp.knn[,c("sum")] <- apply(X = Imp.knn, 1, FUN = sum)
Imp.knn <- Imp.knn[order(Imp.knn$sum, decreasing = T),]
Imp.knn <- rownames(Imp.knn)[c(1:n_model)]

# minimum n time of duplicated 
all <- c(Imp.rf, Imp.svm, Imp.knn, Imp.pls)
#all <- unique(c(Imp.svm, Imp.rf, Imp.pls)) 
#all <- unique(c(Imp.svm,Imp.rf)) 
n_tuple <- 2
all1 <- all[which(tuplicated(all, n_tuple), T)]
all1 <- unique(all1)
ds_d <- cbind(ds[,1], ds[,all1])
#ds_d <- cbind(ds[,1], ds[,all]) 

####################################################################################
# 4th step: Recursive Feature Selection (RFS)
####################################################################################

set.seed(1234)
subsets <- c(1:(ncol(ds_d)-1))
#ctrl_rfe <- rfeControl(functions = rfFuncs,method = "repeatedcv",number = 5, repeats = 5,verbose = FALSE)
ctrl_rfe <- rfeControl(functions = nbFuncs,method = "repeatedcv",number = 5, repeats = 3,verbose = FALSE)
rfe <- rfe(ds_d[,-1],ds_d[,1], sizes = subsets, metric = "Accuracy", rfeControl = ctrl_rfe)
rfe_vi <- rfe[["optVariables"]]
ds_rfe <- data.frame(cbind(as.character(ds[,1]), ds[, rfe_vi]))

####################################################################################
# Supp. step 1: Area Under Receiver Operating Characteristic (AUROC) Calculations
####################################################################################

# raw data (ds_raw)
ds_raw <-as.data.frame(fread(input = "ds_raw.csv", header=T))
rownames(ds_raw) <- ds_raw[,1]
ds_raw <- ds_raw[,-1]
colnames(ds_raw)[1] <- c("Label")

# ds_raw after EigenMS (dsr)
dsr <-as.data.frame(fread(input = "dsr.csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
colnames(dsr)[1] <- c("Label")

# dsr after UVF+MVI (ds)
ds <-as.data.frame(fread(input = "ds.csv", header=T))
rownames(ds) <- ds[,1]
ds <- ds[,-1]
colnames(ds)[1] <- c("Label")

# ds after ML+SFE (ds_d)
ds_d <-as.data.frame(fread(input = "ds_d.csv", header=T))
rownames(ds_d) <- ds_d[,1]
ds_d <- ds_d[,-1]
colnames(ds_d)[1] <- c("Label")

# ds_d after RFE (ds_rfe)
ds_rfe <-as.data.frame(fread(input = "ds_rfe.csv", header=T))
rownames(ds_rfe) <- ds_rfe[,1]
ds_rfe <- ds_rfe[,-1]
colnames(ds_rfe)[1] <- c("Label")

# mean AUROC for ds_raw
set.seed(1234)
auroc_ds_raw <- filterVarImp(ds_raw[,-1], as.factor(ds_raw[,1]))
auroc_ds_raw[,c("sum")] <- apply(X = auroc_ds_raw, 1, FUN = sum)
auroc_ds_raw_mean <- round(mean(auroc_ds_raw$sum/(ncol(auroc_ds_raw)-1)),3)

# mean AUROC for dsr
set.seed(1234)
auroc_dsr <- filterVarImp(dsr[,-1], as.factor(dsr[,1]))
auroc_dsr[,c("sum")] <- apply(X = auroc_dsr, 1, FUN = sum)
auroc_dsr_mean <- round(mean(auroc_dsr$sum/(ncol(auroc_dsr)-1)),3)

# mean AUROC for ds
set.seed(1234)
auroc_ds <- filterVarImp(ds[,-1], as.factor(ds[,1]))
auroc_ds[,c("sum")] <- apply(X = auroc_ds, 1, FUN = sum)
auroc_ds_mean <- round(mean(auroc_ds$sum/(ncol(auroc_ds)-1)),3)

# mean AUROC for ds_d
set.seed(1234)
auroc_ds_d <- filterVarImp(ds_d[,-1], as.factor(ds_d[,1]))
auroc_ds_d[,c("sum")] <- apply(X = auroc_ds_d, 1, FUN = sum)
auroc_ds_d_mean <- round(mean(auroc_ds_d$sum/(ncol(auroc_ds_d)-1)),3)

# mean AUROC for ds_rfe
set.seed(1234)
auroc_ds_rfe <- filterVarImp(ds_rfe[,-1], as.factor(ds_rfe[,1]))
auroc_ds_rfe[,c("sum")] <- apply(X = auroc_ds_rfe, 1, FUN = sum)
auroc_ds_rfe_mean <- round(mean(auroc_ds_rfe$sum/(ncol(auroc_ds_rfe)-1)),3)

# final table with AUROC calculations & # of features
auroc <- data.frame(rbind(auroc_ds_raw_mean, auroc_dsr_mean, auroc_ds_mean, auroc_ds_d_mean, auroc_ds_rfe_mean))
vars_roc <- data.frame(cbind(c((ncol(ds_raw)-1), (ncol(dsr)-1),(ncol(ds)-1),(ncol(ds_d)-1),(ncol(ds_rfe)-1)), auroc))
row.names(vars_roc) <- c("ds_raw", "dsr", "ds", "ds_d", "ds_rfe")
colnames(vars_roc) <- c("# of features","mean AUROC")
vars_roc

# about datasets
# ds_d after RFE (ds_rfe)
ds_rfe <-as.data.frame(fread(input = "ds_rfe.csv", header=T))
rownames(ds_rfe) <- ds_rfe[,1]
ds_rfe <- ds_rfe[,-1]
colnames(ds_rfe)[1] <- c("Label")

# final table with AUROC calculations & # of features (by ds_rfe)
sam <- data.frame(rbind(c(nrow(ds_rfe))))
sam_gr <- data.frame(cbind(sam,c(length(unique(ds_rfe$Label)))))
row.names(sam_gr) <- c("ds_rfe")
colnames(sam_gr) <- c("# samples","# gropus")
sam_gr


####################################################################################
# 5th step: Unsupervised Learning & Projections
####################################################################################

# raw data+MVI+UVF
ds_raw <-as.data.frame(fread(input = "ds_raw_ds.csv", header=T))
rownames(ds_raw) <- ds_raw[,1]
ds_raw <- ds_raw[,-1]
colnames(ds_raw)[1] <- c("Label")

# ds after EigenMS+MVI+UVF
dsr <-as.data.frame(fread(input = "ds.csv", header=T))
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
colnames(dsr)[1] <- c("Label")

# ds_d after RFE
ds_rfe <-as.data.frame(fread(input = "ds_rfe.csv", header=T))
rownames(ds_rfe) <- ds_rfe[,1]
ds_rfe <- ds_rfe[,-1]
colnames(ds_rfe)[1] <- c("Label")

# for plotting
# number of groups
k <- length(unique(ds_rfe$Label)) # groups in HC

# ds_raw
base1 <- ds_raw
mtrx1 <- base1[,-1] 
grp1 <- base1[,1]

# dsr
base2 <- dsr
mtrx2 <- base2[,-1] 
grp2 <- base2[,1] 

# ds_rfe
base3 <- ds_rfe
mtrx3 <- base3[,-1] 
grp3 <- base3[,1] 

# PCA (a, b)
palette_pca <- "lancet" # color
# palette_pca <- JGRAY(length(unique(grp1))) # grey

# a
pca.ds1 <- PCA(mtrx1, scale.unit = T, graph = FALSE)
a <- fviz_pca_ind(pca.ds1,
             title = "",
             geom.ind = "point", # show points only 
             col.ind = grp1, # color by groups
             palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp1)))
             addEllipses = T, # Concentration ellipses
             legend.title = "Groups")

# b
pca.ds2 <- PCA(mtrx2, scale.unit = T, graph = FALSE)
b <- fviz_pca_ind(pca.ds2,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp2, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp2)))
                  addEllipses = T, # Concentration ellipses
                  legend.title = "",
                  )

# c
pca.ds3 <- PCA(mtrx3, scale.unit = T, graph = FALSE)
c <- fviz_pca_ind(pca.ds3,
                  title = "",
                  geom.ind = "point", # show points only 
                  col.ind = grp3, # color by groups
                  palette = palette_pca, # color "jco" gray JGRAY(length(unique(grp3)))
                  addEllipses = T, # Concentration ellipses
                  legend.title = "Groups")

# HCA (c,d)
# supplementary function
# color
Cols = function(vec, ord){
  cols = pal_lancet(palette = c("lanonc"), alpha = 1)(length(unique(vec)))
  return(cols[as.fumeric(vec)[ord]])}

# grey
#Cols = function(vec, ord){
 # cols = JGRAY(length(unique(vec)))
# return(cols[as.fumeric(vec)[ord]])}

# d 
mtrx1_1 <- mtrx1
#mtrx1_1 <- data.frame(scale(mtrx1, center = T, scale = T))
rownames(mtrx1_1) = make.names(grp1, unique=TRUE)
res.dist1 <- dist(mtrx1_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc1 <- hclust(d = res.dist1, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
d <- fviz_dend(res.hc1, k = k, # Cut in k groups
          cex = 0.3, # label size
          k_colors = unique(Cols(grp1,res.hc1$order)), # "lancet" color "jco" gray JGRAY(k_hc)
          color_labels_by_k = F, # color labels by groups
          label_cols = Cols(grp1,res.hc1$order),#Cols(ds_rfe[,1])[res.hc1$order], #as.fumeric(ds_rfe[,1])[res.hc1$order],
          rect = T, # Add rectangle around groups
          rect_fill = T,
          rect_border = unique(Cols(grp1,res.hc1$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
          horiz = F,
          lwd = 0.7, # lines size
          show_labels = T,
          main = "",
          ylab = "")

# e
mtrx2_1 <- mtrx2
#mtrx2_1 <- data.frame(scale(mtrx2, center = T, scale = T))
rownames(mtrx2_1) = make.names(grp2, unique=TRUE)
res.dist2 <- dist(mtrx2_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc2 <- hclust(d = res.dist2, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
e <- fviz_dend(res.hc2, k = k, # Cut in k groups
               cex = 0.3, # label size
               k_colors = unique(Cols(grp2,res.hc2$order)), # "lancet" color "jco" gray JGRAY(k_hc)
               color_labels_by_k = F, # color labels by groups
               label_cols = Cols(grp2,res.hc2$order),#Cols(ds_rfe[,1])[res.hc2$order], #as.fumeric(ds_rfe[,1])[res.hc2$order],
               rect = T, # Add rectangle around groups
               rect_fill = T,
               rect_border = unique(Cols(grp2,res.hc2$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
               horiz = F,
               lwd = 0.7, # lines size
               show_labels = T,
               main = "",
               ylab = "")

# f
mtrx3_1 <- mtrx3
#mtrx3_1 <- data.frame(scale(mtrx3, center = T, scale = T))
rownames(mtrx3_1) = make.names(grp3, unique=TRUE)
res.dist3 <- dist(mtrx3_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc3 <- hclust(d = res.dist3, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
f <- fviz_dend(res.hc3, k = k, # Cut in k groups
               cex = 0.3, # label size
               k_colors = unique(Cols(grp3,res.hc3$order)), # "lancet" color "jco" gray JGRAY(k_hc)
               color_labels_by_k = F, # color labels by groups
               label_cols = Cols(grp3,res.hc3$order),#Cols(ds_rfe[,1])[res.hc3$order], #as.fumeric(ds_rfe[,1])[res.hc3$order],
               rect = T, # Add rectangle around groups
               rect_fill = T,
               rect_border = unique(Cols(grp3,res.hc3$order)), #"lancet"# color "jco" gray JGRAY(k_hc)
               horiz = F,
               lwd = 0.7, # lines size
               show_labels = T,
               main = "",
               ylab = "")

# creating & saving multiplot
# PCA
pca_plot <- plot_grid(a,b,c,  labels = c("A", "B", "C"), ncol = 3, nrow = 1)

# HCA
hca_plot <- plot_grid(d,e,f,  labels = c("A", "B", "C"), ncol = 3, nrow = 1)

# saving plots
ggsave("hca.png", hca_plot, dpi = 600, height = 10, width = 20, limitsize = F, units = "in")

ggsave("pca.png", pca_plot, dpi = 600, height = 10, width = 20, limitsize = F, units = "in")

####################################################################################
# Supp. step 2: Save Outputs
####################################################################################

fwrite(ds, "ds.csv", row.names = T)

fwrite(ds, "ds_raw_ds.csv", row.names = T)

fwrite(ds_d, "ds_d.csv", row.names = T)

fwrite(ds_rfe, "ds_rfe.csv", row.names = T)