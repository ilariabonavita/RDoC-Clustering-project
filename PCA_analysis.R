################################################################
# Rdoc mapping project 
# step 1: PCA on each variable
# subset used: scl, tpq, el variables, complete case only
################################################################

#--- libraries
library("FactoMineR")
library("factoextra")
library("corrplot")
library("nsprcomp")
library("elasticnet")

#--- load data
# choose between: file1: rdoc_scl_tpq.txt     590 individuals
#                 file2: rdoc_scl_tpq_el.txt  292 individuals
file.name <- './rdoc_scl_tpq_el.txt'
pheno <- read.table(file.name,h=T)
exp <- sub(".*rdoc_*(.*?) *.txt.*","\\1",file.name) 


pheno.data <- pheno[,-which(colnames(pheno) %in% c("nid","age","k_sex"))]
rownames(pheno.data) <- make.names(pheno$nid,unique=TRUE)


#=================  PCA  ============================================

pheno.pca <- PCA(pheno.data, ncp=ncol(pheno.data))
print(pheno.pca)

#--- get the eigenvalues
eig.val <- get_eigenvalue(pheno.pca)
eig.val

#--- plot istograms of explained variance per PC
pdf(paste0('./var_explained_',exp,'.pdf'))
fviz_eig(pheno.pca, addlabels = TRUE, ncp=20, ylim = c(0, 50))
dev.off()

#--- graph of variables
#------------------------------------
var <- get_pca_var(pheno.pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
# Coordinates of variables
head(var$coord, 4)

# The quality of representation of the variables on factor map is called cos2 (square cosine, squared coordinates) 
# You can visualize the cos2 of variables on all the dimensions using the corrplot package:

pdf(paste0('./var_cor_',exp,'.pdf'))
corrplot(var$cos2, is.corr=FALSE, title = "variables on factor map cos2",mar = c(0,0,1,0), tl.cex = 0.7)
dev.off()

# The correlation between a variable and a principal component (PC) is used as the coordinates of the variable on the PC. The representation of variables differs from the plot of the observations: The observations are represented by their projections, but the variables are represented by their correlations (Abdi and Williams 2010).


# Color by cos2 values: quality on the factor map
pdf(paste0('./var_pca_0.2cos2',exp ,'.pdf'))
fviz_pca_var(pheno.pca, col.var = "cos2",
             select.var = list(cos2 = 0.2),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             labelsize=3,
             geom.var = c('point','text'),
             title='Variables - PCA, 0.2 cos2 threshold',
             repel = TRUE # Avoid text overlapping
)
dev.off()

#### NOTES on the interpreatations of the two plots above:
# - Positively correlated variables are grouped together.
# - Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).
# - The distance between variables and the origin measures the quality of the variables on the factor map. Variables that are away from the origin are well represented on the factor map.
# - A high cos2 indicates a good representation of the variable on the principal component. 
# - A low cos2 indicates that the variable is not perfectly represented by the PCs. 


#--- Contributions of variables to PCs
#The contributions of variables in accounting for the variability in a given principal component are expressed in percentage.
# - Variables that are correlated with PC1 (i.e., Dim.1) and PC2 (i.e., Dim.2) are the most important in explaining the variability in the data set.
# - Variables that do not correlated with any PC or correlated with the last dimensions are variables with low contribution and might be removed to simplify the overall analysis.

# Use  corrplot() [corrplot package] to highlight the most contributing variables for each dimension:
pdf(paste0('./var_contribution_60',exp,'.pdf'))
corrplot(var$contrib[,1:60], is.corr=FALSE, title = "variables on factor map contribution",mar = c(0,0,1,0))
dev.off()

# barplot of the variable contribution to each PC
# The red dashed line on the graph above indicates the expected average contribution.

# Contributions of variables to PC1
pdf(paste0('./contribution_var_to_dim1',exp,'.pdf'))
fviz_contrib(pheno.pca, choice = "var", axes = 1, top = 50)
dev.off()
# Contributions of variables to PC2
pdf(paste0('./contribution_var_to_dim2_',exp,'.pdf'))
fviz_contrib(pheno.pca, choice = "var", axes = 2, top = 50)
dev.off()
# Contributions of variables to PC7
pdf(paste0('./contribution_var_to_dim7_',exp,'.pdf'))
fviz_contrib(pheno.pca, choice = "var", axes =8 , top = 50)
dev.off()


# The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
pdf(paste0('./var_pca_contr_top20_',exp,'.pdf'))
fviz_pca_var(pheno.pca, col.var = "contrib",
             select.var = list(contrib = 20),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),alpha.var = "contrib",
             labelsize=3,
             geom.var = c('point','text'),
             title='Variables - PCA, 0.2 contrib threshold',
             
)
dev.off()

#--- graph of individuals
#------------------------------------

# Plot of individuals in the first two PCs colored by sex, with concentration ellipses
pdf(paste0('./individuls_',exp,'.pdf'))
fviz_pca_ind(pheno.pca,
             axes=c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(pheno$k_sex), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()

#--- Add loadings of the variables on the PC (FactoMine doesn't provide them)
pheno.pca$loadings <- sweep(
  pheno.pca$var$coord, 2,
  sqrt(pheno.pca$eig[1:ncol(pheno.pca$var$coord),1]),
  FUN="/"
)

#================= sparse PCA  ============================================

#--- elasticnet
# k indicates the number of PCs to compute
# if sparse="varnum", in para= you have to choose how many non-zero values you want for each PC (i.e. you choose the amount of sparseness)
k <- 10
pheno.scpa <- spca(pheno.data, K=k, type = 'predictor',sparse = 'varnum', para = rep(20,k) )
print.spca(pheno.scpa)

# We can check which variables were selected by elasticnet (i.e. which variables are more important), for example by selecting the variables with at least one non-zero values in the first 2 dimensions
var.to.keep <- colnames(pheno.data)[rowSums(abs(pheno.scpa$loadings)) > 0]

# or we can check the loadings (higher loadings -> higher "importance")
pc1load <- abs(pheno.scpa$loadings[,1])
# save the variables with non-zero loadings 
pc1load.nonzero <- rownames(pheno.scpa$loadings)[pheno.scpa$loadings[order(pc1load, decreasing = T),1]!=0]


#--- compare the most important variables in PCA and sPCA
pc1load.top20 <- names(pheno.pca$loadings[order(abs(pheno.pca$loadings[,1]),decreasing = T)[1:20],1])

pc1load.nonzero[which(pc1load.nonzero %in% pc1load.top20)]

# or we can analysis the contribution of each variable...

oo <- fviz_pca_var(pheno.pca, select.var = list(contribut = 20))
var.to.keep[which(var.to.keep %in% rownames(oo$data))]

# normal PCA does not retrieve the el_ variables...something to think about


#================= rdoc domains and variables  =====================================
# We can check if variables who are assumed to map to the same rdoc domain are similarly represented in the PC space


varname <- list(name = tolower( strsplit(
"
scl0_19
scl0_05
scl0_32
tpq_ns3
tpq002
",'\n')[[1]]))
varname

pheno.pca$var$contrib[unlist(varname),1:5]

fviz_pca_var(pheno.pca, col.var = "cos2",axes = c(1,2),
             select.var = varname,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),alpha.var = "contrib",
             labelsize=4,
             title='Variables - PCA'
  )
