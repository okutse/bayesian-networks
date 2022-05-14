## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------
options(java.parameters = "-Xmx15g")

knitr::opts_chunk$set(warning=FALSE,
                      message=FALSE,
                      echo=T,
                      dpi=300,
                      error=FALSE)
 
# Use the table counter that the htmlTable() provides
options(table_counter = TRUE)

# function to install missing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.rstudio.com/')
  sapply(pkg, require, character.only = TRUE)
}

packages =c( "tidyverse","knitr", "kableExtra","skimr", "MatchIt", "RItools","optmatch", "ggplot2", "tufte", "tufterhandout", "plotly", "snowfall", "rstan")

ipak(packages)


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------
bytes <- file.size("README.Rmd")
words <- bytes/10
minutes <- words/200


## ----echo=FALSE------------------------------------------------------------------------------------------------------------------
#the image of the neural network algorithm trained via Bayesian
url<-"blnn_algorithm.jpg"


## ----include=FALSE---------------------------------------------------------------------------------------------------------------
#load the required libraries for the exploration and DEG
library(GEOquery)
library(limma)
library(clusterProfiler)
library(Biobase)


## ----cache=TRUE, include=TRUE----------------------------------------------------------------------------------------------------
#load the data from the GEO
gset <- getGEO("GSE13507", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6102", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("0000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXX22222222222222222222222222222222",
               "22222222222222222222222222222222222222222222222222",
               "22222222222222222222222222222222222222222222222222",
               "222222222222222222222222222222222XXXXXXXXXXXXXXXXX",
               "XXXXXX")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
exprs(gset) <- log2(exprs(gset))

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G2-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf) 
#we can then save the list of 1000 DEGs
saveRDS(tT[1:1000,], "deg.RDS")


## ----table1, include=TRUE, cache=TRUE--------------------------------------------------------------------------------------------
#top regulated genes
upregulated<-tT[which(tT$logFC>0),][1:10,]
#upregulated
upregulated<-subset(upregulated, select=c("Gene.ID", "Gene.symbol","logFC","AveExpr","adj.P.Val","B"))
##getting top 10 downregulated genes
downreg<-tT[which(tT$logFC<0),][1:10,]
downreg<-subset(downreg, select=c("Gene.ID","Gene.symbol","logFC","AveExpr","adj.P.Val","B"))
deg<-rbind(upregulated, downreg) 
rownames(deg)<-NULL
deg%>% kable(format = "html", 
             caption = "Top 10 up and down regulated genes in primary bladder cancer. The first 10 rows represent upregulated genes",
             col.names = c("Gene ID", "Gene Symbol", "logFC", "Average Expression", "Adjusted P-value", "B")) %>%
  kable_styling(full_width = FALSE, latex_options = c("HOLD_position", "stripped", "scale_down"), position = "left")


## ----include=FALSE---------------------------------------------------------------------------------------------------------------
# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)
#volcano plots for all the DEGs
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest


## ----mean-difference, fig.cap="A mean-difference plot showing the statistically significantly up and down regulated genes in primary bladder cancer relative to normal bladder cells.", fig.margin=TRUE, echo=FALSE, fig.pos="hold"----
# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=T, pch=20, cex=1, main = "")
abline(h=0)


## ----include=FALSE, cache=TRUE---------------------------------------------------------------------------------------------------
#Gene Enrichment analysis
library(org.Hs.eg.db)  #to get the GENE IDs/annotation
library(clusterProfiler)
library(enrichplot)
#load the saved gene sample (n=300) to convert gene symbols to ENTREZ IDs
gene_symbols<-read.csv("gene.csv", header=FALSE)
gene_symbols[]<-lapply(gene_symbols, as.character)
#head(gene_symbols)
gene_symbols<-gene_symbols[, 1]
#doing the conversion of IDs to symbols
converted_IDS <- mapIds(org.Hs.eg.db, gene_symbols, 'ENTREZID','SYMBOL')
##remove any missing cases
converted_IDS <- na.omit(converted_IDS)
##GO enrichment analysis
edo<-enrichGO(converted_IDS, org.Hs.eg.db, keyType = "ENTREZID", ont= "ALL", pvalueCutoff = 0.01)
edo<-setReadable(edo, OrgDb = org.Hs.eg.db)


## ----enrich, fig.cap="Sample Enrichment analysis results using the Gene Ontology (GO) enrichment analysis.", fig.pos="hold", fig.margin = TRUE----
#barplots of the similar results showing only 20 enrichment categories
barplot(edo, showCategory=10, cex.names=3)


## ----cache=TRUE, include=FALSE---------------------------------------------------------------------------------------
## data preprocessing for analysis
## load series and platform data from GEO
  gset <- getGEO('GSE13507', GSEMatrix =TRUE, getGPL=FALSE)
## get the expression profiles of the genes in the dataset
  x <- exprs(gset[[1]])
##extract the phenotype information from the dataset
  pData(gset[[1]])
## get the names of the columns in the phenotype dataset
  colnames(pData(gset[[1]]))
## transform the expression data to Z scores
  x <- t(scale(t(x)))
## removing extra rows from the x data take only columns with associated survival information
  x <- x[, 69:233]
## extract information of interest from the phenotype data (pdata) [Info not used in demo analysis]
  idx <- which(colnames(pData(gset[[1]])) %in%
                c('age:ch1', 'grade:ch1', 'survival month:ch1',
                  'SEX:ch1','progression:ch1',
                  'invasiveness:ch1', 'overall survival:ch1'))
  metadata <- data.frame(pData(gset[[1]])[,idx],
                        row.names = rownames(pData(gset[[1]])))
## get metadata for only samples in the dataset used in analysis
  metadata<-metadata[69:233,]
## filter the Z-scores expression data to match the samples in our pdata
  x <- x[, which(colnames(x) %in% rownames(metadata))]
## check that sample names match exactly between pdata and Z-scores
  all((colnames(x) == rownames(metadata)) == TRUE)
## # create a merged pdata and Z-scores object
  coxdata <- data.frame(metadata, t(x))
##load the downloaded separated dataset with the age associated with the samples
  age<-read.csv("age.csv",header=TRUE)
  coxdata$age.ch1<-age$age
## 
## ##convert variables to factors and numeric as necessary
 coxdata$age.ch1<-as.numeric(coxdata$age.ch1)
 coxdata$grade.ch1<-as.factor(coxdata$grade.ch1)
 coxdata$invasiveness.ch1<-as.factor(coxdata$invasiveness.ch1)
 coxdata$survival.month.ch1<-as.numeric(coxdata$survival.month.ch1)
 coxdata$SEX.ch1<-as.factor(coxdata$SEX.ch1)
 coxdata$progression.ch1<-as.factor(coxdata$progression.ch1)
 coxdata$overall.survival.ch1<-as.factor(coxdata$overall.survival.ch1)

## tidy column names
  colnames(coxdata)[1:7] <- c('Age', 'Grade', 'Invasiveness', 'Overall_survival',
                            'Progression', 'Sex', 'Survival_time')

## create a vector of the response variable for analysis
  overall_survival<-coxdata$Overall_survival
  surv<-ifelse(overall_survival=="death", 1, 0)
## coxdata$overall_survival<-as.numeric(ifelse(coxdata$overall_survival == 'survival', 0, 1))


## ----cache=TRUE, include=FALSE, eval=TRUE---------------------------------------------------------------------------------------
## install the package BLNN
  install.packages("remotes")
  remotes::install_github("BLNNdevs/BLNN")
  library(BLNN)
## load the probes to be used in analysis
  probes<-read.csv("probs.csv")
  coxdata<-coxdata[, 8:ncol(coxdata)]
## filter the data to only selected probes from the probes data file
  names<-names(coxdata) %in% probes$x
  coxdata<-coxdata[, names]
## add the survival data column to this subset data
  coxdata$survival<-surv
  trainx<-data.frame(coxdata[, -1376])
  trainy<-coxdata$survival
## save these clean files for use in analysis 
  saveRDS(trainx[, 1:150], "trainxx.RDS")
  saveRDS(trainx[, 1:150], "trainy.RDS")
## remove all the files saved in environment
rm(list = ls())
gc()


## --------------------------------------------------------------------------------------------------------------------------------
#Begin analysis to show on demo and post load the bayesian learning package
library(BLNN)
#load the saved train and test data files
trainy<-readRDS("trainy.RDS")
trainx<-readRDS("trainxx.RDS")


## ----include=FALSE---------------------------------------------------------------------------------------------------------------
#build the BLNN object
survObj<-BLNN_Build(ncov = 150, nout = 1, hlayer_size = 75,
                    actF = "tanh", outF = "sigmoid", costF = "crossEntropy",
                    hp.Err = 10, hp.W1 = rep(0.5, ncol(trainx)), hp.W2 = 0.5,
                    hp.B1 = 0.5, hp.B2 = 0.5)
class(survObj)


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## set the hyperparameters;  change this to evaluate
## network weights
n.par<-length(BLNN_GetWts(survObj))
## number of desired chains
chains<-2
## initials weight values
initials<-lapply(1:chains, function(i) rnorm(n.par, 0, 1/sqrt(n.par)))
## using nnet as baseline for comparisons
library(nnet)
nnetBasesline<-nnet(trainx, trainy, size= 75, maxit = 1000, MaxNWts = 15000)
nnetPredictions<-predict(nnetBasesline)
NET.abserror<-sum(abs(trainy-nnetPredictions)) #0.0736
NET.error<-nnetBasesline$value #7.97425e-05


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## variance for the moments
m1<-rep(1/2, n.par)
## training the model
trainx<-scale(trainx)
survNUTS<-BLNN_Train(survObj, x=trainx, y=trainy,
                      iter = 2000, thin = 10, warmup = 400,
                      init = initials, chains = chains,
                      parallel = TRUE, cores = 2,
                      algorithm = "NUTS", evidence = FALSE,
                      control = list(adapt_delta=0.99, momentum_mass=m1,
                                     stepsize= 1, gamma=0.05, to=100, useDA=TRUE,
                                     max_treedepth=20))


## ----include=FALSE, eval=TRUE---------------------------------------------------------------------------------------------------
## checks
Rhat<-mean(survNUTS$Rhat)
Rhat #1.00
## sample size
ESS<-mean(survNUTS$ess)
ESS #324


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## Bayesian learning with evidence used in re-estimating hyper-parameters
survNUTS.ev<-BLNN_Train(survObj, x=trainx, y=trainy,
                      iter = 5000, thin = 10, warmup = 400,
                      init = initials, chains = chains,
                      parallel = TRUE, cores = 2,
                      algorithm = "NUTS", evidence = TRUE,
                      control = list(adapt_delta=0.99, momentum_mass=m1,
                                     stepsize = 5, gamma = 0.05, to=100, useDA=TRUE,
                                     max_treedepth=20))


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## updating the parameters after the learning process
## update the no evidence neural network
survNUTS<-BLNN_Update(survObj, survNUTS)
## update the evidence neural network
survNUTS.ev<-BLNN_Update(survObj, survNUTS.ev)


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## making predictions using these models
## predictions using no evidence
survpred<-BLNN_Predict(survNUTS, trainx, trainy)
## predictions using bayesian learning
survpred.ev<-BLNN_Predict(survNUTS.ev, trainx, trainy)


## ----eval=TRUE------------------------------------------------------------------------------------------------------------------
## Model evaluations
## extract the errors in the classification
errors<-c(survpred$Errors$Total, survpred.ev$Errors$Total, nnetBasesline$value)

## print out the model evaluations
OutTab<-data.frame(errors)
rownames(OutTab)<-c("NUTS (no evidence)", "NUTS (with evidence)", "NNET")
saveRDS(OutTab, "outTab.RDS")
write.csv(OutTab, "OutTab.csv")
OutTab %>% kable(format = "html", caption = "Artificial Neural Network Model Comparisons") %>% kable_styling(full_width = FALSE, latex_options = c("HOLD_position", "stripped", "scale_down"), position = "left")


## ----bb, echo=FALSE--------------------------------------------------------------------------------------------------------------
options(knitr.kable.NA="")
OutTab<-read.csv("outTab.csv", header = F)
names(OutTab) <- NULL
OutTab %>% kable(format = "html", caption = "Artificial Neural Network Model Comparisons") %>% kable_styling(full_width = FALSE, latex_options = c("HOLD_position", "stripped", "scale_down"), position = "left")

