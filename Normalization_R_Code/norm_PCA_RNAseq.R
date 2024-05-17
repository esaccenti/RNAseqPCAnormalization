# https://www.ebi.ac.uk/training/online/courses/ensembl-browsing-genomes/navigating-ensembl/investigating-a-gene/
# https://www.bioconductor.org/help/course-materials/2015/UseBioconductorFeb2015/A01.5_Annotation.html
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/05_Annotation_and_Visualisation.nb.html

.libPaths("M:/My Documents/R/win-library/4.2") # .libPaths("M:/My Documents/R/win-library/3.6") 

setwd('C:/Users/linge006/OneDrive - Wageningen University & Research/Documenten/WUR_HOME/Diagonal/try_ARCHS4/'); 
#setwd('C:/Users/Gebruiker/OneDrive - Wageningen University & Research/Documenten/WUR_HOME/Diagonal/try_ARCHS4/')

#BiocManager::install("EDASeq") OR "DelayedArray" OR "clusterProfiler"

#library("DESeq2")  # for differential expression
library("ggplot2") # for plotting.
library("org.Hs.eg.db") # for annotation of Human sapiens (use Mm for mouse), etc...
#library(AnnotationDbi)
library(clusterProfiler) # https://youtu.be/JPwdqdo_tRg; Other enrichment packages: topGO, GOfuncR (https://rpubs.com/Akhil_Velluva/GOfuncR), GOrilla
library(edgeR)
library(gridExtra)
library(scone)
#library(qsmooth)
#library(dplyr)
library(cqn) # Does not work in R 3.6
library(EDASeq) # For GC content based normalization
#source('../pqn_FN.R') # load PQN (Probabilistic Quotient Normalization) function
library(KODAMA)
#source('../Maza2016/additional_file2.R') # load MRN (Median Ratio Normalization) function
source('norm_fun.R')
library(pheatmap)
library(enrichplot)

packageVersion("ggplot2")
packageVersion("org.Hs.eg.db")
packageVersion("clusterProfiler")
packageVersion("edgeR")
packageVersion("gridExtra")
packageVersion("scone")
packageVersion("cqn")
packageVersion("EDAseq")
packageVersion("KODAMA")
packageVersion("pheatmap")
packageVersion("enrichplot")

packageVersion("factoextra")
packageVersion("VennDiagram")

packageVersion("ggrepel")
packageVersion("ggvenn")
packageVersion("vegan")

# Load human Genome library details, viz. GeneLength, gene, GeneID and C.G content
genome_lib <- read.csv("details_human_genome_CDS.csv"); head(genome_lib)

#dataset <- "simulated"; xMonthyyyy <- "24April2024"; tumor_data <- rownames(read.table("tumors_full_Raw.csv",sep=",",row.names=1))
#dataset <- "correl_sim"
#dataset <- "tumors" # Gene names; human tissue
#dataset <- 'GSE30017' # Gene names; human+mouse tissue
#dataset <- 'GSE208218' # ENSG numbers/ensemble ID; human tissue
#dataset <- 'GSE218399' # ENSG numbers/ensemble ID; human tissue
dataset <- 'GSE216274' # Gene names
#dataset <- "LungCells" # Gene names; human tissue

if (dataset=="correl_sim") {rawdata <- read.delim("SimulatedCounts_Correlated_26April2024.csv", sep=",", header=F); rownames(rawdata) <- rownames(read.delim("tumors_full_Raw.csv", sep=",",row.names=1, header=F))}
if (dataset=="simulated") {rawdata <- read.delim(paste0("SimulatedCounts_Uncorrelated_",xMonthyyyy,".csv"), sep=",", header=F); rownames(rawdata) <- rownames(read.delim("tumors_full_Raw.csv", sep=",", row.names=1, header=F))}
if (dataset=="tumors") rawdata <- read.delim("../Tuch2010/TumorSeqs.txt", check.names=FALSE, stringsAsFactors=FALSE)
if (dataset=="GSE30017") rawdata <- read.delim("GSE30017_expression_matrix.tsv", row.names=1, sep="\t")
if (dataset=="GSE208218") rawdata <- read.delim("GSE208218_FLUGAZA_Counts.txt", row.names=1, sep=" ")
if (dataset=="GSE218399") rawdata <- read.delim("GSE218399_read_count.txt", row.names=1)
if (dataset=="GSE216274") rawdata <- read.delim("GSE216274_BE2C_rawcounts.csv", row.names=1, header=T, sep=",")
if (dataset=="LungCells") rawdata <- read.delim("../Lung_RNAseq/lung_counts.txt", row.names=1, header=T, sep="\t")

if (dataset!="tumors") rawdata <- rawdata[rowSums(rawdata)!=0,]
str(rawdata); head(rawdata[,1:4]); 
dim(rawdata); typeof(rawdata); class(rawdata)

if (dataset=="GSE218399") rawdata <- rawdata[-which(apply(rawdata, 1, function(x) length(which(x==0))) > 2),] # rm rows/genes with multiple zeros
#if (dataset=="LungCells") rawdata <- rawdata[,-grep("PC", colnames(rawdata))]

##### ##### ##### ##### ##### ##### ##### ##### Start of data filtering section ##### ##### ##### ##### ##### ##### ##### ##### ##### 

if (dataset=="tumors"){
  if (dataset=="tumors") y <- list(counts=rawdata[,4:9], genes=rawdata[,1:3])
  #if (dataset=="GSE216274") y <- list(counts=rawdata[,1:6], genes=rownames(rawdata))
  
  ### Data filtering; only retain genes already in org.Hs.egREFSEQ
  
  # Only retain genes available thru org.Hs.egREFSEQ
  idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
  y <- sapply(y, function(x) x[idfound,])
  
  # Extract Entrez Gene IDs
  egREFSEQ <- toTable(org.Hs.egREFSEQ); head(egREFSEQ)
  
  # Add Entrez Gene IDs to the annotation:
  m <- match(y$genes$RefSeqID, egREFSEQ$accession)
  y$genes$EntrezGene <- egREFSEQ$gene_id[m]
  
  # Extract gene symbols:
  egSYMBOL <- toTable(org.Hs.egSYMBOL); head(egSYMBOL) # gene symbols
  
  # Add gene symbols using the selected Entrez Gene IDs 
  m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
  y$genes$Symbol <- egSYMBOL$symbol[m]
  head(y$genes)
  
  ### Filtering
  o <- order(rowSums(y$counts), decreasing=TRUE)
  y <- sapply(y, function(x) x[o,]) # order genes to number of counts per gene
  d <- duplicated(y$genes$Symbol)
  y <- sapply(y, function(x) x[!d,]) # remove duplicated genes
  nrow(y$genes)
  
  # Use Entrez Gene IDs as row names:
  rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
  #y$genes$EntrezGene <- NULL; y$genes$Symbol <- NULL
  
  # DGEList creates 1) a table of counts (rows=features, columns=samples), 2) group indicator for each column, 3) library size (optional) and 
  # 4) a table of feature annotation (optional).
  y <- DGEList(counts=y$counts, genes=y$genes, samples=NULL)
} # End if statement


if (dataset=="GSE216274"){
  
  y <- list(counts=rawdata[,1:ncol(rawdata)], genes=rownames(rawdata))
  print(dim(y$counts)) # 1st print
  # Extract gene symbols from library:
  egSYMBOL <- toTable(org.Hs.egSYMBOL); head(egSYMBOL) # gene symbols
  head(egSYMBOL) # gene symbols
  
  ### Data filtering; only retain genes already in egSYMBOL
  idfound <- y$genes %in% egSYMBOL$symbol
  
  # Only retain genes available thru egSYMBOL
  y$counts <- y$counts[idfound,]; y$genes <- y$genes[idfound]
  print(dim(y$counts)) # 2nd print
  # Add gene symbols using the selected Entrez Gene IDs 
  m <- match(y$genes, egSYMBOL$symbol)
  y$genes <- data.frame(Symbol=y$genes, EntrezGene=egSYMBOL$gene_id[m])
  head(y$genes)
  
  
  ### Filtering
  o <- order(rowSums(y$counts), decreasing=TRUE)
  y <- sapply(y, function(x) x[o,]) # order genes to number of counts per gene
  d <- duplicated(y$genes$Symbol)
  y <- sapply(y, function(x) x[!d,]) # remove duplicated genes
  nrow(y$genes)
  
  # Use Entrez Gene IDs as row names:
  rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
  #y$genes$EntrezGene <- NULL
  print(dim(y$counts)) # 3nd print
  # DGEList creates 1) a table of counts (rows=features, columns=samples), 2) group indicator for each column, 3) library size (optional) and 
  # 4) a table of feature annotation (optional).
  y <- DGEList(counts=y$counts, genes=y$genes, samples=NULL)
  
  
  ### Discard genes with many zero's
  
  # Extract number of zero's and non-zero's per row or gene
  zero_or_not <- unlist(apply(y$counts[,], 1, function(x) table(x==0)))
  
  # Retain genes with zero's and the number of zero's
  zeros <- zero_or_not[grep(".TRUE",names(zero_or_not))] 
  
  # Retain genes with more than 8 zero's
  many_zeros <- zeros[zeros > 8]
  
  # Remove '.TRUE' from gene names
  rows_with_many_zeros <- as.integer(sub(".TRUE","",names(many_zeros) ))
  
  # Retain records for gene's with more than 8 zero's
  y$counts <- y$counts[-c(which(rownames(y$counts) %in% rows_with_many_zeros)),]
  y$genes <- y$genes[-c(which(rownames(y$genes) %in% rows_with_many_zeros)),] # check line of code
  print(dim(y$counts)) # 4th print
} # End if statement


if (dataset!="tumors" & dataset!="GSE216274") y <- DGEList(counts=rawdata, genes=rownames(rawdata), samples=NULL)

### description of the experiments 
filenameSample <- paste0("sample_description_",dataset,".csv")
if(dataset!="simulated" & dataset!="correl_sim") colData <- read.csv(filenameSample, row.names=1, stringsAsFac=TRUE)
if (!dataset%in%c('GSE30017','simulated')) colData$condition <- as.factor(colData$condition)
str(colData)

#if (dataset=="LungCells") colData <- colData[colData$condition!="PC",]

# count matrix and column data need to be  consistent (same names in same order)
colnames(y$counts)==rownames(colData)  #check that they are in the same order.Everything should be "TRUE"

if (dataset=="tumors") model <- "tissue+condition"
if (dataset=="GSE30017") model <- "tissue"
if (dataset=="GSE208218") model <- "Sample+Patient_no+Sample_ID"
if (dataset=="GSE218399") model <- "condition"
if (dataset=="GSE216274") {rownames(colData) <- colnames(y$counts); model <- "condition"}
if (dataset=="LungCells") model <- "condition"

########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## 

if (dataset=="tumors" | dataset=="simulated" | dataset=="correl_sim" | dataset=="GSE216274"){
  y$counts <- y$counts[rownames(y$counts) %in% genome_lib$GeneID,] # Retain genes available in library
  y$genes <- y$genes[rownames(y$genes) %in% genome_lib$GeneID,] # Retain genes available in library; doesn't work for simulated!!!???
  
  if (dataset=="tumors" | dataset=="GSE216274") {
    #genome_lib <- genome_lib[genome_lib$GeneID %in% rownames(y$genes),] # Retain genes in dataset; doesn't work for simulated!!!???
    #genome_lib <- genome_lib[match(rownames(y$genes),genome_lib$GeneID),] # Order genes in line with dataset
    genome_lib <- genome_lib[genome_lib$GeneID %in% rownames(y$counts),] # Retain genes in dataset; doesn't work for simulated!!!???
    genome_lib <- genome_lib[match(rownames(y$counts),genome_lib$GeneID),] # Order genes in line with dataset
  } # End if()
  
  if (dataset=="simulated" | dataset=="correl_sim") {
    y$genes[rownames(y$genes) %in% genome_lib$GeneID]
    y$genes[rownames(y$genes) %in% as.character(genome_lib$GeneID)]
    genome_lib <- genome_lib[genome_lib$GeneID %in% y$genes,] # Retain genes in dataset; doesn't work for simulated!!!???
    genome_lib <- genome_lib[match(y$genes,genome_lib$GeneID),] # Order genes in line with dataset
  } # End if()
  
} # End if()

########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## 

if (dataset=="GSE30017"){
  y$counts <- y$counts[rownames(y$counts) %in% genome_lib$gene,] # Retain genes available in library
  y$genes <- y$genes[rownames(y$genes) %in% genome_lib$gene,] # Retain genes available in library
  
  genome_lib <- genome_lib[genome_lib$gene %in% rownames(y$counts),] # Retain genes in dataset
  genome_lib <- genome_lib[match(rownames(y$counts),genome_lib$gene),] # Order genes in line with dataset
  
  colData$tissue <- as.character(colData$tissue)
  colData$tissue[grepl(" ",colData$tissue)] <- unlist(strsplit(colData$tissue[grepl(" ",colData$tissue)], " "))[c(1,3)]
  
  colnames(y$counts) <- colData$tissue
} # End if()


########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## <'))>< ########## 

if (dataset %in% c("GSE208218","GSE218399","LungCells")){
  ENSG_IDs <- read.table("HumanGeneIDs_ENSG.txt", header=T, sep="\t")[c("From","To")]
  
  # Remove non-unique records (ENSG that appear twice)
  retain <- which(table(ENSG_IDs$To)==1); retain <- names(retain)
  ENSG_IDs <- ENSG_IDs[ENSG_IDs$To %in% retain,]
  
  # Retain records with ENSG that appear in both data objects
  y$counts <- y$counts[rownames(y$counts) %in% ENSG_IDs$To,]
  y$genes <- y$genes[rownames(y$genes) %in% ENSG_IDs$To,]
  ENSG_IDs <- ENSG_IDs[ENSG_IDs$To %in% rownames(y$counts),]
  
  ENSG_IDs <- ENSG_IDs[order(ENSG_IDs$To),]
  y$counts <- y$counts[order(rownames(y$counts)),]
  
  rownames(y$counts) <- ENSG_IDs$From
  
  y$counts <- y$counts[rownames(y$counts) %in% genome_lib$GeneID,] # Retain genes available in library
  y$counts <- y$counts[order(as.integer(rownames(y$counts))),]# Order genes
  y$genes <- rownames(y$counts) # Retain genes available in library
  
  genome_lib <- genome_lib[genome_lib$GeneID %in% rownames(y$counts),] # Retain genes in dataset
  genome_lib <- genome_lib[match(rownames(y$counts),genome_lib$GeneID),] # Order genes in line with dataset
  ENSG_IDs <- ENSG_IDs[order(ENSG_IDs$From),]
  genome_lib <- data.frame(genome_lib,ENSG=ENSG_IDs$To)
  
  ### Remove duplicate genes on GeneID
  boolean_seq <- rownames(y$counts) %in% names(which(table(rownames(y$counts))==1))
  y$counts <- y$counts[boolean_seq,]
  y$genes <- length(y$genes[boolean_seq])
  ENSG_IDs <- ENSG_IDs[boolean_seq,]
  genome_lib <- genome_lib[boolean_seq,]
} # End if()

##### ##### ##### ##### ##### ##### ##### ##### End of data filtering section ##### ##### ##### ##### ##### ##### ##### ##### 

#Patient <- factor(c(8,8,33,33,51,51)); Tissue <- factor(c("N","T","N","T","N","T"))
#data.frame(Sample=colnames(y), Patient,Tissue); design <- model.matrix(~Patient+Tissue-1)

#### Create DESeq2 object. Note: DESEQ_FN() normalization is the same!
#dds <- DESeqDataSetFromMatrix(countData = y$counts,
 #                             colData = colData,
  #                            design = as.formula(paste("~",model)))

### Normalizing+fitting the data
#dds <- DESeq(dds)

#vsd <- vst(dds)
#class(vsd)
#plotPCA(vsd, "dex")
#pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
 # geom_point(size =3) +
 # xlab(paste0("PC1: ", percentVar[1], "% variance")) +
 # ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 # coord_fixed()

# Create DGE object; normFactors; ...
#y <- DGEList(counts=counts,genes=data.frame(Length=GeneLength)); y <- calcNormFactors(y); estimateSizeFactorsForMatrix(y$counts)

genes_6091 <- FALSE
if (genes_6091 == TRUE){
  common_genes <- rownames(read.csv(paste0(dataset,"_raw.csv")))
  y$counts <- y$counts[rownames(y$counts) %in% common_genes,]
  genome_lib <- genome_lib[genome_lib$GeneID %in% common_genes,]
} # 
###

# TMM normalization is applied to this dataset to account for compositional difference between the libraries.
#y_none <- calcNormFactors(y,"none")
y_TMM <- calcNormFactors(y,"TMM")
#y_TMMwsp <- calcNormFactors(y,"TMMwsp") # Takes long for simulated...?!
#y_RLE <- calcNormFactors(y,"RLE") # Takes long for simulated...?!
#y_UpQ <- calcNormFactors(y,"upperquartile")
#y_none <- y_none$counts


y_raw_norm <- list(Raw=y$counts,
                   TotCounts=SUM_FN(y$counts),
                   UpQuart=UQ_FN(y$counts),
                   #TMM=TMM_FN(y$counts),
                   FullQuant=FQ_FN(y$counts),
                   TMM=cpm(y_TMM), # edgeR TMM
                   #edgeR_TMMwsp=cpm(y_TMMwsp)
                   #MRN=mrnFactors(y$counts,rep(1,ncol(y$counts))),
                   #DESeq2_norm=counts(dds, normalized=T), # RLE normalization
                   RelLogExpr=DESEQ_FN(y$counts), 
                   CenterLogRatio=CLR_FN(y$counts),
                   RPKM=rpkm_2(y$counts,genome_lib$GeneLength), # rpkm(y$counts,genome_lib$GeneLength) also works, function from edgeR
                   TPM=tpm_2(y$counts,genome_lib$GeneLength),
                   GC_w=edaseq(y$counts, genome_lib$C.G, method="With"),
                   GC_wb=edaseq(y$counts, genome_lib$C.G),
                   #QSM=qsmooth(y$counts, group_factor=1:6)@qsmoothData,
                   #PQN=pqn(y$counts), 
                   #PQN_EDO=pqn_Edo(y$counts),
                   PQN=scaling(normalization(1+y$counts, method="pqn")$newXtrain)$newXtrain, # " KODAMA: https://search.r-project.org/CRAN/refmans/KODAMA/html/normalization.html
                   CQN=cqn(y$counts,genome_lib$C.G,genome_lib$GeneLength,verbose=T)) # Takes long...!!!

#y_raw_norm$MRN <- sapply(1:ncol(y_raw_norm$Raw), function(x) y_raw_norm$Raw[,x]/y_raw_norm$MRN$normFactors[x])
#colnames(y_raw_norm$MRN) <- colnames(y_raw_norm$Raw)
#y_raw_norm$PQN <- scaling(y_raw_norm$PQN)$newXtrain # 
y_raw_norm$CQN <- y_raw_norm$CQN$y + y_raw_norm$CQN$offset


sapply(1:length(y_raw_norm), function(x) write.table(round(y_raw_norm[[x]],1), paste0(dataset,"_full_",names(y_raw_norm)[x],".csv"), col.names=F, sep=','))

setwd("C:/Users/linge006/OneDrive - Wageningen University & Research/Documenten/WUR_HOME/Diagonal")


### <'))>< ### Fish ### <'))>< ### Fish ### <'))>< ### Fish ### <'))>< ### Fish ### <'))>< ### Fish ### <'))>< ### Fish ### <'))>< ### Fish

# Remove zero variance genes

sapply(y_raw_norm, function(x) print(dim(x)))

# Remove zero variance genes per normalization method
discard <- lapply(y_raw_norm, function(x){
  discard <- apply(x, 1, sd)
  names(which(discard==0))
  #x <- x[!rownames(x)%in%names(which(discard==0)),]
} )

discard <- unique(unlist(discard))

for (i in 1:length(y_raw_norm)) y_raw_norm[[i]] <- y_raw_norm[[i]][!rownames(y_raw_norm[[i]]) %in% discard,]

#

PCA_ggplot <- function(x, title, no_tr=no_t, s_impose_x=FALSE, s_impose_y=FALSE, leg_pos='none'){
  xx <- prcomp(t(as.matrix(x)), scale=T)
  
  ## more elegant plot using the ggplot library.
  z <- summary(xx)
  print(z$importance[3,]) #print(z$importance[3,which(z$importance[3,] > 0.75)])
  if (title=='Raw') retain_data <- FALSE else retain_data <- TRUE
  PCA_var <- z$importance[3,]; PCA_var[-1] <- diff(as.numeric(PCA_var)); 
  PCA_var <- data.frame(t(PCA_var[1:3]), CumSum=sum(PCA_var[1:3]), row.names=paste0("\ttextbf{",title,"}"))
  write.table(round(100*PCA_var,1), paste0("PCA_var_",dataset,".txt"), sep=" & ", col.names=F, append=retain_data, row.names=T)
  if (s_impose_x==TRUE) xx$x[,1] <- -1*as.numeric(xx$x[,1]); if (s_impose_y==TRUE) xx$x[,2] <- -1*as.numeric(xx$x[,2])
  df <- data.frame(PC1=as.numeric(xx$x[,1]), PC2=as.numeric(xx$x[,2]), Treatment=rownames(xx$x))
  df$PC1 <- df$PC1/max(abs(df$PC1)); df$PC2 <- df$PC2/max(abs(df$PC2))
  #print(df)
  # data object adjustment
  if (dataset == 'tumors'){
    n_c <- nchar(df$Treatment)
    df <- data.frame(df, Patient=substr(df$Treatment,1,n_c-1), Tissue=substr(df$Treatment,n_c,n_c))
    df$Tissue <- sub('T','Tumor',df$Tissue); df$Tissue <- sub('N','Normal',df$Tissue)
    df$Patient <- as.character(rep(rank(unique(df$Patient)),each=2))
    geom_p <- geom_point(aes(colour=Patient, shape=Tissue), size=2)
  } # End tumors if()
  if (dataset=='GSE216274'){
    df$Treatment <- matrix(unlist(strsplit(df$Treatment,"_")),ncol=2,byrow=T)[,1]
    geom_p <- geom_point(aes(colour=Treatment, shape=Treatment), size=2) 
    #geom_p <- geom_p + scale_shape_manual(values=c(15:17,19)) 
  } # End GSE216274 if()
  if (dataset=='GSE218399'){
    df$Treatment[grep('shNC',df$Treatment)] <- 'Control'; df$Treatment[grep('shMIB',df$Treatment)] <- 'Knockout'
    geom_p <- geom_point(aes(colour=Treatment, shape=Treatment), size=2) 
  } # End GSE218399 if()
  if (dataset == 'LungCells'){
    df <- data.frame(df, Treat=as.factor(colData$treatment), Conc=as.factor(colData$concentration))
    geom_p <- geom_point(aes(colour=Treat, shape=Conc), size=2)
  } # End tumors if()
  
  myplot <- ggplot(df,aes(PC1,PC2), color='black') + geom_p + 
    scale_color_manual(values=c('magenta','limegreen','dodgerblue','black','orangered','gold')) + scale_shape_manual(values=c(16,15,17:18)) +
    xlab(paste0("PC1 (",format(100*z$importance[2,1],digits=3),"%)")) +
    ylab(paste0("PC2 (",format(100*z$importance[2,2],digits=3),"%)")) +
    theme_bw() + xlim(-1.25,1.25) + ylim(-1.15,1.15) +
    theme(axis.text.x=element_text(size=11), legend.position=leg_pos, 
          legend.text=element_text(size=11), legend.spacing.x=unit(0,'cm'), legend.spacing.y=unit(0,'cm'), legend.margin=margin(0.1,0,0,0,'cm'), # trbl 
          legend.title=element_text(size=11)) + ggtitle(title)
  ###
  
  #d_matrix <- matrix(ncol=length(PC1),nrow=length(PC2), dimnames=list(rownames(xx$x),rownames(xx$x))) # Initialize empty matrix for PC1+2 distances
  PC_numbers <- 1:ncol(xx$x) # sequence for number of samples
  
  df$PC1 <- df$PC1*z$importance["Proportion of Variance",1]; df$PC2 <- df$PC2*z$importance["Proportion of Variance",2]
  
  d_ij <- function(PC1,PC2,span){ # Compute distances
    dist_matrix <- matrix(ncol=length(PC1),nrow=length(PC2)) # Initialize matrix for distances
    
    while (length(span) > 1){ # Compute distances while distances remained to be calculated
      
      for (p in (min(span)+1):max(span)){ # Compute (remaining) distances point in for loop
        dist_matrix[p,min(span)] <- sqrt((PC1[min(span)] - PC1[p])**2 + (PC2[min(span)] - PC2[p])**2)
      } # End for loop p
      
      span <- span[-1] # Remove used span element while retaining elements for distances to be calculated yet
    } # End while loop
    
    return(dist_matrix/max(dist_matrix[lower.tri(dist_matrix)])) # Return standardized distances
  } # End d_ij()
  
  d_matrix <- d_ij(df$PC1,df$PC2,PC_numbers); rownames(d_matrix) <- df[,ncol(df)]
  
  centroids <- data.frame(matrix(ncol=2,nrow=no_tr,dimnames=list(NULL,c('PC1','PC2')))) # Initialize centroid distance matrix
  for (z in 1:no_tr){ # Compute centroid coordinates
    span <- (1+no_tr*(z-1)):(no_tr*z)
    centroids[z,] <- c(mean(df$PC1[span]), mean(df$PC2[span]))
  } # End for loop z
  c_d_matrix <- d_ij(centroids$PC1,centroids$PC2,1:no_tr)
  
  return(list(myplot,d_matrix,c_d_matrix))
} # End PCA_ggplot()


if (dataset=="GSE218399"){
  leg_xy <- c(0.37,0.45)
} else if (dataset=="GSE216274") {
  leg_xy <- c(0.85,0.625)
} else if (dataset=="LungCells") {
  leg_xy <- c(-0.2,0.5)
} else { # tumors
  leg_xy <- c(0.55,0.5)
} # End if statements


no_t <- length(levels(colData$condition)); reps <- nrow(colData) / no_t
if (dataset=="tumors") no_t <- length(levels(colData$tissue)); reps <- nrow(colData) / no_t

## more elegant plot using the ggplot library.

method_check <- FALSE

if (method_check==TRUE){
  p1 <- PCA_ggplot(y_raw_norm$TMM,"TMM", leg_pos=leg_xy)
  p11 <- PCA_ggplot(y_raw_norm$edgeR_TMMwsp,"edgeR TMMwsp", multip_f1[1]*xy_lim$xlim, multip_f1[2]*xy_lim$ylim)
  p7 <- PCA_ggplot(y_raw_norm$DESEQ_FN,"DESeq_FN")
  p8 <- PCA_ggplot(y_raw_norm$MRN,"MRN")
  p9 <- PCA_ggplot(RLE_Maza(1+y_raw_norm$Raw),"RLE Maza")
  p10 <- PCA_ggplot(MRN_Maza(1+y_raw_norm$Raw),"MRN Maza")
  
  pdf(paste0('methods_norm_PCA_',dataset,'_RNA02.pdf'), height=8, width=12)
  grid.arrange(p1[[1]], p10[[1]], p11[[1]], p7[[1]], p9[[1]], p8[[1]], 
               ncol=3, bottom="")
  dev.off()
} # End method_check

#png(paste0('scaled_RLE_PCA_',dataset,'_RNA.png'), height=4, width=4, res=600, units='in')
#pdf(paste0('scaled_RLE_PCA_',dataset,'_RNA.pdf'), height=4, width=4)
#PCA_ggplot(DESEQ_FN(rawdata[which(apply(rawdata, 1, min) >= 50),]),"", leg_pos=c(0.125,0.575))
#dev.off()

leg_xy <- c(0.5,0.5)
p_leg <- PCA_ggplot(y_raw_norm$Raw,"Raw", s_impose_x=F, s_impose_y=F, leg_pos=leg_xy)
p01 <- PCA_ggplot(y_raw_norm$Raw,"Raw", s_impose_x=F, s_impose_y=F)
p02 <- PCA_ggplot(y_raw_norm$TotCounts,"TC", s_impose_y=F) 
p03 <- PCA_ggplot(y_raw_norm$UpQuart,"UQ", s_impose_y=F)
p04 <- PCA_ggplot(y_raw_norm$FullQuant,"FQ", s_impose_y=F)
p05 <- PCA_ggplot(y_raw_norm$TMM,"TMM") # 
p06 <- PCA_ggplot(y_raw_norm$RelLogExpr,"RLE") ## Same...
p07 <- PCA_ggplot(y_raw_norm$CenterLogRatio,"CLR") ## ...
p08 <- PCA_ggplot(y_raw_norm$RPKM,"RPKM", s_impose_y=F) ##  * 
p09 <- PCA_ggplot(y_raw_norm$TPM,"TPM", s_impose_y=F) ## 
p10 <- PCA_ggplot(y_raw_norm$GC_w,"GCw", s_impose_x=F, s_impose_y=T) ## *
p11 <- PCA_ggplot(y_raw_norm$GC_wb,"GCwb", s_impose_x=F, s_impose_y=F) ## 
p12 <- PCA_ggplot(y_raw_norm$PQN,"PQN", s_impose_x=F, s_impose_y=F) # 
p13 <- PCA_ggplot(y_raw_norm$CQN,"CQN", s_impose_y=F) ## 
legend <- cowplot::get_legend(p_leg[[1]])


#png(paste0('norm_PCA_',dataset,'_RNA.png'), height=12, width=8, units="in", res=600)
pdf(paste0('scaled_norm_PCA_',dataset,'_RNA.pdf'), height=12, width=8)
par(oma=c(0.15,0,0,0))
grid.arrange(p01[[1]], p02[[1]], p03[[1]], 
             p04[[1]], p05[[1]], p06[[1]],
             p07[[1]], p08[[1]], p09[[1]], 
             p10[[1]], p11[[1]], p12[[1]], p13[[1]], 
             legend, ncol=3, bottom="")
dev.off()

# Clustering

d_M <- list(p01[[2]],p02[[2]],p03[[2]],p04[[2]],p05[[2]],p06[[2]],p07[[2]],p08[[2]],
            p09[[2]],p10[[2]],p11[[2]],p12[[2]],p13[[2]]) # Distances
c_d_M <- list(p01[[3]],p02[[3]],p03[[3]],p04[[3]],p05[[3]],p06[[3]],p07[[3]],p08[[3]],
              p09[[3]],p10[[3]],p11[[3]],p12[[3]],p13[[3]]) # Centroid distances


# Silhouette plots
library(cluster); packageVersion("cluster")
library(factoextra); packageVersion("factoextra")

# Design cluster per dataset
if (dataset=="tumors") cluster_design <- rep(1:no_t,reps)
if (dataset=="GSE218399") cluster_design <- rep(1:no_t,each=reps)
if (dataset=="GSE216274") cluster_design <- rep(1:no_t,each=reps)

# Perform silhouette analysis on the raw dataset
sil <- silhouette(cluster_design,d_M[[1]]) # Compute Silhouette scores
print(sil)
gg_main <- paste("Raw: ",round(mean(sil[,3]),2)) # Generate raw data Silhouette plot title

if (dataset=="GSE216274") custom_leg <- scale_fill_discrete(name="Treatment", 
                                  breaks=levels(sil[,"cluster"]), 
                                  labels=c("DMSO", "RA", "PB", "RA+PB")) 
if (dataset=="GSE218399") custom_leg <- scale_fill_discrete(name="Condition", 
                                  breaks=levels(sil[,"cluster"]), 
                                  labels=levels(colData$condition)) 
if (dataset=="tumors") custom_leg <- scale_fill_discrete(name="Tissue", 
                                  breaks=levels(sil[,"cluster"]), 
                                  labels=levels(colData$tissue)) 

s_plot <- list(fviz_silhouette(sil)+ ggtitle(gg_main) + custom_leg + ylab("") + 
                 scale_fill_hue(labels = custom_leg$labels) +
                 guides(color='none', fill=guide_legend(custom_leg$name)) + theme(legend.position="none")) # axis.text.x = element_text(size=15),

s_leg <- list(fviz_silhouette(sil)+ ggtitle(gg_main) + custom_leg + ylab("") + 
                 scale_fill_hue(labels = custom_leg$labels) +
                 guides(color='none', fill=guide_legend(custom_leg$name)) + theme(legend.key.size = unit(0.45, 'cm'), legend.position=c(0.5,0.5),
                                                                                  legend.text=element_text(size=10), legend.title=element_text(size=12),
                                                                                  legend.spacing.x=unit(0.1,'cm'), legend.spacing.y=unit(0.05,'cm'), 
                                                                                  legend.margin=margin(0.05,0.05,0.05,0.05,'cm') # trbl
                 )) #

# Perform silhouette analysis on the normalized datasets
norm_names <- c("Raw", "Total-Count", "Upper-Quartile", "Full-Quantile", "TMM", "Relative Log Expr", "Center Log Ratio",
                "RPKM", "TPM", "Within-GC", "Within+Betw-GC", "Probabilistic Quotient", "Conditional Quantile")
norm_names <- c("Raw", "TC", "UQ", "FQ", "TMM", "RLE", "CLR", "RPKM", "TPM", "GCw", "GCwb", "PQN", "CQN")

for (s in 2:13){
  sil <- silhouette(cluster_design,d_M[[s]]) # Compute Silhouette scores
  #if (s==2|s==3) print(sil)
  gg_add <- fviz_silhouette(sil) + ggtitle(paste0(norm_names[s],": ",round(mean(sil[,3]),2))) # Generate Silhouette plot title given the normalization method
  if (s!=7) gg_add <- gg_add + ylab("") else gg_add <- gg_add + ylab("Silhouette width") + theme(axis.title.y=element_text(size=14))
  if (s==14) gg_add <- gg_add + xlab("PC1 + PC2 points") + theme(axis.title.x=element_text(size=14))
  s_plot[[length(s_plot)+1]]  <- gg_add + theme(legend.position='none') # Prepare Silhouette plot for all normalized data
} # End of for loop s

legend <- cowplot::get_legend(s_leg[[1]])

#names(y_raw_norm)
pdf(paste0("scaled_silhouette_plots_",dataset,".pdf"), width=8.5,height=11)
grid.arrange(s_plot[[1]],s_plot[[2]],s_plot[[3]],s_plot[[4]],s_plot[[5]],s_plot[[6]],
             s_plot[[7]],s_plot[[8]],s_plot[[9]],s_plot[[10]],s_plot[[11]],s_plot[[12]],
             s_plot[[13]],
             legend, ncol=3)
dev.off()

# cluster_distances function
cluster_distances <- function(d_data,no_treats=no_t,no_reps=reps){
  
  distances <- mean(d_data[lower.tri(d_data)]) # Compute mean distance of non-zero matrix elements
  
  for (i in 1:no_treats){
    rows_cols <- (i*no_reps-no_reps+1):(i*no_reps) # Select matrix rows+columns
    block_M <- d_data[rows_cols,rows_cols] # Select matrix block using rows+columns
    distances <- c(distances,mean(block_M[lower.tri(block_M)])) # Compute mean distance of non-zero matrix block elements
  } # End for loop i
  
  return(c(distances[-1],mean(distances[-1]),distances[1]))
} # End cluster_distances(function)

d_M <- sapply(d_M, cluster_distances)
#d_M <- rbind(d_M,d_M[5,]/d_M[6,],
 #            sapply(c_d_M, function(x) mean(x[lower.tri(x)]))
  #           )

#group_names <- unique(matrix(unlist(strsplit(rownames(p1[[2]]),"_")),ncol=2,byrow=T)[,1]) # Extract group names
#group_names <- unique(rownames(p1[[2]])) # Extract group names
#rownames(d_M) <- c(paste0("mean_",group_names),"Means within blocks","All means","Block/all means","Centroid means"); colnames(d_M) <- names(y_raw_norm) # Give block names + Norm names

#write.csv(round(d_M,3), paste0("cluster_dist_",dataset,".csv"))

pdf(paste0("cluster_standdist_",dataset,".pdf"), width=8.5,height=11)
pheatmap::pheatmap(t(apply(d_M, 1, function(x)x/max(x))), 
                   cex=0.95, legend=T,  cluster_rows=F, cluster_cols=F, show_rownames=T, row_names_side=3)
dev.off()

p01 <- p01[[1]]; p02 <- p02[[1]]; p03 <- p03[[1]]; p04 <- p04[[1]]; p05 <- p05[[1]]; p06 <- p06[[1]]; p07 <- p07[[1]]; 
p08 <- p08[[1]]; p09 <- p09[[1]]; p10 <- p10[[1]]; p11 <- p11[[1]]; p12 <- p12[[1]]; p13 <- p13[[1]]


# loads_as_bars: R function to plot gene loadings as bar charts
#
# Inputs:
# y_matrix    - a matrix of gene expression data (rows are genes, columns are samples)
# count_scale - optional, scale factor for the gene expression data (default=5e4)
# y_lim       - optional, range for the y-axis (default=c(0,1.0))
# shorten     - optional, logical indicating whether to shorten the plot (default=TRUE)
#
# Outputs:
# ggplot object with two bar charts of gene loadings for PC1 and PC2

# Start of loads_as_bars function definition
loads_as_bars <- function(y_matrix, normethod, y_lim=c(0,0.025), shorten=TRUE){
  
  # Perform PCA on the gene expression data
  pca <- prcomp(t(y_matrix), scale=T)
  
  # If shorten is TRUE, select the top 25 genes based on the variance
  if (shorten==T){
    varExp <- (pca$sdev^2)/sum(pca$sdev^2)*100 # Extract variance explained per principal component
    loadings1 <- pca$rotation[,1]; loadings2 <- pca$rotation[,2]; loadings3 <- pca$rotation[,3]
    
    varImp <- (loadings1^2)*varExp[1] + (loadings2^2)*varExp[2] #+ (loadings3^2)*varExp[3] # Quantify contribution/prominence of genes
    varImp <- rev(sort(varImp))

    sel_1000 <- data.frame(t(names(head(varImp, n=1000)))) # Assign 1000 most prominent genes to vector/data.frame
    row.names(sel_1000) <- normethod
    
    if (normethod=="Raw") appnd <- FALSE else appnd <- TRUE
    write.table(sel_1000, paste0("sel_genes_",dataset,".csv"), append=appnd, row.names=T, col.names=F, sep=",") # Write 1000 genes to file
    
    sel_genes <- names(head(varImp, n=25)) # Assign 25 most prominent genes to vector/data.frame
    pca$rotation <- pca$rotation[rownames(pca$rotation)%in%sel_genes,]
    pca$rotation <- pca$rotation[match(sel_genes,rownames(pca$rotation)),]
  } # End if()
  
  # Plot 1st PC
  l_xy1 <- data.frame(PC=rownames(pca$rotation),
                      loads=abs(pca$rotation[,1]), var_Imp=head(varImp,n=25)) %>%
    ggplot(aes(x=reorder(PC,-var_Imp),y=loads)) + labs(title="Loadings PC1") + xlab("") + ylim(y_lim) +
    theme(axis.text.x=element_text(angle=75, vjust=1, hjust=1), panel.grid.minor=element_blank(),
          plot.margin = margin(0.125,0.125,0.005,0.005, "cm")) + geom_col() 
  
  # Plot 2nd PC 
  l_xy2 <- data.frame(PC=rownames(pca$rotation),
                      loads=abs(pca$rotation[,2]), var_Imp=head(varImp,n=25)) %>%
    ggplot(aes(x=reorder(PC,-var_Imp),y=loads)) + labs(title="Loadings PC2") + xlab("Genes")  + ylim(y_lim) +
    theme(axis.text.x=element_text(angle=75, vjust=1, hjust=1), panel.grid.minor=element_blank(),
          plot.margin = margin(0.125,0.125,0.005,0.005, "cm")) + geom_col() 
  
  return(grid.arrange(l_xy1, l_xy2, ncol=1)) 
} # loads_as_bars()

###

l1 <- loads_as_bars(y_raw_norm$Raw,"Raw")
l2 <- loads_as_bars(y_raw_norm$TotCounts,"TotalCount") 
l3 <- loads_as_bars(y_raw_norm$UpQuart,"Upper Quartile")
l4 <- loads_as_bars(y_raw_norm$FullQuant,"Full Quantile")
l5 <- loads_as_bars(y_raw_norm$TMM,"TMM")
l6 <- loads_as_bars(y_raw_norm$RelLogExpr,"Rel Log Expr") ## Same...
l7 <- loads_as_bars(y_raw_norm$CenterLogRatio,"Center Log Ratio") ## ...
l8 <- loads_as_bars(y_raw_norm$RPKM,"RPKM")
l9 <- loads_as_bars(y_raw_norm$TPM,"TPM")
l10 <- loads_as_bars(y_raw_norm$GC_w,"Within-Lane GC")
l11 <- loads_as_bars(y_raw_norm$GC_wb,"Within+Between-Lane GC")
l12 <- loads_as_bars(y_raw_norm$PQN,"Probabilistic Quotient")
l13 <- loads_as_bars(y_raw_norm$CQN, "Conditional Quantile")

pdf(paste0("PCA+Loadings_PC12_",dataset,".pdf"), width=6, height=9)
grid.arrange(p01,l1, p02,l2, p03,l3, ncol=2)
grid.arrange(p04,l4, p05,l5, p06,l6, ncol=2)
grid.arrange(p07,l7, p08,l8, p09,l9, ncol=2)
grid.arrange(p10,l10, p11,l11, p12,l12, ncol=2) 
grid.arrange(p11,l11, p12,l12, p13,l13, ncol=2) 
dev.off()

# Print loadings for presenting

presentation <- FALSE

if (presentation==TRUE){
  pdf(paste0("Loadings_PC12_",dataset,"_raw.pdf"), width=6.5, height=4)
  #png(paste0("Loadings_PC12_",dataset,"_raw.png"), width=6.5, height=4, units="in", res=600)
  grid.arrange(l1, ncol=1)
  dev.off()
  pdf(paste0("Loadings_PC12_",dataset,"_TotalCount.pdf"), width=6.5, height=4)
  grid.arrange(l2, ncol=1)
  dev.off()
} # End presentation if()

#

sel_1000 <- read.table(paste0("sel_genes_",dataset,".csv"), row.names=1, header=F, sep=",")

# Select genes that appear influential for all normalization methods
commonground <- sel_1000[1,]
for (q in 2:nrow(sel_1000)) {commonground <- commonground[commonground %in% sel_1000[q,]]; print(length(commonground))}
length(commonground)
# commonground representing genes has length 170 (tumor) vs 14/487 (GSE218399)


gen_1000 <- sel_1000 # 12 x 1000 matrix


# Convert GeneIDs to GeneNames (tumor+)
if (dataset=="tumors" | dataset=="GSE216274"){
  for (k in 1:nrow(gen_1000)) {
    gen_1000[k,] <- y$genes$Symbol[which(rownames(y$counts) %in% unlist(sel_1000[k,]) )]
    gen_1000[k,] <- y$genes$Symbol[match(unlist(sel_1000[k,]),rownames(y$counts) )]
  } # End for loop k
} # End if()

# Convert GeneIDs to GeneNames (HCT116)
if (dataset=="GSE218399"){
  for (k in 1:nrow(gen_1000)) {
    gen_1000[k,] <- ENSG_IDs$To[which(rownames(y$counts) %in% unlist(sel_1000[k,]) )]
    gen_1000[k,] <- ENSG_IDs$To[match(unlist(sel_1000[k,]),rownames(y$counts) )]
  } # End for loop k
} # End if()

# Compute gene overlap per row/normalization method
gene_overlap <- matrix(ncol=nrow(sel_1000),nrow=nrow(sel_1000))
diag(gene_overlap) <- 1e3; colnames(gene_overlap) <- rownames(gen_1000); rownames(gene_overlap) <- rownames(gen_1000)
for (i in 1:(nrow(sel_1000)-1)) for (j in (i+1):nrow(sel_1000)) gene_overlap[i,j] <- length(which(sel_1000[i,] %in% sel_1000[j,])) #/ gene_overlap[i,i]
for (i in 1:(nrow(sel_1000)-1)) for (j in (i+1):nrow(sel_1000)) gene_overlap[j,i] <- length(which(sel_1000[i,] %in% sel_1000[j,])) #/ gene_overlap[j,j]


# Perform clustering and extract order
gene_heat <- pheatmap::pheatmap(gene_overlap,cex=1.25, legend=F, cluster_rows=T, cluster_cols=F)
gene_overlap <- gene_overlap[,gene_heat$tree_row[["order"]]] # reorder according to clustering of pheatmap

dev.off()
pdf(paste0("Heat_gene_overlap_",dataset,".pdf"), height=7.5, width=11) # dimensions...
pheatmap::pheatmap(gene_overlap, cex=1.2, legend=T,  cluster_rows=T, cluster_cols=F)
dev.off()


# VennDiagram of Gene overlap for selected normalization methods
library(VennDiagram)

rownames(sel_1000) <- norm_names

# GGVENN
library(ggvenn)
xyz <- list(
  Raw = sel_1000["Raw",], 
  TMM = sel_1000["TMM",], 
  CondQuant = sel_1000["Conditional Quantile",],
  TPM = sel_1000["TPM",]
)

pdf(paste0("Venn_Genes_4Norms_",dataset,".pdf"), width=8, height=6)
ggvenn(xyz,
       fill_color=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#FAEBD7", "#7FFFD4", "#FF4040", "#98F5FF"),
       stroke_size=1.5, set_name_size=4)
dev.off()


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Perform enrichment analysis   ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

### KEGG enrichment ###

KEGG_res <- list()
for (p in 1:nrow(gen_1000)){
  cat("Processing normalization method",p,"-",rownames(gen_1000)[p],"\n")
  
  if (dataset=="GSE216274") {
    genes2KEGG <- bitr(gen_1000[p,], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    KEGG_res[[length(KEGG_res)+1]] <- enrichKEGG(genes2KEGG$ENTREZID, organism='hsa', pvalueCutoff=0.10)
  } # End if()
  
  if (dataset=="tumors") {
    genes2KEGG <- bitr(gen_1000[p,], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    KEGG_res[[length(KEGG_res)+1]] <- enrichKEGG(genes2KEGG$ENTREZID, organism='hsa', pvalueCutoff=0.10)
  } # End if()
  
  if (dataset=="GSE218399") {
    genes2KEGG <- bitr(gen_1000[p,], fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    KEGG_res[[length(KEGG_res)+1]] <- enrichKEGG(genes2KEGG$ENTREZID, organism='hsa', pvalueCutoff=0.10)
  } # End if() 
} # End for loop p

#

names(KEGG_res) <- norm_names

pqr <- list(Raw=KEGG_res[["Raw"]]$Description,
            TMM=KEGG_res[["TMM"]]$Description,
            CQN=KEGG_res[["CQN"]]$Description,
            TPM=KEGG_res[["TPM"]]$Description)

pdf(paste0("Venn_KEGGpath_4Norms_",dataset,".pdf"), width=8, height=6)
ggvenn(pqr,
       fill_color=c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#FAEBD7", "#7FFFD4", "#FF4040", "#98F5FF"),
       stroke_size=1.5, set_name_size=4)
dev.off()

# Remove first list element because no enrichments were obtained for the colon cancer raw data
if (dataset=="GSE218399") {KEGG_res <- KEGG_res[-1]; sel_1000 <- sel_1000[-1,]; gen_1000 <- gen_1000[-1,]}

# GeneRatio or Counts barplots
pdf(paste0("KEGG_barplots_",dataset,".pdf"),height=11,width=8.5)
lapply(1:length(KEGG_res), function(x) barplot(KEGG_res[[x]], showCategory=20, title=names(KEGG_res)[x], x="Count")) # "Counts" may replace "GeneRatio"
dev.off()


# Gene-Concept Network

pdf(paste0("KEGG_cnetplots_",dataset,".pdf"),height=10,width=15)
lapply(1:length(KEGG_res), function(x) {
  p1 <- cnetplot(KEGG_res[[x]], node_label="category", showCategory=10, cex.params=list(category_label=1.2), 
                 color.params=list(category='firebrick',gene='steelblue')) # foldChange=y_raw_norm[[x]], circular=T, colorEdge=T, 
  # Note: node_label=c("gene","all"), cex.params=list(gene_label=0.8)) 
  cowplot::plot_grid(p1, labels=names(KEGG_res)[x])
})
dev.off()

#

p01 <- heatplot(KEGG_res[[1]], showCategory=5)
p02 <- heatplot(KEGG_res[[2]], showCategory=5)

cowplot::plot_grid(p01, p02, ncol=1, labels=LETTERS[1:2])

#

K01 <- pairwise_termsim(KEGG_res[[1]])
p1 <- treeplot(K01)
p2 <- treeplot(K01, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

#

pdf(paste0("KEGG_EnrichMap_",dataset,".pdf"),height=7,width=12)
lapply(1:length(KEGG_res), function(x) {
  p1 <- emapplot(pairwise_termsim(KEGG_res[[x]]), cex.params=list(category_node=1.5))
  cowplot::plot_grid(p1, labels=names(KEGG_res)[x])
  })
dev.off()

#KEGGs <- c("Estrogen signaling pathway","Cellular senescence","Protein processing in","endoplasmic reticulum","Human cytomegalovirus","infection","Regulation of actin","cytoskeleton","Fluid shear stress and atherosclerosis","Viral myocarditis","Amoebiasis","Adherens junction","Bacterial invasion of epithelial cells","Protein digestion and absorption","Small cell lung cancer","Dilated cardiomyopathy","Arrhythmogenic right","ventricular cardiomyopathy","Hypertrophic cardiomyopathy","PI3K−Akt signaling pathway","Human papillomavirus infection","Proteoglycans in cancer","ECM−receptor interaction","Focal adhesion")
# Biological theme comparison

#pdf(paste0("KEGG_BioTheme_",dataset,".pdf"),height=7,width=12)
#lapply(1:length(KEGG_res), function(x) {
 # if (dataset=="GSE218399") genes2KEGG <- bitr(gen_1000[x,], fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#  if (dataset=="tumors") genes2KEGG <- bitr(gen_1000[x,], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#  if (dataset=="GSE216274") genes2KEGG <- bitr(gen_1000[x,], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#  plot_label <- names(KEGG_res)[x]
#  xx <- compareCluster(genes2KEGG, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05); xx <- pairwise_termsim(xx)                     
#  p4 <- emapplot(xx, pie.params=list(pie="count", legend_n=2), cex.params=list(category_node=1.5), layout.params=list(layout="kk"))
#  cowplot::plot_grid(emapplot(xx, cex.params=list(category_node=0.5)), labels=plot_label)
#  })
#dev.off()

# https://yulab-smu.top/biomedical-knowledge-mining-book/
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/11_FA_functional_class_scoring.html
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# UpSet plot

pdf(paste0("KEGG_Upset_",dataset,".pdf"),height=7,width=12)
lapply(1:length(KEGG_res), function(q) cowplot::plot_grid(upsetplot(KEGG_res[[q]]), labels=names(KEGG_res)[q]))
dev.off()


KEGG_res[[1]]@result$Description %in% "" # Check for empty KEGG Is
for (z in 1:3) browseKEGG(KEGG_res[[1]][,-8], KEGG_res[[1]]$ID[z]) # Visualize KEGG network/Cell cycle online


# Compute total no of KEGG IDs across normalization methods
tot_KEGG <- unique(unlist(sapply(KEGG_res, function(x) x$ID))) # 108 KEGG IDs in total;

# Compute size of the core of the KEGG IDs that appear for every normalization method
common_KEGG <- KEGG_res[[1]]$ID
for (q in 1:length(KEGG_res)) common_KEGG <- common_KEGG[common_KEGG %in% KEGG_res[[q]]$ID] # Compute core KEGG IDs
length(common_KEGG) # 1 common KEGG IDs

# Generate matrix for binary indicating if a KEGG ID was enriched per normalization method
KEGG_binary <- matrix(nrow=length(KEGG_res),ncol=length(tot_KEGG)) # Initialize matrix
rownames(KEGG_binary) <- names(KEGG_res); colnames(KEGG_binary) <- tot_KEGG # Assign rownames (norm. methods) and colnames (KEGG IDs)

for (z in 1:nrow(KEGG_binary)) KEGG_binary[z,] <- as.integer(tot_KEGG %in% KEGG_res[[z]]$ID) # Assign presence of KEGG IDs per norm method (either 0 or 1)

write.csv(KEGG_binary,paste0("KEGG_binary_",dataset,".csv"))

#
PCA_bin_ggplot <- function(x, title, leg_pos='none', s_impose=FALSE){
  xx <- prcomp(t(as.matrix(x)), scale=F)
  
  ## more elegant plot using the ggplot library.
  z <- summary(xx); print(z$importance[3,which(z$importance[3,] > 0.75)])
  PC1 <- as.numeric(xx$x[,1]); PC2 <- as.numeric(xx$x[,2]); if (s_impose==TRUE) PC2 <- -1*PC2
  df <- data.frame(PC1, PC2, Treatment=rownames(xx$x))
  df$PC1 <- df$PC1/max(abs(df$PC1)); df$PC2 <- df$PC2/max(abs(df$PC2))
  print(df)
  myplot <- ggplot(df,aes(PC1,PC2), color='black') +
    geom_text(aes(label=rownames(xx$x))) +
    scale_color_manual(values=c('magenta','limegreen','dodgerblue','black')) + scale_shape_manual(values=c(15:17,19)) +
    xlab(paste0("PC1 (",format(100*z$importance[2,1],digits=3),"%)")) +
    ylab(paste0("PC2 (",format(100*z$importance[2,2],digits=3),"%)")) +
    theme_bw() + xlim(-1.25,1.25) + ylim(-1.25,1.25) +
    theme(axis.text.x=element_text(size=11), legend.position=leg_pos, 
          legend.text=element_text(size=6), legend.title=element_text(size=7)) + ggtitle(title)
  ###
  
  return(myplot)
} # End of PCA_bin_ggplot()

pdf(paste0("PCA_KEGG_binary_",dataset,".pdf"),width=6,height=6)
#png(paste0("PCA_KEGG_binary_",dataset,".png"),width=6,height=6, units="in", res=600)
if (dataset=="tumors") PCA_bin_ggplot(t(KEGG_binary),"") # PCA on KEGG binary matrix
if (dataset=="GSE218399") PCA_bin_ggplot(t(KEGG_binary),"") # Does not work for 218399
if (dataset=="GSE216274") PCA_bin_ggplot(t(KEGG_binary),"") # Does not work for 216274
dev.off()

PCA_bin_ggplot(t(KEGG_binary),"") # Does not work for 216274

PCA_KEGG_bin_tumors <- read.csv("KEGG_binary_tumors.csv",sep=',', row.names=1)
PCA_KEGG_bin_colon <- read.csv("KEGG_binary_GSE218399.csv", row.names=1)
PCA_KEGG_bin_neuroblastoma <- read.csv("KEGG_binary_GSE216274.csv", row.names=1)

pca_bin_01 <- PCA_bin_ggplot(t(PCA_KEGG_bin_tumors),"Tumors") # PCA on KEGG binary matrix
pca_bin_02 <- PCA_bin_ggplot(t(PCA_KEGG_bin_colon),"Colon cancer") # Does not work for 218399
pca_bin_03 <- PCA_bin_ggplot(t(PCA_KEGG_bin_neuroblastoma),"Neuroblastoma") # Does not work for 216274

#png(paste0("PCA_KEGG_binary_",dataset,".png"),width=15,height=6, units="in", res=600)
pdf(paste0("PCA_KEGG_binary_alldata.pdf"),width=15/1.25,height=5/1.25)
grid.arrange(pca_bin_03,pca_bin_01,pca_bin_02,
             ncol=3)
dev.off()


### Perform non-metric multidimensional scaling

# Check: https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/nmds

library(vegan) 

NMDS_ggplot <- function(x, title, x_lim, y_lim, dista="manhattan", kay){
  # Dissimilarity index, partial match to "manhattan", "euclidean", "canberra",
  #  "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita",
  #  "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis",
  #  "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison".
  
  xx <- metaMDS(t(x), distance=dista, k=kay, trymax=500)
  strss <- round(xx$stress,3); print(strss)
  
  #stressplot(xx) # visualise the Shepard stress plot. 
  #ordiplot(xx); plot(zz, type='t')
  
  ## more elegant plot using the ggplot library.
  species.scores <- as.data.frame(scores(xx, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
  species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
  head(species.scores)  #look at the data
  
  species.scores$NMDS1 <- species.scores$NMDS1/max(abs(species.scores$NMDS1))
  species.scores$NMDS2 <- species.scores$NMDS2/max(abs(species.scores$NMDS2))
  
  myplot <- ggplot() + geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=1.5) + 
    theme_bw() + ggtitle(paste0(title," - stress=",strss)) + xlim(-1.5,1.5) + ylim(-1.5,1.5) # add the species labels
  
  return(myplot)
} # NMDS_ggplot()


if (dataset=="tumors") xlims <- 1.5*c(-1,1); ylims <- 1.5*c(-1,1); k_n <- 3
if (dataset=="GSE216274") xlims <- 2.5*c(-5,5); ylims <- 1.5*c(-4,4); k_n <- 3
if (dataset=="GSE218399") xlims <- 3*c(-6,6); ylims <- 1*c(-2,6); k_n <- 4
pdf(paste0("NMDS_KEGG_binary_",dataset,"_Manhattan.pdf"),height=6,width=6)
NMDS_ggplot(KEGG_binary,title="Manhattan",xlims,ylims,dista="manhattan",kay=k_n)
dev.off()

if (dataset=="tumors") xlims <- 2*c(-1,0.75); ylims <- 1.5*c(-1,1.5)
if (dataset=="GSE216274") xlims <- c(-5,5); ylims <- 1.25*c(-1.5,1.5)
if (dataset=="GSE218399") xlims <- c(-5,5); ylims <- c(-2,1)
pdf(paste0("NMDS_KEGG_binary_",dataset,"_Euclidean.pdf"),height=6,width=6)
NMDS_ggplot(KEGG_binary,title="Euclidean",xlims,ylims,dista="euclidean",kay=k_n)
dev.off()

# Jaccard and Bray don't seem to work
if (k_n == 5){
  if (dataset=="tumors") xlims <- c(-2,1); ylims <- c(-1,1.5)
  NMDS_ggplot(KEGG_binary,title="Jaccard",xlims,ylims,dista="jaccard",kay=k_n)
  NMDS_ggplot(KEGG_binary,title="Bray",xlims,ylims,dista="bray",kay=k_n)
}

#

y_max <- 1.05*max(sapply(KEGG_res, nrow))

pdf(paste0("No_KEGG_",dataset,".pdf"),height=8.0, width=11)
data.frame(KEGG_names=names(KEGG_res),
           no_KEGG=sapply(KEGG_res, nrow), KEGG_order=1:length(KEGG_res)) %>%
  ggplot(aes(x=reorder(KEGG_names,KEGG_order),y=no_KEGG)) + labs(title="KEGG IDs per normalization method and core size") + ylab("KEGG ID #") + xlab("") + ylim(c(0,y_max)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), panel.grid.minor=element_blank(),
        plot.margin=margin(0.1,0.0,0.00,0.00, "cm"), text = element_text(size=20)) + geom_col() +
  geom_hline(yintercept=length(common_KEGG), linetype="dashed", color="magenta", linewidth=1.5)
dev.off()

#

KEGG_overlap <- matrix(ncol=nrow(sel_1000),nrow=nrow(sel_1000))
diag(KEGG_overlap) <- sapply(KEGG_res, nrow); colnames(KEGG_overlap) <- names(KEGG_res); rownames(KEGG_overlap) <- names(KEGG_res)
for (i in 1:(nrow(sel_1000)-1)) for (j in (i+1):nrow(sel_1000)) KEGG_overlap[i,j] <- length(which(KEGG_res[[i]]$ID %in% KEGG_res[[j]]$ID)) / KEGG_overlap[j,j]
for (i in 1:(nrow(sel_1000)-1)) for (j in (i+1):nrow(sel_1000)) KEGG_overlap[j,i] <- length(which(KEGG_res[[i]]$ID %in% KEGG_res[[j]]$ID)) / KEGG_overlap[i,i]
diag(KEGG_overlap) <- 1
KEGG_overlap[1:3,1:3]


# Perform clustering and extract order
KEGG_heat <- pheatmap::pheatmap(KEGG_overlap,cex=1.25, legend=F, cluster_rows=T, cluster_cols=F)
KEGG_overlap <- KEGG_overlap[,KEGG_heat$tree_row[["order"]]] # reorder according to clustering of pheatmap

dev.off()
pdf(paste0("Heat_KEGG_overlap_",dataset,".pdf"),height=7.5,width=11)
pheatmap::pheatmap(KEGG_overlap,cex=1.25, legend=T,  cluster_rows=T, cluster_cols=F)
dev.off()


# Evaluate KEGG units/IDs per normalization method

KEGG_res <- lapply(KEGG_res, function(x) x[x$p.adjust < 0.05,]) # Retain KEGG pathways with p.adjust < 0.05
descriptions <- unique(unlist(sapply(KEGG_res, function(x) x$Description))) # Extract all KEGG IDs

KEGG_per_Norm <- sapply(KEGG_res, function(x) as.integer(descriptions %in% x$Description)) # Indicate if KEGG ID
rownames(KEGG_per_Norm) <- descriptions; colnames(KEGG_per_Norm) <- names(KEGG_res) # Assign rownames (KEGG units) and colnames (normalization methods)

pdf(paste0("KEGG_per_Norm_scaled_",dataset,".pdf"), width=8.5, height=11)
pheatmap::pheatmap(KEGG_per_Norm, cex=0.95, legend=F,  cluster_rows=F, cluster_cols=F, color=c("red","limegreen"), show_rownames=T, row_names_side=2)


# Shortened

select_norms <- c(1,5:6,9,12:13)
dscrptns <- unique(unlist(sapply(KEGG_res[select_norms], function(x) x$Description))) # Extract all KEGG IDs

KEGG_per_Norm_short <- sapply(KEGG_res[select_norms], function(x) as.integer(dscrptns %in% x$Description)) # Indicate if KEGG ID
rownames(KEGG_per_Norm_short) <- dscrptns; colnames(KEGG_per_Norm_short) <- rownames(gen_1000)[select_norms] # Assign rownames (KEGG units) and colnames (normalization methods)

dev.off()
pdf(paste0("KEGG_per_Norm_scaled_",dataset,"_short.pdf"), width=8.5, height=11)
pheatmap::pheatmap(KEGG_per_Norm_short, 
                   cex=0.95, legend=F,  cluster_rows=F, cluster_cols=F, color=c("red","limegreen"), show_rownames=T, row_names_side=2)
dev.off()

###

for (p in 1:ncol(KEGG_per_Norm)){
  names4norm <- which(rownames(KEGG_per_Norm) %in% KEGG_res[[p]]$Description)
  KEGG_per_Norm[names4norm,p] <- KEGG_res[[p]]$Count
}

pdf(paste0("KEGG_per_Norm_grad_scaled_",dataset,".pdf"), width=8.5, height=11)
pheatmap::pheatmap(KEGG_per_Norm, cex=0.95, legend=T,  cluster_rows=F, cluster_cols=F, show_rownames=T, row_names_side=2)
pheatmap::pheatmap(t(apply(KEGG_per_Norm, 1, function(x) x/max(x))), 
                   cex=0.95, legend=T,  cluster_rows=F, cluster_cols=F, show_rownames=T, row_names_side=2)
dev.off()

###

print("KEGG pathway enrichment analysis ends here!")


# Then evaluate if KEGG is present in top 20?!


####################################
### Plotting code by Dr Saccenti ###
####################################

library(ggrepel)

#

PCA_bin_ggplot_edo <- function(x, title, leg_pos='none', s_impose=FALSE){
  xx <- prcomp(as.matrix(x), scale=F); z <- summary(xx)
  print(z$importance[3,which(z$importance[3,] > 0.75)])
  
  ## more elegant plot using the ggplot library.
  df <- data.frame(PC1=as.numeric(xx$x[,1]), PC2=as.numeric(xx$x[,2]), Treatment=rownames(xx$x)); 
  df$PC1 <- df$PC1/max(abs(df$PC1)); df$PC2 <- df$PC2/max(abs(df$PC2)); if (s_impose==TRUE) df$PC2 <- -1*df$PC2
  m_cols <- c('Black','grey','limegreen','hotpink','dodgerblue','cyan','bisque3','red','gold','salmon','navy','aquamarine2','darkorchid2'); if (nrow(df)==12) m_cols <- m_cols[-1]
  
  myplot <- ggplot(df, aes(x=PC1, y=PC2)) + geom_point(colour=m_cols[1:nrow(df)], size=4) + theme_bw(base_size=20) +
    theme(aspect.ratio=1, axis.title.x=element_blank(), axis.title.y=element_blank(), plot.title=element_text(size=15, margin=margin(0,0,0,0)),
          axis.text=element_text(size=15), plot.margin=margin(0.15,0.1,0.1,0.15, "cm"),) + ggtitle(title) +
    xlim(-1.1,1.1) + ylim(-1.1,1.1) +
    geom_label_repel(aes(label = Treatment),
                     size=4, color=m_cols[1:nrow(df)],
                     box.padding=0.75, point.padding=0.1,
                     segment.color='grey',
                     max.overlaps=50, label.size=0.3)
  ###
  return(myplot)
} # End of PCA_bin_ggplot_edo()

rot_overlap <- read.csv("OverlapGenes_May2024_4Figure.csv", row.names=1)

pdf("Overlap_PCA_rotation.pdf", height=6, width=6)
PCA_bin_ggplot_edo(rot_overlap[rot_overlap$Difference >= 0.1,],"Rotation overlap")
dev.off()

#

p_Henk_n <- PCA_bin_ggplot_edo(PCA_KEGG_bin_neuroblastoma,"A) Neuroblastoma"); p_Henk_t <- PCA_bin_ggplot_edo(PCA_KEGG_bin_tumors, "B) Tumors")
p_Henk_c <- PCA_bin_ggplot_edo(PCA_KEGG_bin_colon, "C) Colon cancer")

pdf(paste0("PCA_KEGG_bin_alldata_Edo.pdf"),width=15/1.25,height=5/1.25)
grid.arrange(p_Henk_n, p_Henk_t, p_Henk_c,
             ncol=3)
dev.off()

###
