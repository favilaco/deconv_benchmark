args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=9){

  	print("Please check that all required parameters are indicated or are correct")
  	print("Example usage for bulk deconvolution methods: 'Rscript Master_deconvolution.R baron none bulk TMM all nnls 100 none 1'")
  	print("Example usage for single-cell deconvolution methods: 'Rscript Master_deconvolution.R baron none sc TMM TMM MuSiC 100 none 1'")
  	stop()
} 

### arguments
dataset = args[1]
transformation = args[2]
deconv_type = args[3]

if(deconv_type == "bulk"){
	normalization = args[4]
	marker_strategy = args[5]
} else if (deconv_type == "sc") {
	normalization_scC = args[4]
	normalization_scT = args[5]
} else {
	print("Please enter a valid deconvolution framework")
	stop()
}

method = args[6]
number_cells = round(as.numeric(args[7]), digits = -2) #has to be multiple of 100
to_remove = args[8]
num_cores = min(as.numeric(args[9]),parallel::detectCores()-1)

#-------------------------------------------------------
### Helper functions + CIBERSORT external code + CT proportions
source('./helper_functions.R')
# source('CIBERSORT.R')

p2 <- as.matrix(read.table("./p2.txt"))
p3 <- as.matrix(read.table("./p3.txt"))
p4 <- as.matrix(read.table("./p4.txt"))
p5 <- as.matrix(read.table("./p5.txt"))

#-------------------------------------------------------
### Read data and metadata
data = readRDS(list.files(path = dataset, pattern = "rds", full.names = TRUE))
full_phenoData = read.table(list.files(path = dataset, pattern = "phenoData", full.names = TRUE), header=TRUE)

#-------------------------------------------------------
### QC 
require(dplyr); require(Matrix)

# First: cells with library size, mitochondrial or ribosomal content further than three MAD away were discarded
filterCells <- function(filterParam){
	cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
	cellsToRemove
}

libSizes <- colSums(data)
gene_names <- rownames(data)

mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mtPercent <- colSums(data[mtID, ])/libSizes
rbPercent <- colSums(data[rbID, ])/libSizes

lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
	unlist() %>% 
	unique() -> cellsToRemove

if(length(cellsToRemove) != 0){
	data <- data[,-cellsToRemove]
	full_phenoData <- full_phenoData[-cellsToRemove,]
}

# Keep only "detectable" genes: at least 30% of cells within a CT have a read/UMI count different from 0
cellType <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])
keep <- sapply(unique(cellType), function(x) {
        CT_hits = which(cellType %in% x)
        size = ceiling(0.30*length(CT_hits))
        Matrix::rowSums(data[,CT_hits,drop=FALSE] != 0) >= size
    })

data = data[Matrix::rowSums(keep) > 0,]

#-------------------------------------------------------
### Data split into training/test 
set.seed(24)
require(limma); require(dplyr); require(pheatmap)

original_cell_names = colnames(data)
colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])

# Keep CTs with >= 50 cells after QC
cell_counts = table(colnames(data))
to_keep = names(cell_counts)[cell_counts >= 50]
pData <- full_phenoData[full_phenoData$cellType %in% to_keep,]
to_keep = which(colnames(data) %in% to_keep)   
data <- data[,to_keep]
original_cell_names <- original_cell_names[to_keep]


# Data split into train & test  
training <- as.numeric(unlist(sapply(unique(colnames(data)), function(x) {
            sample(which(colnames(data) %in% x), cell_counts[x]/2) })))
test <- which(!1:ncol(data) %in% training)

# Generate phenodata for reference matrix C
pDataC = pData[training,]

train <- data[,training]
test <- data[,test]

# "write.table" & "saveRDS" statements are optional, for users willing to avoid generation of matrix C every time:    
# write.table(pDataC, file = paste(dataset,"phenoDataC",sep="_"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

train_cellID = train
colnames(train_cellID) = original_cell_names[training]
# saveRDS(object = train_cellID, file = paste(dataset,"qc_filtered_train.rds",sep="_")) #It has to contain cellID as colnames, not cellType (for scRNA-seq methods)
# saveRDS(object = test, file = paste(dataset,"qc_filtered_test.rds",sep="_"))

# reference matrix (C) + refProfiles.var from TRAINING dataset
cellType <- colnames(train)
group = list()
for(i in unique(cellType)){ 
	group[[i]] <- which(cellType %in% i)
}
C = lapply(group,function(x) Matrix::rowMeans(train[,x])) #C should be made with MEAN, not sum! (to agree with the way markers were selected)
C = round(do.call(cbind.data.frame, C))
# write.table(C, file = "C",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE,)

refProfiles.var = lapply(group,function(x) train[,x])
refProfiles.var = lapply(refProfiles.var, function(x) matrixStats::rowSds(Matrix::as.matrix(x)))
refProfiles.var = round(do.call(cbind.data.frame, refProfiles.var))
rownames(refProfiles.var) <- rownames(train)
# write.table(refProfiles.var, "refProfiles.var", quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

#-------------------------------------------------------
#Normalization of "train" followed by marker selection 
train2 = Normalization(train)

# INITIAL CONTRASTS for marker selection WITHOUT taking correlated CT into account 
#[compare 1 grp with avg expression of all other groups (COLUMN-WISE)]
annotation = factor(colnames(train2))
design <- model.matrix(~0+annotation)
colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
colnames(cont.matrix) <- colnames(design)
diag(cont.matrix) <- (ncol(design)-1)/ncol(design)

v <- limma::voom(train2, design=design, plot=FALSE) 
fit <- limma::lmFit(v, design)
fit2 <- limma::contrasts.fit(fit, cont.matrix)
fit2 <- limma::eBayes(fit2, trend=TRUE)

topTable_RESULTS = limma::topTable(fit2,coef=1:ncol(cont.matrix),number=Inf, adjust.method="BH", p.value=0.05, lfc=round(log2(1.5),2))
AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]

markers <- apply(topTable_RESULTS,1,function(x){
         temp = sort(x)
         ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= 0.3) | (abs(temp[1] - temp[2]) >= 0.3)
     })

topTable_RESULTS = topTable_RESULTS[markers,]

topTable_RESULTS <- cbind.data.frame(rownames(topTable_RESULTS),
                                t(apply(topTable_RESULTS, 1, function(x){
                                    temp = max(x)
                                    if(temp < round(log2(1.5),2)){
                                        temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                    } else {
                                        temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                    } 
                                    temp
                                })))
                                    
colnames(topTable_RESULTS) <- c("gene","log2FC","CT")
topTable_RESULTS$log2FC = as.numeric(as.character(topTable_RESULTS$log2FC))
topTable_RESULTS <- topTable_RESULTS %>% dplyr::arrange(CT,desc(log2FC))

if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }
topTable_RESULTS$AveExpr <- AveExpr_pval$AveExpr[match(topTable_RESULTS$gene,rownames(AveExpr_pval))]
#write.table(topTable_RESULTS, file = "markers",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE,)

#-------------------------------------------------------
### Generation of 1000 pseudo-bulk mixtures (T) (on test data)
cellType <- colnames(test)
generator <- Generator(test,full_phenoData,number_cells,dataset)
T <- generator[["T"]]
P <- generator[["P"]]

#-------------------------------------------------------
### Transformation, scaling/normalization, marker selection for bulk deconvolution methods and deconvolution:

if(deconv_type == "bulk"){

	T = Transformation(T, transformation)
	C = Transformation(C, transformation)

	T = Scaling(T, normalization)
	C = Scaling(C, normalization)

	# marker selection (on training data) 
	marker_distrib = marker_strategies(topTable_RESULTS, marker_strategy, C)

	#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
	if(to_remove != "none"){

		T <- T[,P[to_remove,] != 0]
		C <- C[, colnames(C) %in% rownames(P) & (!colnames(C) %in% to_remove)]
		P <- P[!rownames(P) %in% to_remove, colnames(T)]
		refProfiles.var = refProfiles.var[,colnames(refProfiles.var) %in% rownames(P) & (!colnames(refProfiles.var) %in% to_remove)]
		marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(P) & (marker_distrib$CT != to_remove),]

	}

	RESULTS = Deconvolution(T = T, C = C, method = method, P = P, elem = to_remove, marker_distrib = marker_distrib, refProfiles.var = refProfiles.var) 

} else if (deconv_type == "sc"){

	T = Transformation(T, transformation)
	C = Transformation(train_cellID, transformation)

	T = Scaling(T, normalization_scT)
	C = Scaling(C, normalization_scC)

	#If a cell type is removed, only meaningful mixtures where that CT was present (proportion < 0) are kept:
	if(to_remove != "none"){

		T <- T[,P[to_remove,] != 0]
		C <- C[,pDataC$cellType != to_remove]
		P <- P[!rownames(P) %in% to_remove, colnames(T)]
		pDataC <- pDataC[pDataC$cellType != to_remove,]

	}

	RESULTS = Deconvolution(T = T, C = C, method = method, phenoDataC = pDataC, P = P, elem = to_remove, refProfiles.var = refProfiles.var) 
		
}

RESULTS = RESULTS %>% dplyr::summarise(RMSE=sqrt(mean((observed_values-expected_values)^2)), Pearson=cor(observed_values,expected_values))
print(RESULTS)