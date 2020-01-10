#source('./CIBERSORT.R')

Normalization <- function(data){

    data <- edgeR::DGEList(data)
    CtrlGenes <- grep("ERCC-",rownames(data))

    if(length(CtrlGenes)>1){

        spikes <- data[CtrlGenes,]
        spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
        data$samples$norm.factors <- spikes$samples$norm.factors

    } else {

        data <- edgeR::calcNormFactors(data, method = "TMM")  

    }

    return(data)

}



Generator <- function(test, full_phenoData, ncells, dataset){

    require(plyr); require(data.table); require(Matrix); require(matrixStats); require(dplyr); require(gtools)

    for(i1 in ncells){
        
        print(i1)
        random_string = as.character(i1) 

        ################################################################################################################
        #INPUT generator
        cellType <- colnames(test)
        cellTypes <- unique(cellType)
        group = sapply(cellTypes, function(x) which(cellType %in% x))
        cell_numbers <- table(cellType)

        CT_composition = list()

        count = 1
        while(count <= 1000){ #this one ensures proportions and number of cells are feasible with the dataset

            p_num_CTs <- sample(2:min(length(cellTypes),5),1)
            candidate_comp <- cbind.data.frame(sample(cellTypes,p_num_CTs),get(paste("p",p_num_CTs,sep=""))[sample(1:nrow(get(paste("p",p_num_CTs,sep=""))),1,replace = FALSE),])
            colnames(candidate_comp) <- c("cell_type","proportion")

            if(sum(!sapply(1:nrow(candidate_comp), function(x) cell_numbers[as.character(candidate_comp$cell_type[x])] >= candidate_comp$proportion[x] * i1))==0){
            CT_composition[[paste("tissue",count,sep="_")]] <- candidate_comp
            count = count + 1
            }

        }
        
        #composition in matrix form:
        COMP = suppressMessages(as.data.frame(data.table::dcast(data.table::melt(CT_composition),cell_type ~ L1, fill = 0)))
        rownames(COMP) <- COMP$cell_type
        COMP$cell_type <- NULL
        P <- COMP
        #P: add "0" for those present in C but not randomly assigned a proportion
        not_present <- cellTypes[!cellTypes %in% rownames(COMP)]
        not_present <- matrix(0,nrow=length(not_present), ncol = ncol(P),dimnames = list(not_present,colnames(P)))
        
        if(nrow(not_present) > 0){P <- rbind.data.frame(P, not_present)}

        P = P[,gtools::mixedsort(colnames(P))]
        P = P[gtools::mixedsort(rownames(P)),]
        #write.table(P, paste("P",random_string,sep="_"), quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")

        cell_composition <- P * i1
        cell_composition$CT <- rownames(cell_composition)
        cell_composition <- data.table::melt(cell_composition)
        cell_composition <- cell_composition[cell_composition$value != 0,]

        chosen_cells <- sapply(names(CT_composition), function(x){
            temp = cell_composition[cell_composition$variable == x,]
            unlist(apply(temp,1, function(y){sample(which(cellType==y[1]),size=y[3])}))
        })

        colnames(chosen_cells) <- colnames(P)

        T <- as.data.frame(apply(chosen_cells,2, function(x){
            Matrix::rowSums(test[,x])
        }))
        
       return(list(T=T,P=P))

    }

}



marker_strategies <- function(marker_distrib, marker_strategy, C){

    set.seed(4)

    if(marker_strategy == "all"){
        
        #using all markers that were found
        markers = marker_distrib

    } else if (marker_strategy == "pos_fc"){

        # using only markers with positive FC (=over-expressed in cell type of interest)
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% as.data.frame()

    } else if (marker_strategy == "top_50p_logFC"){

        # top 50% of markers (per CT) based on logFC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "bottom_50p_logFC"){

        # bottom 50% of markers based on logFC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "top_50p_AveExpr"){

        # top 50% of markers based on average gene expression (baseline expression)
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(ceiling(n()*0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "bottom_50p_AveExpr"){

        # low 50% based on average gene expression.
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(AveExpr)) %>% group_by(CT) %>% dplyr::top_n(floor(n()*-0.5), wt=log2FC) %>% as.data.frame()

    } else if (marker_strategy == "top_n2"){

        # using the top 2 genes/CT with highest log2FC
        markers = marker_distrib %>% dplyr::filter(log2FC > 0) %>% dplyr::arrange(CT, desc(log2FC)) %>% group_by(CT) %>% dplyr::top_n(2, wt=log2FC) %>% as.data.frame()
    
    } else if (marker_strategy == "random5"){

        # using 5 random markers for each different cell types
        markers = marker_distrib[1:(ncol(C)*5),]
        markers$CT = rep(colnames(C),5) #labelling purposes: important for semi-supervised
        markers$gene = sample(rownames(C), nrow(markers), replace = FALSE)

    }

    return(markers) 

}



Transformation <- function(matrix, option){
    
    #############################################################
    ##########    DATA TRANSFORMATION (on full data)   ##########
    if(option=="none"){
        
        matrix = matrix

    }

    if(option=="log"){
        
        matrix = log1p(matrix)

    }

    if(option=="sqrt"){
        
        matrix = sqrt(matrix)

    }

    if(option=="vst"){
        
        matrix = DESeq2::varianceStabilizingTransformation(as.matrix(matrix))

    }

    return(matrix)

}




Scaling <- function(matrix, option, phenoDataC=NULL){

    ##########    Remove rows & columns full of zeroes   ##########
    #Error: NMF::nmf - Input matrix x contains at least one null or NA-filled row.
    matrix = matrix[rowSums(matrix)!=0,]
    #OLS with error if all elements within a row are equal (e.g. all 0, or all a common value after log/sqrt/vst transformation)
    matrix = matrix[!apply(matrix, 1, function(x) var(x) == 0),]
    matrix = matrix[,colSums(matrix)!=0]

    if(option=="column"){
        
        matrix = apply(matrix,2,function(x) x/sum(x)) 

    } else if(option=="row"){ 
        
        matrix = t(apply(matrix,1,function(x) x/sum(x))) 

    } else if(option=="mean"){ 
        
        matrix = apply(matrix,2,function(x) x - mean(x)) 

    } else if(option=="column_z-score"){ 
        
        matrix = scale(matrix, center = TRUE, scale = TRUE)

    } else if(option=="global_z-score"){
        
        matrix = (matrix - mean(as.matrix(matrix))) / sd(as.matrix(matrix))

    } else if(option=="column_min-max"){
        
        matrix = apply(matrix, 2, function(x) (x - min(x))/(max(x) - min(x)))

    } else if(option=="global_min-max"){
        
        matrix = (matrix - min(matrix))/(max(matrix) - min(matrix))

    } else if (option=="LogNormalize"){
        
        matrix = as.matrix(expm1(Seurat::LogNormalize(matrix, display.progress = FALSE)))

    } else if (option=="QN"){

        matrix_rownames <- rownames(matrix); matrix_colnames <- colnames(matrix)

        matrix = preprocessCore::normalize.quantiles(as.matrix(matrix))

        rownames(matrix) <- matrix_rownames; colnames(matrix) <- matrix_colnames

    } else if (option=="TMM"){# CPM counts coming from TMM-normalized library sizes; https://support.bioconductor.org/p/114798/

        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(data))
  
        if(length(CtrlGenes)>1){
            
            spikes <- data[CtrlGenes,]
            spikes <- edgeR::calcNormFactors(spikes, method = "TMM") 
            matrix$samples$norm.factors <- spikes$samples$norm.factors
        
        } else {
        
            matrix <- edgeR::calcNormFactors(matrix, method = "TMM")  
        
        }

        matrix <- edgeR::cpm(matrix)

    } else if (option=="UQ"){

        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        matrix <- edgeR::DGEList(counts=matrix, group=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(data))
  
        if(length(CtrlGenes)>1){
            
            spikes <- data[CtrlGenes,]
            spikes <- edgeR::calcNormFactors(spikes, method = "upperquartile") 
            matrix$samples$norm.factors <- spikes$samples$norm.factors
        
        } else {
        
            matrix <- edgeR::calcNormFactors(matrix, method = "upperquartile")  
        
        }

        matrix <- edgeR::cpm(matrix)

    } else if (option=="median_ratios"){#requires integer values
        
        if(!is.null(phenoDataC)){#use CT info for scRNA-seq

            Celltype = as.character(phenoDataC$cellType[phenoDataC$cellID %in% colnames(matrix)])

        } else {

            Celltype = colnames(matrix)

        }

        metadata <- data.frame(Celltype=Celltype)
        CtrlGenes <- grep("ERCC-",rownames(matrix))
        matrix = DESeq2::DESeqDataSetFromMatrix(matrix, colData = metadata, design = ~ Celltype)

        if(length(CtrlGenes)>1){

            dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio", controlGenes = CtrlGenes)

        } else {

            dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio")
            
        }

        matrix <- DESeq2::counts(dds, normalized=TRUE)
    
    } else if (option=="TPM"){

        require(SingleR)
        data(human_lengths)
        matrix = SingleR::TPM(counts = matrix, lengths = human_lengths)
        rownames(matrix) = toupper(rownames(matrix))
        detach(package:SingleR, unload=TRUE)
    
    ####################################################################################
    ## scRNA-seq specific  

    } else if (option=="RNBR"){

        #Model formula is y ~ log_umi
        
        ##following line needed to solve "Wrong R type for mapped matrix"
        #https://github.com/ChristophH/sctransform/issues/24
        matrix = as(matrix, "dgCMatrix")
        matrix = sctransform::vst(matrix, return_corrected_umi=TRUE, show_progress = FALSE)$umi_corrected
        matrix = as(matrix, "matrix")

    } else if (option=="scran"){
        
        sf = scran::computeSumFactors(as.matrix(matrix), clusters=NULL) 
        sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts=as.matrix(matrix)))
        sizeFactors(sce) <- sf

        sce = scater::normalize(sce,exprs_values = "counts", return_log = FALSE) 
        matrix = normcounts(sce)
        
    } else if (option=="scater"){  

        size_factors = scater::librarySizeFactors(matrix)
        matrix <- scater::normalizeCounts(as.matrix(matrix), size_factors = size_factors, return_log = FALSE)

    } else if (option=="Linnorm"){#It is not compatible with log transformed datasets. 

        matrix = expm1(Linnorm::Linnorm(as.matrix(matrix))) #Main function contains log1p(datamatrix)

    }

    return(matrix)

}




#################################################
##########    DECONVOLUTION METHODS    ##########
Deconvolution <- function(T, C, method, phenoDataC, P=NULL, elem = NULL, STRING = NULL, marker_distrib, refProfiles.var){ #marker subsetting performed right before for bulk methods!

    bulk_methods = c("CIBERSORT","DeconRNASeq","OLS","nnls","FARDEEP","RLR","DCQ","elastic_net","lasso","ridge","EPIC","DSA","ssKL","ssFrobenius","dtangle")
    sc_methods = c("MuSiC","BisqueRNA","DWLS","deconvSeq","SCDC")

    ########## Using marker information for bulk_methods
    if(method %in% bulk_methods){

        C = C[rownames(C) %in% marker_distrib$gene,]
        T = T[rownames(T) %in% marker_distrib$gene,]
        refProfiles.var = refProfiles.var[rownames(refProfiles.var) %in% marker_distrib$gene,]

    } else { ### For scRNA-seq methods 

        #BisqueRNA requires "SubjectName" in phenoDataC
        sample_column = grep("[S-s]ample|[S-s]ubject",colnames(phenoDataC))
        colnames(phenoDataC)[sample_column] = "SubjectName"
        rownames(phenoDataC) = phenoDataC$cellID

        require(xbioc)
        C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C),phenoData = Biobase::AnnotatedDataFrame(phenoDataC))
        T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T))

    }

    ##########    MATRIX DIMENSION APPROPRIATENESS    ##########
    keep = intersect(rownames(C),rownames(T)) 
    C = C[keep,]
    T = T[keep,]

    ###################################
    if(method=="CIBERSORT"){ #without QN. By default, CIBERSORT performed QN (only) on the mixture.

        RESULTS = CIBERSORT(sig_matrix = C, mixture_file = T, QN = FALSE) 
        RESULTS = t(RESULTS[,1:(ncol(RESULTS)-3)])

    } else if(method=="DeconRNASeq"){ #nonnegative quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
        #datasets and reference matrix: signatures, need to be non-negative. 
        #"use.scale": whether the data should be centered or scaled, default = TRUE
        unloadNamespace("Seurat") #needed for PCA step
        library(pcaMethods) #needed for DeconRNASeq to work
        RESULTS = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(T), signatures = as.data.frame(C), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
        colnames(RESULTS) = colnames(T)
        require(Seurat)

    } else if (method=="OLS"){

        RESULTS = apply(T,2,function(x) lm(x ~ as.matrix(C))$coefficients[-1])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))

    } else if (method=="nnls"){

        require(nnls)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- colnames(C)

    } else if (method=="FARDEEP"){

        require(FARDEEP)
        RESULTS = t(FARDEEP::fardeep(C, T, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method=="RLR"){ #RLR = robust linear regression

        require(MASS)
        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) MASS::rlm(x ~ as.matrix(C), maxit=100)), function(y) y$coefficients[-1]))
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))

    } else if (method=="DCQ"){#default: alpha = 0.05, lambda = 0.2. glmnet with standardize = TRUE by default

        require(ComICS)
        RESULTS = t(ComICS::dcq(reference_data = C, mix_data = T, marker_set = as.data.frame(row.names(C)) , alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 10)$average)
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method=="elastic_net"){#standardize = TRUE by default. Me: alpha = 0.2 (pg73, ESL book); lambda=NULL by default 

        require(glmnet)# gaussian is the default family option in the function glmnet. https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0.2, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method=="ridge"){ #alpha=0

        require(glmnet)
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 0, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method=="lasso"){ #alpha=1; shrinking some coefficients to 0. Deal with correlation/collinearity

        require(glmnet)
        RESULTS = apply(T, 2, function(z) coef(glmnet::glmnet(x = as.matrix(C), y = z, alpha = 1, standardize = TRUE, lambda = glmnet::cv.glmnet(as.matrix(C), z)$lambda.1se))[1:ncol(C)+1,])
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

        RESULTS[is.na(RESULTS)] <- 0 #Needed for models where glmnet drops all terms of a model and fit an intercept-only model (very unlikely but possible).

    } else if (method=="EPIC"){

        require(EPIC)
        marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),]
        markers = as.character(marker_distrib$gene)
        C_EPIC <- list()

        common_CTs <- intersect(colnames(C),colnames(refProfiles.var))

        C_EPIC[["sigGenes"]] <- rownames(C[markers,common_CTs])
        C_EPIC[["refProfiles"]] <- as.matrix(C[markers,common_CTs])
        C_EPIC[["refProfiles.var"]] <- refProfiles.var[markers,common_CTs]

        RESULTS <- t(EPIC::EPIC(bulk=as.matrix(T), reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices
        RESULTS = RESULTS[!rownames(RESULTS) %in% "otherCells",]

    } else if (method=="DSA"){ #DSA algorithm assumes that the input mixed data are in linear scale; If log = FALSE the data is left unchanged

        require(CellMix)
        md = marker_distrib
        ML = CellMix::MarkerList()
        ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)
        RESULTS = CellMix::ged(as.matrix(T), ML, method = "DSA", log = FALSE)@fit@H
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method=="ssKL"){ 

        require(CellMix)
        md = marker_distrib #Full version, irrespective of C
        ML = CellMix::MarkerList()
        ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

        RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssKL", sscale = FALSE, maxIter=500, log = FALSE)@fit@H 

    } else if (method=="ssFrobenius"){

        require(CellMix)
        md = marker_distrib #Full version, irrespective of C
        ML = CellMix::MarkerList()
        ML@.Data <- tapply(as.character(md$gene),as.character(md$CT),list)

        RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssFrobenius", sscale = TRUE, maxIter = 500, log = FALSE)@fit@H #equivalent to coef(CellMix::ged(T,...)

    } else if(method=="dtangle"){#Only works if T & C are log-transformed

        require(dtangle)
        mixture_samples = t(T)
        reference_samples = t(C)

        marker_distrib <- tidyr::separate_rows(marker_distrib,"CT",sep="\\|")    
        marker_distrib = marker_distrib[marker_distrib$gene %in% rownames(C),] 
        MD <- tapply(marker_distrib$gene,marker_distrib$CT,list)
        MD <- lapply(MD,function(x) sapply(x, function(y) which(y==rownames(C))))

        RESULTS = t(dtangle::dtangle(Y=mixture_samples, reference=reference_samples, markers = MD)$estimates)

    ###################################
    ###################################
    } else if (method == "MuSiC"){

        require(MuSiC)
        RESULTS = t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType',
                                            markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                                            verbose = F)$Est.prop.weighted)

    } else if (method == "DWLS"){
        
        require(DWLS)
        path=paste(getwd(),"/results_",STRING,sep="")

        if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created

            dir.create(path)
            Signature <- DWLS::buildSignatureMatrixMAST(scdata = C, id = as.character(phenoDataC$cellType), path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)

        } else {#re-load signature and remove CT column + its correspondent markers

            load(paste(path,"Sig.RData",sep="/"))
            Signature <- Sig
            
            if(!is.null(elem)){#to be able to deal with full C and with removed CT
                
                Signature = Signature[,!colnames(Signature) %in% elem]
                CT_to_read <- dir(path) %>% grep(paste(elem,".*RData",sep=""),.,value=TRUE)
                load(paste(path,CT_to_read,sep="/"))
            
                Signature <- Signature[!rownames(Signature) %in% cluster_lrTest.table$Gene,]

            }
            
        }
        
        RESULTS <- apply(T,2, function(x){
            b = setNames(x, rownames(T))
            tr <- DWLS::trimData(Signature, b)
            RES <- t(DWLS::solveDampenedWLS(tr$sig, tr$bulk))
        })

        rownames(RESULTS) <- as.character(unique(phenoDataC$cellType))
        RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
        RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint

    } else if (method == "BisqueRNA"){#By default, Bisque uses all genes for decomposition. However, you may supply a list of genes (such as marker genes) to be used with the markers parameter

        require(BisqueRNA)
        RESULTS <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers=NULL, use.overlap=FALSE)$bulk.props #use.overlap is when there's both bulk and scRNA-seq for the same set of samples

    } else if (method == "deconvSeq"){
      
        singlecelldata = C.eset 
        celltypes.sc = as.character(phenoDataC$cellType) #To avoid "Design matrix not of full rank" when removing 1 CT 
        tissuedata = T.eset 

        design.singlecell = model.matrix(~ -1 + as.factor(celltypes.sc))
        colnames(design.singlecell) = levels(as.factor(celltypes.sc))
        rownames(design.singlecell) = colnames(singlecelldata)

        dge.singlecell = deconvSeq::getdge(singlecelldata,design.singlecell, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")
        b0.singlecell = deconvSeq::getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 4)
        dge.tissue = deconvSeq::getdge(tissuedata, NULL, ncpm.min = 1, nsamp.min = 4, method = "bin.loess")

        RESULTS = t(deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0.singlecell, dge.tissue)$x1) #genes with adjusted p-values <0.05 after FDR correction

    } else if (method == "SCDC"){ ##Proportion estimation with traditional deconvolution + >1 subject

        require(SCDC)
        RESULTS <- t(SCDC::SCDC_prop(bulk.eset = T.eset, sc.eset = C.eset, ct.varname = "cellType", sample = "SubjectName", ct.sub = unique(as.character(phenoDataC$cellType)), iter.max = 200)$prop.est.mvw)

    }

    RESULTS = RESULTS[gtools::mixedsort(rownames(RESULTS)),]
    RESULTS = data.table::melt(RESULTS)

    if(!is.null(P)){

        P = P[gtools::mixedsort(rownames(P)),]
        P$ID = rownames(P)
        P = data.table::melt(P, id.vars="ID")

        RESULTS$expected_values <- P$value
        colnames(RESULTS)[1:3] <- c("CT","tissue","observed_values")
        RESULTS$observed_values <- round(RESULTS$observed_values,3)

    }

    return(RESULTS) 

}