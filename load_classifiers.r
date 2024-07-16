# Format the classifiers of the study to be handled in adverSCarial functions
# This file is loaded in the run_attacks.ipynb notebook

CHETAH_classifier <- function(expr, clusters, target){
    library(CHETAH)
    if ( !exists("reference_3k")){
        train_3k <- readRDS("repr_data/sce_pbmc_train.rds")
    }
    input <- CHETAHclassifier(input = expr, ref_cells = train_3k, ref_ct = 'cell_type')
    input <- Classify(input = input, 0.01)
    final_predictions = input$celltype_CHETAH[clusters == target]
    ratio <- as.numeric(sort(table(final_predictions), decreasing = TRUE)[1]) /
        sum(as.numeric(sort(table(final_predictions), decreasing = TRUE)))
    predicted_class <- names(sort(table(final_predictions), decreasing = TRUE)[1])
    resCHETAH <- list(prediction=predicted_class, odd=ratio)
    return(resCHETAH)
}

scMLP_classifier <- function(expr, clusters, target){
    expr = as.matrix(expr)
    library(reticulate)
    use_python("/usr/bin/python3", required = TRUE)
    library(keras)
    mlpModel <- load_model_hdf5("repr_data/classifiers/scMLP//dl_model.h5")
    newColnames <- read.table("repr_data/classifiers/scMLP/new_colnames.txt")$V1
    predictions <- predict(mlpModel, expr)
    colnames(predictions) <- newColnames
    rownames(predictions) <- rownames(expr)
    predictions <- as.data.frame(predictions)
    if (sum(clusters == target) == 0 ){
        return( c("UNDETERMINED",1))
    }
    cell_types <- apply(predictions[clusters == target,], 1, function(x){
        names(x[x == max(x)])[1]
    })
    table_cell_type <<- table(cell_types)
    str_class <- names(table_cell_type[order(table_cell_type, decreasing=T)][1])
    resSCMLP <- list(prediction=str_class,
                     odd=1,
                     typePredictions=as.data.frame(t(predictions)),
                     cellTypes=cell_types)
    return(resSCMLP)
}

scType_classifier = function(expr, clusters, target){
    expr = t(expr)
    if ( !exists("sctype_score")){
        library(HGNChelper)
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
        sctype_score <<- sctype_score
    }
    if (!exists("gs_list")){
        library(HGNChelper)
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
        db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
        tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
        # prepare gene sets
        gs_list <<- gene_sets_prepare(db_, tissue)
    }
    es.max = sctype_score(scRNAseqData = expr, scaled = T, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    if (sum(clusters == target) == 0 ){
        return( c("UNDETERMINED",1))
    }
    cell_types <- apply(t(es.max[, clusters == target]), 1, function(x){
        names(x[x == max(x)])[1]
    })
    table_cell_type <<- table(cell_types)
    str_class <- names(table_cell_type[order(table_cell_type, decreasing=T)][1])
    resSCtype <- list(prediction=str_class,
                      odd=1,
                      typePredictions=es.max,
                      cellTypes=cell_types)
    return(resSCtype)
}

scAnnotatR_classifier = function(expr, clusters, target){
    if (!"scAnnotatR" %in% loadedNamespaces()){
        library(scAnnotatR)
    }
    if ( !exists("sca_default_models")){
        sca_default_models <<- load_models("repr_data/classifiers/scAnnotatR/trainedModels/")
        sca_cell_types <<- names(sca_default_models)
    }
    seurat.obj <- classify_cells(classify_obj = expr, 
                             assay = 'RNA', slot = 'scale.data',
                             cell_types = sca_cell_types, 
                             path_to_models = "repr_data/classifiers/scAnnotatR/trainedModels/")
    typePredictions <- seurat.obj@meta.data[, stringr::str_replace_all( paste0(sca_cell_types, "_p"), " ", "_")]
    pred_table <- table(seurat.obj@meta.data[clusters == target, "most_probable_cell_type"])
    cluster_pred <- names(pred_table[order(pred_table, decreasing = TRUE)])[1]
    odd <- unname(pred_table[cluster_pred]/sum(pred_table))
    colnames(typePredictions) <- unlist(lapply(colnames(typePredictions), function(x){
        unlist(strsplit(x,"_p$"))[1]
    }))
    typePredictions <- as.data.frame(t(typePredictions))
    result <- list(prediction=stringr::str_replace_all(cluster_pred," ","_"),
                   odd=odd, 
                   typePredictions=typePredictions,
                   cellTypes=seurat.obj@meta.data$most_probable_cell_type)
    return(result)
}

scRF_classifier <- function(expr, clusters, target){
    if ( !exists("rfModel") || 0==0){
        message("load model")
        library(randomForest)
        rfModel <- readRDS("repr_data/classifiers/scRF/random_forest_model.rds")
    }
    predictions <- predict(rfModel, expr, type="prob")
    if (sum(clusters == target) == 0 ){
        return( c("UNDETERMINED",1))
    }
    cell_types <- apply(predictions[clusters == target,], 1, function(x){
        names(x[x == max(x)])[1]
    })
    table_cell_type <<- table(cell_types)
    str_class <- names(table_cell_type[order(table_cell_type, decreasing=T)][1])
    resSCtype <- list(prediction=str_class,
                      odd=1,
                      typePredictions=as.data.frame(t(predictions)),
                      cellTypes=cell_types)
    return(resSCtype)
}

# We store each classifier in a list with associated parameters:
# - the input type needed => SingleCellExperiment, data.frame or seurat
# - which slot of data to use => "raw", "scale.data", NULL
classifList <- list(CHETAH=list(inputType="SingleCellExperiment", fct=CHETAH_classifier),
                    scMLP=list(inputType="data.frame", fct=MLP_classifier),
                    scRF=list(inputType="data.frame", fct=RF_classifier),
                    scType=list(inputType="data.frame", fct=scType_classifier),
                    scAnnotatR=list(inputType="seurat", slot="scale.data", fct=scAnnotatR_classifier)
                    )






