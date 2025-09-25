## Input Seurat object with PCA already run

# packages necessary: Seurat, utils, enrichR, tidyr, dplyr


## sc_varipath() ##

#' sc_varipath()
#' @description Identify novel pathway shifts associated with perturbations in single-cell RNA-seq data using a varimax rotation of scRNA-seq PCA.
#' @import Seurat
#' @import utils
#' @import dplyr
#' @import enrichR
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data with PCA already run.
#' @param perturbation_column The name of the metadata column in the Seurat object containing cell perturbation/guide information.
#' @param enrichR_dbs A character vector of enrichR database names to use for pathway enrichment analysis. Default is 'GO_Biological_Process_2021'.
#' @param loading_cutoff A numeric value indicating the loading cutoff to use for selecting genes associated with each rotated component. Default is 0.1.
#' @param FDR A numeric value indicating the FDR threshold to use for selecting significant components based on guide association p-values. Default is 0.25.
#'
#'
#' It is advised that you save the output of this function as an RDS file using saveRDS() immediately after running it for easy loading in future sessions and because memory failure can occur when running downstream plotting functions on large datasets.
#' @returns A list with the following components:
#' \itemize{
#'  \item{vari_weights}{Varimax weights of genes for each rotated component.}
#'  \item{fitted_models}{Linear models fitting perturbation to each rotated component.}
#'  \item{fit_pvals_adj}{Adjusted p-values of the association between perturbation and each rotated component to identify best fit for each perturbation.}
#'  \item{perturbation_RC_df}{Adjusted p-values of the association between perturbation and each rotated component to plot using plot_sc_varipath().}
#'  \item{weights_dataframe}{Rotated component gene loadings to be used for plot_sc_varipath().}
#'  \item{adjacency_table}{Contains 3 most associated pathways for each rotated component to plot using plot_sc_varipath().}
#' }
#'
#' @examples
#' out = sc_varipath(seurat_obj = adata.R, perturbation_column = 'miR.family')
#' @export





sc_varipath <- function(seurat_obj, perturbation_column = NULL, enrichR_dbs = 'GO_Biological_Process_2021', loading_cutoff = 0.1, FDR = 0.25) {
    
    if(is.null(perturbation_column)){
        stop('Error: Please provide the name of the metadata column containing perturbation/guide information using the perturbation_column parameter.')
    }
    message('Please adjust loading_cutoff parameter as needed after plotting result of this function using sc_varipth_plot(type = "RC_to_gene"). Default is 0.1.')

    ## check if PCA has been run and embeddings exist ##
    result <- tryCatch({Embeddings(seurat_obj, reduction = 'pca')}, error = function(e) {TRUE})
    # This function is executed if an error occurs and returns TRUE if an error occurred

    if(is.logical(result)){
        if(result == TRUE){
            stop('Error: PCA embeddings not found in Seurat object. Please run RunPCA() on the Seurat object before running this function.')
        }
    }
    ## done ##

    meta = seurat_obj@meta.data
    meta$guide = meta[[perturbation_column]]

    scaled_scores <- scale(Embeddings(seurat_obj, reduction = "pca"))
    var_rot <- varimax(Loadings(seurat_obj, reduction = "pca"))$rotmat
    rotated_varimax_scores <- scaled_scores %*% var_rot


    vari_weights = varimax(Loadings(seurat_obj, reduction = "pca"))
    vari_weights = vari_weights$loadings

    ###### guide/perturbation to RC ######
    ## build linear model between rotated scores and cell guide identity, pick best fit model for each RC from log likelihood ratio test p-value

    uniq = unique(meta$guide)
    fitted_models = list()
    fit_pvals = list()
    pb = utils::txtProgressBar(min = 0, max = ncol(rotated_varimax_scores), style = 3)
    for (j in 1:ncol(rotated_varimax_scores)){
        fitted_models[[j]] <- suppressWarnings(lapply(as.list(uniq), function(x) glm(identity ~ scores, data = data.frame(scores = rotated_varimax_scores[meta$guide == x | meta$guide == 'TuD_NC',j], identity = factor(meta$guide[meta$guide == x | meta$guide == 'TuD_NC'] == x)), family = 'binomial')))
        fit_pvals[[j]] = lapply(fitted_models[[j]], function(x) with(x, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)))
        utils::setTxtProgressBar(pb, j)
    }
    close(pb)

    fit_pvals_adj = lapply(fit_pvals, function(x) p.adjust(as.numeric(unlist(x)), method = 'BH'))
    fit_pvals_adj = lapply(fit_pvals_adj, function(x) -log10(x))

    ## guide to RC association pvals
    uniq = unique(meta$guide)
    df = data.frame(matrix(unlist(fit_pvals_adj), nrow = length(fit_pvals_adj), byrow = T))
    colnames(df) = uniq
    df$y = paste0('RC', 1:ncol(rotated_varimax_scores))
    df_long <- tidyr::pivot_longer(df, cols = -c(y), names_to = 'x', values_to = 'value')
    d <- as.data.frame(df_long)
    perturbation_RC_df = d

    ### return perturb_RC_df (d2) for plot A (perturbation_to_RC) ----------------------------------------

    ###### RC to gene weights ######

    ## replace with zero all below cutoffs
    zero_below_thresh <- function(df, cutoffs){
        df = as.matrix(df)
        for (i in 1:ncol(df)){
        df[,i][abs(df[,i]) < cutoffs[i]] = 0
        colnames(df)[i] = paste0('RC', i)
        }
        return(df)
    }

    vari_weights_abv_threshold = zero_below_thresh(vari_weights, rep(loading_cutoff, ncol(vari_weights))) 

    ## genes with loadings above loading_cutoff kept for each RC for pathway analysis
    vari_weights_abv_threshold = data.frame(matrix(as.numeric(vari_weights_abv_threshold), attributes(vari_weights_abv_threshold)$dim, dimnames=attributes(vari_weights_abv_threshold)$dimnames))
    vari_weights_abv_threshold$y = rownames(vari_weights_abv_threshold)
    weights_dataframe <- tidyr::pivot_longer(vari_weights_abv_threshold, cols = -c(y), names_to = 'x', values_to = 'value')
    weights_dataframe <- as.data.frame(weights_dataframe)

    ### return weights_dataframe (weights_df2) for plot B (RC_to_gene) ----------------------------------------

    ###### gene weights to enrichR pathways ######
    gene_cluster = perturbation_RC_df
    colnames(gene_cluster) = c('cluster', 'Gene', 'count')
    # make data square to calculate euclidean distance
    mat <- gene_cluster %>% # drop unused columns to faciliate widening
    tidyr::pivot_wider(names_from = cluster, values_from = count) %>%
    data.frame() # make df as tibbles -> matrix annoying
    row.names(mat) <- mat$Gene  # put gene in `row`
    mat <- mat[,-1] #drop gene column as now in rows
    clust <- hclust(dist(mat %>% as.matrix())) 
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix

    significant_RCs = gene_cluster %>% 
    mutate(Gene = factor(Gene, levels = clust$labels[clust$order]),
            cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%
            filter(count > (-log10(FDR))) %>% dplyr::select(cluster) %>% unique()

    significant_RCs = as.character(significant_RCs$cluster)

    ## get significant RCs from guide association analysis
    gene_cluster = weights_dataframe
    colnames(gene_cluster) = c('Gene', 'cluster', 'count')
    # make data square to calculate euclidean distance
    mat <- gene_cluster %>% # drop unused columns to faciliate widening
    tidyr::pivot_wider(names_from = cluster, values_from = count) %>%
    data.frame() # make df as tibbles -> matrix annoying
    row.names(mat) <- mat$Gene  # put gene in `row`
    mat <- mat[,-1] #drop gene column as now in rows
    clust <- hclust(dist(mat %>% as.matrix())) 
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix

    weights_dataframe_sig_genes = gene_cluster %>% 
    dplyr::mutate(Gene = factor(Gene, levels = clust$labels[clust$order]),
            cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%
            dplyr::filter(cluster %in% significant_RCs) %>% dplyr::filter(abs(count) > 0) %>% dplyr::arrange(cluster) %>% as.data.frame()

    setEnrichrSite("Enrichr")
    dbs = enrichR_dbs


    list_convert <- function(n, df) {
        clusters = unique(df$cluster)
        filtered_df = df %>% dplyr::filter(cluster == clusters[n])
        return(as.character(filtered_df$Gene))
    }

    significant_genes_abv_load_cutoff = lapply(1:length(unique(weights_dataframe_sig_genes$cluster)), list_convert, weights_dataframe_sig_genes)
    names(significant_genes_abv_load_cutoff) = unique(weights_dataframe_sig_genes$cluster)

    ## there is a problem with the enrichr package, it requries a delay between queries to prevent repeat results for earlier inputs
    ## make sure none are identical -- make sure all false at bottom of output of this cell
    enriched_list = list()

    for (i in 1:length(significant_genes_abv_load_cutoff)) {
        enriched_list[[i]] = enrichr(as.character(significant_genes_abv_load_cutoff[[i]]), dbs)
        Sys.sleep(1)
    }
    enriched_list = lapply(enriched_list, '[[', 1)
    fail_condition = sapply(enriched_list[-1], function(x) identical(x,enriched_list[[1]]))

    if(all(fail_condition)){
        stop('Error: enrichr results are identical for multiple inputs, please rerun the code block. This is a known issue with the enrichR package and requires a delay between queries to prevent repeat results for earlier inputs. Please increase the Sys.sleep time if this error persists.')
    }


    ## get top 3 pathways for each RC
    Terms = lapply(enriched_list, function(x) x %>% dplyr::select(Term, Combined.Score) %>% dplyr::arrange(desc(Combined.Score)) %>% dplyr::pull(Term) %>% head(3))
    Scores = lapply(enriched_list, function(x) x %>% dplyr::select(Term, Combined.Score) %>% dplyr::arrange(desc(Combined.Score)) %>% dplyr::pull(Combined.Score) %>% head(3) %>% as.numeric())
    adjacency_table = data.frame(pathways = unlist(Terms), components = as.character(sapply(names(significant_genes_abv_load_cutoff), rep, each = 3)), scores = unlist(Scores))

    ### adjacency_table return for plot C (RC_to_pathway) ----------------------------------------

    return(list(vari_weights = vari_weights, fitted_models = fitted_models, fit_pvals_adj = fit_pvals_adj, perturbation_RC_df = perturbation_RC_df, weights_dataframe = weights_dataframe, adjacency_table = adjacency_table))
    
}