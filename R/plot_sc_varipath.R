## plot_sc_varipath() ##

# packages necessary: ggplot2, ggtree, cowplot, patchwork, aplot, igraph, dplyr, tidyr, viridis

#' plot_sc_varipath()
#' @description Plot results from sc_varipath() function to identify pathways associated with perturbations.
#' @import dplyr
#' @import ggplot2
#' @import ggtree
#' @import cowplot
#' @import patchwork
#' @import aplot
#'
#' @param sc_varipath_out The output list from the sc_varipath() function.
#' @param type Type of plot to generate. Options are 'perturbation_to_RC', 'RC_to_gene', or 'RC_to_pathway'.
#'
#' @returns A dotplot with dendrograms showing associations between perturbations and rotated components, rotated components and genes, or rotated components and pathways.
#' @details This function visualizes the relationships identified by the sc_varipath() function. Depending on the specified type, it generates a dotplot with dendrograms to illustrate the associations between perturbations and rotated components, rotated components and genes, or rotated components and pathways. The size and color of the dots represent the strength of the associations.
#' @examples
#' out = sc_varipath(seurat_obj = adata.R, perturbation_column = 'miR.family')
#'
#' options(repr.plot.width=20, repr.plot.height=10)
#' plot_sc_varipath(sc_varipath_out = out, type = 'perturbation_to_RC')
#' @export

plot_sc_varipath <- function(sc_varipath_out, type = NULL){
    if (is.null(type)){
        stop("Please specify a type of plot to generate: 'perturbation_to_RC', 'RC_to_gene', or 'RC_to_pathway'")
    }
    if (type == 'perturbation_to_RC'){
        FDR = sc_varipath_out$FDR
        gene_cluster = sc_varipath_out$perturbation_RC_df
        colnames(gene_cluster) = c('cluster', 'Gene', 'count')

        # make data square to calculate euclidean distance
        mat <- gene_cluster %>% # drop unused columns to faciliate widening
        tidyr::pivot_wider(names_from = cluster, values_from = count) %>%
        data.frame() # make df as tibbles -> matrix annoying
        row.names(mat) <- mat$Gene  # put gene in `row`
        mat <- mat[,-1] #drop gene column as now in rows
        clust <- hclust(dist(mat %>% as.matrix())) 
        v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
        ############ NOTICE THE t() above)

        ddgram_col <- as.dendrogram(v_clust)
        ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

        dotplot <- gene_cluster %>% 
        mutate(Gene = factor(Gene, levels = clust$labels[clust$order]),
                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%
                filter(count > (-log10(FDR))) %>%
        ggplot(aes(x=cluster, y = Gene, color = count, size = count)) + 
        geom_point() + 
        cowplot::theme_cowplot() + 
        theme(axis.line  = element_blank()) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab('') +
        theme(axis.ticks = element_blank()) +
        scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = '-Log10(P-value)') +
        scale_size_continuous(name = '-Log10(P-value)', range = c(4,10)) +
        scale_y_discrete(position = "right")
        #################################################
        ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)


        ddgram <- as.dendrogram(clust) # create dendrogram
        ggtree_plot <- ggtree(ddgram)
        ggtree_plot <- ggtree_plot + ylim2(dotplot)

        labels <- ggplot(gene_cluster %>% 
                        mutate(`Cell Type` = rep('A', nrow(gene_cluster)),
                                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
                        aes(x = cluster, y = 1, fill = `Cell Type`)) + 
        geom_tile() + 
        scale_fill_brewer(palette = 'Set1') + 
        theme_nothing() +
        xlim2(dotplot)

        legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))

        plot = plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + labels +
        plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + dotplot + 
        plot_spacer() + plot_spacer() + legend + 
        plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1))

        return(plot)
        
    }


    if (type == 'RC_to_gene'){
        significant_RCs = sc_varipath_out$significant_RCs
        gene_cluster = sc_varipath_out$weights_dataframe
        colnames(gene_cluster) = c('Gene', 'cluster', 'count')
        # make data square to calculate euclidean distance
        mat <- gene_cluster %>% # drop unused columns to faciliate widening
        tidyr::pivot_wider(names_from = cluster, values_from = count) %>%
        data.frame() # make df as tibbles -> matrix annoying
        row.names(mat) <- mat$Gene  # put gene in `row`
        mat <- mat[,-1] #drop gene column as now in rows
        clust <- hclust(dist(mat %>% as.matrix())) 
        v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
        ############ NOTICE THE t() above)

        ddgram_col <- as.dendrogram(v_clust)
        ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

        dotplot <- gene_cluster %>% 
        mutate(Gene = factor(Gene, levels = clust$labels[clust$order]),
                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%
                filter(cluster %in% significant_RCs) %>% filter(count > 0 | count < 0) %>%
        ggplot(aes(x=cluster, y = Gene, color = count, size = count)) + 
        geom_point() + 
        cowplot::theme_cowplot() + 
        theme(axis.line  = element_blank()) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 6)) +
        ylab('') +
        theme(axis.ticks = element_blank()) +
        scale_color_gradientn(guide = 'legend', name = 'Loading', colours = viridis::viridis(20), limits = c(round(min(gene_cluster$count), 1), round(max(gene_cluster$count), 1)), oob = scales::squish, breaks = round(seq(min(gene_cluster$count), max(gene_cluster$count), length.out = 4), 1)) +
        scale_size_continuous(guide = 'legend', name = 'Loading', range = c(0,5), limits = c(round(min(gene_cluster$count), 1), round(max(gene_cluster$count), 1)), breaks = round(seq(min(gene_cluster$count), max(gene_cluster$count), length.out = 4), 1)) +
        scale_y_discrete(position = "right")
        #################################################
        ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)


        ddgram <- as.dendrogram(clust) # create dendrogram
        ggtree_plot <- ggtree::ggtree(ddgram)
        ggtree_plot <- ggtree_plot + ylim2(dotplot)

        labels <- ggplot(gene_cluster %>% 
                        mutate(`Cell Type` = rep('A', nrow(gene_cluster)),
                                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
                        aes(x = cluster, y = 1, fill = `Cell Type`)) + 
        geom_tile() + 
        scale_fill_brewer(palette = 'Set1') + 
        theme_nothing() +
        xlim2(dotplot)

        legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))

        plot = plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + labels +
        plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + dotplot + 
        plot_spacer() + plot_spacer() + legend + 
        plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1))

        return(plot)
        
    }

    if (type == 'RC_to_pathway'){
        significant_RCs = sc_varipath_out$significant_RCs
        gene_cluster = sc_varipath_out$adjacency_table
        colnames(gene_cluster) = c('Gene', 'cluster', 'count')
        mygraph <- igraph::graph_from_data_frame(gene_cluster)
        mat <- as.matrix(igraph::get.adjacency(mygraph, sparse = FALSE, attr='count'))
        mat <- select(as.data.frame(mat), starts_with("RC")) %>% filter(!startsWith(row.names(mat), "RC"))

        clust <- hclust(dist(mat %>% as.matrix())) 
        v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
        ############ NOTICE THE t() above)

        ddgram_col <- as.dendrogram(v_clust)
        ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

        dotplot <- gene_cluster %>% 
        mutate(Gene = factor(Gene, levels = clust$labels[clust$order]),
                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])) %>%
                filter(cluster %in% significant_RCs) %>% filter(count > 0) %>%
        ggplot(aes(x=cluster, y = Gene, color = count, size = count)) + 
        geom_point() + 
        cowplot::theme_cowplot() + 
        theme(axis.line  = element_blank()) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14), axis.text.y = element_text(size = 14)) +
        ylab('') +
        theme(axis.ticks = element_blank()) +
        scale_color_gradientn(colours = viridis::viridis(20), limits = c(0.0001, max(gene_cluster$count)), oob = scales::squish, name = 'Enrichr Score') +
        scale_size_continuous(name = 'Enrichr Score', range = c(0.2,5), limits = c(0.0001, max(gene_cluster$count))) +
        scale_y_discrete(position = "right")
        #################################################
        ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)


        ddgram <- as.dendrogram(clust) # create dendrogram
        ggtree_plot <- ggtree::ggtree(ddgram)
        ggtree_plot <- ggtree_plot + ylim2(dotplot)

        labels <- ggplot(gene_cluster %>% 
                        mutate(`Cell Type` = rep('A', nrow(gene_cluster)),
                                cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
                        aes(x = cluster, y = 1, fill = `Cell Type`)) + 
        geom_tile() + 
        scale_fill_brewer(palette = 'Set1') + 
        theme_nothing() +
        xlim2(dotplot)

        legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))

        plot = plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + labels +
        plot_spacer() + plot_spacer() + plot_spacer() +
        plot_spacer() + plot_spacer() + dotplot + 
        plot_spacer() + plot_spacer() + legend + 
        plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1))

        return(plot)

    }


}