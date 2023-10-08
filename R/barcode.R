#' use the compare_entities method to assess 'closeness' of embedded elements
#' @param simbaref instance of python module for simba
#' @param refemb AnnData instance with 'reference' embedding
#' @param queryemb AnnData instance with 'query' embedding
#' @param n_top integer(1) number of 'reference' elements to
#' be examined for 'max' output summarizing closeness (maximum of normed dot product)
#' @param temperature numeric(1), a tuning parameter whose reciprocal scales softmax
#' transformations of dot products between reference and query entities, defaulting
#' to 1.0 for this application 
#' @examples
#'   simbaref = simba_ref()
#'   cemb = simbaref$read_h5ad( get_3k_cell_emb() )
#'   gemb = simbaref$read_h5ad( get_3k_gene_emb() )
#'   cc = compare_entities( simbaref, cemb, gemb )
#'   cc
#'   head(cc$var)
#' @export
compare_entities = function( simbaref, refemb, queryemb,
   n_top = 50L, temperature = 1.0 ) {
   simbaref$tl$compare_entities(cemb, gemb, n_top = as.integer(n_top),
      T = temperature)
}

#' 'barcode' data for figure 2d
#' @param symbol character(1) gene symbol for barcode visualization
#' @param origad AnnData reference to original RNA-seq data
#' @param clabel character(1) column name for entity type in origad
#' @param simbaref instance of python module for simba
#' @param cembpath character(1) path to h5ad file for cell embedding
#' @param gembpath character(1) path to h5ad file for gene embedding
#' @param layer character(1) either 'softmax' (default) or 'norm', used
#' as source of scoring of strength of cell-gene association
#' @examples
#' sref = simba_ref()
#' p3k = sref$read_h5ad(get_10x3kpbmc_path(overwrite=TRUE))
#' myd = simba_barcode_data( symbol = "MS4A1",
#'    origad = p3k, clabel = "celltype",
#'    simbaref = sref, cembpath = get_3k_cell_emb(),
#'    gembpath = get_3k_gene_emb() )
#' head(myd)
#' @export
simba_barcode_data = function(symbol, origad, clabel, simbaref,
   cembpath, gembpath, layer="softmax") {
   cemb = simbaref$read_h5ad( cembpath )
   gemb = simbaref$read_h5ad( gembpath )
   cemb$obs$tmp = origad[cemb$obs_names]$obs[[clabel]]
   names(cemb$obs) = clabel
   cmp <- simbaref$tl$compare_entities(cemb, gemb )
   scos = cmp$layers[layer]
   ind = which(rownames(cmp$var)==symbol)
   if (length(ind)==0) stop(sprintf("%s not present in compare_entities $var result",
               symbol))
   oog = order(scos[,ind], decreasing=TRUE)
   data.frame(rank=seq_len(nrow(scos)), score=sort(scos[,ind], decreasing=TRUE),
      attr=cmp$obs[[clabel]][oog])
   #plot(sort(sm[,456]), pch=19, col=factor(cmp$obs$celltype[oog]))
}

