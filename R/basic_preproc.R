#' Transform an AnnData input to one that has had gene filtering,
#' normalization, log transformation, low-variance gene exclusion,
#' and discretization applied.  
#' @param ad AnnData instance (i.e., inherits from `anndata._core.anndata.AnnData`)
#' @param workdir character(1) defines working directory for graph serialization; set to NULL to take simba default
#' @param min_n_cells integer(1) used to filter genes, excluding those expressed in fewer than `min_n_cells` cells
#' @param norm_method character(1) defaults to 'lib_size', see simba doc
#' @param do_log logical(1) if TRUE, conducts log normalization
#' @param n_top_genes NULL or numeric(1) if non-null, genes are excluded if
#' variance is smaller than that obtained for the "nth" top
#' @param n_bins integer(1) number of bins for expression discretization
#' @param simba_ref instance of python.builtin.module, checked to have component 'tl'
#' @return instance of AnnData with layers raw and simba, the latter consisting
#' of a sparse matrix, and other components as determined by the operations
#' selected through argument bindings.
#' @examples
#' p3k = get_10x3kpbmc_path(overwrite=TRUE)
#' ref = simba_ref()
#' pp = ref$read_h5ad(p3k)
#' pp
#' bb = basic_preproc(pp, simba_ref=ref)
#' bb
#' @export
basic_preproc = function(ad, workdir=tempdir(),
   min_n_cells = 3L, norm_method = 'lib_size',
   do_log = TRUE, n_top_genes = NULL, n_bins=6L, simba_ref ) {
 stopifnot(inherits(simba_ref, "python.builtin.module"),
    "tl" %in% names(simba_ref))
 stopifnot(inherits(ad, "anndata._core.anndata.AnnData"))
 if (!is.null(workdir)) simba_ref$settings$set_workdir(workdir)
 if (!is.null(min_n_cells)) simba_ref$pp$filter_genes(ad, min_n_cells=
      as.integer(min_n_cells))
 simba_ref$pp$cal_qc_rna(ad)
 simba_ref$pp$normalize(ad, method=norm_method)
 if (do_log) simba_ref$pp$log_transform(ad)
 if (!is.null(n_top_genes)) {
    simba_ref$pp$select_variable_genes(ad, n_top_genes = as.integer(n_top_genes))
    }
 simba_ref$tl$discretize(ad, n_bins = as.integer(n_bins) )
 ad
}
 
 
 
