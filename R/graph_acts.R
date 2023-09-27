#' build and train SIMBA PBG graph
#' @param adSimb an instance of AnnData with simba layer as
#' produced by basic_preproc
#' @param copy logical(1) passed to gen_graph in simba.tl module, defaults
#' to FALSE
#' @param use_hvg logical(1) passed as `use_highly_variable` to gen_graph, defaults to FALSE
#' @param dirname character(1) defaults to 'graph0'
#' @param auto_wd logical(1) defaults to TRUE, for automatic setting of weight decay
#' @param save_wd logical(1) defaults to TRUE
#' @param output character(1) defaults to 'model'
#' @param simba_ref instance of python.builtin.module, checked to have component 'tl'
#' @return NULL
#' @examples
#' p15 = get_paul15(overwrite=TRUE) # allow repetition
#' ref = simba_ref()
#' pp = ref$read_h5ad(p15)
#' bb = basic_preproc(pp, simba_ref=ref)
#' gg = build_and_train_pbg( bb, simba_ref=ref )
#' dict = ref$read_embedding()
#' dict 
#' @export
build_and_train_pbg = function( adSimb, copy=FALSE,
  use_hvg = FALSE, dirname = 'graph0', auto_wd =  TRUE,
  save_wd = TRUE, output = 'model', simba_ref ) {
 stopifnot(inherits(simba_ref, "python.builtin.module"),
    "tl" %in% names(simba_ref))
 stopifnot(inherits(adSimb, "anndata._core.anndata.AnnData"))
 stopifnot(inherits(adSimb$layers['simba'], "dgRMatrix"))
 simba_ref$tl$gen_graph(list_CG = list(adSimb), copy=copy,
    use_highly_variable = use_hvg, dirname = dirname)
 simba_ref$tl$pbg_train(auto_wd=auto_wd, save_wd=save_wd,
    output=output)
}
