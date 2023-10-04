#' build and train SIMBA PBG graph
#' @param adSimb an instance of AnnData with simba layer as
#' produced by basic_preproc
#' @param copy logical(1) passed to gen_graph in simba.tl module, defaults
#' to FALSE
#' @param use_hvg logical(1) passed as `use_highly_variable` to gen_graph, defaults to FALSE
#' @param dirname character(1) defaults to 'graph0'; training stats will be saved
#' to `paste0(simba_ref$settings$workdir, "/", dirname, "/", model, "/", "training_stats.json")`
#' @param auto_wd logical(1) defaults to TRUE, for automatic setting of weight decay
#' @param save_wd logical(1) defaults to TRUE
#' @param output character(1) defaults to 'model'
#' @param simba_ref instance of python.builtin.module, checked to have component 'tl'
#' @return NULL; components of `simba_ref` are updated to permit retrieval
#' of embedding.
#' @examples
#' p3k = get_10x3kpbmc_path(overwrite=TRUE)
#' ref = simba_ref()
#' pp = ref$read_h5ad(p3k)
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

#' Helper function to recursively flatten nested lists
#' produced by Chat-GPT 3
#' @param x list as produced by fromJSON
#' @param prefix character(1)
flatten_json <- function(x, prefix = "") {
  out <- list()
  for (key in names(x)) {
    new_key <- if (prefix == "") key else paste(prefix, key, sep = "_")
    if (is.list(x[[key]])) {
      out <- c(out, flatten_json(x[[key]], new_key))
    } else {
      out[[new_key]] <- x[[key]]
    }
  }
  return(out)
}

#' get pbg_metrics data frames
#' @importFrom jsonlite fromJSON
#' @param simba_ref reference to simba module
#' @param gg_dirname character(1) the dirname specified for build_and_train_pbg,
#' defaults to 'graph0'
#' @param tr_output character(1) the output specified for build_and_train_pbg
#' defaults to 'model'
#' @param jsonpath character(1) path to training_stats.json
#' @return a list with two data.frames of training statistics
#' @examples
#' # full run-based, commented out as too long
#' # p3k = get_10x3kpbmc_path(overwrite=TRUE) # allow repetition
#' # ref = simba_ref()
#' # pp = ref$read_h5ad(p3k)
#' # bb = basic_preproc(pp, simba_ref=ref)
#' # gg = build_and_train_pbg( bb, simba_ref=ref )
#' # ts = ingest_training_stats(ref)
#' # head(ts$df1)
#' #
#' # use archived pbg output
#' #
#' tpath = system.file(file.path("extdata", "pbg3k.tar.xz"), package="Simba4Bioc")
#' untar(tpath, exdir = tempdir())
#' jsonpath = paste0(tempdir(), "/pbg/graph0/model/training_stats.json")
#' ts = ingest_training_stats(jsonpath=jsonpath)
#' head(ts$df1)
#' @export
ingest_training_stats = function(simba_ref,
  gg_dirname = 'graph0', tr_output = 'model', jsonpath) {
  if (missing(jsonpath)) jsonpath = sprintf(
     paste0(simba_ref$settings$workdir, "/pbg/%s/%s/training_stats.json"),
     gg_dirname, tr_output)
  txt = readLines(jsonpath)
  parsed_stats = lapply(txt, jsonlite::fromJSON)
  # seems to be alternating document schemata
  ntxt = length(parsed_stats)
  odd = seq(1, ntxt, 2)
  ev = seq(2, ntxt, 2)
  df1 = do.call(rbind, lapply(parsed_stats[odd], 
     function(x) data.frame(flatten_json(x)))) 
  df2 = do.call(rbind, lapply(parsed_stats[ev], 
     function(x) data.frame(flatten_json(x)))) 
  list(df1=df1, df2=df2)
}
 
retrieve_embedding = function(simba_ref) {
 simba_ref$read_embedding()
}

# for scRNA-seq, the embedding "cell" element
# needs to be reordered if there are metadata
# characteristics (such as type or cluster labels)
# to be propagated for graph investigation

# the embedding $C component will have an empty
# data frame, and will also have an obs_names
# element that is an index.  this can be used
# to re-order the original AnnData object's 'obs'
# element, from which a label can be extracted

#' propagation -- updates a post-train cell embedding, in
#' which embedding "rows" corresponding to cells have
#' been reordered, adding a column to the `obs` component
#' @examples
#' p3k = get_10x3kpbmc_path(overwrite=TRUE) # allow repetition
#' ref = simba_ref()
#' pp = ref$read_h5ad(p3k)
#' # bb = basic_preproc(pp, simba_ref=ref)
#' # gg = build_and_train_pbg( bb, simba_ref=ref )
#' # labc = propagate_label( ref, pp, "celltype" )
#' cemb = ref$read_h5ad( get_3k_cell_emb() )
#' cemb$obs$celltype = pp[cemb$obs_names]$obs[["celltype"]] # what propagate_label does
#' um = uwot::umap( cemb['X'] )
#' plot(um, pch=19, col=factor(cemb$obs$celltype) )
#' @export
propagate_label = function( simba_ref, origad, label_name ) {
  cglist = retrieve_embedding( simba_ref )
  cglist$C$obs$tmp =
    origad[cglist$C$obs_names]$obs[[label_name]] # reorder + assign
  names(cglist$C$obs) = label_name   # rename 'tmp'
  cglist$C
}

