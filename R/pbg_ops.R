#' retrieve and HDF5 representation of an embedding produced by pbg in simba
#' @importFrom rhdf5 h5read
#' @param ctarpath character(1) path to a compressed tar of pbg folder
#' @param type character(1) 'cell' or 'gene'
#' @param numtag character(1) defaults to "_0"
#' @param gname character(1) used as dirname in build_and_train_pbg
#' @param mname character(1) used as output in build_and_train_pbg
#' @param spec character(1) defaults to "v10"
#' @param exdir character(1) path where tar operation outputs will be placed,
#' defaults to tempdir()
#' @return result of rhdf5::h5read on the embedding values
#' @examples
#' tpath = system.file(file.path("extdata", "pbg3k.tar.xz"), package="Simba4Bioc")
#' cemb = get_emb_h5(tpath)
#' dim(cemb)
#' @export
get_emb_h5 = function (ctarpath, type = "cell", numtag = "_0", gname = "graph0", 
    mname = "model", spec = "v10", exdir = tempdir()) 
{
    tag = switch(type, cell = "C", gene = "G")
    if (is.null(tag)) 
        stop("type not in ('cell', 'gene')")
    entity = paste0("embeddings_", tag, numtag, ".", spec, ".h5")
    file_in_tar = paste0("pbg/", gname, "/", mname, "/", entity)
    untar(ctarpath, files = file_in_tar, exdir = exdir)
    rhdf5::h5read(file.path(exdir, file_in_tar), "/embeddings")
}

#get_emb_h5 = function(ctarpath, type="cell",  numtag="_0", gname="graph0", 
#   mname="model", spec="v10", exdir=tempdir()) {
#  tag = switch(type, cell="C", gene="G")
#  if (is.null(tag)) stop("type not in ('cell', 'gene')")
#  entity = paste0("embeddings_", tag, numtag, ".", spec, ".h5")
#  file_in_tar = paste0("pbg/", gname, "/", mname, "/", entity)
#  untar(ctarpath, files=file_in_tar, exdir=exdir)
#  rhdf5::h5read(file.path(exdir, file_in_tar), "/embeddings")
#}
