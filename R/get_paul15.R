#' retrieve path to Paul 2015 (Cell 163:7, 1663) AnnData instance
#' @importFrom R.utils decompressFile
#' @param \dots passed to R.utils::decompressFile, note that ext, FUN,
#' temporary and remove are set.
#' @export
get_paul15_path = function (...) 
{
    gzpath = system.file("extdata", "paul15.h5ad.gz", package = "Simba4Bioc")
    R.utils::decompressFile(gzpath, ext = "gz", FUN = gzfile, 
        temporary = TRUE, remove = FALSE, ...)
}

#' retrieve path to 10x 3k PBMC h5ad (compressed)
#' @param \dots passed to R.utils::decompressFile, note that ext, FUN,
#' temporary and remove are set.
#' @export
get_10x3kpbmc_path = function (...) 
{
    gzpath = system.file("extdata", "tenx3k.h5ad.gz", package = "Simba4Bioc")
    R.utils::decompressFile(gzpath, ext = "gz", FUN = gzfile, 
        temporary = TRUE, remove = FALSE, ...)
}
