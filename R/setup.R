


#' hand sfaira from specified basilisk environment
#' @importFrom reticulate import
#' @note simba 1.2
#' @export
simba_ref = function() {
 proc = basilisk::basiliskStart(bsklenv, testload="simba") # avoid package-specific import
 basilisk::basiliskRun(proc, function() {
     reticulate::import("simba")
   })
}

#' hand scanpy from specified basilisk environment
#' @note assumed implicitly installed with simba 1.2
#' @export
scanpy_ref = function() {
 proc = basilisk::basiliskStart(bsklenv, testload="scanpy") # avoid package-specific import
 #on.exit(basilisk::basiliskStop(proc))
 basilisk::basiliskRun(proc, function() {
     reticulate::import("scanpy")
   })
}
