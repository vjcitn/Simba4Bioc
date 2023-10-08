
test_that("barcode data for MS4A1 exhibits presence in B cells", {
 sref = simba_ref()
 p3k = sref$read_h5ad(get_10x3kpbmc_path())
 myd = simba_barcode_data( symbol = "MS4A1",
    origad = p3k, clabel = "celltype",
    simbaref = sref, cembpath = get_3k_cell_emb(),
    gembpath = get_3k_gene_emb() )
 expect_true(myd$attr[1] == "B")
})

