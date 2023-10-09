library(shiny)
library(Simba4Bioc)
library(uwot)
library(ggplot2)
library(plotly)
sref = simba_ref()

ui = fluidPage(
 sidebarLayout(
  sidebarPanel(
   helpText("SIMBA single cell analysis"),
   fileInput("h5ad", "h5ad upload"),
   uiOutput("varbox"),
   textOutput("stuff"),
   uiOutput("preproc"),
   uiOutput("step2"),
   width=3
   ),
  mainPanel(
   tabsetPanel(
    tabPanel("Cell embedding", plotlyOutput("cellum")),
    tabPanel("Cell and gene embedding", plotlyOutput("jointum"))
   )
  )
 )
)

suffix = function(x) sub(".*\\.(.*)", "\\1", x)

admess = function (x) 
{
    sprintf("Analyzing AnnData with %d cells, %d genes.", 
        nrow(x["X"]), ncol(x["X"]))
}


server = function(input, output, session) {
  options(shiny.maxRequestSize=40*1024^2)
  get_input = reactive({
   sref$read_h5ad(input$h5ad$datapath)
   })
  output$stuff = renderText({
   shiny::validate(need(!is.null(input$h5ad), "waiting for h5ad file selection"))
   shiny::validate(need(suffix(input$h5ad)=="h5ad", "file suffix must be 'h5ad'"))
   bb = get_input()
   output$varbox = renderUI({
     checkboxGroupInput("obsvar", "cell vbls", 
        choices=names(bb$obs), selected=names(bb$obs)[1])
     })
   output$preproc = renderUI({
    actionButton("dopreproc", "Basic preprocessing")
    })
   admess(bb)
   })
  observeEvent( input$dopreproc, {
   bb = get_input()
   post <<- basic_preproc(bb, simba_ref = sref)
   output$step2 = renderUI({
    actionButton("dopbg", "Build graph")
    })
   })
  make_graph = reactive({
    bb = get_input() # original object
    build_and_train_pbg( bb, simba_ref = sref)  # assigns info into sref
    propagate_label(sref, post, names(bb$obs)[1] )  # NEED SELECTOR FIXME
    })
#
# cell embedding
#
  observeEvent( input$dopbg, {
   output$cellum = renderPlotly({
    bb = get_input() # original object
    post2 <<- make_graph()
    um = uwot::umap(post2['X'])
    post2$obs$cllab = post[post2$obs_names]$obs[[names(bb$obs)[1]]] # FIXME 'cllab' reserved?
    myd = data.frame(x=um[,1], y=um[,2], ty=post2$obs$cllab)
    g1 = ggplot(myd, aes(x=x,y=y,colour=ty)) + geom_point()
    ggplotly(g1)
    })
   })
#
# joint embedding
#
   output$jointum = renderPlotly({
    bb = get_input() # original object
# phase 3 - joint embedding by querying gene embedding to cell embedding
    emb = sref$read_embedding()
    cemb = emb$C
    gemb = emb$G
    # bivariate embedding
    bvout = sref$tl$embed(adata_ref = cemb, list_adata_query = list(gemb) )
    cemb$obs$cllab = post2$obs$cllab
    eanno = c(as.character(post2$obs$cllab), rep("gene", nrow(gemb['X'])))
    bvout$obs$entity_anno = eanno
    uuu = umap(bvout['X'])
    fullanno = c(as.character(post2$obs$cllab), rownames(bb$var))
    bvdf = data.frame(x=uuu[,1], y=uuu[,2], ent=bvout$obs$entity_anno,
       fullanno=fullanno)
    g1 = ggplot(bvdf, aes(x=x, y=y, colour=ent, text=fullanno)) + geom_point()
    ggplotly(g1)
  })
}

runApp(list(ui=ui, server=server))
