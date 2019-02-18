constructGeneRegionTracks <- function(myTxDb, BioneerV5) {
  gr <- GeneRegionTrack(myTxDb)
  gr <- gr[which(grepl("^ENSRNA", group(gr))==F)]
  displayPars(gr) <- list(shape=c("smallArrow", "arrow"), showId=T, transcriptAnnotation="gene", fill="gray", just.group="above")
  grGroups <- group(gr)
  grGroups <- sub("\\.\\d+$", "", grGroups)
  grBioneer <- gr[grGroups %in% BioneerV5$Systematic.ID]
  grOther <- gr[grGroups %in% BioneerV5$Systematic.ID == F]
  names(grBioneer) <- "Bioneer V5\nTranscripts"
  names(grOther) <- "Other\nTranscripts"
  return(list(grBioneer, grOther))
}
constructAlignedReadTrack <- function(tags, name, ymax=200) {
  
  tagsGRanges <- GRanges(Rle(tags$chr), IRanges(tags$start, tags$end), Rle(tags$strand), col=tags$col)
  
  tagsGRanges <-  split(tagsGRanges, seqnames(tagsGRanges))
  
  tagsGRanges <- lapply(tagsGRanges, function(i) {split(i, i$col)})
  
  tagsTracks <- lapply(tagsGRanges, function(i) {lapply(i, function(j) {AlignedReadTrack(j, lwd=0, fill=j$col[1], alpha.title=1, alpha=0.8, name=name, ylim=c(0, ymax))})})
  
  tagsTracks <- lapply(tagsTracks, OverlayTrack)
  
  return(tagsTracks)
  
}
getGenesSharingBarcodes <- function(tags, gene, TagsGeneFreq, minTagGeneFreq) {
  
  tagsSub <- tags[tags$gene %in% gene,]
  
  tagsFreq <- table(as.character(tagsSub$tag_BestMatch))
  
  tagsFreq <- tagsFreq[tagsFreq>=minTagGeneFreq]
  
  tagsAllGenes <- TagsGeneFreq[TagsGeneFreq$barcode %in% names(tagsFreq),]
  
  return(tagsAllGenes)
  
}
plotBarcodeGeneNetwork <- function(gene, tag, tags, TagsGeneFreq, totalTagsTable, minTagGeneFreq, scale=1) {
  
  geneNodes <- getGenesSharingBarcodes(tags, gene, TagsGeneFreq, minTagGeneFreq)
  geneNodes <- as.character(unique(geneNodes$gene[geneNodes$freq >= minTagGeneFreq]))
  geneNodes <- unique(c(gene, geneNodes))
  
  links <- getGenesSharingBarcodes(tags, geneNodes, TagsGeneFreq, 1)
  links <- subset(links, select=c("gene", "barcode", "freq"))
  links <- links[links$gene %in% geneNodes,]
  links$gene <- as.character(links$gene)
  links$barcode <- as.character(links$barcode)
  
  if(nrow(links)==0) {
    return()
  }
  
  totalTagsTableSub <- totalTagsTable[totalTagsTable$tag==tag,]
  links2plot <- links[which(paste(links$gene, links$barcode, sep=":") %in% paste(totalTagsTableSub$gene, totalTagsTableSub$barcode, sep=":")),]
  if(nrow(totalTagsTableSub)==0) {
    links2merge <- links
  } else {
    links2merge <- links[which(paste(links$gene, links$barcode, sep=":") %in% paste(totalTagsTableSub$gene, totalTagsTableSub$barcode, sep=":")==F),]
  }
  links2merge <- plyr::ddply(links2merge, plyr::.(gene), plyr::summarize, barcode=paste0(gene[1], "_Others"), freq=sum(freq))
  links2plot <- rbind(links2plot, links2merge)
  
  barcodeNodes <- unique(links2plot$barcode)
  
  g <- graph.empty()
  
  g <- add.vertices(g, length(geneNodes), attr=list(name=geneNodes, type=rep(T, length(geneNodes))))
  g <- add.vertices(g, length(barcodeNodes), attr=list(name=barcodeNodes), type=rep(F, length(barcodeNodes)))
  
  vertexLinks <- subset(links2plot, select=c("gene", "barcode"))
  vertexLinks <- as.vector(t(as.matrix(vertexLinks)))
  
  linkCols <- character(0)
  for(i in 1:nrow(links2plot)) {
    if(grepl("_Others$", links2plot$barcode[i])) {
      linkCols[i] <- "gray"
    } else {
      linkCols[i] <- totalTagsTableSub$col[totalTagsTableSub$gene==links2plot$gene[i] & totalTagsTableSub$barcode==links2plot$barcode[i]]
    }
  }
  
  g <- add.edges(g, vertexLinks, weight=c(links2plot$freq), color=linkCols)
  g <- as.undirected(g, "each")
  
  nodeSizes <- integer(0)
  counter <- 0
  for(i in geneNodes) {
    counter <- counter + 1
    nodeSizes[counter] <- sum(links2plot$freq[links2plot$gene==i])
  }
  for(i in barcodeNodes) {
    counter <- counter + 1
    nodeSizes[counter] <- sum(links2plot$freq[links2plot$barcode==i])
  }
  
  cols <- c("black", "white", "gray", totalTagsTable$col)
  
  pie.cols <- list()
  counter <- 0
  for(i in geneNodes) {
    counter <- counter + 1
    if(i==gene) {
      pie.cols[[counter]] <- c(1, rep(0, length(cols)-1))
    } else {
      pie.cols[[counter]] <- c(0, 1, rep(0, length(cols)-2))
    }
  }
  for(i in barcodeNodes) {
    counter <- counter + 1
    if(grepl("_Others$", i)) {
      pie.cols[[counter]] <- c(0, 0, 1, rep(0, length(cols)-3))
    } else {
      pie.cols[[counter]] <- c(rep(0, 3), sapply(1:nrow(totalTagsTable), function(j) {
        ifelse(i == totalTagsTable$barcode[j], totalTagsTable$freq[j], 0)
      }))
    }
  }
  
  layout <- layout_as_bipartite(g)
  layout <- layout.norm(layout[,2:1])
  plot(g, layout=layout, xlim=c(0, 1), vertex.size=nodeSizes*0.4*scale, vertex.label=NA, edge.color=E(g)$color, edge.width=E(g)$weight*0.1*scale, vertex.shape="pie", vertex.pie=pie.cols, vertex.pie.color=list(cols))
  
  myText <- ifelse(grepl("_Others$", names(V(g))), "Other Barcodes", names(V(g)))
  for(i in 1:nrow(layout)) {
    text(x=layout[i, 1], y=layout[i, 2], myText[i], pos=c(4, 2)[V(g)$type[i]+1], offset=2.5)
  }
  
}
createTotalTagsTable <- function(upTags, dnTags, upTagsAgg, dnTagsAgg, gene, minTagGeneFreq) {
  dnTags$col <- "gray"
  upTags$col <- "gray"
  
  dnTagsShared <- getGenesSharingBarcodes(dnTags, gene, dnTagsAgg, minTagGeneFreq)
  upTagsShared <- getGenesSharingBarcodes(upTags, gene, upTagsAgg, minTagGeneFreq)
  
  nDnTags <- sum(dnTagsShared$freq>=minTagGeneFreq)
  nUpTags <- sum(upTagsShared$freq>=minTagGeneFreq)
  
  totalTags <- sum(nDnTags, nUpTags)
  
  if(totalTags>0) {
    
    tagsCol <- sample(colorspace::rainbow_hcl(totalTags))
    
    dnTagsShared$tag <- if(nrow(dnTagsShared)>0) {dnTagsShared$tag <- "dnTags"} else {dnTagsShared$tag <- character(0)}
    upTagsShared$tag <- if(nrow(upTagsShared)>0) {upTagsShared$tag <- "upTags"} else {upTagsShared$tag <- character(0)}
    totalTagsTable <- rbind(dnTagsShared, upTagsShared)
    totalTagsTable <- totalTagsTable[totalTagsTable$freq>=minTagGeneFreq,]
    totalTagsTable$col <- tagsCol
    
    for(i in 1:nrow(totalTagsTable)) {
      tags <- get(totalTagsTable$tag[i])
      tags$col[tags$gene==as.character(totalTagsTable$gene[i]) & tags$tag_BestMatch==as.character(totalTagsTable$barcode[i])] <- totalTagsTable$col[i]
      assign(totalTagsTable$tag[i], tags)
    }
    
  } else {
    
    totalTagsTable <- data.frame(gene=character(0), tag_BestMatch=character(0), freq=integer(0), tag=character(0), col=character(0))
    
  }
  
  return(list(dnTags, upTags, totalTagsTable))
  
}
createAllFeatures <- function(upTags, dnTags, myTxDb, gene, totalTagsTable, minDist) {
  allFeatures <- rbind(dnTags[dnTags$col!="gray",], upTags[upTags$col!="gray",])
  allFeatures <- GRanges(Rle(allFeatures$chr), IRanges(allFeatures$start, allFeatures$end))
  allFeatures <- suppressWarnings(c(allFeatures, transcripts(myTxDb, columns=NULL, filter=list(gene_id=c(gene, as.character(totalTagsTable$gene))))))
  allFeatures <- reduce(allFeatures, ignore.strand=T, min.gapwidth=minDist)
  return(allFeatures)
}
plotAllTracks <- function(i, allFeatures, gene, totalTagsTable, minDist, ax, gr, upTagTracks, dnTagTracks, iTrack) {
  chr <- as.character(allFeatures[i]@seqnames)
  from <- allFeatures[i]@ranges@start - minDist
  to <- allFeatures[i]@ranges@start + allFeatures[i]@ranges@width + minDist
  
  grBioneer <- gr[[1]]
  grOther <- gr[[2]]
  feature(grBioneer) <- rep("gray", length(feature(grBioneer)))
  feature(grBioneer)[group(grBioneer) %in% paste0(totalTagsTable$gene, ".1")] <- "white"
  feature(grBioneer)[group(grBioneer)==paste0(gene, ".1")] <- "black"
  colArgs <- unique(feature(grBioneer))
  colArgs <- colArgs[-which(colArgs=="gray")]
  colArgs <- as.list(colArgs)
  names(colArgs) <- colArgs
  
  do.call(plotTracks, c(list(trackList=list(iTrack[[chr]], ax, grOther, grBioneer, upTagTracks[[chr]], dnTagTracks[[chr]]), from=from, to=to, chromosome=chr), colArgs))
}
createTagChoice <- function(tag, totalTagsTable, gene, n, tagsAutoSafe, tagsManualSafe) {
  
  if(nrow(totalTagsTable)==0) {
    return(radioButtons(inputId = paste0(tag, "Choice", n), label = paste0("Select ", sub("s$", "", tag)), choiceNames = list("NA"), choiceValues = list("NA")))
  }
  choices <- as.list(as.character(totalTagsTable$barcode))
  choiceNames <- list()
  for(i in 1:nrow(totalTagsTable)) {
    choiceNames[[i]] <- shiny::tags$span(style=paste0("color:", totalTagsTable$col[i]), totalTagsTable$barcode[i])
  }
  choiceValues <- choices
  choiceNames <- c(list(shiny::tags$span(style="color: black", "NA")), choiceNames)
  choiceValues <- c(list("NA"), choiceValues)
  safeTags <- read.csv(tagsAutoSafe, stringsAsFactors = F)
  default <- safeTags$barcode[safeTags$gene==gene]
  if(length(default)==0) {
    default <- tryCatch({
      manualTags <- read.csv(tagsManualSafe, stringsAsFactors = F)
      manualTags$barcode[manualTags$gene==gene]
    }, warning=function(w) {"NA"}, error=function(e) {"NA"})
  }
  if(length(default)==0) {
    default <- "NA"
  }
  radioButtons(inputId = paste0(tag, "Choice", n), label = paste0("Select ", sub("s$", "", tag)), choiceNames = choiceNames, choiceValues = choiceValues, selected = default)
  
}
getMissingGenes <- function(TagsAutoSafe, TagsManualSafe) {
  
  TagsAutoSafeTable <- read.csv(TagsAutoSafe)
  TagsSafe <- as.character(TagsAutoSafeTable$gene)
  TagsManualSafeTable <- tryCatch({read.csv(TagsManualSafe)}, error=function(e) {return(NULL)}, warning=function(w) {return(NULL)})
  if(is.null(TagsManualSafeTable)==F) {
    TagsSafe <- c(TagsSafe, as.character(TagsManualSafeTable$gene))
  }
  missing <- BioneerV5$Systematic.ID[BioneerV5$Systematic.ID %in% TagsSafe==F]
  return(missing)
  
}
runCheckSave <- function(tag, tagsAutoSafe, tagsManualSafe) {
  currentTags <- read.csv(tagsAutoSafe, stringsAsFactors = F)$barcode
  currentTags <- tryCatch({
    c(currentTags, read.csv(tagsManualSafe, stringsAsFactors = F)$barcode)
  }, warning=function(w) {return(currentTags)}, error=function(e) {return(currentTags)})
  currentTags <- currentTags[is.na(currentTags)==F]
  if(tag %in% currentTags) {
    return(T)
  } else {
    return(F)
  }
}
runManualSave <- function(gene, upTag, dnTag, upInSafe, dnInSafe, upTagsManualSafe, dnTagsManualSafe) {
  if(upInSafe==F) {
    colNamesUp <- ifelse(file.exists(upTagsManualSafe), F, T)
    suppressWarnings({write.table(data.frame(barcode=upTag, gene=gene, method="manual"), upTagsManualSafe, append=T, sep=",", col.names = colNamesUp, row.names = F)})
  }
  if(dnInSafe==F) {
    colNamesDn <- ifelse(file.exists(dnTagsManualSafe), F, T)
    suppressWarnings({write.table(data.frame(barcode=dnTag, gene=gene, method="manual"), dnTagsManualSafe, append=T, sep=",", col.names = colNamesDn, row.names = F)})
  }
}

#'@export
BahlerBarcodeBrowser <- function(upTags, dnTags, upTagsAgg, dnTagsAgg, upTagsAutoSafe, dnTagsAutoSafe, upTagsManualSafe, dnTagsManualSafe, myTxDb) {
  
  #Data input
  upTags <- read.csv(upTags)
  dnTags <- read.csv(dnTags)
  upTagsAgg <- read.csv(upTagsAgg)
  dnTagsAgg <- read.csv(dnTagsAgg)
  
  #Set up axis track
  ax <- GenomeAxisTrack()
  
  #Set up gene region tracks for Bioneer and other transcripts - need to colour later
  gr <- constructGeneRegionTracks(myTxDb, BioneerV5)
  
  ui <- fluidPage(
    shinyjs::useShinyjs(),
    h1("Bahler Barcode Browser", style = "text-align:center"),
    br(),
    fluidRow(
      column(
        3,
        actionButton(inputId = "newGene", label = "New gene!"),
        br(),
        br(),
        uiOutput("geneInput"),
        numericInput(inputId = "minTagGeneFreq", label = "Minimum co-occurence frequency between barcode and gene to highlight barcode", value = 100),
        numericInput(inputId = "minDist", label = "Minimum Distance between features to plot on the same track / bp", value = 10000),
        numericInput(inputId = "ymax", label = "ymax for aligned read track", value = 200),
        numericInput(inputId = "scale", label = "Scaling node and edge sizes for network visualisation", value = 0.05),
        br(),
        shinyjs::disabled(actionButton("saveChoice", "Save & new gene!"))
      ),
      column(
        9,
        fluidRow(
          uiOutput("iUI"),
          plotOutput("GvizPlots"),
          br()
        ),
        fluidRow(
          column(
            6,
            uiOutput("upTagChoiceUI"),
            plotOutput("upTagNetwork")
          ),
          column(
            6,
            uiOutput("dnTagChoiceUI"),
            plotOutput("dnTagNetwork")
          )
        )
      )
    )
  )
  
  server <- function(input, output) {
    
    myReactiveValues <- reactiveValues()
    myReactiveValues$gene <- ""
    myReactiveValues$geneCounter <- 0
    
    observeEvent(input$newGene, {
      myReactiveValues$missingUpGenes <- getMissingGenes(upTagsAutoSafe, upTagsManualSafe)
      myReactiveValues$missingDnGenes <- getMissingGenes(dnTagsAutoSafe, dnTagsManualSafe)
      myReactiveValues$gene <- sample(unique(c(myReactiveValues$missingUpGenes, myReactiveValues$missingDnGenes)), 1)
    })
    
    output$geneInput <- renderUI({
      textInput(inputId = "gene", label = "Gene", value = myReactiveValues$gene)
    })
    
    gene <- reactive({
      req(input$gene)
      validate(need(paste0(input$gene) %in% BioneerV5$Systematic.ID, "Please enter the systematic name of a Bioneer V5 mutant."))
      return(input$gene)
    })
    
    minTagGeneFreq <- reactive({
      msg <- "Please make sure that the minimum co-occurence frequency between barcode and gene is a positive integer."
      validate(need(is.integer(input$minTagGeneFreq), msg))
      validate(need(input$minTagGeneFreq>0, msg))
      return(input$minTagGeneFreq)
    })
    
    minDist <- reactive({
      msg <- "Please make sure that the minimum distance is an integer above 100."
      validate(need(is.integer(input$minDist), msg))
      validate(need(input$minDist > 100, msg))
      return(input$minDist)
    })
    
    totalTagsTableOutput <- reactive({
      totalTagsTableOutput <- createTotalTagsTable(upTags, dnTags, upTagsAgg, dnTagsAgg, gene(), minTagGeneFreq())
    })
    
    dnTagsCol <- reactive({
      totalTagsTableOutput()[[1]]
    })
    
    upTagsCol <- reactive({
      totalTagsTableOutput()[[2]]
    })
    
    totalTagsTable <- reactive({
      totalTagsTableOutput()[[3]]
    })
    
    upTagTracks <- reactive({
      constructAlignedReadTrack(upTagsCol(), "upTags", input$ymax)
    })
    
    dnTagTracks <- reactive({
      constructAlignedReadTrack(dnTagsCol(), "dnTags", input$ymax)
    })
    
    allFeatures <- reactive({
      createAllFeatures(upTagsCol(), dnTagsCol(), myTxDb, gene(), totalTagsTable(), minDist())
    })
    
    output$nPanels <- reactive({
      length(allFeatures())
    })
    
    outputOptions(output, "nPanels", suspendWhenHidden=F)
    
    output$iUI <- renderUI({
      if(length(allFeatures())==1) {return()}
      sliderInput(inputId = "i", min=1, max=length(allFeatures()), label="The barcodes are scattered across multiple genomic regions. Scroll to view different regions", value=1, step=1, width="40%")
    })
    
    i <- reactive({
      if(length(allFeatures())==1) {return(1)} else {return(input$i)}
    })
    
    output$GvizPlots <- renderPlot({
      req(i())
      plotAllTracks(i(), allFeatures(), gene(), totalTagsTable(), minDist(), ax, gr, upTagTracks(), dnTagTracks(), iTrack)
    })
    
    output$upTagNetwork <- renderPlot({
      plotBarcodeGeneNetwork(gene(), "upTags", upTags, upTagsAgg, totalTagsTable(), minTagGeneFreq(), input$scale)
    })
    
    output$dnTagNetwork <- renderPlot({
      plotBarcodeGeneNetwork(gene(), "dnTags", dnTags, dnTagsAgg, totalTagsTable(), minTagGeneFreq(), input$scale)
    })
    
    output$upTagChoiceUI <- renderUI({
      totalTagsTableSub <- totalTagsTable()[totalTagsTable()$tag=="upTags" & totalTagsTable()$gene==gene(),]
      createTagChoice("upTags", totalTagsTableSub, gene(), myReactiveValues$geneCounter, upTagsAutoSafe, upTagsManualSafe)
    })
    
    output$dnTagChoiceUI <- renderUI({
      totalTagsTableSub <- totalTagsTable()[totalTagsTable()$tag=="dnTags" & totalTagsTable()$gene==gene(),]
      createTagChoice("dnTags", totalTagsTableSub, gene(), myReactiveValues$geneCounter, dnTagsAutoSafe, dnTagsManualSafe)
    })
    
    observeEvent(totalTagsTable(), {
      myReactiveValues$geneCounter <- myReactiveValues$geneCounter + 1
    })
    
    observeEvent(input[[paste0("upTagsChoice", myReactiveValues$geneCounter)]], {
      if(gene() %in% myReactiveValues$missingUpGenes) {
        myReactiveValues$upInSafe <- F
        shinyjs::enable(paste0("upTagsChoice", myReactiveValues$geneCounter))
      } else {
        myReactiveValues$upInSafe <- T
        shinyjs::disable(paste0("upTagsChoice", myReactiveValues$geneCounter))
      }
    })
    
    observeEvent(input[[paste0("dnTagsChoice", myReactiveValues$geneCounter)]], {
      if(gene() %in% myReactiveValues$missingDnGenes) {
        myReactiveValues$dnInSafe <- F
        shinyjs::enable(paste0("dnTagsChoice", myReactiveValues$geneCounter))
      } else {
        myReactiveValues$dnInSafe <- T
        shinyjs::disable(paste0("dnTagsChoice", myReactiveValues$geneCounter))
      }
    })
    
    observeEvent(totalTagsTable(), {
      shinyjs::enable("saveChoice")
    })
    
    observeEvent(input$saveChoice, {
      if(myReactiveValues$upInSafe == F) {
        problemUp <- runCheckSave(input[[paste0("upTagsChoice", myReactiveValues$geneCounter)]], upTagsAutoSafe, upTagsManualSafe)
      } else {
        problemUp <- F
      }
      if(problemUp) {
        showNotification("Error: The chosen UpTag has already been assigned to another gene. Please make another selection.", type="error")
      } else {
        if(myReactiveValues$dnInSafe == F) {
          problemDn <- runCheckSave(input[[paste0("dnTagsChoice", myReactiveValues$geneCounter)]], dnTagsAutoSafe, dnTagsManualSafe)
        } else {
          problemDn <- F
        }
        if(problemDn) {
          showNotification("Error: The chosen DnTag has already been assigned to another gene. Please make another selection.", type="error")
        } else {
          shinyjs::disable("saveChoice")
          runManualSave(gene(), input[[paste0("upTagsChoice", myReactiveValues$geneCounter)]], input[[paste0("dnTagsChoice", myReactiveValues$geneCounter)]], myReactiveValues$upInSafe, myReactiveValues$dnInSafe, upTagsManualSafe, dnTagsManualSafe)
          shinyjs::click("newGene")
          invalidateLater(100)
        }
      }
    })
    
  }
  
  shinyApp(ui, server)
  
}
