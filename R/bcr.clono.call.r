#' BCR 10x clonotype calling with hamming distance
#'
#' Makes clonotype calls on 10x data using the hamming distance calls 
#' @param contig.list  normalized log counts
#' @param seq type of sequence to use for clonotype calling (nt or aa) default is nt
#' @param V.gene uses V gene to make clonotype calls (default is TRUE)
#' @param J.gene uses J gene to make clonotype calls (default is FALSE)
#' @param CDR3 uses CDR3 sequence to make clonotype calls (default is TRUE)
#' @param chain use light, heavy or both chains to make clonotype (default is both)
#' @param hammingthreshold threshold for which to group cells/sequences in the same clonotype (default =0.7)
#' @param cluster.plot output hamming distance matrix plot 
#' @param graph.plot output graph plot with hamming distance 
#' @param results_folder directory to output information to (default is getwd()), verbose=TRUE whether to save figures in PDF and PNG file format
#' @param verbose print messages (default is TRUE)
#' @export
#' @import e1071
#' @import stringr
#' @import dplyr
#' @import data.table
#' @import ComplexHeatmap
#' @import circlize
#' @import igraph
#' @examples
#' contig.list <- CallClonoHD(contig.listFromCellRanger)
#' 

BCR.CallClono.HD <- function(contig.list, seq="nt", V.gene=TRUE, CDR3=TRUE, J.gene=FALSE, chain="both",
     hammingthreshold=0.7, cluster.plot=TRUE, graph.plot=FALSE, results_folder=getwd(), verbose=FALSE) {

    calculate_hamming_distance <-function(x) {#, all_combinationstmp, totaldatafinal, seqofinterest ) {
        seq1 <-unlist(strsplit(x[[1]], split=""))
        seq2 <-unlist(strsplit(x[[2]], split=""))

        value <- e1071::hamming.distance(seq1, seq2)
        value <- value/length(seq1)
        return(c(1-value))
    }

    #parse through options 
    if (!(seq %in% c("nt", "aa"))) {
        stop("ERROR: seq value needs to be aa and nt")
    }
    if (!(chain %in% c("heavy", "light", "both"))) {
        stop("ERROR: chain value needs to be heavy, light, or both")
    }
    cdrseqs <- c()
    seqofinterest <- c()
    geneseqs <- c()
    if (seq=="nt") {
        if (chain=="heavy" || chain=="both") {
            if (isTRUE(V.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt"))
                geneseqs <- c(geneseqs, c("fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt"))
            }
            if (isTRUE(J.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr4_nt"))
                geneseqs <- c(geneseqs, c("fwr4_nt"))
            }
            if (isTRUE(CDR3)) {
                seqofinterest <- c(seqofinterest, c("cdr3_nt"))
                cdrseqs <- c(cdrseqs,  c("cdr3_nt"))
            }
        }
        if((chain=="light" || chain=="both")) {
            if (isTRUE(V.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr1_nt.light", "cdr1_nt.light", "fwr2_nt.light", "cdr2_nt.light", "fwr3_nt.light"))
                geneseqs <- c(geneseqs, c("fwr1_nt.light", "cdr1_nt.light", "fwr2_nt.light", "cdr2_nt.light", "fwr3_nt.light"))
            }
            if (isTRUE(J.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr4_nt.light"))
                geneseqs <- c(geneseqs, c("fwr4_nt.light"))

            }
            if (isTRUE(CDR3)) {
                seqofinterest <- c(seqofinterest, c("cdr3_nt.light"))
                cdrseqs <- c(cdrseqs,  c("cdr3_nt.light"))

           }
        }
    } else {
        if (chain=="heavy" || chain=="both") {
            if (isTRUE(V.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3"))
                geneseqs <- c(geneseqs, c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3"))
            }
            if (isTRUE(J.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr4"))
                geneseqs <- c(geneseqs, c("fwr4"))

            }
            if (isTRUE(CDR3)) {
                seqofinterest <- c(seqofinterest, c("cdr3"))
                cdrseqs <- c(cdrseqs,  c("cdr3"))

            }
        }
        if((chain=="light" || chain=="both")) {
            if (isTRUE(V.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr1.light", "cdr1.light", "fwr2.light", "cdr2.light", "fwr3.light"))
                geneseqs <- c(geneseqs, c("fwr1.light", "cdr1.light", "fwr2.light", "cdr2.light", "fwr3.light"))

            }
            if (isTRUE(J.gene)) {
                seqofinterest <- c(seqofinterest, c("fwr4.light"))
                geneseqs <- c(geneseqs, c("fwr4.light"))

            }
            if (isTRUE(CDR3)) {
                seqofinterest <- c(seqofinterest, c("cdr3.light"))
                cdrseqs <- c(cdrseqs,  c("cdr3"))
            }
        }
    }

    totaldata <- c()
    for (n in names(contig.list)) {
        contig.list[[n]]$barcodev2 <- paste(n, contig.list[[n]]$barcode, sep="_")
        contig.list[[n]]$sample_individual <- rep(n, nrow(contig.list[[n]]))
        totaldata <- rbind(totaldata, contig.list[[n]])
    }
    totaldata <- as.data.frame(totaldata)
    df <- totaldata %>% 
        group_by(barcodev2) %>% # Variable to be transformed
        dplyr::count()

    df <- as.data.frame(df)
    removebarcodes <- df[df$n==1, "barcodev2"]
    totaldata <- totaldata[!(totaldata$barcodev2 %in% removebarcodes), ]
    if (isTRUE(verbose)) {
        message("STATUS: build table of VDJ information")
    }
    write.csv(totaldata, "totaltable.csv", quote=F, row.names=F)

    #construct table 
    totaldatafinal <- build_data_table(totaldata)
    #get sequence lengths 
    if (isTRUE(verbose)) {
        message("STATUS: get sequence length")
    }
    x <- c()
    x2 <- c()
    x3 <- c()
    x4 <- c()
    x5 <- c()
    for (i in rownames(totaldatafinal)) {
        x <-c(x, nchar(paste(totaldatafinal[i,seqofinterest], collapse = "")))
        x2 <-c(x2, paste(totaldatafinal[i,seqofinterest], collapse = ""))
        if (length(cdrseqs)>0) {
            x3 <- c(x3, nchar(paste(totaldatafinal[i,cdrseqs], collapse = "")))
            x4 <-c(x4, paste(totaldatafinal[i,cdrseqs], collapse = ""))
        }
        if (length(geneseqs)>0) {
            x5 <- c(x5, nchar(paste(totaldatafinal[i,geneseqs], collapse = "")))
        }
    }
    totaldatafinal$length4alignment <- x
    totaldatafinal$seqofinterest <- x2
    if (length(cdrseqs)>0) {
        totaldatafinal$cdr3length4alignment <- x3
        totaldatafinal$cdr3seqofinterest <- x4
    }
    if (length(geneseqs)>0) {
        totaldatafinal$genelength4alignment <- x5
    }

    if (isTRUE(verbose)) {
        message("STATUS: get combinations of barcodes to calculate hamming distance")
    }
    # all_combinations <- comboGrid(totaldatafinal$barcodev2, totaldatafinal$barcodev2, repetition = F)
    all_combinations <- expand.grid(totaldatafinal$barcodev2, totaldatafinal$barcodev2)
    all_combinations <- as.data.frame(all_combinations)

    if (isTRUE(verbose)) {
        message("STATUS: reducing hamming distance calculations based on equal lengths")
    }
    
    if (isTRUE(J.gene)) {
        totaldatafinal$totalvgene <- paste(totaldatafinal$v_gene,totaldatafinal$v_gene.light,totaldatafinal$j_gene,totaldatafinal$j_gene.light, sep="_")
    } else {
        totaldatafinal$totalvgene <- paste(totaldatafinal$v_gene,totaldatafinal$v_gene.light, sep="_")
    }

    # map data on to all combinations 
    all_combinations <- all_combinations %>%
        mutate(seqlengthinterestbarcode1 = mapvalues(Var1,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$length4alignment)))
    all_combinations <- all_combinations %>%
        mutate(seqinterestbarcode1 = mapvalues(Var1,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$seqofinterest)))
    all_combinations <- all_combinations %>%
        mutate(totalvgenebarcode1 = mapvalues(Var1,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$totalvgene)))
    if (length(cdrseqs)>0) {
        all_combinations <- all_combinations %>%
            mutate(cdr3lengthbarcode1 = mapvalues(Var1,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$cdr3length4alignment)))
        all_combinations <- all_combinations %>%
            mutate(cdr3seqbarcode1 = mapvalues(Var1,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$cdr3seqofinterest)))
    }
    if (length(geneseqs)>0) {
        all_combinations <- all_combinations %>%
            mutate(genelengthbarcode1 = mapvalues(Var1,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$genelength4alignment)))
    }
    all_combinations <- all_combinations %>%
        mutate(seqlengthinterestbarcode2 = mapvalues(Var2,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$length4alignment)))
    all_combinations <- all_combinations %>%
        mutate(seqinterestbarcode2 = mapvalues(Var2,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$seqofinterest)))
    all_combinations <- all_combinations %>%
        mutate(totalvgenebarcode2 = mapvalues(Var2,
            from = as.character(totaldatafinal$barcodev2),
            to = as.character(totaldatafinal$totalvgene)))
    if (length(cdrseqs)>0) {
        all_combinations <- all_combinations %>%
            mutate(cdr3lengthbarcode2 = mapvalues(Var2,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$cdr3length4alignment)))
        all_combinations <- all_combinations %>%
            mutate(cdr3seqbarcode2 = mapvalues(Var2,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$cdr3seqofinterest)))
    } 
    if (length(geneseqs)>0) {
        all_combinations <- all_combinations %>%
            mutate(genelengthbarcode2 = mapvalues(Var2,
                from = as.character(totaldatafinal$barcodev2),
                to = as.character(totaldatafinal$genelength4alignment)))
    }
    # filter out combinations of cells that are not the same length or the same gene call
    if(length(cdrseqs)==0) { 
        all_combinations$length.equal <- all_combinations$seqlengthinterestbarcode1 == all_combinations$seqlengthinterestbarcode2
        all_combinations$genecall <- all_combinations$totalvgenebarcode1 == all_combinations$totalvgenebarcode2 
        all_combinations4hd <- all_combinations[(all_combinations$length.equal==TRUE  & all_combinations$genecall==TRUE),]

        all_combinations4hd_NA <- all_combinations[(all_combinations$length.equal==FALSE | all_combinations$genecall==FALSE),]
        all_combinations4hd_NA$hd <- rep(NA, nrow(all_combinations4hd_NA))

    } else if ((length(cdrseqs)>0 & length(geneseqs)>0)) {
        all_combinations$genelength.equal <- all_combinations$genelengthbarcode1 == all_combinations$genelengthbarcode2
        all_combinations$CDR3length.equal <- all_combinations$cdr3lengthbarcode1 == all_combinations$cdr3lengthbarcode2
        all_combinations$genecall <- all_combinations$totalvgenebarcode1 == all_combinations$totalvgenebarcode2

        all_combinations4hd <- all_combinations[(all_combinations$genelength.equal==TRUE & all_combinations$CDR3length.equal==TRUE & all_combinations$genecall==TRUE ),]
        all_combinations4hd_NA <- all_combinations[(all_combinations$genelength.equal!=TRUE | all_combinations$CDR3length.equal!=TRUE |all_combinations$genecall!=TRUE),]
        all_combinations4hd_NA$hd <- rep(NA, nrow(all_combinations4hd_NA))
    }
    #calculate subsetted hamming distances for sequences that are the same length and have the same gene calls
    if (isTRUE(verbose)) {
        message(paste0("STATUS: calculating hamming distance for ", nrow(all_combinations4hd), " combinations"))
    }
    start.time <- Sys.time()
    totalhd = apply(all_combinations4hd[c("seqinterestbarcode1", "seqinterestbarcode2")], 1, calculate_hamming_distance)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    if (isTRUE(verbose)) {
        message("STATUS: ", round(time.taken, 2), " time taken for hamming distances")
    }
    all_combinations4hd$hd <- totalhd
    write.csv(all_combinations4hd, file.path(results_folder, "hammingdistances.csv"), quote=FALSE)
    
    all_combinations4hdtotal <- rbind(all_combinations4hd, all_combinations4hd_NA)
    all_combinations4hdtotal$genecall <- all_combinations4hdtotal$totalvgenebarcode1 == all_combinations4hdtotal$totalvgenebarcode2
    all_combinations4hdtotal[all_combinations4hdtotal$genecall==FALSE, "hd"] <- 0
    dffinal <- dcast(data.table(all_combinations4hdtotal), Var1 ~ Var2, value.var="hd")
    dffinal <- as.data.frame(dffinal)
    rownames(dffinal) <- dffinal$Var1
    dffinal$Var1 <- NULL

    # output complete hamming distancematrix 
    dffinal <- dffinal[order(rownames(dffinal)), order(colnames(dffinal))]
    write.csv(dffinal, file.path(results_folder, "hammingdistancesmatrix.csv"), quote=FALSE)

    #generate hamming distance distance matrix plot 
    if (isTRUE(cluster.plot)) {
        if(isTRUE(verbose)) {
            message("STATUS: Generating hamming distance cluster plot")
        }
        hamming_dist_cluster_plot(dffinal = dffinal, results_folder = results_folder)
    }

    if (isTRUE(verbose)) {
        message ("STATUS: get clonotype numbers")
    }

    # get rid of max values that compare cells to themselves 
    for (b in colnames(dffinal)) {
        dffinal[b,b] <- NA
    }

    # generate graph
    dffinal[is.na(dffinal)] <- 0
    dffinalthreshmatrix <- as.matrix(as.data.frame(dffinal))
    dffinalthreshmatrix[dffinalthreshmatrix<hammingthreshold] <- 0
    adj_matrix <- dffinalthreshmatrix
    diag(adj_matrix) <- 0 # Ensure no self-loops
    g <- graph_from_adjacency_matrix(adj_matrix, diag=FALSE, mode = "undirected", weighted = TRUE)
    if(isTRUE(graph.plot)) {
        png(file.path(results_folder, "graph.png"),width=8, height=8, units="in", res=500)
        plot(g,
                vertex.label = NA, #V(g)$name, # Use entity names as labels
                vertex.size = 1, 
                vertex.color = "skyblue",
                edge.width = E(g)$weight,
                edge.color = "black")
        dev.off()
    }

    # call clonotypes by grouping all connnected nodes (edges are kept with all hamming distances above hamming thresholds)
    clusters <- components(g)
    cluster_membership <- igraph::membership(clusters)

    # Print the edge IDs
   hammingdistclonotypes <- as.data.frame(cluster_membership)

   colnames(hammingdistclonotypes) <- c("hd.clonotypes")
   finalresults <- merge(totaldatafinal, hammingdistclonotypes, by.x="barcodev2", by.y="row.names")
   print(paste0("STATUS: number of unique clonotypes is ", length(unique(finalresults$hd.clonotypes))))
   contig.list_final <- split( finalresults , f = finalresults$sample_individual )
   write.csv(finalresults, file.path(results_folder, "finalresults.csv"), quote=FALSE, row.names=FALSE)

    dfvdjcounts <- finalresults  %>% group_by(sample_individual, hd.clonotypes) %>%
    dplyr::summarise(count = n())
    finalresults$CDR3.heavy_light.seq <- paste(finalresults$cdr3, finalresults$cdr3.light, sep="_")
    clonotypeswithseqs <- finalresults  %>% group_by(hd.clonotypes, CDR3.heavy_light.seq) %>%
    dplyr::summarise(count = n())
    clonotypeswithseqs <- as.data.frame(clonotypeswithseqs)
    dfvdjcounts <- as.data.table(dfvdjcounts)
    dfvdjcountspivot <- dfvdjcounts %>% dcast(sample_individual ~ hd.clonotypes, value.var="count")
    dfvdjcountspivot <- as.data.frame(dfvdjcountspivot)
    dfvdjcountspivot[is.na(dfvdjcountspivot)] <- 0
    rownames(dfvdjcountspivot) <- dfvdjcountspivot$sample_individual
    dfvdjcountspivot$sample_individual <- NULL
    tmppivot <- t(dfvdjcountspivot)
    if (ncol(tmppivot)>1) {
        tmppivot <- tmppivot[,order(colnames(tmppivot))]
        }
    A <- rowSums(tmppivot)
    tmppivot <- cbind(tmppivot,A)
    colnames(tmppivot)[length(colnames(tmppivot))] <- c("sum")
    seq <- c()
    for (c in rownames(tmppivot)) {
        maxvalue <- max(clonotypeswithseqs[clonotypeswithseqs$hd.clonotypes==c,"count"])
        aaseq <- clonotypeswithseqs[(clonotypeswithseqs$hd.clonotypes==c & clonotypeswithseqs$count==maxvalue),"CDR3.heavy_light.seq"] 
        seq <- c(seq, aaseq[1])
    }
    tmppivot <- as.data.frame(tmppivot)
    tmppivot$aa.CDR3.heavy_light.seq <-seq
    tmppivot <- tmppivot[rev(order(tmppivot$sum)),]
    write.csv(tmppivot, file.path(results_folder, "clonotypecount.csv"), quote=FALSE)
    write.csv(clonotypeswithseqs,file.path(results_folder, "clonotypeCDR3sequences.csv"), quote=FALSE,row.names=FALSE)
   return(list("graph"=g, "BCRclonotypes"=contig.list_final))
}


hamming_dist_cluster_plot <- function(dffinal, results_folder) {
    dffinal[is.na(dffinal)] <- 0 
    dffinal <- as.matrix(data.frame(dffinal, check.names = FALSE))
    class(dffinal) <- "numeric"
    row_hclust <- hclust(as.dist(dffinal), method = "complete") # "average" is a common method
    col_hclust <- hclust(as.dist(dffinal), method = "complete") # "average" is a common method
    col_fun <- colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("white", "orange", "red", "#bb0202", "darkred","#500000" ))

    hmC <- Heatmap(dffinal,
        name = "Hamming\ndist.",
        col = col_fun,
        show_row_names = FALSE,
        show_column_names = FALSE,
        cluster_columns = col_hclust,
        cluster_rows = row_hclust,
        cluster_column_slices = FALSE,
        border = F,
        heatmap_legend_param = list(
            legend_height = unit(3, "cm"),
            labels_gp = gpar(fontsize = 10),
            title_gp = gpar(fontsize = 10, fontface = "bold")
        ), use_raster = TRUE,  raster_quality = 5
    )
    png(file.path(results_folder, "hammingdist.png"), width = 8, height = 8, units = "in", res = 400)
    heatmapCOV <- draw(hmC,
        heatmap_legend_side = "left"
    )
    dev.off()
}

build_data_table <- function(totaldata, seqofinterest) {
    totaldatafinal <- c()
    for (b in unique(totaldata$barcodev2)) {
        tmp <- totaldata[totaldata$barcodev2==b, ]
        if(nrow(tmp)==1) {
            print("only one row... skipping .. ")
        }else if (nrow(tmp)==2) {
            chains <- unique(tmp$chain)
            for (c in chains) {
                if(c=="IGH") {
                    IGH <- tmp[tmp$chain==c,  ]
                } else {
                    IGL <- tmp[tmp$chain==c,  ]
                }
            }
            dftmp <- c(IGH, IGL)
            dftmp <- as.data.frame(dftmp)
            totaldatafinal <- rbind(totaldatafinal, dftmp)
        } else if (nrow(tmp) > 2) {
            chains <- unique(tmp$chain)
            if (("IGL" %in% chains) & ("IGK" %in% chains)) {
                d <- tmp[tmp$chain %in% c("IGL","IGK"),  ]
                if(nrow(d) > 1) {
                    x = max(d$umis)
                    v = d[d$umis==x, ]
                    if (nrow(v) >1) {
                        x = max(v$reads)
                        v = v[v$reads==x, ]
                    }
                } else { 
                    v=d
                }
                digl <- v
                d <- tmp[tmp$chain=="IGH",  ]
                if(nrow(d) > 1) {
                    x = max(d$umis)
                    v = d[d$umis==x, ]
                    if (nrow(v) >1) {
                        x = max(v$reads)
                        v = v[v$reads==x, ]
                    }
                } else { 
                    v=d
                }
                digh <- v
            } else {
                for (c in chains) {
                    d <- tmp[tmp$chain==c,  ]
                    if(nrow(d) > 1) {
                        x = max(d$umis)
                        v = d[d$umis==x, ]
                        if (nrow(v) >1) {
                            x = max(v$reads)
                            v = v[v$reads==x, ]
                        }
                    } else { 
                        v=d
                    }
                    if(c=="IGH") {
                        digh <- v
                    } else {
                        digl <- v
                    }
                }
            }
            dftmp <- c( digh, digl)
            dftmp <- as.data.frame(dftmp)
            totaldatafinal <- rbind(totaldatafinal, dftmp)
        }
    }

    colnames(totaldatafinal) <- str_replace_all(colnames(totaldatafinal), "\\.1","\\.light")
    totaldatafinal$barcodev2.light <- NULL
    totaldatafinal$sample.light <- NULL
    totaldatafinal$sample_individual.light <- NULL
    return(totaldatafinal)
}