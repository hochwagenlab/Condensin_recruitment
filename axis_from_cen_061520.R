library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)

if(!(dir.exists("axis_pdf"))) {
dir.create("axis_pdf")
}

# prepare summits data
summits <- hwglabr2::get_Red1_summits(genome='SK1Yue')

centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
midpoint <- floor(width(centromeres) / 2)
start(centromeres) <- start(centromeres) + midpoint
end(centromeres) <- start(centromeres)

d <- distanceToNearest(summits,centromeres,ignore.strand=T)
mcols(summits) <- DataFrame(distance=mcols(d),score=summits$score)

centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
midpoint <- floor(width(centromeres) / 2)
start(centromeres) <- start(centromeres) + midpoint - 15000
end(centromeres) <- start(centromeres) + 30000

tmp <- intersect(summits,centromeres)
tmp$class <- "cen"
tmp2 <- setdiff(summits,centromeres)
tmp2$class <- "ncen"
tmp3 <- merge(summits,tmp,all=T)
tmp4 <- merge(tmp3,tmp2,all=T)

summits_sorted <- sort(tmp4, by =~ distance.distance, decreasing = F)
summits_1pos <- summits_sorted
start(summits_1pos) <- end(summits_1pos)

# bedgraph files to analyze
bedgraphs <- list.files(path="~/Desktop/MACS2_FE_reps",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
    #bedgraph_file <- "~/Desktop/YueSK1/AH6408I-144-148-183-185-Reps-SK1Yue-B3W3-MACS2/AH6408I-144-148-183-185-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
    bedgraph_file <- bedgraphs[i]
    

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][6],split="_")[[1]][1]

Cond_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=F)
genAvg <- hwglabr2::average_chr_signal(Cond_bg)$genome_avrg
Cond_bg$score <- Cond_bg$score/genAvg

mat1 <- normalizeToMatrix(Cond_bg, summits_1pos, value_column = "score",
                         extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)

########                         
if (grepl("Cond",bedgraph_file)) { # Condensin files
    if(!(exists("col_fun1"))) {
        col_fun1 <- colorRamp2(quantile(mat1, c( 0, 0.5, 0.95),na.rm=T), c("gold", "white", "blue4"))
    }
    col_fun <- col_fun1
    #	col_fun <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("royalblue", "white", "red"))
    YLIM <- c(0.8,1.5)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
    if(!(exists("col_fun2"))) {
        col_fun2 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("gold", "white", "blue4"))
    }
    col_fun <- col_fun2
    YLIM <- c(0.9,1.8)	
} else { # Red1 files
    if(!(exists("col_fun3"))) {
        col_fun3 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("gold", "white", "blue4"))
    }
    col_fun <- col_fun3
    YLIM <- c(0.9,3.25)
}

pdf(paste0("axis_pdf/",naming,"_around_axis_by_cendist_divgenavg.pdf"))
print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal", row_title_rot = 0,pos_line=FALSE,
        top_annotation_height = unit(2, "cm"), 
    row_order = 1:length(summits_sorted),
       split=summits_sorted$class,
    axis_name = c("-1.5 kb", "Axis sites", "1.5 kb"),
    column_title = naming))
dev.off()
}