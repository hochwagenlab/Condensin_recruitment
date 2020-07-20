library(EnrichedHeatmap)
library(circlize)

if(!(dir.exists("DSB_pdf"))) {
dir.create("DSB_pdf")
}

# prepare summits data
hotspots <- hwglabr2::get_dsb_hotspots(genome='SK1Yue')

hotspots_sorted <- sort(hotspots, by =~ score, decreasing = T)

midpoint <- floor(width(hotspots_sorted) / 2)
start(hotspots_sorted) <- start(hotspots_sorted) + midpoint
end(hotspots_sorted) <- start(hotspots_sorted)

#mcols(hotspots_sorted) <- DataFrame(class=c(rep(1:5, each=length(hotspots_sorted)/5),5))
mcols(hotspots_sorted) <- DataFrame(class=c(rep(1:3, each=length(hotspots_sorted)/3),3,3))

# bedgraph files to analyze
bedgraphs <- list.files(path="~/Desktop/MACS2_FE_reps",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
#bedgraph_file <- "~/Desktop/YueSK1/AH6408I-144-148-183-185-Reps-SK1Yue-B3W3-MACS2/AH6408I-144-148-183-185-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
bedgraph_file <- bedgraphs[i]

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][8],split="_")[[1]][1]

Cond_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
genAvg <- hwglabr2::average_chr_signal(Cond_bg)$genome_avrg
Cond_bg$score <- Cond_bg$score/genAvg

mat1 <- normalizeToMatrix(Cond_bg, hotspots_sorted, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)


if (grepl("Cond",bedgraph_file)) { # Condensin files
		if(!(exists("col_fun1"))) {
			col_fun1 <- colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
		}
		col_fun <- col_fun1
		YLIM <- c(0.85,1.25)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
		if(!(exists("col_fun2"))) {
		col_fun2 <- colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
		}
		col_fun <- col_fun2
		YLIM <- c(0.9,1.5)	
} else { # Red1 files
	if(!(exists("col_fun3"))) {
		col_fun3 <- colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("skyblue", "aliceblue","white", "pink2","deeppink4"))
		}
		col_fun <- col_fun3
		YLIM <- c(0.6,1.3)
}

pdf(paste0("DSB_pdf/",naming,"_around_DSB_by_hotness.pdf"))
# plot sorted by hotspot hotness
partition <- hotspots_sorted$class
print( EnrichedHeatmap(mat1, col = col_fun, name = "Signal", row_title_rot = 0, pos_line=FALSE,
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 1:3),
                                                             show_error = TRUE,ylim=YLIM,pos_line=FALSE,)),
    top_annotation_height = unit(2, "cm"), 
   split=hotspots_sorted$class,
    axis_name = c("-1.5 kb", "DSB hotspots", "1.5 kb"),
    column_title = naming)+
    Heatmap(partition, col = structure(1:5, names = as.character(1:5)), name = "",
              show_row_names = FALSE, width = unit(3, "mm"))
)
    dev.off()
    
}
