library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

if(!(dir.exists("tRNA_pdf"))) {
dir.create("tRNA_pdf")
}

gffYue <- hwglabr2::get_gff(genome="SK1Yue")
tRNA <- gffYue[gffYue$type=="tRNA",]
tRNA <- tRNA[tRNA@seqnames!="chrMT",]

tRNAtotal_starts <- tRNA
start(tRNAtotal_starts) <- ifelse(strand(tRNA) == '+', start(tRNA), end(tRNA))
end(tRNAtotal_starts) <- start(tRNAtotal_starts)

bedgraphs <- list.files(path="~/Desktop/MACS2_FE_reps",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
#bedgraph_file <- "~/Desktop/YueSK1/AH6408I-144-148-183-185-Reps-SK1Yue-B3W3-MACS2/AH6408I-144-148-183-185-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
bedgraph_file <- bedgraphs[i]

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][8],split="_")[[1]][1]
if (is.na(naming)) {
	naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][9],split="_")[[1]][1]
}

Cond <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)

mat1 <- normalizeToMatrix(Cond, tRNAtotal_starts, value_column = "score",
                          extend = 1000, mean_mode = "w0", w = 25,empty_value=NA)
                          

if (grepl("Cond",bedgraph_file)) { # Condensin files
		if(!(exists("col_fun1"))) {
			col_fun1 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("thistle1", "white", "chocolate4"))
		}
		col_fun <- col_fun1
	#	col_fun <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("royalblue", "white", "red"))
		YLIM <- c(0.5,1.6)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
		if(!(exists("col_fun2"))) {
		col_fun2 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("thistle1", "white", "chocolate4"))
		}
		col_fun <- col_fun2
		YLIM <- c(0.9,1.7)	
} else { # Red1 files
	if(!(exists("col_fun3"))) {
		col_fun3 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("thistle1", "white", "chocolate4"))
		}
		col_fun <- col_fun3
		YLIM <- c(0.9,3.5)
}

pdf(paste0("tRNA_pdf/",naming,"_around_tRNA_bysize.pdf"))
print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE,
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                             show_error = TRUE,ylim=c(0.5,1.6),pos_line=FALSE)),
    top_annotation_height = unit(2, "cm"), row_title_rot = 0,	
    axis_name = c("-1 kb", "tRNA starts", "1 kb"),
    row_order = order(width(tRNA)),
    column_title =naming))
dev.off()

########## 
    
data <- hwglabr2::signal_mean_and_ci(signal_data=mat1, ci=0.95, rep_bootstrap=1000, na_rm=TRUE)
data2<- data.frame(Position=seq(-1000,1000-1,by=25),data)

pdf(paste0("tRNA_pdf/",naming,"_around_tRNA_CI.pdf"))
p <- ggplot(data2, aes(x=Position,y=Mean))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
theme_classic() +
ylim(YLIM) +
labs( x = "Distance from tRNA translational start site (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()


}
