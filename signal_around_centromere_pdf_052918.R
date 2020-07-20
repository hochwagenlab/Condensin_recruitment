library(EnrichedHeatmap)
library(circlize)
library(GenomicRanges)
library(ggplot2)

if(!(dir.exists("centromere_pdf"))) {
dir.create("centromere_pdf")
}

# prepare centromere data
centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)

midpoint <- floor(width(centromeres) / 2)
start(centromeres) <- start(centromeres) + midpoint
end(centromeres) <- start(centromeres)

# bedgraph files to analyze
bedgraphs <- list.files(path="~/heatmaps/",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)

for (i in 1:length(bedgraphs)) {
bedgraph_file <- bedgraphs[i]

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][8],split="_")[[1]][1]

Cond_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=F)
genAvg <- hwglabr2::average_chr_signal(Cond_bg)$genome_avrg
Cond_bg$score <- Cond_bg$score/genAvg

j=2
extend <- 30000
windowL <- 100

mat1 <- normalizeToMatrix(Cond_bg, centromeres, value_column = "score",
                          extend = extend, mean_mode = "weighted", w = windowL,empty_value=NA)

if (grepl("Cond",bedgraph_file)) { # Condensin files
if(j==1) {
		if(!(exists("col_funA"))) {
			col_funA <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun1 <- col_funA
		} else {
			if(!(exists("col_funB"))) {
			col_funB <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun2 <- col_funB
		} 
		YLIM <- c(0.25,9)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
		if(j==1) {
		if(!(exists("col_funC"))) {
		col_funC <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun1 <- col_funC
		} else {
		if(!(exists("col_funD"))) {
		col_funD <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun2 <- col_funD
		}
		YLIM <- c(0.25,9)	
} else { # Red1 files
	if(j==1) {
	if(!(exists("col_funE"))) {
		col_funE <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun1 <- col_funE
		} else {
		if(!(exists("col_funF"))) {
		col_funF <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("lemonchiffon", "white", "darkgreen"))
		}
		col_fun2 <- col_funF
		}
		YLIM <- c(0,2.6)
}

if (j==1) {
	col_fun <- col_fun1
} else {
	col_fun <- col_fun2
}

pdf(paste0("centromere_pdf/",naming,"_around_centromere_",extend,"_heatmap.pdf"))

print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal", pos_line=FALSE,
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "darkgreen"),
                                                             show_error = TRUE,ylim=YLIM, pos_line=FALSE)),
    #top_annotation_height = unit(2, "cm"), 
    row_title_rot = 0,
    row_order = 1:length(centromeres),
    axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
    column_title = naming))
dev.off()

data <- hwglabr2::signal_mean_and_ci(signal_data=mat1, ci=0.95, rep_bootstrap=1000, na_rm=TRUE)
data2<- data.frame(Position=seq(-(extend),extend-1,by=windowL),data)

pdf(paste0("centromere_pdf/",naming,"_around_centromere_",extend,"_CI.pdf"))
p <- ggplot(data2, aes(x=Position,y=Mean))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
theme_classic() +
ylim(YLIM) +
labs( x = "Distance to centromere (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()

}

