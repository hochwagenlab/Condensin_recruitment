library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

if(!(dir.exists("axis_pdf"))) {
dir.create("axis_pdf")
}

# prepare summits data
summits <- hwglabr2::get_Red1_summits(genome='SK1Yue')

summits_sorted <- sort(summits, by =~ score, decreasing = T)
summits_1pos <- summits_sorted
start(summits_1pos) <- end(summits_1pos)

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
		YLIM2 <- c(0.7,2)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
		if(!(exists("col_fun2"))) {
		col_fun2 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("gold", "white", "blue4"))
		}
		col_fun <- col_fun2
		YLIM <- c(0.9,1.8)	
		YLIM2 <- c(0.7,2.2)
} else { # Red1 files
	if(!(exists("col_fun3"))) {
		col_fun3 <- colorRamp2(quantile(mat1, c( 0,0.5, 0.95),na.rm=T), c("gold", "white", "blue4"))
		}
		col_fun <- col_fun3
		YLIM <- c(0.9,3.25)
		YLIM2 <- c(0.7,6)
}

# plot sorted by hotspot hotness
pdf(paste0("axis_pdf/",naming,"_around_axis_by_hotness_divgenavg.pdf"))
print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal", row_title_rot = 0, pos_line=FALSE,
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = "blue4"),
                                                             show_error = TRUE,ylim=YLIM, pos_line=FALSE)),
    top_annotation_height = unit(2, "cm"), 
    row_order = 1:length(summits_sorted),
    axis_name = c("-1.5 kb", "Axis sites", "1.5 kb"),
    column_title = naming))
    dev.off()
########## 
    
data <- hwglabr2::signal_mean_and_ci(signal_data=mat1, ci=0.95, rep_bootstrap=1000, na_rm=TRUE)
data2<- data.frame(Position=seq(-1500,1500-1,by=50),data)

pdf(paste0("axis_pdf/",naming,"_around_axis_CI.pdf"))
p <- ggplot(data2, aes(x=Position,y=Mean))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
theme_classic() +
ylim(YLIM) +
labs( x = "Distance to axis site (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()
##############
nGroups <- 3

out <- data.frame()
for (j in 1:nGroups) {
	starts <- round(seq(1,nrow(mat1),length.out=nGroups+1))
	tmpmat <- mat1[starts[j]:(starts[j+1]-1),]

	tmp <- hwglabr2::signal_mean_and_ci(signal_data=tmpmat,
                                                  ci=0.95, rep_bootstrap=1000,
                                                  na_rm=TRUE)
maxScore <- round(summits_1pos$score[starts[j]])
minScore <- round(summits_1pos$score[(starts[j+1]-1)])

tmp <- data.frame(Quantile=rep(paste0(maxScore,"-",minScore),nrow(data)),Position=seq(-1500, 1500-1,by=50), tmp)
out <- rbind(out,tmp)
}

out <- dplyr::group_by(out, Quantile)
cbPalette <- c("#009E73", "#CC79A7", "#F0E442", "#E69F00", "#56B4E9", "#0072B2", "#D55E00")

pdf(paste0("axis_pdf/",naming,"_around_axis_by_hotness_CI_3groups.pdf"))
p <- ggplot(out, aes(x=Position,y=Mean,group=Quantile,fill=Quantile))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
scale_fill_manual(name = "Signal", values = cbPalette) +
theme_classic() +
ylim(YLIM2) +
labs( x = "Distance to axis site (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()

}

