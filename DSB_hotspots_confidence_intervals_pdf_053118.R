library(EnrichedHeatmap)
library(ggplot2)

if(!(dir.exists("DSB_pdf"))) {
dir.create("DSB_pdf")
}

nGroups <- 3

# prepare summits data
hotspots <- hwglabr2::get_dsb_hotspots(genome='SK1Yue')

hotspots_sorted <- sort(hotspots, by =~ score, decreasing = T)

midpoint <- floor(width(hotspots_sorted) / 2)
start(hotspots_sorted) <- start(hotspots_sorted) + midpoint
end(hotspots_sorted) <- start(hotspots_sorted)


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
				YLIM <- c(0.8,1.4)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
				YLIM <- c(0.8,1.45)	
} else { # Red1 files
				YLIM <- c(0.6,1.3)
}

data <- hwglabr2::signal_mean_and_ci(signal_data=mat1, ci=0.95, rep_bootstrap=1000, na_rm=TRUE)
data2<- data.frame(Position=seq(-1500,1500-1,by=50),data)

pdf(paste0("DSB_pdf/",naming,"_around_DSB_CI.pdf"))
p <- ggplot(data2, aes(x=Position,y=Mean))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
theme_classic() +
ylim(YLIM) +
labs( x = "Distance to hotspot midpoint (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()

out <- data.frame()
for (j in 1:nGroups) {
	starts <- round(seq(1,nrow(mat1),length.out=nGroups+1))
	tmpmat <- mat1[starts[j]:(starts[j+1]-1),]

	data <- hwglabr2::signal_mean_and_ci(signal_data=tmpmat,
                                                  ci=0.95, rep_bootstrap=1000,
                                                  na_rm=TRUE)
maxScore <- round(hotspots_sorted$score[starts[j]]/1000)
minScore <- round(hotspots_sorted$score[(starts[j+1]-1)]/1000)

data <- data.frame(Quantile=rep(paste0(maxScore,"-",minScore),nrow(data)),Position=seq(-1500, 1500-1,by=50), data)
out <- rbind(out,data)
}

out <- dplyr::group_by(out, Quantile)
cbPalette <- c("#009E73", "#CC79A7", "#F0E442", "#E69F00", "#56B4E9", "#0072B2", "#D55E00")

pdf(paste0("DSB_pdf/",naming,"_around_DSB_by_hotness_CI_",nGroups,"groups.pdf"))
p <- ggplot(out, aes(x=Position,y=Mean,group=Quantile,fill=Quantile))
# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3) +
scale_fill_manual(name = "Hotness (10^3)", values = cbPalette) +
theme_classic() +
ylim(YLIM) +
labs( x = "Distance to hotspot midpoint (bp)", y = "Average ChIP-seq signal") 
# Add signal line
print(p + geom_line())
dev.off()
}
