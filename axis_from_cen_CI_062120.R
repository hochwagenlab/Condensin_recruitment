#library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ggplot2)

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
bedgraphs <- list.files(path="C:/Users/tmark/Documents/heatmaps",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
    #bedgraph_file <- "~/Desktop/YueSK1/AH6408I-144-148-183-185-Reps-SK1Yue-B3W3-MACS2/AH6408I-144-148-183-185-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
    bedgraph_file <- bedgraphs[i]
    

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][7],split="_")[[1]][1]

Cond_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=F)
genAvg <- hwglabr2::average_chr_signal(Cond_bg)$genome_avrg
Cond_bg$score <- Cond_bg$score/genAvg

mat1 <- normalizeToMatrix(Cond_bg, summits_1pos, value_column = "score",
                         extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)

if (grepl("Cond",bedgraph_file)) { # Condensin files
    YLIM2 <- c(0.8,5)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
    YLIM2 <- c(0.8,5)
} else { # Red1 files
#    YLIM2 <- c(0.7,6)
}


out <- data.frame()
for (j in 1:2) {
    if (j == 1) {
        SetName <- "cen"
    } else {
        SetName <- "ncen"
    }
    rows <- which(summits_1pos$class == SetName)
    tmpmat <- mat1[rows,]
    
    data <- hwglabr2::signal_mean_and_ci(signal_data=tmpmat,
                                         ci=0.95, rep_bootstrap=1000,
                                         na_rm=TRUE)
    
    data <- data.frame(Set=rep(SetName,nrow(data)),Position=seq(-1500, 1500-1,by=50), data)
    out <- rbind(out,data)
}

cbPalette <- c("#009E73", "#CC79A7", "#F0E442", "#E69F00", "#56B4E9", "#0072B2", "#D55E00")

pdf(paste0("axis_pdf/",naming,"_around_axis_by_distance_from_centromere_CI.pdf"))
p <- ggplot(out, aes(x=Position,y=Mean,group=Set,fill=Set))
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