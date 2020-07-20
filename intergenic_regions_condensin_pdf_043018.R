library(EnrichedHeatmap)
library(circlize)

if(!(dir.exists("intergen_pdf"))) {
dir.create("intergen_pdf")
}

intergen <- hwglabr2::get_intergenic_regions("SK1Yue",as_gr=T)
div <- intergen[intergen$type=="divergent",]
div2 <- div
midpoint <- floor(width(div2) / 2)
start(div2) <- start(div2) + midpoint
end(div2) <- start(div2)
conv <- intergen[intergen$type=="convergent",]
conv2 <- conv
midpoint <- floor(width(conv2) / 2)
start(conv2) <- start(conv2) + midpoint
end(conv2) <- start(conv2)
negtand <- intergen[which((intergen$type=="tandem")&(intergen$right_gene_strand=="-")),]
negtand2 <- negtand
strand(negtand2) <- "-"

postand <- intergen[which((intergen$type=="tandem")&(intergen$right_gene_strand=="+")),]
postand2 <- c(postand,negtand2)
postand3 <- postand2
midpoint <- floor(width(postand3) / 2)
start(postand3) <- start(postand3) + midpoint
end(postand3) <- start(postand3)

#########################################
# bedgraph files to analyze
bedgraphs <- list.files(path="~/Desktop/MACS2_FE_reps",pattern="MACS2_FE.bdg.gz",recursive=T,full.names=T)
bedgraphs <- bedgraphs[grep("SK1Yue-PM",bedgraphs)]

for (i in 1:length(bedgraphs)) {
bedgraph_file <- bedgraphs[i]

naming <- strsplit(strsplit(bedgraph_file,split="/")[[1]][8],split="_")[[1]][1]

Cond_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=F)
genAvg <- hwglabr2::average_chr_signal(Cond_bg)$genome_avrg
Cond_bg$score <- Cond_bg$score/genAvg

convmat <- normalizeToMatrix(Cond_bg, conv2, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)
divmat <- normalizeToMatrix(Cond_bg, div2, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)
posmat <- normalizeToMatrix(Cond_bg, postand3, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = 50,empty_value=NA)


if (grepl("Cond",bedgraph_file)) { # Condensin files
		if(!(exists("col_fun1"))) {
			col_fun1 <- colorRamp2(quantile(c(convmat,divmat,posmat), c(0,0.5, 0.95),na.rm=T), c("seagreen3", "white", "purple4"))
			#          0%         50%         95% 
			# 0.002444527 0.991892023 1.598845157 
		}
		col_fun <- col_fun1
		YLIM <- c(0.5,1.4)
} else if(grepl("Scc2",bedgraph_file)) { #Scc2 files
		if(!(exists("col_fun2"))) {
		col_fun2 <- colorRamp2(quantile(c(convmat,divmat,posmat), c(0,0.5, 0.95),na.rm=T), c("seagreen3", "white", "purple4"))
		}
		col_fun <- col_fun2
		YLIM <- c(0.8,1.4)	
} else { # Red1 files
	if(!(exists("col_fun3"))) {
		col_fun3 <- colorRamp2(quantile(c(convmat,divmat,posmat), c(0,0.5, 0.95),na.rm=T), c("seagreen3", "white", "purple4"))
		#          0%         50%         95% 
		# 0.001420697 0.782422081 2.825267901 
#		col_fun3 <- colorRamp2(quantile(c(convmat,divmat,posmat), c(0,0.65, 0.95),na.rm=T), c("seagreen3", "white", "purple4"))
		#		         0%         65%         95% 
		# 0.001420697 0.998256389 2.825267901 
		}
		col_fun <- col_fun3
		YLIM <- c(0.5,3)
}
    
pdf(paste0("intergen_pdf/",naming,"_around_convergent_regions_sortedbygenedist.pdf"))
print(EnrichedHeatmap(convmat, col = col_fun, name = "Signal",
 row_title_rot = 0, pos_line=FALSE,
 top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "convergent regions", "1.5 kb"),
    row_order = order(width(conv)),
    column_title = naming))
dev.off()
pdf(paste0("intergen_pdf/",naming,"_around_divergent_regions_sortedbygenedist.pdf"))
print(EnrichedHeatmap(divmat, col = col_fun, name = "Signal",
    row_title_rot = 0, pos_line=FALSE,
    top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "divergent regions", "1.5 kb"),
    row_order = order(width(div)),
    column_title = naming))
dev.off()
pdf(paste0("intergen_pdf/",naming,"_around_all_tandem_regions_sortedbygenedist.pdf"))
print(EnrichedHeatmap(posmat, col = col_fun, name = "Signal",
     row_title_rot = 0, pos_line=FALSE,
     top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "positive tandem regions", "1.5 kb"),
    row_order = order(width(postand2)),
        column_title = naming))
dev.off()
}
