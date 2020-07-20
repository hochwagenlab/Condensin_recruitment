library(GenomicRanges)
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

extend <- 5

conv3S <- conv
end(conv3S) <- start(conv3S) + extend
conv3E <- conv
start(conv3E) <- end(conv3E) - extend
conv3 <- c(conv3E,conv3S)
conv3$score <- 1
div3S <- div
end(div3S) <- start(div3S) + extend
div3E <- div
start(div3E) <- end(div3E) - extend
div3 <- c(div3E,div3S)
div3$score <- 1
postand4S <- postand2
end(postand4S) <- start(postand4S) + extend
postand4E <- postand2
start(postand4E) <- end(postand4E) - extend
postand4 <- c(postand4E,postand4S)
postand4$score <- 1


convmat <- normalizeToMatrix(conv3, conv2, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = extend)
divmat <- normalizeToMatrix(div3, div2, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = extend)
posmat <- normalizeToMatrix(postand4, postand3, value_column = "score",
                          extend = 1500, mean_mode = "weighted", w = extend)

col_fun <- colorRamp2(c(0,1), c("white","black"))

pdf(paste0("intergen_pdf/convergent_regions_sortedbygenedist_",extend,".pdf"))
print(EnrichedHeatmap(convmat, col = col_fun, name = "Signal",
row_title_rot = 0, pos_line=FALSE,
 top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "convergent regions", "1.5 kb"),
    row_order = order(width(conv)),
        column_title = "convergent"))
dev.off()
pdf(paste0("intergen_pdf/divergent_regions_sortedbygenedist_",extend,".pdf"))
print(EnrichedHeatmap(divmat, col = col_fun, name = "Signal",
 row_title_rot = 0, pos_line=FALSE,
  top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "divergent regions", "1.5 kb"),
    row_order = order(width(div)),
            column_title = "divergent"))
dev.off()
pdf(paste0("intergen_pdf/alltandem_regions_sortedbygenedist_",extend,".pdf"))
print(EnrichedHeatmap(posmat, col = col_fun, name = "Signal",
row_title_rot = 0, pos_line=FALSE,
 top_annotation_height = unit(0, "cm"), 
    axis_name = c("-1.5 kb", "positive tandem regions", "1.5 kb"),
    row_order = order(width(postand2)),
            column_title = "positive tandem"))
dev.off()

## removing all white values in illustrator