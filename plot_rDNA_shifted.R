setwd("rDNA_shifted/")

plotPrep <- function(rDNA, tilesize) {
  rDNA <- sort(GenomeInfoDb::sortSeqlevels(rDNA))
  
  seqlens <- 9242
  names(seqlens) <- "rDNA"
  
  bins <- GenomicRanges::tileGenome(seqlens, tilewidth=tilesize,
                                    cut.last.tile.in.chrom=TRUE)
  
  rDNAscore <- GenomicRanges::coverage(rDNA, weight="score")
  bins <- GenomicRanges::binnedAverage(bins, rDNAscore, "rDNA_score")
  positions <- bins@ranges@start + floor(bins@ranges@width / 2)
  out <- data.frame(position=positions, rDNA=bins$rDNA_score)
  return(out)
}

AH6408I <- rtracklayer::import.bedGraph("AH6408I-144-183-rDNA_shift_W3_MACS2_FE.bdg.gz")
AH6408Ip <- plotPrep(AH6408I,tilesize=2)

AH7797I <- rtracklayer::import.bedGraph("AH7797I-V5-185-148-rDNA_shift_W3_MACS2_FE.bdg.gz") 
AH7797Ip <- plotPrep(AH7797I,tilesize=2)
YLIM <- c(0,max(c(AH6408Ip$rDNA,AH7797I$rDNA)))

pdf("AH6408I-rDNA_shifted.pdf")
plot(AH7797Ip$position,AH7797Ip$rDNA, xlab='rDNA', ylab='Average signal', type='l', lwd=2,col="darkgrey",ylim=YLIM,main="AH6408I")
points(AH6408Ip$position,AH6408Ip$rDNA, type='l', lwd=2,col="darkblue")
abline(h=1,col="grey60")
hwglabr::plot_gene_arrow(geneEnd=5508,geneLength=8903-5508,orientation=-1,yPos=0) #RDN25
hwglabr::plot_gene_arrow(geneEnd=562,geneLength=2361-562,orientation=-1,yPos=0) #RDN18
hwglabr::plot_gene_arrow(geneEnd=4304,geneLength=4424-4304,orientation=-1,yPos=0) #RDN5
dev.off()

AH6408B <- rtracklayer::import.bedGraph("AH6408B-28-37-77-rDNA_shift_W3_MACS2_FE.bdg.gz")
AH7797B <- rtracklayer::import.bedGraph("AH7797B-noreson-rDNA_shift_W3_MACS2_FE.bdg.gz")
AH6408Bp <- plotPrep(AH6408B,tilesize=2)
AH7797Bp <- plotPrep(AH7797B,tilesize=2)
YLIM <- c(0,max(c(AH6408Bp$rDNA,AH7797Bp$rDNA)))

pdf("AH6408B-rDNA_shifted.pdf")
plot(AH7797Bp$position,AH7797Bp$rDNA, xlab='rDNA', ylab='Average signal', type='l', lwd=2,col="darkgrey",ylim=YLIM,main="AH6408B")
points(AH6408Bp$position,AH6408Bp$rDNA, type='l', lwd=2,col="darkblue")
abline(h=1,col="grey60")
hwglabr::plot_gene_arrow(geneEnd=5508,geneLength=8903-5508,orientation=-1,yPos=0) #RDN25
hwglabr::plot_gene_arrow(geneEnd=562,geneLength=2361-562,orientation=-1,yPos=0) #RDN18
hwglabr::plot_gene_arrow(geneEnd=4304,geneLength=4424-4304,orientation=-1,yPos=0) #RDN5
dev.off()

AH6408K <- rtracklayer::import.bedGraph("AH6408k-PK9-751-375-rDNA_shift_W3_MACS2_FE.bdg.gz")
AH7797K <- rtracklayer::import.bedGraph("AH7797K-PK9-755-377-rDNA_shift_W3_MACS2_FE.bdg.gz")
AH6408Kp <- plotPrep(AH6408K,tilesize=2)
AH7797Kp <- plotPrep(AH7797K,tilesize=2)
YLIM <- c(0,max(c(AH6408Kp$rDNA,AH7797Kp$rDNA)))

pdf("AH6408K-rDNA_shifted.pdf")
plot(AH7797Kp$position,AH7797Kp$rDNA, xlab='Position on rDNA', ylab='Average signal', type='l', lwd=2,col="darkgrey",ylim=YLIM,main="AH6408K")
points(AH6408Kp$position,AH6408Kp$rDNA, type='l', lwd=2,col="darkblue")
abline(h=1,col="grey60")
hwglabr::plot_gene_arrow(geneEnd=5508,geneLength=8903-5508,orientation=-1,yPos=0) #RDN25
hwglabr::plot_gene_arrow(geneEnd=562,geneLength=2361-562,orientation=-1,yPos=0) #RDN18
hwglabr::plot_gene_arrow(geneEnd=4304,geneLength=4424-4304,orientation=-1,yPos=0) #RDN5
dev.off()
