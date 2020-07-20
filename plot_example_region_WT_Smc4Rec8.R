Rec8 <- hwglabr2::import_bedGraph("Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
Cond <- hwglabr2::import_bedGraph("AH6408I-144-183-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)

# normalize against genome average
CondgenAvg <- hwglabr2::average_chr_signal(Cond)$genome_avrg
Cond$score <- Cond$score/CondgenAvg

Rec8genAvg <- hwglabr2::average_chr_signal(Rec8)$genome_avrg
Rec8$score <- Rec8$score/Rec8genAvg

# incorporate genome info
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

Rec8 <- sort(GenomeInfoDb::sortSeqlevels(Rec8))
Cond <- sort(GenomeInfoDb::sortSeqlevels(Cond))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
GenomeInfoDb::seqlengths(Rec8) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(Cond) <- GenomeInfoDb::seqlengths(genome_info)

# bin data and combine into a mega table
#tilesize =3000
tilesize =750
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(Rec8),
                                  tilewidth=tilesize,
                                  cut.last.tile.in.chrom=TRUE)

Rec8score <- GenomicRanges::coverage(Rec8, weight="score")
Condscore <- GenomicRanges::coverage(Cond, weight="score")
bins <- GenomicRanges::binnedAverage(bins, Rec8score, "Rec8_score")
bins <- GenomicRanges::binnedAverage(bins, Condscore, "Cond_score")

chr <- "chrII"
startPos <- 630
endPos <- 800
bins2 <- GenomeInfoDb::keepSeqlevels(bins, chr, pruning.mode="coarse")
rows <- which((bins2@ranges@start/1000 > startPos) & (bins2@ranges@start/1000 < endPos))
bins3 <- bins2[rows,]
positions <- bins3@ranges@start + floor(bins3@ranges@width / 2)
dfChrII <- data.frame(position=positions / 1000, Rec8=bins3$Rec8_score, Cond=bins3$Cond_score)

pdf(paste0("WT_Rec8_Smc4PK9_chrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$Rec8, xlab='Position on chrII (Kb)', ylab='Rec8 signal', type='l', lwd=2,col="dodgerblue4",ylim=c(-2,5))
points(dfChrII$position,4.5*dfChrII$Cond-3.5, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
axis(4,at=(seq(0.5,2.5,by=0.5)*4.5)-3.5,labels=seq(0.5,2.5,by=0.5))
mtext(side=4,line=2.5,"Condensin signal")
dev.off()

chr <- "chrII"
bins2 <- GenomeInfoDb::keepSeqlevels(bins, chr, pruning.mode="coarse")
positions <- bins2@ranges@start + floor(bins2@ranges@width / 2)
dfChrII <- data.frame(position=positions / 1000, Rec8=bins2$Rec8_score, Cond=bins2$Cond_score)

pdf(paste0("WT_Rec8_Smc4PK9_entirechrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$Rec8, xlab='Position on chrII (Kb)', ylab='Rec8 signal', type='l', lwd=2,col="dodgerblue4",ylim=c(-4,3))
points(dfChrII$position,4.5*dfChrII$Cond-3.5, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
axis(4,at=(seq(0.5,2.5,by=0.5)*4.5)-3.5,labels=seq(0.5,2.5,by=0.5))
mtext(side=4,line=2.5,"Condensin signal")
dev.off()

