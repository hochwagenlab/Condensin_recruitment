nrCond <- hwglabr2::import_bedGraph("AH6408-3h-notreson-PK9-28-37-77-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
rCond <- hwglabr2::import_bedGraph("AH6408I-144-183-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
nrCl <- hwglabr2::import_bedGraph("AH7797B-v5-140401-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
rCl <- hwglabr2::import_bedGraph("AH7797I-148-185-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)

# normalize against genome average
nrCondgenAvg <- hwglabr2::average_chr_signal(nrCond)$genome_avrg
nrCond$score <- nrCond$score/nrCondgenAvg
rCondgenAvg <- hwglabr2::average_chr_signal(rCond)$genome_avrg
rCond$score <- rCond$score/rCondgenAvg

nrClgenAvg <- hwglabr2::average_chr_signal(nrCl)$genome_avrg
nrCl$score <- nrCl$score/nrClgenAvg
rClgenAvg <- hwglabr2::average_chr_signal(rCl)$genome_avr
rCl$score <- rCl$score/rClgenAvg

# incorporate genome info
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

nrCl <- sort(GenomeInfoDb::sortSeqlevels(nrCl))
rCl <- sort(GenomeInfoDb::sortSeqlevels(rCl))
nrCond <- sort(GenomeInfoDb::sortSeqlevels(nrCond))
rCond <- sort(GenomeInfoDb::sortSeqlevels(rCond))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
GenomeInfoDb::seqlengths(nrCl) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(rCl) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(nrCond) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(rCond) <- GenomeInfoDb::seqlengths(genome_info)

# bin data and combine into a mega table
#tilesize =3000
tilesize =750
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(nrCl),
                                  tilewidth=tilesize,
                                  cut.last.tile.in.chrom=TRUE)

nrClscore <- GenomicRanges::coverage(nrCl, weight="score")
rClscore <- GenomicRanges::coverage(rCl, weight="score")
nrCondscore <- GenomicRanges::coverage(nrCond, weight="score")
rCondscore <- GenomicRanges::coverage(rCond, weight="score")
bins <- GenomicRanges::binnedAverage(bins, nrClscore, "nrCl_score")
bins <- GenomicRanges::binnedAverage(bins, rClscore, "rCl_score")
bins <- GenomicRanges::binnedAverage(bins, nrCondscore, "nrCond_score")
bins <- GenomicRanges::binnedAverage(bins, rCondscore, "rCond_score")

chr <- "chrII"
startPos <- 630
endPos <- 800
bins2 <- GenomeInfoDb::keepSeqlevels(bins, chr, pruning.mode="coarse")
rows <- which((bins2@ranges@start/1000 > startPos) & (bins2@ranges@start/1000 < endPos))
bins3 <- bins2[rows,]
positions <- bins3@ranges@start + floor(bins3@ranges@width / 2)
dfChrII <- data.frame(position=positions / 1000, nrCl=bins3$nrCl_score, rCl=bins3$rCl_score, nrCond=bins3$nrCond_score, rCond=bins3$rCond_score)

pdf(paste0("Noreson_Smc4PK9_chrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$nrCl, xlab='Position on chrII (Kb)', ylab='PK9 signal', type='l', lwd=2,col="grey40",ylim=c(0.5,2))
points(dfChrII$position,dfChrII$nrCond, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
dev.off()

pdf(paste0("Reson_Smc4PK9_chrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$rCl, xlab='Position on chrII (Kb)', ylab='PK9 signal', type='l', lwd=2,col="grey40",ylim=c(0.5,2))
points(dfChrII$position,dfChrII$rCond, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
dev.off()

chr <- "chrII"
bins2 <- GenomeInfoDb::keepSeqlevels(bins, chr, pruning.mode="coarse")
positions <- bins2@ranges@start + floor(bins2@ranges@width / 2)
dfChrII <- data.frame(position=positions / 1000, nrCl=bins2$nrCl_score, rCl=bins2$rCl_score, nrCond=bins2$nrCond_score, rCond=bins2$rCond_score)

pdf(paste0("Noreson_Smc4PK9_entirechrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$nrCl, xlab='Position on chrII (Kb)', ylab='PK9 signal', type='l', lwd=2,col="grey40",ylim=c(0.5,2))
points(dfChrII$position,dfChrII$nrCond, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
dev.off()

pdf(paste0("Reson_Smc4PK9_entirechrII_",tilesize,"bpwindow.pdf"))
par(mar=c(5,4,4,4)+0.1)
plot(dfChrII$position,dfChrII$rCl, xlab='Position on chrII (Kb)', ylab='PK9 signal', type='l', lwd=2,col="grey40",ylim=c(0.5,2))
points(dfChrII$position,dfChrII$rCond, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
dev.off()
