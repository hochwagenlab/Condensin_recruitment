rec8 <- hwglabr2::import_bedGraph("~/Desktop/Reps/Condensin/AH7660I-164-330-Reps-SK1Yue-B3W3-MACS2/AH7660I-164-330-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
WT <- hwglabr2::import_bedGraph("~/Desktop/Reps/Condensin/AH6408I-144-183-Reps-SK1Yue-B3W3-MACS2/AH6408I-144-183-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)
asyn <- hwglabr2::import_bedGraph("~/Desktop/Reps/Condensin/AH6408K-375-751-Reps-SK1Yue-B3W3-MACS2/AH6408K-375-751-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz",local_copy=F)

# normalize against genome average
rec8genAvg <- hwglabr2::average_chr_signal(rec8)$genome_avrg
rec8$score <- rec8$score/rec8genAvg

WTgenAvg <- hwglabr2::average_chr_signal(WT)$genome_avrg
WT$score <- WT$score/WTgenAvg

asyngenAvg <- hwglabr2::average_chr_signal(asyn)$genome_avrg
asyn$score <- asyn$score/asyngenAvg

# incorporate genome info
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

asyn <- sort(GenomeInfoDb::sortSeqlevels(asyn))
rec8 <- sort(GenomeInfoDb::sortSeqlevels(rec8))
WT <- sort(GenomeInfoDb::sortSeqlevels(WT))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
GenomeInfoDb::seqlengths(asyn) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(rec8) <- GenomeInfoDb::seqlengths(genome_info)
GenomeInfoDb::seqlengths(WT) <- GenomeInfoDb::seqlengths(genome_info)

# bin data and combine into a mega table
tilesize =5000
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(asyn),
                                  tilewidth=tilesize,
                                  cut.last.tile.in.chrom=TRUE)

asynscore <- GenomicRanges::coverage(asyn, weight="score")
rec8score <- GenomicRanges::coverage(rec8, weight="score")
WTscore <- GenomicRanges::coverage(WT, weight="score")
bins <- GenomicRanges::binnedAverage(bins, asynscore, "asyn_score")
bins <- GenomicRanges::binnedAverage(bins, WTscore, "WT_score")
bins <- GenomicRanges::binnedAverage(bins, rec8score, "rec8_score")

chr <- "chrXII"
bins2 <- GenomeInfoDb::keepSeqlevels(bins, chr)
positions <- bins2@ranges@start + floor(bins2@ranges@width / 2)
dfChrXII <- data.frame(position=positions / 1000, asyn=bins2$asyn_score, WT=bins2$WT_score,rec8D=bins2$rec8_score)

pdf(paste0("Smc4PK9_asynWTrec8D_entirechrXII_",tilesize,"bpwindow.pdf"))
plot(dfChrXII$position,dfChrXII$asyn, xlab='Position on chrXII (Kb)', ylab='signal', type='l', lwd=2,col="darkgreen",ylim=c(0.5,3))
points(dfChrXII$position,dfChrXII$WT, type='l', lwd=2,col="black")
points(dfChrXII$position,dfChrXII$rec8D, type='l', lwd=2,col="salmon")
abline(h=1,col="grey60")
legend(800,3,col=c("darkgreen","black","salmon"),pch=19,legend=c("asynchronous","wild type","rec8D"))
dev.off()
