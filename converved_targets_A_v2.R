setwd("~/Desktop/4-NGS/TargetScan/targetscan_70/")
targetscan_species <- read.table(file="Selected_Sequences.txt", header = FALSE, fill = TRUE, sep =  "\t")

colnames(targetscan_species) <- c("SYMBOL","species", "UTR", "nospace")
targetscan_species$nospace <- toupper(targetscan_species$nospace)
human_UTR <- unique(targetscan_species[which(targetscan_species$species==9606),])

setwd("~/Desktop/RNA-seq Acong/")
require(gdata)
ts27a3p <- read.xls ("TargetScan7.2__miR-27-3p.predicted_targets.xlsx", sheet = "TargetScan7", header = TRUE) #all targets 27a
ts27a3p$SYMBOL <- ts27a3p$Target.gene

#Potential G:U combinations on the seedless
#miR27a <- "TTCACAGTGGCTAAGTTCCGC"
#seedless <- substr(miR27a, nchar(miR27a)-6, nchar(miR27a))
#heptamers_mir27a_seedless <- toupper(c2s(rev(comp(s2c(seedless)))))
#heptamers_mir27a_seedless <- gsub("A","[AG]",gsub("C","[CT]",heptamers_mir27a_seedless)) #degenerate seedless
#heptamers_mir27a_seedless <- paste0("A",heptamers_mir27a_seedless)

#Potential G:U combinations on the seedless
nt_tested <- "A"
heptamers_mir27a_seedless <- expand.grid(n0=c(nt_tested), n1 = c("G"), n2 = c("C", "U"), n3 = c("G"), n4 = c("G"), n5 = c("A", "G"), n6 = c("A", "G"), n7 = c("C", "U"))
heptamers_mir27a_seedless$seedT <-paste0(heptamers_mir27a_seedless$n0,heptamers_mir27a_seedless$n1,heptamers_mir27a_seedless$n2,heptamers_mir27a_seedless$n3,heptamers_mir27a_seedless$n4,heptamers_mir27a_seedless$n5,heptamers_mir27a_seedless$n6, heptamers_mir27a_seedless$n7) 
require(stringdist)
heptamers_mir27a_seedless$dist <- stringdist(substr(heptamers_mir27a_seedless$seedT, 2, 8), "GCGGAAC")
heptamers_mir27a_seedless$GC <- 1-(nchar(gsub("C", "", gsub("G", "", heptamers_mir27a_seedless$seedT)))/nchar(heptamers_mir27a_seedless$seedT))

heptamers_mir27a_seedless <- heptamers_mir27a_seedless[which(heptamers_mir27a_seedless$dist<=3),]

toMatch <- unique(heptamers_mir27a_seedless$seedT)

human_UTR$t7 <- gregexpr(paste(toMatch,collapse="|"), human_UTR$nospace, perl = TRUE)
num_sites <- human_UTR[which(human_UTR$t7!="-1"),]
num_sites$num_sites <- nchar(num_sites$t7)-nchar(gsub(",", "", num_sites$t7))+1
num_sites <- num_sites[,c("SYMBOL","num_sites")]
num_sites<-aggregate(num_sites$num_sites, by=list(num_sites$SYMBOL), FUN=max, na.rm=FALSE)
colnames(num_sites) <- c("SYMBOL","num_sites")

human_candidates <- unique(human_UTR[which(human_UTR$t7!="-1"),1])

seedless27a_species <- targetscan_species[(targetscan_species$SYMBOL %in% human_candidates),]
seedless27a_species$t7 <- gregexpr(paste(toMatch,collapse="|"), seedless27a_species$nospace, perl = TRUE)
seedless27a_species <- unique(seedless27a_species[which(seedless27a_species$t7!="-1"),])
seedless27a_species$count <- 1
seedless27a_species$maxlen <- nchar(as.character(seedless27a_species$nospace))
len_aggregate<-aggregate(seedless27a_species$maxlen, by=list(seedless27a_species$SYMBOL, seedless27a_species$species), FUN=max, na.rm=FALSE)
colnames(len_aggregate) <- c("SYMBOL","species","maxlen")
seedless27a_species <- merge(seedless27a_species,len_aggregate,by=c("SYMBOL","species","maxlen"))
seedless27a_species_aggregate <-aggregate(seedless27a_species$count, by=list(seedless27a_species$SYMBOL), FUN=sum, na.rm=FALSE)
colnames(seedless27a_species_aggregate) <- c("SYMBOL", "Counts")
seedless27a_species_aggregate <- seedless27a_species_aggregate[order(-seedless27a_species_aggregate$Counts),]
seedless27a_species_aggregate <- merge(seedless27a_species_aggregate, num_sites, by="SYMBOL")
over10 <- seedless27a_species_aggregate[which(seedless27a_species_aggregate$Counts>=10),]
write.table(seedless27a_species_aggregate, "seedless27a_species_aggregate.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)
names_over10 <- unique(over10$SYMBOL)

UTMR_27a_only <- over10[ !(over10$SYMBOL %in% ts27a3p$SYMBOL), ]
UTMR_27a_only$group <- "UTMR_only"
UTMR_27a_only <- UTMR_27a_only[,c(1,4)]
both <- over10[ (over10$SYMBOL %in% ts27a3p$SYMBOL), ]
both$group <- "both"
both <- both[,c(1,4)]
TS_27a_only <- ts27a3p[ !(ts27a3p$SYMBOL %in% over10$SYMBOL), ]
TS_27a_only$group <- "TargetScan_only"
TS_27a_only <- TS_27a_only[,c(19,20)]

venn_27a <- rbind(UTMR_27a_only, both, TS_27a_only)
#write.table(venn_27a, "venn_27a2.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)

#Expression cumulative curves analysis

#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationDbi")
#biocLite("org.Hs.eg.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
#install.packages("gtools")
library(gtools)
total4 <- read.table("Acong-RNAseq-expression-table.txt", header = TRUE)
total4 <- na.omit(total4)

total4$logexpr <- log2(total4$WT_OE27a_)
hist(total4$logexpr)

threshold <- 1 #3cpm ~ 1FPKM
total4 <- total4[which(total4$WT_OE27a_>=threshold),]
hist(total4$logexpr)

#install.packages("gdata")
ts243p <- read.xls ("TargetScan7.2__miR-24-3p.predicted_targets.xlsx", sheet = "TargetScan7", header = TRUE) #all targets 27a
ts243p$SYMBOL <- ts243p$Target.gene
ts233p <- read.xls ("TargetScan7.2__miR-23-3p.predicted_targets.xlsx", sheet = "TargetScan7", header = TRUE) #all targets 27a
ts233p$SYMBOL <- ts233p$Target.gene

total4$FC27awt <- log2(total4$WT_OE27a_/total4$WT_OE21_U)
total4$FC27adko <- log2(total4$DKO_OE27a/total4$DKO_OE21_)
total4$FC27atko <- log2(total4$TKO_OE27a/total4$TKO_OE21_)
total4$FC_used <- total4$FC27awt

seedless27a <- merge(seedless27a_species_aggregate, total4, by="SYMBOL")

notargets <- total4[ !(total4$SYMBOL %in% seedless27a$SYMBOL), ]
notargets <- notargets[ !(notargets$SYMBOL %in% ts27a3p$SYMBOL), ]
notargets <- notargets[ !(notargets$SYMBOL %in% ts243p$SYMBOL), ]
notargets <- notargets[ !(notargets$SYMBOL %in% ts233p$SYMBOL), ]

both_sets <- seedless27a[ (seedless27a$SYMBOL %in% ts27a3p$SYMBOL), ]
both_sets <- both_sets[which(both_sets$Counts>=10),]

seedless27a <- seedless27a[ !(seedless27a$SYMBOL %in% ts27a3p$SYMBOL), ]
seedless27a <- seedless27a[ !(seedless27a$SYMBOL %in% ts243p$SYMBOL), ]
seedless27a <- seedless27a[ !(seedless27a$SYMBOL %in% ts233p$SYMBOL), ]

ts27a3p <- ts27a3p[ !(ts27a3p$SYMBOL %in% both_sets$SYMBOL), ]

breaks <-  seq(-5, 5, by=0.05)
#breaks <-  seq(-5, 5, by=0.1) #smooth curve

#Baseline genes TargetScan
par(mfrow=c(1,1))
datLog2FC <-  as.numeric(notargets$FC_used) #canonical FC
notargets1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(notargets))) 

plot(breaks, notargets1, type="n", xlab = "mRNA fold-change Log2", ylab = "Cumulative fraction", ylim=c(0,1), xlim=c(-2,2))
lines(breaks, notargets1)
avg_notarget1 <- format(round(notargets1[c(81,101)], 2), nsmall = 2)

ts27a3p <- merge(ts27a3p, total4, by="SYMBOL")
datLog2FC <-  as.numeric(ts27a3p$FC_used) #canonical FC without targets
targetscan1 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(ts27a3p))) 
lines(breaks, targetscan1, col="green")
avg_targetscan1 <- format(round(targetscan1[c(81,101)], 2), nsmall = 2)
ts27a3p$group <- "TargetScan"

filtered2 <- unique(seedless27a[,c("SYMBOL","WT_OE27a_", "FC27awt", "FC27adko", "FC_used")])
datLog2FC <-  as.numeric(filtered2$FC_used) #canonical FC without targets
allseedless = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(filtered2))) 
lines(breaks, allseedless, col="red")
avg_allseedless <- format(round(allseedless[c(81,101)], 2), nsmall = 2)

both_sets <- unique(both_sets[,c("SYMBOL","WT_OE27a_", "FC27awt", "FC27adko", "FC_used")])
datLog2FC <-  as.numeric(both_sets$FC_used) #canonical FC without targets
both_sets2 = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(both_sets))) 
lines(breaks, both_sets2, col="purple")
avg_both_sets2 <- format(round(both_sets2[c(81,101)], 2), nsmall = 2)


top_conserved <- seedless27a_species_aggregate[which(seedless27a_species_aggregate$Counts>=10),]
#top_conserved <- seedless27a_species_aggregate[which(seedless27a_species_aggregate$num_sites>=3),]
filteredA_con <- seedless27a[ (seedless27a$SYMBOL %in% top_conserved$SYMBOL), ]
datLog2FC <-  as.numeric(filteredA_con$FC_used) #canonical FC without targets
allseedlessAcon = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(filteredA_con))) 
lines(breaks, allseedlessAcon, col="blue")
avg_allseedlessAcon <- format(round(allseedlessAcon[c(81,101)], 2), nsmall = 2)
filteredA_con$group <- "UTMR"



filteredA_sum <- unique(filteredA_con[,c("SYMBOL","WT_OE27a_", "FC27awt", "FC27adko", "group")])
ts27a3p_sum <- unique(ts27a3p[,c("SYMBOL","WT_OE27a_", "FC27awt", "FC27adko", "group")])
sum_sum <- rbind(filteredA_sum,ts27a3p_sum)
sum_sum <- sum_sum[which(sum_sum$FC27awt<0),]
#write.table(sum_sum, "summary_plot_UTMR_targetscan.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)

#Legend
legend(-2, 1, legend=c(paste("No target", nrow(notargets), avg_notarget1[2], collapse="|"), 
                       paste("TS-27a", nrow(ts27a3p), avg_targetscan1[2], collapse="|"), 
                       paste("Seedless-27a", nrow(filtered2), avg_allseedless[2],collapse="|"), 
                       paste("Seedless-Acons-27a", nrow(filteredA_con), avg_allseedlessAcon[2], collapse="|"),
                       paste("Both_TS_SeedlessAcons-27a", nrow(both_sets), avg_both_sets2[2], collapse="|")),
       col=c("black", "green", "red", "blue", "purple"), lty=1, cex=0.7)

summary <- rbind(breaks, notargets1, targetscan1, allseedless,allseedlessAcon,both_sets2)
write.table(summary, "summary_for_prism_dkoA.txt", sep="\t", append = FALSE, row.names = TRUE, col.names = TRUE)

conservation <- ""
i <- 1
for (i in 1:max(seedless27a_species_aggregate$Counts)) {
  top_conserved <- seedless27a_species_aggregate[which(seedless27a_species_aggregate$Counts>=i),]
  filteredA_con <- seedless27a[ (seedless27a$SYMBOL %in% top_conserved$SYMBOL), ]
  datLog2FC <-  as.numeric(filteredA_con$FC_used) #canonical FC without targets
  allseedlessAcon = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(filteredA_con))) 
  avg_allseedlessAcon <- format(round(allseedlessAcon[c(101)], 5), nsmall = 5)
  test <- cbind(i, nrow(filteredA_con), avg_allseedlessAcon)
  conservation <- rbind(conservation,test)
}
write.table(conservation, "summary_conservation_A.txt", sep="\t", append = FALSE, row.names = TRUE, col.names = TRUE)
sites_effect <- ""
cons=1
sites=1

for (cons in 1:11) {
  for (sites in 1:4) {
    top_conserved <- seedless27a_species_aggregate[which(seedless27a_species_aggregate$num_sites>=sites & seedless27a_species_aggregate$Counts>=cons),]
    filteredA_con <- seedless27a[ (seedless27a$SYMBOL %in% top_conserved$SYMBOL), ]
    datLog2FC <-  as.numeric(filteredA_con$FC_used) #canonical FC without targets
    allseedlessAcon = c(0, cumsum(table(cut(datLog2FC, breaks, right=FALSE))/nrow(filteredA_con))) 
    avg_allseedlessAcon <- format(round(allseedlessAcon[c(101)], 5), nsmall = 5)
    test <- cbind(cons, sites, nrow(filteredA_con), avg_allseedlessAcon)
    sites_effect <- rbind(sites_effect,test)
  }
}
write.table(sites_effect, "summary_sites_effect_A.txt", sep="\t", append = FALSE, row.names = TRUE, col.names = TRUE)
