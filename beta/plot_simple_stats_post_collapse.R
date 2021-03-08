# --------------------------------------
library(ggplot2)
library(dplyr)
library(ggthemes)
# --------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input.file1 <- args[1]
input.file2 <- args[2]
PREFIX <- args[3]

x.1 <- read.table(input.file1, sep='\t',header=T)
x.2 <- read.table(input.file2,sep='\t',header=T)

x <- x.1
# plot length histogram
ggplot(x, aes(length)) + geom_histogram(binwidth=200, fill='lightblue', color='black') + xlab("Transcript Length (bp)") + ylab("Number of Unique Transcripts") + labs(title="Mapped Unique Transcript Lengths")  + theme_tufte() + xlim(c(0,10000))
ggsave(paste(PREFIX,".Rplot.histogram.png",sep=''), width=6, height=4.5)

# plot genomic length histogram
ggplot(x, aes(genomic_length/1000)) + geom_histogram(binwidth=100, fill='lightgreen', color='black') + xlab("Transcript Genomic Length (kb)") + ylab("Number of Unique Transcripts") + labs(title="Mapped Unique Transcript Genomic Lengths")  + theme_tufte()
ggsave(paste(PREFIX,".Rplot.genomic_len_histogram.png",sep=''), width=6, height=4.5)

# plot number of exons
x$exon_cat <- "1"
x[x$num_exon>=2,"exon_cat"] <- "2-5"
x[x$num_exon>=10,"exon_cat"] <- "10-20"
x[x$num_exon>20,"exon_cat"] <- ">20"
x$exon_cat <- factor(x$exon_cat, levels=c("1","2-5","10-20",">20"))
ggplot(x, aes(exon_cat)) + geom_bar(fill='darkgreen') + xlab("Number of Exons") + ylab("Number of Transcripts") + theme_tufte() + labs(title="Number of Exons Per Transcript")
ggsave(paste(PREFIX,".Rplot.num_exons.png",sep=''), width=6, height=4.5)

# number oof isoforms per locus
t <- x %>% group_by(locus) %>% summarise(count=n())
t$isocat <- "1"
t[t$count>=2, "isocat"] <- "2-5"
t[t$count>=6, "isocat"] <- "6-10"
t[t$count>10, "isocat"] <- "11-20"
t[t$count>20, "isocat"] <- ">20"
t$isocat <- factor(t$isocat, levels=c("1","2-5","6-10","11-20",">20"))
ggplot(t, aes(isocat)) + geom_bar(fill='pink') + xlab("Number of Isoforms") + ylab("Number of Locus") + theme_tufte() + labs(title="Number of Isoforms Per Locus")
ggsave(paste(PREFIX,".Rplot.num_isoforms.png",sep=''), width=6, height=4.5)

# ---------------
# read exon/intron size
# ---------------
x <- x.2
x$exon_size_cat <- "<100bp"
x[x$exon_size>=100,"exon_size_cat"] <- "100-200"
x[x$exon_size>200,"exon_size_cat"] <- "201-500"
x[x$exon_size>500,"exon_size_cat"] <- "501-1000"
x[x$exon_size>1000,"exon_size_cat"] <- "1001-10000"
x[x$exon_size>10000,"exon_size_cat"] <- ">10000"
x$exon_size_cat <- factor(x$exon_size_cat, levels=c("<100bp","100-200","201-500","501-1000","1001-10000",">10000"))
ggplot(x, aes(exon_size_cat)) + geom_bar() + theme_tufte() + xlab("Exon Lengths (bp)") + ylab("Count") + labs(title="Distribution of Exon Lengths")
ggsave(paste(PREFIX,".Rplot.num_exon_sizes.png",sep=''), width=6, height=4.5)

x1 <- subset(x, !is.na(intron_size))
x1$intron_size_cat <- "<100bp"
x1[x1$intron_size>100, "intron_size_cat"] <- "100-1000"
x1[x1$intron_size>1000, "intron_size_cat"] <- "1000-10000"
x1[x1$intron_size>10000, "intron_size_cat"] <- "10kb-100kb"
x1[x1$intron_size>100000, "intron_size_cat"] <- "100kb-1Mb"
x1[x1$intron_size>1000000, "intron_size_cat"] <- ">1Mb"
x1$intron_size_cat <- factor(x1$intron_size_cat, levels=c("<100bp","100-1000","1000-10000", "10kb-100kb", "100kb-1Mb", ">1Mb")
)
ggplot(x1, aes(intron_size_cat)) + geom_bar() + theme_tufte() + xlab("Intron Lengths (bp)") + ylab("Count") + labs(title="Distribution of Intron Lengths")
ggsave(paste(PREFIX,".Rplot.num_intron_sizes.png",sep=''), width=6, height=4.5)

