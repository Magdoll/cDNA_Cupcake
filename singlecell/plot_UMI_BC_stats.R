# -----------------------
# UMI-BC report
# -----------------------
library(ggplot2)
library(dplyr)
library(ggseqlogo)
library(grid)
library(gridExtra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

input.csv <- args[1]
report.file <- paste(input.csv, '_report.pdf', sep='')


x <- read.table(input.csv, sep=',',header=T)

pdf(file=report.file, width=6.5, height=6.5)
grid.newpage()
cover <- textGrob("UMI-BC report", gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)

fields <- c("# of loci", "# of UMI", "# of BC", "# of FLNC", "# of distinct FLNC (locus,UMI,BC)");
n_loci <- length(unique(x$index))
n_umi <- length(unique(x$UMI))
n_bc <- length(unique(x$BC))
x$n_reads <- str_count(x$members, ",")+1
n_reads <- sum(x$n_reads)
n_molecules <- dim(x)[1]

d <- matrix(c(n_loci, n_umi, n_bc, n_reads, n_molecules), ncol=1)
rownames(d) <- fields

# draw: table of UMI/BC stats
grid.newpage()
grid.draw(tableGrob(as.table(d), cols=NULL))

# draw: histogram of duplication rate
p.umi_hist <- ggplot(x, aes(x=log2(n_reads))) + 
  geom_histogram(binwidth=1) + 
  xlab("log2(# of duplicates)") + 
  ylab("Count") +
  labs(title="Histogram of PCR duplciation rate")
print(p.umi_hist)


# draw: histogram of molecules covering each locus
p.loci_hist <- ggplot(x, aes(x=index)) + 
  geom_histogram(binwidth=1) + 
    xlab("Unique locus (sorted in chromosomal order)") +
	ylab("Count of unique molecules") + 
	labs(title="Histogram of loci coverage by unique molecules")

t.loci <- table(x$locus)
t.loci <- sort(t.loci, decreasing=T)
table.loci_top10 <- tableGrob(t.loci[1:10], rows=rownames(t.loci)[1:10], cols=c("count"))
grid.arrange(p.loci_hist, table.loci_top10, ncol=1)

# draw: UMI and BC motifs
p.logo.umi <- ggseqlogo(substring(as.character(x$UMI), 1, 8), method='prob')
p.logo.umi.text <- textGrob("Frequency of UMI base at each position", 
							gp=gpar(fontface="italic", fontsize=12), 
							vjust=0)
if (any(!is.na(x$BC))) {
  p.logo.bc  <- ggseqlogo(as.character(x$BC), method='prob')
  p.logo.bc.text <- textGrob("Frequency of BC base at each position", 
                             gp=gpar(fontface="italic", fontsize=12), 
                             vjust=0)
  grid.arrange(p.logo.umi.text, p.logo.umi, 
               p.logo.bc.text, p.logo.bc,
               ncol=1)
} else
{
  grid.arrange(p.logo.umi.text, p.logo.umi,ncol=1);
}


dev.off()

