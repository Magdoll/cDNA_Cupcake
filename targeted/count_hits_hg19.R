x<-read.table('isoseq_flnc.fasta.hg19.sam.probe_hit.txt',sep='\t',header=T)
print("# hit align hg19: ")
print(length(unique(x$read_id)))

good <- subset(x, num_probe>0) # reads that hit one or more probes

print("# hit hg19 probe: ")
print(length(unique(good$read_id)))
print("# of hit hg19 genes: ")
print(length(unique(good$genes)))
