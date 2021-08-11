import os, sys

f = open('dedup.5merge.collapsed.read_stat.txt', 'w')
f.write("id\tpbid\tis_fl\tstat\n")
for line in open('dedup.5merge.collapsed.group.txt'):
    a,b=line.strip().split()
    for x in set(b.split(',')): f.write(f"{x}\t{a}\tY\tunique\n")
    
f.close()
