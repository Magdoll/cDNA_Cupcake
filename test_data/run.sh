dataset create --type SubreadSet test.subreadset.xml test.subreads.bam

# -------------------
# run CCS and classify
# -------------------
ccs --numThreads=20 --noPolish --minLength=50 --maxLength=15000 --minPasses=1 --minPredictedAccuracy=0.7 --minZScore=-999 --maxDropFraction=0.8 --minPredictedAccuracy=0.8 --minSnr=3.75 test.subreads.bam test.ccs.bam
bamtools convert -format fasta -in test.ccs.bam > ccs.fasta
pbtranscript classify --flnc isoseq_flnc.fasta --nfl isoseq_nfl.fasta -d output --cpus 2 --min_seq_len 100 ccs.fasta isoseq_draft.fasta

# -------------------
# running ToFU2
# -------------------
run_preCluster.py --cpus=2 --num_seqs_per_batch=10000
generate_batch_cmd_for_preCluster_out.py preCluster.cluster_info.csv preCluster_out/ --cpus=4 --cmd_filename cmds
bash cmds
ls -1d preCluster_out/* > preCluster_out_dirs.list.txt
collect_IceIterative2_result.py --num_chunks 1 preCluster_out_dirs.list.txt collected_final
generate_batch_cmd_for_polishing.py --cpus=12 --cmd_filename=cmds2 collected_final.chunk isoseq_nfl.fasta test.subreadset.xml

### instead of going through qsub (which are the commands in `cmd2`), just directly execute locally
bash collected_final.chunk0/collected_final.chunk0.sh
