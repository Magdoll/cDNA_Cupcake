#!/usr/bin/bash

DIR=$1
for file in $DIR/*/Analysis_Results/*.1.bax.h5
do
	l=`expr ${#file} - 9`  # now l is just the prefix <l>.1.bax.h5
	p=${file:0:$l}
	echo "bax2bam $p.1.bax.h5 $p.2.bax.h5 $p.3.bax.h5"
done
