#!/bin/bash

#FILES=/scratch/xc406/hg19_wgEncodeMapability/wgEncodeCrgMapabilityAlign*mer.bigWig
chroms=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

        for file in "$@"
        do
                echo $file
                pathfname="${file%.*}"
                for chrom in ${chroms[@]}
                do
                        echo $chrom
                        bigWigToBedGraph -chrom=${chrom} ${file} ${pathfname}${chrom}.bedGraph
			bzip2 ${pathfname}${chrom}.bedGraph
                done
        done