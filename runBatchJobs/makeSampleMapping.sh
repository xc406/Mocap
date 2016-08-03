##directory to put pbs scripts in
pbsdir=/data/cgsb/bonneau/xchen/mocap/pbsdir

#Number of samples
for sample in {1..98}
do
        cd $pbsdir
        echo "#PBS -V
#PBS -r n
#PBS -m ae
#PBS -M xc406@nyu.edu
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=48,walltime=48:00:00,mem=1000gb
#PBS -N smbisquare${sample}

module load r/intel/3.2.2

        Rscript /home/xc406/code/mocap/runBatchJobs/sample_mapping.R --xdir \"/data/cgsb/bonneau/xchen/motifs/xpv3/\" \\
								--ydir \"/data/cgsb/bonneau/xchen/motifs/xpv4/\" \\
								--gamma 1.3 \\
								--ind.test $sample \\
								--cores 48 \\
								--exp DNase-Seq \\
								--verbose TRUE

" > smbisquare${sample}.pbs

done
