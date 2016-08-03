tfctv=( ARID3AHepg2 ARID3AK562 ATF3Hepg2 ATF3K562 BHLHE40A549 \
BHLHE40Hepg2 BHLHE40K562 CEBPBA549 CEBPBHepg2 CEBPBK562 \
CEBPDHepg2 CEBPDK562 CHD2Hepg2 CHD2K562 CREB1Hepg2 \
CREB1K562 CTCFA549 CTCFHepg2 CTCFK562 E2F4K562 E2F6A549 \
E2F6K562 ELF1Hepg2 ELF1K562 EP300Hepg2 EP300K562 ETS1K562 \
FOSL2Hepg2 FOXA1Hepg2 FOXA2Hepg2 GABPAHepg2 GABPAK562 \
HDAC2Hepg2 HDAC2K562 HEY1Hepg2 HEY1K562 JUNDHepg2 JUNDK562 \
JUNHepg2 JUNK562 MAFFHepg2 MAFFK562 MAFKHepg2 MAFKK562 MAXA549 \
MAXHepg2 MAXK562 MAZHepg2 MAZK562 MXI1Hepg2 MXI1K562 MYCA549 \
MYCHepg2 MYCK562 NFYAK562 NR2C2Hepg2 NR2F2Hepg2 NR2F2K562 \
NRF1Hepg2 NRF1K562 RAD21A549 RAD21Hepg2 RAD21K562 RESTK562 \
RFX5Hepg2 RFX5K562 SIN3AHepg2 SIN3AK562 SIX5K562 SMC3Hepg2 \
SMC3K562 SP1Hepg2 SP1K562 SP2Hepg2 SP2K562 SRFHepg2 SRFK562 \
TAF1Hepg2 TAF1K562 TBPHepg2 TBPK562 TCF12Hepg2 TEAD4A549 \
TEAD4Hepg2 TEAD4K562 USF1Hepg2 USF1K562 USF2Hepg2 USF2K562 \
YY1Hepg2 YY1K562 ZBTB33Hepg2 ZBTB33K562 ZBTB7AHepg2 ZBTB7AK562 \
ZNF143K562 ZNF274K562 ZNF384K562)

##directory to put X Y files in
xdir=/data/cgsb/bonneau/xchen/mocha/xdir
ydir=/data/cgsb/bonneau/xchen/mocha/ydir
##directory to put pbs scripts in
pbsdir=/data/cgsb/bonneau/xchen/mocha/pbsdir

for i in "$@"
do
        for j in ${tfctv[@]}
        do
		cd $pbsdir
		echo "#PBS -V
#PBS -r n
#PBS -m a
#PBS -M xc406@nyu.edu
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1,walltime=4:00:00,mem=60gb
#PBS -N slrx${i}${j}

module load r/intel/3.2.2

        Rscript genXYmat.R --tf_name_cell_type ${i} --pred_tf_name_cell_type ${j} --xdir ${xdir} --ydir ${ydir}

" > slrx${i}${j}.pbs

	done
done
