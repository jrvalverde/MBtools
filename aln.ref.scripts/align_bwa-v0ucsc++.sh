#!/bin/bash

#et -x
 
me=`basename "${BASH_SOURCE[0]}"`
mydir=`dirname "${BASH_SOURCE[0]}"`
IN_FQ=$1
OUT_PFX=$2
ASSEMBLY=$3
PAIRED=$4

#REFBASE="$WORK_COMMON/ref_genome/$BWA_DIR/base"
REFBASE=~/data/
# name of the installation folders for mosdepth and sambamba
mosdepth=~/contrib/mosdepth/mosdepth-master/
sambamba=~/contrib/sambamba/

# number of simultaneous processors to use in speeding up calculations
nproc=`nproc`
nproc=$(( nproc / 2 ))
export nproc

# show usage if we don't have 4 command-line arguments
#if [ $# -ne 4 ] ; then
if [ "$PAIRED" == "" ]; then
    echo "-----------------------------------------------------------------";
    echo "Align fastq data with bwa, producing a sorted indexed BAM file.";
    echo " ";
    echo "align_bwa.sh in_file out_pfx assembly paired(0|1|2)";
    echo " ";
    echo "  in_file   For single-end alignments, path of the input fastq file.";
    echo "            For paired-end alignemtts, path to the the R1 fastq file"
    echo "            which must contain the string '_R1.' in its name. The";
    echo "            corresponding 'R2' must have the same path except for '_R1'";
    echo "  out_pfx   Desired prefix of output files.";
    echo "  assembly  One of: hg19 hg18 mm10 mm9 sacCer3 sacCer3 ecoli or your own";
    echo "            fasta file.";
    echo "  paired    0 = single end alignment; 1|2 = paired end.";
    echo "            (Use 1 or 2 to indicate which read file has been used)"; 
    echo " ";
    echo "Examples:";
    echo "  align_bwa.sh my.fastq mrna_b1_ln1 hg18 0";
    echo "  align_bwa.sh my_L001_R1.fastq swi6_b2_ln1 sacCer3 1";
    echo "  align_bwa.sh my_L001_R2.fastq swi6_b2_ln1 sacCer3 2";
    exit 1;
fi

# redirect all output to log files (one for stdout and one for stderr)
#log=./$name/log.$myname.outerr
#$appendlog='-a'
#$appendlog=''
#exec >& >(stdbuf -o 0 -e 0 tee $appendlog "$log")
stdout=./log.$me.out
stderr=./log.$me.err
exec > $stdout 2> $stderr


# general function that exits after printing its text argument
#   in a standard format which can be easily grep'd.
err() {
  echo "$1...exiting";
  exit 1; # any non-0 exit code signals an error
}
# function to check return code of programs.
# exits with standard message if code is non-zero;
# otherwise displays completiong message and date.
#   arg 1 is the return code (usually $?)
#   arg2 is text describing what ran
ckRes() {
  if [ "$1" == "0" ]; then
    echo "..Done $2 `date`";
  else
    err "$2 returned non-0 exit code $1";
  fi
}
# function that checks if a file exists
#   arg 1 is the file name
#   arg2 is text describing the file (optional)
ckFile() {
  if [ ! -e "$1" ]; then
    err "$2 File '$1' not found";
  fi
}
# function that checks if a file exists and
#   that it has non-0 length. needed because
#   programs don't always return non-0 return
#   codes, and worse, they also create their
#   output file with 0 length so that just
#   checking for its existence is not enough
#   to ensure the program ran properly
#   arg 1 is the file name
#   arg2 is text describing the file (optional)
ckFileSz() {
  if [ ! -s "$1" ] ; then
      err "$2 file '$1' doesn't exist or is zero length";
  fi
  return
  ckFile $1 $2;  
  SZ=`ls -l $1 | awk '{print $5}'`;
  if [ "$SZ" == "0" ]; then
    err "$2 file '$1' is zero length";
  fi
}
 
# --------------
# Defaulting
# --------------
# If aligning read pairs, find the name of the R2 file
#   based on the R1 filename.
if [ "$PAIRED" == "1" ]; then
    IN_FQ_R1="$IN_FQ";
    IN_FQ_R2=${IN_FQ_R1/_R1/_R2};	# substitute _R1 by _R2
fi
# samewise find R1 filename if we got the R2 file instead
if [ "$PAIRED" == "2" ]; then
    IN_FQ_R2="$IN_FQ";
    IN_FQ_R1=${IN_FQ_R2/_R2/_R1};	# substitute _R1 by _R2
fi

# Find the path of the BWA reference based on the assembly
#   name provided. Assumes a common directory structure
# At TACC the structure is rooted a common BioITeam directory.
#   Set WORK_COMMON appropriately outside this script if
#   not running at TACC.
#: ${WORK_COMMON:=/corral-repl/utexas/BioITeam/}
# At CNB we keep it in ~/data
: ${WORK_COMMON:=~/data}
# Note: reference indexes created with bwa 0.5.9 are *not*
#   compatible with bwa 0.6 or higher. So we go to some
#   trouble to figure out which BWA version we have
BWA_VER=`bwa 2>&1 | grep Version | awk '{print $2}' | awk -F '.' '{print $1 "." $2}'`;
# BWA_DIR is not further used.
#if [ "$BWA_VER" == "0.5" ]; then BWA_DIR=bwa;
#elif [ "$BWA_VER" == "0.6" ]; then BWA_DIR=bwa6;
#elif [ "$BWA_VER" == "0.7" ]; then BWA_DIR=bwa7;
#else BWA_DIR=unknown_bwa_version; fi

# select the assembly
if [ "$ASSEMBLY" == "hg18" ]; then
    REF_PFX="$REFBASE/genomes/$ASSEMBLY/Homo_sapiens_assembly18.fasta";
elif [ "$ASSEMBLY" == "hg19" ]; then
    REF_PFX="$REFBASE/genomes/$ASSEMBLY/ucsc.hg19.fasta";
elif [ "$ASSEMBLY" == "ecoli" ]; then
    REF_PFX="$REFBASE/genomes/$ASSEMBLY/REL606.5.fasta";
else
    REF_PFX="$REFBASE/$ASSEMBLY/${ASSEMBLY}.fa";
    if [ ! -e $REF_PFX ] ; then
        REF_PFX=`realpath "${ASSEMBLY}"`
    fi
fi
 
# Read group information is part of the SAM/BAM header that desribes
#   what is being aligned. When multiple lanes of data are combined
#   from separate BAM files, read groups provide identification of the
#   source of each read. Variant callers such as GATK depend on
#   having well-defined read group information.
# Here we set the RG variable to be the read group line we want
#   inserted in the header.
READ_GRP=$OUT_PFX;
RG='@RG\tID:1\tPL:ILLUMINA\tSM:'$READ_GRP'\tDS:ref='$ASSEMBLY',pfx='$REF_PFX

# We can do better using more info:
#id=1
#bc=$BARCODE_SEQ
#cn=$SEQUENCING_CENTER_NAME
#dt=$DATE_RUN
#sm=$SAMPLE"
#pl="$PLATFORM"		#CAPILLARY,DNBSEQ,ELEMENT,HELICOS,ILLUMINA,IONTORRENT,LS454,ONT,PACBIO,SOLID,ULTIMA
#pm="$PLATFORM_MODEL"
#pu="$PLATFORM_UNIT"	#$FLOCELL_BARCODE.$LANE.$SAMPLE_BARCODE
#lb="$LIBRARY_ID"
#ds="$DESCRIPTION"
#RG="@RG\tID:$rgid\tPU=$pu\tSM=$sm\tPL=$pl\tLB=$ln\tDS:ref=$ASSEMBLY,pfx=$REF_PFX""

 
# Display how the program will be run, including
#   defaulted arguments. Do this before running
#   checks so user can see what went wrong.
echo "=================================================================";
echo "align_bwa.sh - `date`";
if [ "$PAIRED" == "1" ]; then
echo "  fastq read1 file:  $IN_FQ_R1";
echo "  fastq read2 file:  $IN_FQ_R2";
else
echo "  input file:        $IN_FQ";
fi
echo "  output prefix:     $OUT_PFX";
echo "  assembly:          $ASSEMBLY";
echo "  bwa version:       $BWA_VER";
echo "  ref prefix:        $REF_PFX";
echo "  read group line:   $RG";
echo "---------------------------------------------------------";
 
# ------------------
# Error Checks
# ------------------
# Make sure the fastq file(s) exist.
# For paired end data, also make sure we have
#   two different files.
if [ "$PAIRED" == "1" ]; then
    ckFile "$IN_FQ_R1" "Fastq read1";
    ckFile "$IN_FQ_R2" "Fastq read2";
    if [ "$IN_FQ_R1" == "$IN_FQ_R2" ]; then
        err "Fastq read1 and read2 files are the same: '$IN_FQ_R2'";
    fi
    IN_FQ_R1=`realpath "$IN_FQ_R1"`
    IN_FQ_R2=`realpath "$IN_FQ_R2"`
else
    ckFile "$IN_FQ" "Input fastq";
    IN_FQ=`realpath "$IN_FQ"`
fi
 
# Make sure we have found an appropriate reference
#   by checking that one of the standard files exists.
# If it doesn't create it
if [ ! -s "${REF_PFX}.amb" ] ; then
    #printenv | sort
    echo "indexing $REF_PFX"
    bwa index -a is "${REF_PFX}"
else
    echo "using existing ${REF_PFX}.amb"
fi
ckFile "${REF_PFX}.amb" "$ASSEMBLY Reference";
 
# Make sure version information for our two programs is
#   part of our execution record. This is done by
#   calling the programs with no arguments
echo "---------------------------------------------------------";
echo "Program version information";
echo "---------------------------------------------------------";
bwa |& head -n 5
samtools |& head -n 4
 
# ------------------
# The actual work!
# ------------------
 
if [ "$PAIRED" == "0" ]; then
    echo "---------------------------------------------------------";
    echo "Running bwa aln (single end reads)";
    echo "---------------------------------------------------------";
    if [ ! -s "OUT_PFX.sai" ] ; then
        bwa aln $REF_PFX $IN_FQ > $OUT_PFX.sai
        ckRes $? "bwa aln";
        ckFileSz "$OUT_PFX.sai";
    fi
    
    echo "---------------------------------------------------------";
    echo "Running bwa samse";
    echo "---------------------------------------------------------";
    if [ ! -s "$OUT_PFX.bam" ] ; then
        bwa samse -r "$RG" $REF_PFX $OUT_PFX.sai $IN_FQ \
            | samtools view -b -S - \
            > $OUT_PFX.bam;
        ckRes $? "bwa samse";
        ckFileSz "$OUT_PFX.bam";
    fi
else
    echo "---------------------------------------------------------";
    echo "Running bwa aln on read1, read2 paired ends";
    echo "`date`";
    echo "---------------------------------------------------------";
     
    echo "Aligning '$IN_FQ_R1'...";
    if [ ! -s "$OUT_PFX.R1.sai" ] ; then
        bwa aln -t $nproc $REF_PFX $IN_FQ_R1 \
            > $OUT_PFX.R1.sai;
        ckRes $? "bwa aln read1";
        ckFileSz "$OUT_PFX.R1.sai";
    fi
    
    echo "Aligning '$IN_FQ_R2'...";
    if [ ! -s "$OUT_PFX.R2.sai" ] ; then
        bwa aln -t $nproc $REF_PFX $IN_FQ_R2 \
            > $OUT_PFX.R2.sai;
        ckRes $? "bwa aln read2";
        ckFileSz "$OUT_PFX.R2.sai";
    fi
     
    echo "---------------------------------------------------------";
    echo "Running bwa sampe";
    echo "---------------------------------------------------------";
    if [ ! -s "$OUT_PFX.bam" ] ; then
        bwa sampe -r "$RG" $REF_PFX $OUT_PFX.R1.sai $OUT_PFX.R2.sai \
            $IN_FQ_R1 $IN_FQ_R2 \
            | samtools view -@ $nproc -b -S - \
            > $OUT_PFX.bam;
        ckRes $? "bwa sampe";
        ckFileSz "$OUT_PFX.bam";
    fi
fi
 
echo "---------------------------------------------------------";
echo "Creating sorted, indexed bam file";
echo "---------------------------------------------------------";
if [ ! -s "$OUT_PFX.sorted.bam" ] ; then
    echo "Sorting '$OUT_PFX.bam'...";
    #samtools sort -@ $nproc $OUT_PFX.bam -O bam > $OUT_PFX.sorted.bam ;
    #samtools sort -@ $nproc $OUT_PFX.bam -f $OUT_PFX.sorted.bam
    samtools sort -@ $nproc $OUT_PFX.bam -O bam -o $OUT_PFX.sorted.bam
    ckRes $? "samtools sort";
    ckFileSz "$OUT_PFX.sorted.bam";
fi

if [ ! -s  "$OUT_PFX.sorted.bam.bai" ] ; then
    echo "Indexing '$OUT_PFX.sorted.bam'...";
    samtools index $OUT_PFX.sorted.bam;
    ckRes $? "samtools index";
    ckFileSz "$OUT_PFX.sorted.bam.bai";
fi
 
echo "---------------------------------------------------------";
echo "Collecting alignment statistics";
echo "---------------------------------------------------------";
echo "Running flagstat...";
if [ ! -s "$OUT_PFX.flagstat.txt" ] ; then
    samtools flagstat $OUT_PFX.sorted.bam | tee $OUT_PFX.flagstat.txt
    ckRes $? "samtools flagstat";
    ckFileSz "$OUT_PFX.flagstat.txt";
fi
 
 
echo "---------------------------------------------------------";
echo "Colmputing coverage depth";
echo "---------------------------------------------------------";
echo "Running samtools depth...";
if [ ! -e "${OUT_PFX}.sorted.sam_cov.depth" ] ; then
    samtools depth -d 0 ${OUT_PFX}.sorted.bam > ${OUT_PFX}.sorted.sam_cov.depth
    ckFileSz ${OUT_PFX}.sorted.sam_cov.depth
fi

mkdir -p cov.plots
if [ ! -s cov.plots/sam.cov.all.png ] ; then
    echo "$IN_FQ: plotting coverage"
    if [ -s ${OUT_PFX}.sorted.sam_cov.depth ] ; then
	# we have a SAM pre-calculated coverage
	echo "$IN_FQ: plotting all chromosomes together"
	# generate input data
	cut -f 1,2,4 ${OUT_PFX}.sorted.sam_cov.depth \
	| grep '^[0-9]' ${comment# select numbered chromosomes} \
	${comment# feedgnuplot wants 'pos ref depth' not 'ref pos depth'} \
	| awk '{print $2 "\t" $1 "\t" $3 }' ${comment# put columns in order} \
	| feedgnuplot --lines \
	              --hardcopy cov.plots/sam.cov.all.png \
		      --autolegend \
		      --dataid \
		      --domain \
		      --ylabel depth \
		      --xlabel "chr pos" \
		      --title "Chromosome-wide depth" \
		      --exit 
	# this plot may be too overcrowded, it may help to make one
	# separate plot for each reference chromosome
	# first get a list of chromosomes
	cat ${OUT_PFX}.sorted.sam_cov.depth | cut -f1 | sort | uniq \
	| while read chr ; do
	    # only do it if needed
	    if [ -s cov.plots/sam.cov.$chr.png ] ; then continue ; fi
	    echo ">>>$IN_FQ: Plotting coverage for chromosome $chr"
	    # repeat for each chromosome
            cut -f 1,2,4 ${OUT_PFX}.sorted.sam_cov.depth \
	    | grep "^${chr}	" ${comment# select required chromosome} \
	    | awk '{print $2 "\t" $1 "\t" $3 }' ${comment# order columns} \
	    | feedgnuplot --lines \
	          --hardcopy cov.plots/sam.cov.$chr.png \
			  --autolegend \
			  --dataid \
			  --domain \
			  --ylabel depth \
			  --xlabel "chr pos" \
			  --title "Chromosome-wide depth" \
			  --exit 
	done
	#
	# At this point we should consider removing the file
	# ${OUT_PFX}.sorted_sam_cov.depth
	# for it can be very large...
    fi
fi

if [ "$DO_ALL_COVERAGES_PLEASE" == "YES" ] ; then
    echo "Running bedtools genomecov...";
    if [ ! -e "${OUT_PFX}_bed_cov_depth" ] ; then
        bedtools genomecov -d -ibam ${OUT_PFX}.sorted.bam > ${OUT_PFX}.sorted_bed_cov.depth
	plotcov.ch ${OUT_PFX}.sorted_bed_cov.depth all
    fi

    echo "Running mosdepth to calculate coverage..."
    if [ ! -s ${OUT_PFX}.sorted_mosdepth_cov.depth ] ; then
	# to get the same results as samtools and bedtools, we should
	# remove mate-pair overlap detection (effectively double-counting
	# mate-pair overlap regions).
	echo "Computing depth with mosdepth"
	$mosdepth/mosdepth ${OUT_PFX} ${OUT_PFX}.sorted.bam
	zcat ${OUT_PFX}.per-base.bed.gz | cut -f1,3,4 > ${OUT_PFX}_mosdepth_cov.depth
	python3 $mosdepth/scripts/plot-dist.py \
        	${OUT_PFX}.mosdepth.global.dist.txt
	mv dist.html ${OUT_PFX}.mosdepth.cum.cov.plot.html
    fi

    echo "Sunning bb to calculate coverage..."
    if [ ! -s ${OUT_PFX}_sambamba_cov.depth ] ; then
    	# without --fix-mate-overlaps, this gives results similar
    	# to samtools and bedtools
    	# with --fix-mate-overlaps, this gives results similar to
    	# mosdepth
    	echo "Computing depth with sambamba"
    	$sambamba/sambamba depth base \
        	--min-coverage=0 ${comment# same as -c 0 } \
        	--fix-mate-overlaps ${comment# same as -m } \
        	${comment# send to stodout -o ${OUT_PFX}.sorted_sambamba_cov.stats} \
        	${OUT_PFX}.sorted.bam \
        	2> /dev/null	\
    	| while read r p0 d o ; do 
        	echo "$r	$((p0+1))	$d"
    	done > ${OUT_PFX}.sorted_sambamba_cov.depth
    fi
fi

# If we make it here, all went well. Exit with a standard
#   message that can be easily grep'd
echo "---------------------------------------------------------";
echo "All bwa alignment tasks completed successfully!";
echo "`date`";
echo "---------------------------------------------------------";
exit 0;
