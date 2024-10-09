#!/bin/bash

set +x

vcfutils=/usr/share/samtools/vcfutils.pl


function abspath() {
    # generate absolute path from relative path
    # $1     : relative filename
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    else
        #echo "/dev/null"
        echo " ---FILE-DOES-NOT-EXIST---"
    fi
}

function relpath() {
    # both $1 and $2 are absolute paths beginning with /
    # returns relative path to $2($target) from $1($source)
    if [ $# -eq 2 ] ; then
       source="$1"
       target="$2"
    elif [ $# -eq 1 ] ; then
       source=`pwd`
       target="$1"
    else
       echo "---relpath error, usage: relpath [source] target---"
       return
    fi
    
    # example use (ensure a symlink contains the relative path to 
    # a file from our current position instead of its absolute path):
    #    f=`abspath $file`
    #    ln -sf `relpath $f` target.lnk
     
    if [ ! -e "$source" -o ! -e "$target" ] ; then
        echo "---NO-SUCH-FILE---"
        return
    fi

    # example:
    #    f=`abspath $file`
    #    ln -s `relpath `pwd` $f` f.lnk

    common_part=$source # for now
    result="" # for now

    while [[ "${target#$common_part}" == "${target}" ]]; do
        # no match, means that candidate common part is not correct
        # go up one level (reduce common part)
        common_part="$(dirname $common_part)"
        # and record that we went back, with correct / handling
        if [[ -z $result ]]; then
            result=".."
        else
            result="../$result"
        fi
    done

    if [[ $common_part == "/" ]]; then
        # special case for root (no common path)
        result="$result/"
    fi

    # since we now have identified the common part,
    # compute the non-common part
    forward_part="${target#$common_part}"

    # and now stick all parts together
    if [[ -n $result ]] && [[ -n $forward_part ]]; then
        result="$result$forward_part"
    elif [[ -n $forward_part ]]; then
        # extra slash removal
        result="${forward_part:1}"
    fi

    echo $result
}

function posix_relpath () {
    [ $# -ge 1 ] && [ $# -le 2 ] || return 1
    current="${2:+"$1"}"
    target="${2:-"$1"}"
    [ "$target" != . ] || target=/
    target="/${target##/}"
    [ "$current" != . ] || current=/
    current="${current:="/"}"
    current="/${current##/}"
    appendix="${target##/}"
    relative=''
    while appendix="${target#"$current"/}"
        [ "$current" != '/' ] && [ "$appendix" = "$target" ]; do
        if [ "$current" = "$appendix" ]; then
            relative="${relative:-.}"
            echo "${relative#/}"
            return 0
        fi
        current="${current%/*}"
        relative="$relative${relative:+/}.."
    done
    relative="$relative${relative:+${appendix:+/}}${appendix#/}"
    echo "$relative"
}

function analyze_sam_file() {
    # IMPORTANT: FILE MUST BE IN CURRENT DIRECTORY!!!
    #	Usage: analyze_sam_file file.sam reference.fna
    
    # extract filename
    dir=`dirname "$1"`
    base=`basename "$1" .sam`
    cd "$dir"

    ref=`abspath "$2"`

    if [ ! -e "$1" -o ! -e "$2" ] ; then 
    	echo "analyze_sam_file: '$1' and '$2' must exist"
	return
    fi

    echo "Checking SAM file"
    head -n 1 $base.sam | grep '@SQ'
    if [ $? -eq 1 ] ; then
        echo "Sanitizing sam file (old file named $base.sam.orig)"
        # insert sequence information from the reference sequence
        mv $base.sam $base.sam.orig
        sq=`head -1 $ref | cut -d' ' -f1 | tr -d '>'`
	ln=`tail -n+2 $ref | tr -d '\n' | wc -c`
	echo "@SQ	SN:$sq	LN:$ln" | cat - $base.sam.orig > $base.sam
    fi

    echo "Analyzing SAM file"

    # Postprocess to BAM
    if [ ! -e ${base}.sorted.bam.bai ] ; then
	echo "Creating ${base}.sorted.bam.bai"
        samtools view -bS ${base}.sam -o ${base}.bam
        samtools sort ${base}.bam ${base}.sorted
        samtools index ${base}.sorted.bam
    fi

    # Do variant calling
    # variant calling (pre-samtools 1.x)
    #samtools faidx $ref
    #samtools pileup -f $ref ${base}.sorted.bam | \
    #	bcftools view -bvcg - > ${base}.bcf
    #bcftools view ${base}.bcf > ${base}.vcf
    if [ ! -e ${base}.variants.vcf ] ; then
	echo "Doing variant calling"
        samtools mpileup -d 10000 -g -f $ref ${base}.sorted.bam > \
 	    ${base}.variants.bcf
        bcftools view -c -v ${base}.variants.bcf > ${base}.variants.vcf
    fi

    # Do consensus
    # consensus (very old samtools) ***
    if [ ! -e ${base}.cons.fastq ] ; then
	echo "Calculating consensus"
        # old samtools
	#~/src/samtools/samtools-527/samtools pileup -cf $ref ${base}.sorted.bam > \
	#    ${base}.cons.pileup
        #perl /usr/share/samtools/samtools.pl pileup2fq ${base}.cons.pileup > \
	#    ${base}.cons.fastq

	# consensus (samtools 1.x) DOES NOW WORK
	samtools mpileup -u -f $ref ${base}.sorted.bam | \
		bcftools view -gcv - | $vcfutils vcf2fq > ${base}.cons.fastq
	# But has also changes quite a bit over times!
	#	bcftools consensus -f $ref - | vcfutils.pl vcf2fq > ${base}.cons.fastq
	
	# and yet another try
	#samtools mpileup -d 10000 -v -f $ref ${base}.sorted.bam > E.vcf 
	#bcftools consensus -f $ref E.vcf -o E.c.vcf
	#cat E.c.vcf | vcfutils.pl vcf2fq > E.fq
    fi

    # Compute coverage
    if [ ! -e ${base}.cov.stats ] ; then
	echo "Computing coverage"
        samtools mpileup -d 10000 -ABR ${base}.sorted.bam > ${base}.coverage
        # requires processing to compute averages, mean, etc...
        cut -f4 ${base}.coverage | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
              -e 'cat(min(d), max(d), median(d), mean(d), "\n", sep="\t")' \
              -e 'summary(d)' > ${base}.cov.stats
    fi

    ###
    if [ ! -e ${base}.sorted.fastq ] ; then
	echo "Piling up"
        samtools mpileup -u -f $ref ${base}.sorted.bam > ${base}.sorted.bcf
        bcftools view -v -c -g ${base}.sorted.bcf > ${base}.sorted.vcf
        perl /usr/share/samtools/vcfutils.pl vcf2fq ${base}.sorted.vcf > ${base}.sorted.fastq
    fi
}

#
# map reads against a reference.
#	We need the reference in "reference.fna" and the reads (direct
#	and indirect) as "R1.fastq" and "R2.fastq"
#

if [ ! -e ref.prj ] ; then
  echo "Making reference database"
  lastdb ref reference.fna
fi

echo "Mapping individual reads"
echo "R1"
if [ ! -s R1.maf ] ; then
    lastal -Q1 -i1 ref R1.fastq > R1.maf
fi
if [ ! -s R1.split.maf ] ; then
    last-split R1.maf > R1.split.maf
fi

echo "R2"
if [ ! -s R2.maf ] ; then
    lastal -Q1 -i1 ref R2.fastq > R2.maf
fi
if [ ! -s R2.split.maf ] ; then
    last-split R2.maf > R2.split.maf
fi

echo "Mapping paired reads"
if [ ! -e mapped.maf ] ; then
    last-pair-probs R1.maf R2.maf > mapped.maf 2> last-pair-stats
fi

#last-pair-probs R1.split.maf R2.split.maf > split.mapped.maf 2> last-split-stats

if [ ! -e mapped.1.maf ] ; then
    head -n 18 R1.maf | cat - mapped.maf > mapped.H.maf
    last-split mapped.H.maf > mapped.1.maf 
fi


echo "Making plots"
for i in *.maf ; do
    if [ ! -e `basename $i maf`png ] ; then
        echo $i
        last-dotplot $i `basename $i maf`png
    fi
done

echo "Making SAM files"
#maf-convert sam R1.maf > R1.sam
#maf-convert sam R2.maf > R2.sam
if [ ! -e mapped.sam ] ; then
    maf-convert sam mapped.maf > mapped.sam
fi
if [ ! -e mapped.1.sam ] ; then
    maf-convert sam mapped.1.maf > mapped.1.sam
fi

echo "Analyzing SAM files"
if [ ! -e mapped.cons.fastq ] ; then
    analyze_sam_file mapped.sam reference.fna
fi
if [ ! -e mapped.1.cons.fastq ] ; then
    analyze_sam_file mapped.1.sam reference.fna
fi
