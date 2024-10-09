#!/bin/bash
#	make a reference-based genome reconstruction using LAST aligner
#
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
        bcftools call -c -v ${base}.variants.bcf > ${base}.variants.vcf
    fi

    # Do consensus
    # consensus (very old samtools) ***
    if [ ! -e ${base}.cons.fastq ] ; then
	echo "Calculating consensus"
        ~/src/samtools/samtools-527/samtools pileup -cf $ref ${base}.sorted.bam > \
	    ${base}.cons.pileup
        perl /usr/share/samtools/samtools.pl pileup2fq ${base}.cons.pileup > \
	    ${base}.cons.fastq
    fi
    # consensus (samtools 1.x) DOES NOW WORK
    #samtools mpileup -u -f $ref ${base}.sorted.bam | \
    #	bcftools consensus -f $ref - | vcfutils.pl vcf2fq > ${base}.cons.fastq
    #	bcftools view -gc - | vcfutils.pl vcf2fq > ${base}.cons.fastq
    #samtools mpileup -d 10000 -v -f $ref ${base}.sorted.bam > E.vcf 
    #bcftools consensus -f $ref E.vcf -o E.c.vcf
    #cat E.c.vcf | vcfutils.pl vcf2fq > E.fq

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
last_reference_assemble() {
    # assemble usinf a reference and R1 + R2 reads
    ref=`abspath "$1"`
    R1=`abspath "$2"`
    R2=`abspath "$3"`

    if [ ! -d last ] ; then
        echo "Creating last"
        mkdir last
    fi
    echo "Entering last"
    cd last
    if [ ! -e reference.fna ] ; then
        if [ -e "$ref" ] ; then
	    # this is safe because we know reference.fna does not exist
	    #	but it could be a dangling symlink, so we force removal
            #   of any existing symlink before creating the new one.
	    ln -sf `relpath "$ref"` reference.fna
            #ln -sf $ref reference.fna
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a reference sequence for LAST"
	    return
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
	    ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run LAST"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
	    ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run LAST"
            return
        fi
    fi

    # We require two mate-pair reads files in fastq format
    # and a reference file in fasta format in the current
    # directory
    #
    #	Further, they must be called reference.fna, R1.fastq and R2.fastq
    #
    # Check we have everything we need
    if [ ! -e reference.fna -o ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need a reference and two mate pairs files to run LAST"
        return
    fi
    
    ##############
    # DO THE WORK
    ##############
    
    echo "Making reference database"
    lastdb ref reference.fna

    echo "Mapping individual reads"
    if [ ! -s R1.maf ] ; then
	lastal -Q1 -e120 -i1 ref R1.fastq > R1.maf
	last-split R1.maf > R1.split.maf
    fi
    if [ ! -s R2.maf ] ; then
	lastal -Q1 -e120 -i1 ref R2.fastq > R2.maf
	last-split R2.maf > R2.split.maf
    fi

    echo "Mapping paired reads"
    if [ ! -e mapped.maf ] ; then
	last-pair-probs R1.maf R2.maf > paired.maf 2> last-pair-stats
    fi

    #last-pair-probs R1.split.maf R2.split.maf > split.mapped.maf 2> last-split-stats

    echo "Making SAM file"
    #maf-convert R1.maf > R1.sam
    #maf-convert R2.maf > R2.sam
    maf-convert sam paired.maf > paired.sam

    echo "Making plots"
    for i in *.maf ; do
	if [ ! -e `basename $i maf`png ] ; then
            last-dotplot $i `basename $i maf`png
	fi
    done

    cd ..
}


last_reference_assemble reference.fna R1.fastq R2.fastq
cd last
echo "Analyzing SAM file"
#analyze_sam_file R1.sam reference.fna
#analyze_sam_file R2.sam reference.fna
analyze_sam_file paired.sam reference.fna
cd ..
