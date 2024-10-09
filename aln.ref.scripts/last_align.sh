# Assemble reads
#
# We assume that reads are mated pairs provided as two *_R[12]_001.fastq files
#
# We take care not to repeat any work, hence, re-running is safe.
#
# The graphics routines may fail on servers where no X support is available.
# If so, in order to re-run them, we need to remove the generated last, 
# mummer3 and nucmer directories and re-run in a different computer with
# appropriate support.
#
# If we want to genome-align against different reference genomes we can, either
# repeat the full process or only align by renaming reference.fna to something
# else (e.g. reference.fna.orig), renaming the last/mummer3/nucmer directories
# to allow re-running these programs against the new reference and 
# installing/symlinking a new reference. This is confusing, I know, and I
# should certainly separate assembly and comparison into separate programs
# one of these days, eventually.
#
# (C) José R. Valverde, CNB/CSIC. 2015-2016
#
set +x

# Default values
MUMMER="no"			# -m --mummer3
LAST="no"			# -l --last
NUCMER="no"			# -n --nucmer
VERBOSE=0                       # -v --verbose



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

function fastasize() {
    awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' \
        "$1"
}

function fastasizes() {
    # adds total summaries at the end
    awk '/^>/ {if (seqlen){print seqlen};print;seqtotal+=seqlen;seqlen=0;seq+=1;next;}
        {seqlen=seqlen+length($0)}
        END{print seqlen;print seq" sequences, total length " seqtotal+seqlen}' \
        "$1"
}

function fq2f() {
    # Convert fastq file to fasta format
    # we'll simply strip the quality codes, without using them
    # output goes to stdout!!
    #set +x
    fq=`abspath "$1"`
    #echo "#fq2f: Converting $fq to fasta format"
    if [ ! -e "$fq" ] ; then
        echo "fq2f: error, nothing to do"
        return
    fi
    count=0
    deleted=0
    n=0
    # we have a problem when reading lines: by default read interprets 
    # backslashes, therefore, a line ending in \ will include the next 
    # line as well
    # We must use read -r to avoid this problem.
    cat "$fq" | \
    while read -r line ; do
        #n=$((n + 1))
        if [ "$line" = "+" ] ; then
            count=$((count - 1))	# discount the title line
            #echo "XXX $count $line"
            #while [ $count -gt 0 ] ; do
            while [ $deleted -lt $count ] ; do
                deleted=$((deleted + 1))
                read -r line
                #n=$((n + 1))
                #echo ">>>" $count $deleted $n $line
            done
            count=0	# reset counter
            deleted=0
        else
            #echo $count: $line
            echo $line
            count=$((count + 1))
         fi
    done | \
    tr '@' '>' 
    #set -x
}


function mummer3_align() {
    if [ "$MUMMER" != "yes" ] ; then return ; fi
    ref=`abspath "$1"`
    qontigs=`abspath "$2"`

    f="${qontigs##*/}"
    ext="${qontigs##*.}"
    name="${f%.*}"

    p="${qontigs##*/}" ; p="${p%.*}"
    r="${ref##*/}" ; r="${r%.*}"
    wd=mummer3_${p}_x_${r}
    #echo "mummer3 $ref $qontigs"    
    if [ ! -d $wd ] ; then
        echo "Creating $wd"
        mkdir $wd
    fi
    echo "Entering $wd"
    cd $wd
    if [ ! -e reference.fna ] ; then
        echo linking `relpath "$ref"` as reference.fna
	ln -sf `relpath "$ref"` reference.fna
        #ln -sf $ref reference.fna
        #symlinks -c .
    fi

    if [ ! -e contigs.fna ] ; then
        if [ "$ext" = "fastq" ] ; then
    	    # convert to fasta
            fq2f $qontigs > contigs.fna
        else
            ln -sf `relpath $qontigs` contigs.fna
            #ln -sf $qontigs contigs.fna
            #symlinks -c .
        fi
    fi
    out="$name"

    #set -x
    if [ ! -e $out.out ] ; then
        echo "Comparing $out contigs.fna against reference.fna"
        run-mummer3 reference.fna contigs.fna $out
	# plot
	# --oplor	color plot lines by percent similarity
	# -c		generate coverage plot
	# -f		only display best hit
	# -l            layout a .delta multiplot intelligibly
	# -R -Q
	# -t png
        mummerplot --color -c -t png -p ${out}-cov $out.out
	mummerplot            -t png -p ${out}-dp  $out.out
	
    fi
    
    if [ ! -e ${out}-concat.mums ] ; then
        echo "Concatenating contigs and aligning against refereece"
	echo ">concatenated_fragments" > ${out}-concat.fna
	grep -v "^>" contigs.fna >> ${out}-concat.fna
	mummer -mum -b -c reference.fna ${out}-concat.fna > ${out}-concat.mums
	mummerplot --color -c -t png -p ${out}-concat-cov ${out}-concat.mums
	mummerplot --color    -t png -p ${out}-concat-dp  ${out}-concat.mums
	mummer -mum -b -c ${out}-concat.fna reference.fna > ${out}-tacnoc.mums
	mummerplot --color -c -t png -p ${out}-tacnoc-cov ${out}-tacnoc.mums
	mummerplot --color    -t png -p ${out}-tacnoc-dp  ${out}-tacnoc.mums
    fi
    cd ..
    #set +x
}

function nucmer_align() {
    if [ "$NUCMER" != "yes" ] ; then return ; fi
    ref=`abspath $1`
    sqs=`abspath $2`
    
    if [ ! -e $1 -o ! -e $2 ] ; then
        echo "nucmer: $1 and $2 must exist"
        return
    fi

    p="${qontigs##*/}" ; p="${p%.*}"
    r="${ref##*/}" ; r="${r%.*}"
    wd=nucmer_${p}_x_${r}
    if [ ! -d $wd ] ; then
        echo "Creating $wd"
        mkdir $wd
    fi
    echo "Entering $wd"
    cd $wd
    if [ ! -e reference.fna ] ; then
        ln -sf `relpath $ref` reference.fna
        #ln -sf $ref reference.fna
        #symlinks -c .
    fi
    
    f="${sqs##*/}"
    ext="${sqs##*.}"
    name="${f%.*}"
    echo "nucmer_align: $sqs - $f - $name - $ext"
    if [ ! -e $name.fna ] ; then
        if [ "$ext" = "fastq" ] ; then
    	    # convert to fasta and use fasta file instead
	    echo "Converting $sqs to fasta format"
            fq2f $sqs > $name.fna
            sqs=$name.fna
        else
            ln -sf `relpath $sqs` $name.fna
            #ln -sf $sqs $name.fna
            #symlinks -c .
        fi
    fi
    out=$name

    #set -x
    #echo `awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' reference.fna`
    #reflen=`awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' \
    #        reference.fna | tail -1`
    reflen=`fastasize reference.fna | tail -1`
    echo "Referece" > SIZES
    fastasizes reference.fna >> SIZES
    echo "Probes" >> SIZES 
    fastasizes $name.fna >> SIZES

    if [ ! -s $out.delta ] ; then
        nucmer -maxmatch -p $out $ref $name.fna
	mummerplot --color -c -t png -p ${out}-cov $out.delta
	mummerplot --color    -t png -p ${out}-dp  $out.delta
    fi
    if [ ! -s $out.coords ] ; then
        show-coords -L 100 -I 50 -r -c -l $out.delta > $out.coords
	#mapview -n 1 -p ${out}-mapview $out.coords
	mapview -n 1 -f png -p ${out}-mapview $out.coords
	#mapview -n 1 -f pdf -p ${out}-mapview $out.coords
    fi
    
    if [ ! -s $out.snps ] ; then
        show-snps -C $out.delta > $out.snps
    fi
    
    if [ ! -s $out.aligns ] ; then
        show-aligns -r $out.delta $ref $sqs > $out.aligns
    fi
    
    if [ ! -s $out.tiling ] ; then
        # assume circular reference (-c)
        show-tiling -c -p $out.pseudo $out.delta > $out.tiling
	mummerplot --color -c -t png -p ${out}-tile-cov ${out}.tiling
    fi

    # compute coverage by summing-up inter-contig gaps
    if [ ! -s $out.tiling ] ; then
        echo "0" > $out.coverage
    else
        #cat $out.tiling | tail -n +2 | cut -d'	' -f3
	echo -n "coverage: "
        cat $out.tiling | tail -n +2 | cut -d'	' -f3 | \
	    awk "{s+=\$1} END {print (1-(s/$reflen))*100 }" 
        cat $out.tiling | tail -n +2 | cut -d'	' -f3 | \
	    awk "{s+=\$1} END {print (1-(s/$reflen))*100 }" > $out.coverage
    fi

    # Now, repeat for the concatenated contigs
    #
    if [ ! -e ${out}-concat.delta ] ; then
        echo "nucmer: Concatenating $sqs and comparing against reference"
        echo ">Concatenated_contigs" > ${out}-concat.fna
	grep -v "^>" $out.fna >> $out-concat.fna
	# compare contigs to ref
	nucmer -maxmatch -c 100 -p ${out}-concat $ref ${out}-concat.fna
	mummerplot --color -c -t png -p ${out}-concat-cov ${out}-concat.delta
	mummerplot --color    -t png -p ${out}-concat-dp  ${out}-concat.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-concat.delta > ${out}-concat.coords
	mapview -n 1 -f png -p ${out}-concat-map ${out}-concat.coords
	
	show-tiling ${out}-concat.delta > ${out}-concat.tiling
	mummerplot --color -c -t png -p ${out}-concat-tiling-cov ${out}-concat.tiling
	
	# compare ref to contigs
	nucmer -maxmatch -c 100 -p ${out}-tacnoc ${out}-concat.fna $ref
	mummerplot --color -c -t png -p ${out}-tacnoc-cov ${out}-tacnoc.delta
	mummerplot --color    -t png -p ${out}-tacnoc-dp  ${out}-tacnoc.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-tacnoc.delta > ${out}-tacnoc.coords
	mapview -n 1 -f png -p ${out}-tacnoc ${out}-tacnoc.coords

	show-tiling ${out}-tacnoc.delta > ${out}-tacnoc.tiling
	mummerplot --color -c -t png -p ${out}-tacnoc-tiling-cov ${out}-tacnoc.tiling
    fi

    cd ..
    #set +x
}

function last_reference_align() {
    if [ "$LAST" != "yes" ] ; then return ; fi
    #set -x
    probe=`abspath $1`
    ref=`abspath $2`

    f="${probe##*/}"
    ext="${probe##*.}"
    out="${f%.*}"

    if [ ! -e $ref -o ! -e $probe ] ; then
        echo "Error, $ref and $probe must exist"
	return
    fi

    p="${probe##*/}" ; p="${p%.*}"
    r="${ref##*/}" ; r="${r%.*}"
    wd=last_${p}_x_${r}
    if [ ! -d $wd ] ; then
        echo "Creating $wd"
        mkdir $wd
    fi
    echo "Entering $wd"
    cd $wd
    
    if [ ! -e reference.fna ] ; then
    	ln -sf `relpath $ref` reference.fna
    	#ln -sf $ref reference.fna
        #symlinks -c .
    fi
    if [ ! -e $out.fastq ] ; then
        if [ "$ext" = "fastq" ] ; then
    	    ln -sf `relpath $probe` $out.fastq
    	    #ln -sf $probe probe.fastq
            #symlinks -c .
            fq2f $out.fastq > $out.fna
        else
	    ln -s $probe $out.fna
	fi
    fi
    if [ ! -e refdb.prj ] ; then
    	echo "Making default db"
    	lastdb -c refdb reference.fna
    fi
    if [ ! -e $out.maf ] ; then
    	echo "aligning"
	# we use here the less stringent db and fasta without quality codes
	#	to get a minimum restriction (maximal information) alignment
    	lastal -Q0 refdb $out.fna > $out.maf
    	maf-convert html $out.maf > $out.html
    	maf-convert tab  $out.maf > $out.tab
    fi
    if [ ! -e $out.png ] ; then
        last-dotplot $out.maf $out.png
    fi
    
    # to require strong alignments:
    if [ ! -e refdbs.prj ] ; then
    	echo "Making strong db"
    	lastdb -m1111110 refdbs reference.fna
    fi
    if [ ! -e $out.strong.maf ] ; then
        # to indicate fastq-sanger format use -Q1, to indicate fasta -Q0
        echo "aligning with strong DB"
        lastal -Q1 refdbs $out.fastq > $out.strong.maf
    fi
    if [ ! -e $out.strong.png ] ; then
        last-dotplot $out.strong.maf $out.strong.png
    fi
    if [ ! -e $out.e1000.maf ] ; then
        # to ask for score > 1000
        echo "aligning score > 1000 with strong DB"
        lastal -Q1 -e1000 refdbs $out.fastq > $out.e1000.maf
    fi
    if [ ! -e $out.e1000.png ] ; then
        last-dotplot $out.e1000.maf $out.e1000.png
    fi
    if [ ! -e $out.1.maf ] ; then
	#set -x

        # to look for unique best alignments
        echo "selecting unique alignments"
	#
	# for some reason the fastq-based alignments fail to compute
        #last-split -g refdbs $out.strong.maf > $out.strong.1.maf
        last-split $out.strong.maf > $out.strong.1.maf
        last-dotplot $out.strong.1.maf $out.strong.1.png

        #last-split -g refdbs $out.e1000.maf > $out.e1000.1.maf
        last-split $out.e1000.maf > $out.e1000.1.maf
        last-dotplot $out.e1000.1.maf $out.e1000.1.png

	# but the fasta-based alignment has no problem 
	#last-split -g refdb  $out.maf > $out.1.maf
	last-split $out.maf > $out.1.maf
        last-dotplot $out.1.maf $out.1.png

        # sort by seqname, strand, start, end, of the top sequence
        maf-sort $out.1.maf > $out.1.srt.maf
	maf-convert tab $out.1.srt.maf > $out.1.srt.maf.tab
	
	#set +x
    fi
    for i in *.maf ; do
        echo -n $i
        if [ ! -e `basename $i maf`png ] ; then
            echo " ..plotting"# `pwd`/$i
            last-dotplot $i `basename $i maf`png
        else
            echo "."
        fi
    done
    
    cd ..
}


function usage()
{
    echo "Usage: $0 -m -l -n -v -h"
    echo ""
    echo "       $0 --mummer3 --last --nucmer contigs.fna reference.fna"
    echo ""
    echo "    -m --mummer3     Compare fasta contigs against reference"
    echo "                     using MUMMER3"
    echo "    -l --last        Compare fasta contigs against reference"
    echo "                     using LAST aligner"
    echo "    -n --nucmer      Compare fasta contigs against reference"
    echo "                     using NUCMER"
    echo "                     The comparison options require a reference"
    echo "                     genome in fasta format, available as "
    echo "                     refs/reference.fna"
    echo "    -h --help        Pint this help and exit"
    echo "    -v --verbose     Output additional information"
    echo ""
    exit
}

# --------------------------------------------------------------------------
#                              ASSEMBLE MATES
# --------------------------------------------------------------------------

# Parse the command line
# ----------------------
# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o hvribBMAmln \
     --long help,verbose,roche454,illumina,bwa,bowtie2,maq,a5,mummer3,last,nucmer \
     -n "$0" -- "$@"`

# an invalid option was given
if [ $? != 0 ] ; then usage >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
        case "$1" in
                -h|--help) 
                    usage ; shift ;;

                -v|--verbose) 
                    VERBOSE=$(( $VERBOSE + 1 )) ; shift ;;
                    
                -m|--mummer3)
                    MUMMER="yes" ; shift ;;
                
                -l|--last)
                    LAST="yes" ; shift ;;
                
                -n|--nucmer)
                    NUCMER="yes" ; shift ;;

                --) shift ; break ;;
                *) echo "Internal error!" >&2 ; usage ; exit 1 ;;
        esac
done

if [ $VERBOSE -gt 0 ] ; then
cat <<END
Selected values:

MUMMER="$MUMMER"			# -m --mummer3
LAST="$LAST"				# -l --last
NUCMER="$NUCMER"			# -n --nucmer
VERBOSE=$VERBOSE                       	# -v --verbose

END
fi


# Prepare input
# =============

# Analyze sequences
#
echo "Aligning with MUMMER3"
mummer3_align $1 $2
echo "Aligning with NUCMER"
nucmer_align $1 $2
echo "Aligning with LAST"
last_reference_align $1 $2

