#!/bin/bash
#
# Compare a set of (one or more) sequences against a reference genome
#
# We take care not to repeat any work, hence, re-running is safe.
#
# The graphics routines may fail on servers where no X support is available.
# If so, in order to re-run them, we need to remove the generated last, 
# mummer3 and nucmer directories and re-run in a different computer with
# appropriate support.
#
#
# (C) José R. Valverde, CNB/CSIC. 2015-2016
#
set +x

# Default values
MUMMER="no"			# -m --mummer3
LAST="no"			# -l --last
NUCMER="no"			# -n --nucmer
VERBOSE=0                       # -v --verbose
QUERY="aligned.cons.fasta"	# -q --query $cons
REF="reference.fna"		# -r --reference $ref
TYPE="png"			# -t --type x11 postscript png
SIZE="small"			# -s --size small medium large



function abspath() {
    # generate absolute path from relative path
    # $1     : relative filename
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    #elif [ -f "$1" ]; then
    elif [ -e "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    else
        #echo "/dev/null"
        echo " ---FILE-$1-DOES-NOT-EXIST--- "
    fi
}

function relpath() {
    # both $1 and $2 are absolute paths beginning with /
    # returns relative path to $2($target) from $1($source)
    if [ $# -eq 2 ] ; then
       local source="$1"
       local target="$2"
    elif [ $# -eq 1 ] ; then
       local source=`pwd`
       local target="$1"
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

    local common_part=$source # for now
    local result="" # for now

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
    local current="${2:+"$1"}"
    local target="${2:-"$1"}"
    [ "$target" != . ] || target=/
    target="/${target##/}"
    [ "$current" != . ] || current=/
    current="${current:="/"}"
    current="/${current##/}"
    local appendix="${target##/}"
    local relative=''
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
    local fq=`abspath "$1"`
    #echo "#fq2f: Converting $fq to fasta format"
    if [ ! -e "$fq" ] ; then
        echo "fq2f: error, nothing to do"
        return
    fi
    local count=0
    local deleted=0
    local n=0
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
    local ref=`abspath "$1"`
    local qontigs=`abspath "$2"`
    local type=$TYPE	# externally defined!
    local size=$SIZE	# externally defined!

    local f="${qontigs##*/}"
    local ext="${qontigs##*.}"
    local name="${f%.*}"

    local p="${qontigs##*/}" ; p="${p%.*}"
    local r="${ref##*/}" ; r="${r%.*}"
    local wd=mummer3_${p}_x_${r}
    #echo "mummer3 $ref $qontigs $wd"    
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
    local out="$name"

    #set -x
    if [ ! -e $out.out ] ; then
        echo "Comparing $out contigs.fna against reference.fna"
        run-mummer3 reference.fna contigs.fna $out
	# plot
	# --color	color plot lines by percent similarity
	# -c		generate coverage plot
	# -f		only display best hit
	# -l            layout a .delta multiplot intelligibly
	# -R -Q
	# -t png
        mummerplot --color -c -t $type -s $size -p ${out}-cov $out.out
	mummerplot            -t $type -s $size -p ${out}-dp  $out.out
    fi	
    
    if [ ! -e ${out}-concat.mums ] ; then
        echo "Concatenating contigs and aligning against refereece"
	echo "  by the way, we'll also remove all 'n's!"
	echo ">concatenated_fragments" > ${out}-concat.fna
	#grep -v "^>" contigs.fna >> ${out}-concat.fna
	grep -v "^>" contigs.fna | tr -d n >> ${out}-concat.fna
	mummer -mum -b -c reference.fna ${out}-concat.fna > ${out}-concat.mums
	mummerplot --color -c -t $type -s $size -p ${out}-concat-cov ${out}-concat.mums
	mummerplot --color    -t $type -s $size -p ${out}-concat-dp  ${out}-concat.mums
	mummer -mum -b -c ${out}-concat.fna reference.fna > ${out}-tacnoc.mums
	mummerplot --color -c -t $type -s $size -p ${out}-tacnoc-cov ${out}-tacnoc.mums
	mummerplot --color    -t $type -s $size -p ${out}-tacnoc-dp  ${out}-tacnoc.mums
    fi
    cd ..
    #set +x
}

function nucmer_align() {
    if [ "$NUCMER" != "yes" ] ; then return ; fi
    local ref=`abspath $1`
    local sqs=`abspath $2`
    local type=$TYPE	# externally defined!
    local size=$SIZE	# externally defined!
    if [ "$type" == "postscript" ] ; then
        local mvtyp="ps"
    else
    	local mvtyp=$type
    fi

    
    if [ ! -e $1 -o ! -e $2 ] ; then
        echo "nucmer: $1 and $2 must exist"
        return
    fi

    local p="${sqs##*/}" ; p="${p%.*}"
    local r="${ref##*/}" ; r="${r%.*}"
    local wd=nucmer_${p}_x_${r}
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
    
    local f="${sqs##*/}"
    local ext="${sqs##*.}"
    local name="${f%.*}"
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
    local out=$name

    #set -x
    #echo `awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' reference.fna`
    #reflen=`awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' \
    #        reference.fna | tail -1`
    local reflen=`fastasize reference.fna | tail -1`
    echo "Referece" > SIZES
    fastasizes reference.fna >> SIZES
    echo "Probes" >> SIZES 
    fastasizes $name.fna >> SIZES

    if [ ! -s $out.delta ] ; then
        nucmer -maxmatch -p $out $ref $name.fna
	mummerplot --color -c -t $type -s $size -p ${out}-cov $out.delta
	mummerplot --color    -t $type -s $size -p ${out}-dp  $out.delta
    fi
    if [ ! -s $out.coords ] ; then
        show-coords -L 100 -I 50 -r -c -l $out.delta > $out.coords
	#mapview -n 1 -p ${out}-mapview $out.coords
	mapview -n 1 -f $mvtyp -p ${out}-mapview $out.coords
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
	mummerplot --color -c -t $type -s $size -p ${out}-tile-cov ${out}.tiling
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
        echo "nucmer: Concatenating contigs and aligning against refereece"
	echo "  by the way, we'll also remove all 'n's!"
	echo ">concatenated_fragments" > ${out}-concat.fna
	#grep -v "^>" contigs.fna >> ${out}-concat.fna
	grep -v "^>" contigs.fna | tr -d n >> ${out}-concat.fna
	# compare contigs to ref
	nucmer -maxmatch -c 100 -p ${out}-concat $ref ${out}-concat.fna
	mummerplot --color -c -t $type -s $size -p ${out}-concat-cov ${out}-concat.delta
	mummerplot --color    -t $type -s $size -p ${out}-concat-dp  ${out}-concat.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-concat.delta > ${out}-concat.coords
	mapview -n 1 -f $mvtyp -p ${out}-concat-map ${out}-concat.coords
	
	show-tiling ${out}-concat.delta > ${out}-concat.tiling
	mummerplot --color -c -t $type -s $size -p ${out}-concat-tiling-cov ${out}-concat.tiling
	
	# compare ref to contigs
	nucmer -maxmatch -c 100 -p ${out}-tacnoc ${out}-concat.fna $ref
	mummerplot --color -c -t $type -s $size -p ${out}-tacnoc-cov ${out}-tacnoc.delta
	mummerplot --color    -t $type -s $size -p ${out}-tacnoc-dp  ${out}-tacnoc.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-tacnoc.delta > ${out}-tacnoc.coords
	mapview -n 1 -f $mvtyp -p ${out}-tacnoc ${out}-tacnoc.coords

	show-tiling ${out}-tacnoc.delta > ${out}-tacnoc.tiling
	mummerplot --color -c -t $type -s $size -p ${out}-tacnoc-tiling-cov ${out}-tacnoc.tiling
    fi

    cd ..
    #set +x
}


function maf_to_tab_ends() {
    # take a MAF file and convert it into something more useful that
    # can be imported into a spreadsheet. It will contain the initial
    # coordinate, fragment length, computed termination coordinate, etc...
    # and only a limited number of nt on each end so they can be used
    # to locate the ends on the sequence.
    
    local file=$1

    local new='n'
    # CAUTION: THIS NEEDS TO BE MADE MORE INTELLIGENT!!!
    local out=`basename $1 maf`tab

    # print out the file header
    echo "name	start	length	end	strand	tot	5'up	3'down" > $out

    cat $1 \
    | grep -v '^#' \
    | while read line ; do
	#echo $line
	if [[ ${line:0:1} == "a" ]] ; then
	#if [[ $line == a* ]] ; then
	#if [[ $line =~ ^/ ]]; then
            new='y'
            #echo "NEW!"
    	    continue
	fi
	if [[ $new == 'y' ]] ; then
            # sanity check
            if [[ ${line:0:1} == "a" ]] ; then
        	echo "maf_to_tab_ends: Error: unexpected line after 'a'lignment line"
        	echo "$line"
        	exit
            fi
            # split line
            read s name st len strand tot seq <<< $line 
            seqbeg=`grep -o '^.\{80\}' <<< $seq`
            seqend=`grep -o '.\{80\}$' <<< $seq`
	    end=$((st+len))
            echo "$name	$st	$len	$end	$strand	$tot	$seqbeg	$seqend" \
	          >> $out
            new='n'
	else
            continue
	fi
    done
}




### NOTE: PROBE MUST BE IN FASTQ FORMAT!!!
function last_reference_align() {
    if [ "$LAST" != "yes" ] ; then return ; fi
    #set -x
    #echo $1: `abspath $1`
    local ref=`abspath $1`
    local probe=`abspath $2`

    local f="${probe##*/}"
    local ext="${probe##*.}"
    local out="${f%.*}"

    #echo "REF: $ref ($1)"
    #echo "PROBE: $probe ($2)"
	
    if [ ! -e $ref -o ! -e $probe ] ; then
        echo "Error, $ref and $probe must exist"
	return
    fi
    
    if [ "$ext" == "fastq" -o "$ext" == "fq" ] ; then
        local FASTQ='y'
	echo "Last reference align: using FastQ"
    else
    	local FASTQ='n'
    fi

    local p="${probe##*/}" ; p="${p%.*}"
    local r="${ref##*/}" ; r="${r%.*}"
    local wd=last_${p}_x_${r}
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
        if [ "$FASTQ" == "y" ] ; then
	    # use probe file ensuring it has a ".fastq" extension
    	    ln -sf `relpath $probe` $out.fastq
    	    #ln -sf $probe probe.fastq
            #symlinks -c .
	    # generate the fasta file from fastq ourselves
            fq2f $out.fastq > $out.fna
	else
	    # use probe file ensuring it has a ".fna" extension
	    ln -s $probe $out.fna
	fi
    fi

    ###
    #	FIRST STEP: GENERAL ALIGNMENT
    ###

    if [ ! -e refdb.prj ] ; then
    	echo "Making default db"
        # default is better at long, weak alignments
        lastdb refdb reference.fna
        # mask lowercase letters
        #lastdb -c refdb reference.fna
    fi
    if [ ! -e $out.maf ] ; then
    	echo "    aligning with minimal requirements"
	# we use here the less stringent db and fasta without quality codes
	#	to get a minimum restriction (maximal information) alignment
    	lastal -Q0 refdb $out.fna > $out.maf
    	maf-convert html $out.maf > $out.html
    	maf-convert tab  $out.maf > $out.tab
    fi
    if [ ! -e $out.png ] ; then
        last-dotplot $out.maf $out.png
    fi
    if [ ! -e $out.1.srt.maf ] ; then
	echo "    singling out alignment"
	#last-split -g refdb  $out.maf > $out.1.maf
	last-split $out.maf > $out.1.maf
        last-dotplot $out.1.maf $out.1.png

        # sort by seqname, strand, start, end, of the top sequence
        maf-sort $out.1.maf > $out.1.srt.maf
	maf-convert tab $out.1.srt.maf > $out.1.srt.maf.tab
	maf_to_tab_ends $out.1.srt.maf
    fi


    ###
    #	SECOND STEP: STRONG ALIGNMENTS
    ###
    # We now want to require short, strong alignments:
    if [ ! -e refdbs.prj ] ; then
    	echo "Making strong db"
    	lastdb -m1111110 refdbs reference.fna
    fi
    if [ "$FASTQ" == 'y' ] ; then
	# align with the short/strong reference database using FastQ
	if [ ! -e $out.strong.fq.maf ] ; then
            # to indicate fastq-sanger format use -Q1, to indicate fasta -Q0
            echo "    aligning with strong DB using FastQ"
            lastal -Q1 refdbs $out.fastq > $out.strong.fq.maf
	fi
	if [ ! -e $out.strong.fq.png ] ; then
            last-dotplot $out.strong.fq.maf $out.strong.fq.png
	fi
	if [ ! -e $out.strong.fq.1.srt.maf ] ; then
	    echo "    singling out strong alignment using FastQ"
	    last-split -g refdbs $out.strong.fq.maf > $out.strong.fq.1.maf
            last-dotplot $out.strong.fq.1.maf $out.strong.fq.1.png

            # sort by seqname, strand, start, end, of the top sequence
            maf-sort $out.strong.fq.1.maf > $out.strong.fq.1.srt.maf
	    maf-convert tab $out.strong.fq.1.srt.maf > $out.strong.fq.1.srt.maf.tab
	fi
    fi
    
    # Do the alignment with FastN sequence in any case
    if [ ! -e $out.strong.1.srt.maf ] ; then
	echo "    singling out strong alignment using FastN"
	# for some reason the fastq-based alignments fail to compute
        # but the fasta ones do work
        # so, we will recompute the alignments using fasta format

	lastal refdbs $out.fna > $out.strong.maf
        last-split $out.strong.maf > $out.strong.1.maf
        last-dotplot $out.strong.1.maf $out.strong.1.png

        # sort by seqname, strand, start, end, of the top sequence
        maf-sort $out.strong.1.maf > $out.strong.1.srt.maf
	maf-convert tab $out.strong.1.srt.maf > $out.strong.1.srt.maf.tab
	maf_to_tab_ends $out.strong.1.srt.maf
    fi

    # third step: additionally ask for a score > 1000
    if [ "$FASTQ" == 'y' ] ; then
	if [ ! -e $out.e1000.fq.maf ] ; then
            # to ask for score > 1000
            echo "    aligning score > 1000 with strong DB using FastQ"
            lastal -Q1 -e1000 refdbs $out.fastq > $out.e1000.fq.maf
	fi
	if [ ! -e $out.e1000.fq.png ] ; then
            last-dotplot $out.e1000.fq.maf $out.e1000.fq.png
	fi
	if [ ! -e $out.e1000.fq.1.srt.maf ] ; then
	    #set -x
	    echo "    singling out alignment with score > 1000 using fastQ"
	    last-split -g refdbs $out.e1000.fq.maf > $out.e1000.fq.1.maf
            last-dotplot $out.e1000.fq.1.maf $out.e1000.fq.1.png

            # sort by seqname, strand, start, end, of the top sequence
            maf-sort $out.e1000.fq.1.maf > $out.e1000.fq.1.srt.maf
	    maf-convert tab $out.e1000.fq.1.srt.maf > $out.e1000.fq.1.srt.maf.tab
	    #set +x
	fi
    fi
    # try with score > 1000 and FastN sequence in any case
    if [ ! -e $out.e1000.1.srt.maf ] ; then
	#set -x
	echo "    singling out alignment with score > 1000"

	lastal -e1000 refdbs $out.fna > $out.e1000.maf
        last-split $out.e1000.maf > $out.e1000.1.maf
        last-dotplot $out.e1000.1.maf $out.e1000.1.png

        # sort by seqname, strand, start, end, of the top sequence
        maf-sort $out.e1000.1.maf > $out.e1000.1.srt.maf
	maf-convert tab $out.e1000.1.srt.maf > $out.e1000.1.srt.maf.tab
	#set +x
    fi
    
    echo "last_reference_align: making useful TAB files"
    # Plot any unplotted maf files
    for i in *.maf ; do
        echo -n $i
        if [ ! -e `basename $i maf`tab ] ; then
            echo " ..converting" 
            maf_to_tab_ends $i
        else
            echo "."
        fi
    done

    echo "last_reference_align: plotting MAF files"
    # Plot any unplotted maf files
    for i in *.maf ; do
        echo -n $i
        if [ ! -e `basename $i maf`png ] ; then
            echo " ..plotting" 
            last-dotplot $i `basename $i maf`png
        else
            echo "."
        fi
    done
    
    cd ..
}


function usage()
{
    echo "Usage: $0 -m -l -n -v -h -p probes -r reference"
    echo ""
    echo "       $0 --mummer3 --last --nucmer "
    echo "            --probes contigs.fna --reference reference.fna"
    echo ""
    echo "    -m --mummer3     Compare fasta contigs against reference"
    echo "                     using MUMMER3"
    echo "    -l --last        Compare fasta contigs against reference"
    echo "                     using LAST aligner"
    echo "    -n --nucmer      Compare fasta contigs against reference"
    echo "                     using NUCMER"
    echo "    -q --query       File containing one or more sequeces"
    echo "    -r --reference   File containing one reference genome"
    echo "    -s --size        small or medium or large"
    echo "    -t --type        x11 or png or postscript"
    echo "    -h --help        Pint this help and exit"
    echo "    -v --verbose     Output additional information"
    echo ""
    exit
}

# --------------------------------------------------------------------------
#                              COMPARE GENOMES
# --------------------------------------------------------------------------

# Parse the command line
# ----------------------
# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o q:r:s:t:hvmln \
     --long query:,reference:,size:,type:,help,verbose,mummer3,last,nucmer \
     -n "$0" -- "$@"`

# an invalid option was given
if [ $? != 0 ] ; then usage >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
#  the set command takes any arguments after the options (-- signals the end
#  of the options) and assigns them to positional parameters $0..$n.
#  We pass set through eval so that bash will honor embedded quotes in the
#  string, i.e.
#	whithout eval, "ab cd" is split at the space into "ab and cd"
#	with eval, bash honors the quotes and gives "ab cd" as a single value
#
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
                
                -q|--query)
                    QUERY="$2" ; shift 2 ;;
                    
                -r|--reference)
                    REFERENCE="$2" ; shift 2 ;;
                    
		-s|--size)
		    SIZE="$2" ; shift 2 ;;
		    
		-t|--type)
		    TYPE="$2" ; shift 2 ;;
                --) shift ; break ;;
                *) echo "$0: Internal error!" >&2 ; usage ; exit 1 ;;
        esac
done

# test size and type
function is_valid () {
    local opt=$1 ; shift
    local values=( $@ )
    #echo $opt ${values[@]}
    for (( i = 0; i < ${#values[@]}; i++ )) ; do
       #echo "$i ${values[$i]} $opt"
       if [ "${values[$i]}" == "$opt" ]; then
           return 0
       fi
    done
    #echo no match
    return 1

}

if ! is_valid $SIZE 'small' 'medium' 'large' ; then
    usage
    exit
fi
if ! is_valid $TYPE 'x11' 'png' 'postscript' ; then
    usage
    exit
fi

#exit
if [ $VERBOSE -gt 0 ] ; then
    cat <<END
Selected values:

MUMMER="$MUMMER"
LAST="$LAST"
NUCMER="$NUCMER"
VERBOSE="$VERBOSE"
QUERY="$QUERY"
REFERENCE="$REFERENCE"
SIZE="$SIZE"
TYPE="$TYPE"

END

fi

QRYD=`dirname "$QUERY"`	# dir	(we use this to ensure we get at least '.')
QRYF="${QUERY##*/}"	# file
QRYE="${QRYF##*.}"	# ext
QRY="${QRYF%.*}"	# nam

REFD=`dirname "$REFERENCE"`
REFF="${REFERENCE##*/}"
REFE="${REFF##*.}"
REF="${REFF%.*}"

if [ $VERBOSE -gt 0 ] ; then
    # redirect all output to log files (one for stdout and one for stderr)
    #log=./LOG.$myname.$QRY.$REF.out
    #err=./LOG.$myname.$QRY.$REF.err
    #exec > $log 2> $err
    #
    # use a single log file, backing it up to an increasing count if one
    # already exists.
    #
    log=./log.$QRY.$REF.std
    if [ -e $log ] ; then
	cnt=1
	while [ -e $log.`printf %03d $cnt` ] ; do
            cnt=$((cnt + 1))
	done
	mv $log   $log.`printf %03d $cnt`
    fi
    exec > $log 2>&1
    #exec |& tee $log
fi



# Prepare input
# =============
# if either is a fastq file and there is no fasta version to use
# then make a fasta counterpart
if [ "$QRYE" == "fastq" -a "$QRYE" == "fq" ] ; then
    if [ ! -e "$QRYD/$QRY".fasta -o ! -e "$QRYD/$QRY".fna ] ; then
        fq2f "$QUERY" > "$QRYD/$QRY".fna
        probe="$QRYD/$QRY.fna"
    elif [ -e "$QRYD/$QRY.fna" ] ; then
        probe="$QRYD/$QRY.fna"
    else
        probe="$QRYD/$QRY.fna"
    fi
else
    # assume it is in fasta format
    probe="$QUERY"
fi
if [ "$REFE" == "fastq" -o "$REFE" == "fq" ] ; then
    if [ ! -e "$REFD/$REF".fasta -o ! -e "$REFD/$REF".fna ] ; then
        fq2f "$REFERENCE" > "$REFD/$REF".fna
        ref="$REFD/$REF.fna"
    elif [ -e "$REFD/$REF.fna" ] ; then
        ref="$REFD/$REF.fna"
    else
        ref="$REFD/$REF.fna"
    fi
else
    # assume it is in fasta format
    ref="$REFERENCE"
fi

# Analyze sequences
#
if [ "$MUMMER" == "yes" ] ; then
    echo "Aligning with MUMMER3"
    mummer3_align $ref $probe
fi
if [ "$NUCMER" == "yes" ] ; then
    echo "Aligning with NUCMER"
    nucmer_align $ref $probe
fi
if [ "$LAST" == "yes" ] ; then
    #echo "REF: $ref"
    #echo "QUERY: $QUERY"
    #echo "PROBE: $probe"
    #echo "Aligning with LAST"
    if [ "$QRYE" == "fastq" -a "$QRYE" == "fq" ] ; then
	# prefer FASTQ whenever possible
	last_reference_align $ref $QUERY
    else
        last_reference_align $ref $probe
    fi
fi
