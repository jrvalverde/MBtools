
function plotcov.ch() {
    # by allowing use as a filter we can avoid generating space-demanding
    # temporary files
    if [ $# -ge 1 ] ; then
        local depth=$1
    else
        local depth='-'
    fi
    # check there is an input
    if [ ! -s "$depth" -a "$depth" != "-" ] ; then 
    	echo "Error: file $depth does not exist"
        return 1 ; 
    fi
    
    if [ $# -gt 1 ] ; then local chr="$2" ; else chr="all" ; fi
    if [ "$chr" == "all" ] ; then pat="[[:alnum:]]*" ; else pat="$chr" ; fi

    echo "plotting coverage for chromosome(s) $chr ($pat) from $depth"
    
    # we assume $depth is in REF-CHR POS DEPTH format
    # feedgenuplot requires it in POS REF_CHR DEPTH format, hence the awk
    cat "$depth" \
    | grep "^$pat	" \
    | awk '{print $2 "\t" $1 "\t" $3 }' \
    | feedgnuplot --lines \
        	  --hardcopy "$depth.$chr.png" \
        	  --autolegend \
        	  --dataid \
        	  --domain \
        	  --ylabel depth \
        	  --xlabel "chr pos" \
        	  --title "Chromosome-wide depth" \
        	  --exit 

    display -alpha remove "$depth.$chr.png" &
}

