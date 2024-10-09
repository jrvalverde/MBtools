#!/bin/bash

file=$1

new='n'

echo "name	start	length	end	strand	tot	5'up	3'down"

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
            echo "Error: unexpected line after 'a'lignment line"
            echo "$line"
            exit
        fi
        # split line
        read s name st len strand tot seq <<< $line 
        seqbeg=`grep -o '^.\{80\}' <<< $seq`
        seqend=`grep -o '.\{80\}$' <<< $seq`
	end=$((st+len))
        echo "$name	$st	$len	$end	$strand	$tot	$seqbeg	$seqend"
        new='n'
    else
        continue
    fi
done
