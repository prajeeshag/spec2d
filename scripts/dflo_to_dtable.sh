#!/bin/bash

dflo=$1

dflo="${dflo:-diag_field_log.out}"

dtable=diag_table

tail -n +2 "$dflo" | sed 's/|/ /g' | awk '{print $1", "$2", "$2", "$1"dstamp, all, .true., none, 2"}' >> ${dtable}_temp
sed -i 's/dstamp/%4yr%2mo%2dy%2hr%2mi%2sc/g' ${dtable}_temp

echo AMFI > $dtable
echo 1900 01 01 00 00 00 >> $dtable
echo " " >> $dtable
echo "#----------Files-----------------------------------#" >> $dtable
awk '{print $4}' ${dtable}_temp | sort | uniq | awk '{print $1" 1, months, 1, days, time, 1, months"}' >> $dtable
echo " " >> $dtable
echo " " >> $dtable
echo "#----------Fields-----------------------------------#" >> $dtable
echo " " >> $dtable
cat ${dtable}_temp >> ${dtable}
echo " " >> $dtable
echo "#----------End diag_table-----------------------------------#" >> $dtable

rm -f ${dtable}_temp
