#!/bin/bash

dflo=diag_field_log.out
dtable=diag_table_dflo

echo AMFI > $dtable
echo 1900 01 01 00 00 00 >> $dtable
echo " " >> $dtable
echo "#----------Files-----------------------------------#" >> $dtable
echo atm_out, 1, hours, 1, hours, time >> $dtable
echo atm_rad, 1, hours, 1, hours, time >> $dtable
echo " " >> $dtable
echo " " >> $dtable
echo "#----------Fields-----------------------------------#" >> $dtable
echo " " >> $dtable
tail -n +2 "$dflo" | sed 's/|/ /g' | awk '$1!="am_rad" {print $1", "$2", "$2", atm_out, all, .false., none, 2"}' >> $dtable
tail -n +2 "$dflo" | sed 's/|/ /g' | awk '$1=="am_rad" {print $1", "$2", "$2", atm_rad, all, .false., none, 2"}' >> $dtable
echo " " >> $dtable
echo "#----------End diag_table-----------------------------------#" >> $dtable
