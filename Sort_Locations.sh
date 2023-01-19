#!/bin/bash
# date: 19.01.23
# author: Mirjam Karlsson-Müller
# name: Sort_Locations.sh
#The location table outputs from vcf_location_table.py are not sorted. This script fixes that for easier genotype assignment.
#Note that the files are given per chromosome, so that entry can be ignored.
# To be run in folder where location tables are.

for file in *.tsv; do echo $file; sort $file -n -t "_" -k2 > ../Sorted_Location_Table/sorted_$file; done