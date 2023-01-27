#!/home/mi3258mu-se/miniconda3/envs/NewBash/bin/bash

#Starts variant script then sorts output and then starts genotype script, separately for each chromosome.

#I want it to run maximal 10 processes at a time, cause i dont want it to stop randomly or slow down.
max_children=10

#I do not know how to do this, so the function parallel is from stackoverflow
#https://stackoverflow.com/questions/3004811/how-do-you-run-multiple-programs-in-parallel-from-a-bash-script

#echo ${BASH_VERSION}

function parallel {
  local time1=$(date +"%H:%M:%S")
  local time2=""

  # 
  echo "starting $2 ($time1)..."
  "$@" && time2=$(date +"%H:%M:%S") && echo "finishing $2 ($time1 -- $time2)..." &

  local my_pid=$$
  local children=$(ps -eo ppid | grep -w $my_pid | wc -w)
  children=$((children-1))
  if [[ $children -ge $max_children ]]; then
    wait -n
  fi
}


#Now comes the part where I need to add to it. So I want it to parallel run per chromosome.
#But for each chromosome it needs to do 3 steps, each waiting for the previous to finish

function steps {
    #run variants
    #python variants_in_AS_Pipeline/vcf_location_table.py -c "chr$1:1-250000000" -s ../Sample_Data/ -o Trial_Run/chr$1_location_table.tsv >> step1.txt;
    #Exclude all lines with only ND and -
    #wait && cat Trial_Run/chr$1_location_table.tsv | grep "#\|0/1\|1/1\|1/2\|2/2\|1/3\|2/3\|3/3\|0/2\|0/3" > Trial_Run/filtered_chr$1_location_table.tsv;
    #Sort after location
    #wait && sort Trial_Run/filtered_chr$1_location_table.tsv -n -t "_" -k2 > Trial_Run/sorted_chr$1_location_table.tsv;
    #wait && python variants_in_AS_Pipeline/genotype.py -s ../Sample_Data/ -i Trial_Run/sorted_chr$1_location_table.tsv -r gene_ranges.tsv -c "chr$1:1-250000000" -o Trial_Run/chr$1_genotype_table.tsv >> step3.txt;
    #Takes input a chromosome marker. $2 cause right after function name
    #first variant script
    python vcf_location_table.py -c "chr$1:1-250000000" -s /raidset/mi3258mu-se/mirjam/Sample_Data/ -o /raidset/mi3258mu-se/mirjam/New_Variant_Run/Location_Table/chr$1_location_table.tsv;
    wait && cat /raidset/mi3258mu-se/mirjam/New_Variant_Run/Location_Table/chr$1_location_table.tsv | grep "#\|0/1\|1/1\|1/2\|2/2\|1/3\|2/3\|3/3\|0/2\|0/3" > /raidset/mi3258mu-se/mirjam/New_Variant_Run/Location_Table/filtered_chr$1_location_table.tsv;
    wait && sort /raidset/mi3258mu-se/mirjam/New_Variant_Run/Location_Table/chr$1_location_table.tsv -n -t "_" -k2 > /raidset/mi3258mu-se/mirjam/New_Variant_Run/Sorted_Location_Table/sorted_chr$1_location_table.tsv;
    wait && python genotype.py -s /raidset/mi3258mu-se/mirjam/Sample_Data/ -i /raidset/mi3258mu-se/mirjam/New_Variant_Run/Sorted_Location_Table/sorted_chr$1_location_table.tsv -r gene_ranges.tsv -c "chr$1:1-250000000" -o /raidset/mi3258mu-se/mirjam/New_Variant_Run/Sorted_Genotype_Table/chr$1_genotype_table.tsv;
}

#For loop, going through the chromosomes
chromosomes=("M" "Y" "X" 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1)

for c in ${chromosomes[@]};
    do
    echo $c
    parallel steps $c;
    done
wait
