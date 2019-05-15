#!/bin/bash

#### ensure you're inside the folder containing the sub-folders with the sequence files #### 
#1. Unzip the fastq.compressed
for f in */*.gz; do echo "$f"; done
for f in */*.gz; do gunzip "${f}"; done


#2. Assemble the R1 and R2 files with PEAR
R1=($(ls -d */*R1_001.fastq))
R2=($(ls -d */*R2_001.fastq))

output=(${R1[@]/R1_001.fastq/PEAR})
length=${#R1[@]}
length=$(($length-1))

for ((i=0;i<=$length;i++)); do echo $i; echo ${R1[$i]}; echo ${R2[$i]}; echo ${output[$i]}; done
for ((num=0;num<=$length;num++)); do echo $num; pear -n 86 -m 86 -f ${R1[$num]} -r ${R2[$num]} -o ${output[$num]} > ${output[$num]}_log; done


#3. Remove every other file from the PEAR assembly analysis except for the assembled file
for f in */*001.fastq; do echo "$f"; done
for f in */*001.fastq; do rm "$f"; done

for f in */*PEAR.unassembled*.fastq; do echo "$f"; done
for f in */*PEAR.unassembled*.fastq; do rm "$f"; done

for f in */*discarded.fastq; do echo "$f"; done
for f in */*discarded.fastq; do rm "$f"; done

#4. Remove exact reads using rmdup from dedupe.sh

SCRIPT_PATH="/home/ucbtomi/bbmap/dedupe.sh"

for f in */*assembled.fastq; do echo "$f"; done
for f in */*assembled.fastq; do source "$SCRIPT_PATH" "$f" out="${f}_rmdup.fastq"; done

#5. Run barcount on dntag
for f in */*dn*rmdup.fastq; do barcount --fastq "$f" --flanking_left AGTATC --flanking_right TTTAAA --max_distance_flanks 1 --max_distance_barcode 3 --barcode_table /home/ucbtomi/dntag_reference.csv --debug --save_extracted_barcodes --verbose --umiA_position " 4:8" --umiB_position " -8:-4" --out "${f}_Barcode_filter"; done


#6. Run barcount on uptag
for f in */*up*rmdup.fastq; do barcount --fastq "$f" --flanking_left GATATC --flanking_right TTTAAA --max_distance_flanks 1 --max_distance_barcode 3 --barcode_table /home/ucbtomi/uptag_reference.csv --debug --save_extracted_barcodes --verbose --umiA_position " 4:8" --umiB_position " -8:-4" --out "${f}_Barcode_filter"; done
