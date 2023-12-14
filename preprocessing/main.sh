#!/bin/bash

## Defining some colors to echo pipeline
BBlue='\033[1;94m'
BRed='\033[1;91m'
NC='\033[0m' # No Color

BASE_DIR=${1}
TYPE=${2}

[ -d ${BASE_DIR}/hisat2 ] || mkdir ${BASE_DIR}/hisat2
[ -d ${BASE_DIR}/rsem ] || mkdir ${BASE_DIR}/rsem

## setting hisat2 and rsem references paths for fresh or frozen samples
REF_HISAT2_FROZEN=/home/jjl78/genome/human/hg19_hisat2_tran_index_premrna_clean/hisat2_tran_index_premrna_hg19
REF_RSEM_FROZEN=/home/jjl78/genome/human/hg19_rsem_tran_index_premrna_clean/rsem_tran_index_premrna_hg19
REF_HISAT2_FRESH=/home/jjl78/genome/human/hg19_hisat2_tran_index_clean/hisat2_tran_index_hg19
REF_RSEM_FRESH=/home/jjl78/genome/human/hg19_rsem_tran_index_clean/rsem_tran_index_hg19

echo -e "${BBlue}STEP 01: Aligning and Mapping reads (hisat2 and rsem)${NC}"

a=$(ls -1 ${BASE_DIR}/Fastq | cut -f1 -d"_" | sed 's/-/_/g' | uniq)


ls -1 ${BASE_DIR}/Fastq | cut -f1 -d"_" | uniq | while read SAMPLE; do 
  FQ1=${SAMPLE}_R1.fastq.gz
  FQ2=${SAMPLE}_R2.fastq.gz

  if [[ $TYPE == 'frozen' ]]
  then
    hisat2 -t \
      -x ${REF_HISAT2_FROZEN} \
      -1 ${BASE_DIR}/Fastq/${FQ1} \
      -2 ${BASE_DIR}/Fastq/${FQ2} \
      --rg-id=${SAMPLE} --rg SM:${SAMPLE} --rg LB:${SAMPLE} \
      --rg PL:ILLUMINA --rg PU:${SAMPLE} \
      --new-summary --summary-file ${BASE_DIR}/hisat2/${SAMPLE}.log \
      --met-file ${BASE_DIR}/hisat2/${SAMPLE}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p 3 -S ${BASE_DIR}/hisat2/${SAMPLE}.sam

    rsem-calculate-expression --estimate-rspd --paired-end -sam -p 8 ${BASE_DIR}/hisat2/${SAMPLE}.sam ${REF_RSEM_FROZEN} ${BASE_DIR}/rsem/${SAMPLE}
  elif [[ $TYPE == 'fresh' ]]
  then
    hisat2 -t \
      -x ${REF_HISAT2_FRESH} \
      -1 ${BASE_DIR}/Fastq/${FQ1} \
      -2 ${BASE_DIR}/Fastq/${FQ2} \
      --rg-id=${SAMPLE} --rg SM:${SAMPLE} --rg LB:${SAMPLE} \
      --rg PL:ILLUMINA --rg PU:${SAMPLE} \
      --new-summary --summary-file ${BASE_DIR}/hisat2/${SAMPLE}.log \
      --met-file ${BASE_DIR}/hisat2/${SAMPLE}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p 3 -S ${BASE_DIR}/hisat2/${SAMPLE}.sam

    rsem-calculate-expression --estimate-rspd --paired-end -sam -p 8 ${BASE_DIR}/hisat2/${SAMPLE}.sam ${REF_RSEM_FRESH} ${BASE_DIR}/rsem/${SAMPLE}
  fi
  
  mv ${BASE_DIR}/rsem/${SAMPLE}.stat/${SAMPLE}.cnt ${BASE_DIR}/rsem
  rm -rfv ${BASE_DIR}/rsem/${SAMPLE}.transcript.bam
  rm -rfv ${BASE_DIR}/hisat2/${SAMPLE}.sam
done

echo -e "${BBlue}\nSTEP 02: Aggregating results (counts and tpm)${NC}"
Rscript /n/scratch/users/s/sad167/EPN/scripts/preprocessing/aggregate.R ${BASE_DIR}/rsem

echo -e "${BBlue}STEP 03: Creating QC files${NC}"
Rscript n/scratch/users/s/sad167/EPN/scripts/preprocessing/qc_rsem.R ${BASE_DIR}/rsem

echo -e "${BBlue}STEP 04: Cleaning files/directories${NC}"
rm -rfv ${BASE_DIR}/rsem/*.stat
ls -d ${BASE_DIR}/rsem/*/ | xargs rm -rf
gzip ${BASE_DIR}/rsem/*