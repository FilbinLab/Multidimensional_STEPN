while IFS=, read -r SampleName SampleDeidName; do
  sh /n/scratch/users/s/sad167/EPN/scRNAseq/scripts/7b-NMF_rank-batch.sh $SampleName
done </n/scratch/users/s/sad167/EPN/scRNAseq/scripts/metadata_NMF.csv