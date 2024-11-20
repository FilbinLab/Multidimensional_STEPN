#!/bin/bash

# Loop through each line in the input file
while IFS=, read -r SampleName SampleID; do
  # Remove leading and trailing whitespace (including carriage return characters)
  SampleName=$(echo "$SampleName" | tr -d '\r')

  # Execute your shell script with the sanitized FolderName
  sh /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/3c_NicheAnalysis.sh $SampleName
done < /n/scratch/users/s/sad167/EPN/Xenium/scripts_revisions/SampleIdentifier.csv