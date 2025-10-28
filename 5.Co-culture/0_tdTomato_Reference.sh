########################################################################################
#> References per sample
########################################################################################
- 1LP, 2LP, 4LP, 6LP, 8LP, 9LP, EP1NSmono: /n/scratch/users/c/cao385/refs/tdTomato/GRCh38
- 3LP, 5LP, EP1NSrat, VBT242rat: /n/scratch/users/c/cao385/refs/tdTomato/GRCh38_and_mRatBN7-2
- 7LP: /n/scratch/users/c/cao385/refs/tdTomato/mRatBN7-2



########################################################################################
#> tdTomato GTF
########################################################################################
tdTomato_1	unknown	exon	1	2374	.	+	.	gene_id "tdTomato_1"; transcript_id "tdTomato_1"; gene_name "tdTomato_1"; gene_biotype "protein_coding";



########################################################################################
#> tdTomato FASTA
########################################################################################
>tdTomato_1
ATGGTGAGCAAGGGCGAGGAGGTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCATGGAG
GGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAG
GGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGAC
ATCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCCGACATC
CCCGATTACAAGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTC
GAGGACGGCGGTCTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCACGCTGATC
TACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAGAAGAAG
ACCATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTGAAGGGC
GAGATCCACCAGGCCCTGAAGCTGAAGGACGGCGGCCACTACCTGGTGGAGTTCAAGACC
ATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTACTACGTGGACACCAAGCTG
GACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCTCCGAGGGC
CGCCACCACCTGTTCCTGGGGCATGGCACCGGCAGCACCGGCAGCGGCAGCTCCGGCACC
GCCTCCTCCGAGGACAACAACATGGCCGTCATCAAAGAGTTCATGCGCTTCAAGGTGCGC
ATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCC
TACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCC
TGGGACATCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCC
GACATCCCCGATTACAAGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATG
AACTTCGAGGACGGCGGTCTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCACG
CTGATCTACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAG
AAGAAGACCATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTG
AAGGGCGAGATCCACCAGGCCCTGAAGCTGAAGGACGGCGGCCACTACCTGGTGGAGTTC
AAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTACTACGTGGACACC
AAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCTCC
GAGGGCCGCCACCACCTGTTCCTGTACGGCATGGACGAGCTGTACAAGTAGCAACTTTAT
TATACATAGTTGATCAATTCCGATAATCAACCTCTGGATTACAAAATTTGTGAAAGATTG
ACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCT
TTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGG
TTGCTGTCTCTTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACT
GTGTTTGCTGACGCAACCCCCACTGGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCC
GGGACTTTCGCTTTCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCC
CGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAG
CTGACGTCCTTTCCATGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGGGACGTCC
TTCTGCTACGTCCCTTCGGCCCTCAATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCG
GCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGG
GCCGCCTCCCCGCATCGGGAATTCCCGCGGTTCGCTTTAAGACCAATGACTTACAAGGCA
GCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCC
CAACGAAGACAAGATCTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGA
GCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCT
TGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA



########################################################################################
#> Human + tdTomato reference
########################################################################################
human_genome="GRCh38"
version="2024-A"

build="GRCh38"
mkdir -p "$build"

# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"

# Using release 109 for GRCh38 instead of release 110 -- release 110 moved from GRCh38.p13 to GRCh38.p14,
# which unmasked the pseudo-autosomal region. This causes ambiguous mappings to PAR locus genes.
human_fasta_url="http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
human_fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
human_gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
human_gtf_in="${source}/gencode.v44.primary_assembly.annotation.gtf"

if [ ! -f "$human_fasta_in" ]; then
    curl -sS "$human_fasta_url" | zcat > "$human_fasta_in"
fi
if [ ! -f "$human_gtf_in" ]; then
    curl -sS "$human_gtf_url" | zcat > "$human_gtf_in"
fi

# String patterns used for both genomes
ID="(ENS(RNOG)?[GTE][0-9]+)\.([0-9]+)"

MAIN_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${MAIN_PATTERN}\""
TX_PATTERN="transcript_type \"${MAIN_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Process FASTA -- translate chromosome names
human_fasta_modified="$build/$(basename "$human_fasta_in").modified"
cat "$human_fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$human_fasta_modified"

# Process GTF -- split Ensembl IDs from version suffixes
human_gtf_modified="$build/$(basename "$human_gtf_in").modified"
cat "$human_gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$human_gtf_modified"

# Process GTF -- filter based on gene/transcript tags
cat "$human_gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

human_gtf_filtered="${build}/$(basename "$human_gtf_in").filtered"
grep -E "^#" "$human_gtf_modified" > "$human_gtf_filtered" 
grep -Ff "${build}/gene_allowlist" "$human_gtf_modified" \
    | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
    >> "$human_gtf_filtered"

## Add tdTomato custom sequences 
cp /n/shared_db/GRCh38/uk/cellranger/7.0.0/7.0.0/refdata-gex-GRCh38-2020-A/fasta/genome.fa /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genome_custom.fa
cat /n/scratch/users/c/cao385/refs/tdTomato/genome.tdTomato.fa >> /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genome_custom.fa

cp /n/shared_db/GRCh38/uk/cellranger/7.0.0/7.0.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genes_custom.gtf
cat /n/scratch/users/c/cao385/refs/tdTomato/genes.tdTomato.gtf >> /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genes_custom.gtf

## Run cellranger
module load cellranger/7.1.0
cellranger mkref \
    		 --genome GRCh38 --fasta /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genome_custom.fa --genes /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genes_custom.gtf \
    		 --nthreads 12



########################################################################################
#> Rat + tdTomato reference
########################################################################################
curl -sS http://ftp.ensembl.org/pub/release-109/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz | zcat > /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
curl -sS http://ftp.ensembl.org/pub/release-109/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.109.gtf.gz | zcat > /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf
    
# Process FASTA -- translate chromosome names
cat /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.modified

cellranger mkgtf /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered_tmp  \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:antisense \
    --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene

awk 'BEGIN {OFS="\t"} {
    if ($1 !~ /^#/) {
        if ($1 ~ /^[0-9]+$/ || $1 == "X" || $1 == "Y") {
            $1 = "chr" $1
        } else if ($1 == "MT") {
            $1 = "chrM"
        }
    }
    # Join the first 8 fields with tabs
    header = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8
    # Join the attributes with spaces
    attributes = $9
    for (i=10; i<=NF; i++) {
        attributes = attributes " " $i
    }
    print header, attributes
}' /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered_tmp > /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered

cat /n/scratch/users/c/cao385/refs/tdTomato/genome.tdTomato.fa >> /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.modified
cat /n/scratch/users/c/cao385/refs/tdTomato/genes.tdTomato.gtf >> /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered

## Run cellranger
module load cellranger/7.1.0
cd /n/scratch/users/c/cao385/refs/tdTomato
cellranger mkref \
    		 --genome mRatBN7-2 --fasta /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.modified --genes /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered \
    		 --nthreads 10



########################################################################################
#> Human + Rat + tdTomato reference
########################################################################################
cellranger mkref \
    		 --genome GRCh38 --fasta /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genome_custom.fa --genes /n/scratch/users/c/cao385/refs/tdTomato/files/GRCh38-genes_custom.gtf \
    		 --genome mRatBN7-2 --fasta /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.modified --genes /n/scratch/users/c/cao385/refs/tdTomato/files/Rattus_norvegicus.mRatBN7.2.109.gtf.filtered \
    		 --nthreads 10
