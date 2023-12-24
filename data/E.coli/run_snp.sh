outdir='/data/HWK/DeepGene/data/E.coli/outputs'
ref='/data/HWK/DeepGene/data/E.coli/REF_genomic.fna.gz'
g1='/data/HWK/DeepGene/data/E.coli/ERR434268_1.fastq.gz'
g2='/data/HWK/DeepGene/data/E.coli/ERR434268_2.fastq.gz'
# !snippy --cpus 16 --outdir /data/HWK/DeepGene/data/E.coli/outputs --ref /data/HWK/DeepGene/data/E.coli/REF_GCF_904425475.1_MG1655_genomic.fna --R1 /data/HWK/DeepGene/data/E.coli/ERR434268_1.fastq.gz --R2 /data/HWK/DeepGene/data/E.coli/ERR434268_2.fastq.gz
# !snippy --cpus 16 --outdir /data/HWK/DeepGene/data/E.coli/outputs2 --ref /data/HWK/DeepGene/data/E.coli/REF_GCF_904425475.1_MG1655_genomic.fna --ctgs /data/HWK/DeepGene/data/E.coli/contig.fa
!snippy-multi input.tab --cpus 32 --ref /data/HWK/DeepGene/data/E.coli/REF_GCF_904425475.1_MG1655_genomic.fna > RunSNP.sh