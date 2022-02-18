
ml bcftools

### Extract allele frequency at each position
bcftools query -f '%CHROM %POS %AF\n' Chr09_2000001-3000000.vcf.gz > largedata/frq.txt
bcftools query -f '%CHROM %POS %REF %ALT [\t%GT]\n' Chr09_2000001-3000000.vcf.gz > largedata/geno.txt
