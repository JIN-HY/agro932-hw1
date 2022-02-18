
ml bcftools

### Extract allele frequency at each position
bcftools query -f '%CHROM %POS %AF\n' Chr09_2000001-3000000.vcf.gz > frq.txt
bcftools query -f '%CHROM %POS %REF %ALT [\t%GT]\n' Chr09_2000001-3000000.vcf.gz > geno.txt
