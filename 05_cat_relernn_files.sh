grep "#" Scaffold_1.recode.vcf > nuthatch_relernn.vcf

for i in $(seq 39); do grep -v "#" Scaffold_${i}.recode.vcf >> nuthatch_relernn.vcf; done
