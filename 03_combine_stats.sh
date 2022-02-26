# add header for heterozygosity and fst stats files
grep 'pop1' Scaffold_10__10000001__10100000__stats.txt > ../window_heterozygosity.txt
grep 'pop1' Scaffold_10__10000001__10100000__stats.txt > ../window_fst.txt

# pull out only the relevant stat for each file
for i in $( ls ); do grep 'heterozygosity' $i >> ../window_heterozygosity.txt; done
for i in $( ls ); do grep 'Fst' $i >> ../window_fst.txt; done
