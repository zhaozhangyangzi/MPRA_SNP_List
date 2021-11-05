#!/bin/bash

ml legacy tbb python R

cd /oak/stanford/groups/khavari/users/yzhao8/David/MPRA/skin-disease

source ./script/ExtractSnpLD.sh

ExtractSnpLD /oak/stanford/groups/khavari/users/yzhao8/Jeet_DAPID/data/LD_AFR.tsv ./data/indexSNP.txt 1 ./result/LD_AFR.txt
ExtractSnpLD /oak/stanford/groups/khavari/users/yzhao8/Jeet_DAPID/data/LD_AMR.tsv ./data/indexSNP.txt 1 ./result/LD_AMR.txt
ExtractSnpLD /oak/stanford/groups/khavari/users/yzhao8/Jeet_DAPID/data/LD_ASN.tsv ./data/indexSNP.txt 1 ./result/LD_ASN.txt
ExtractSnpLD /oak/stanford/groups/khavari/users/yzhao8/Jeet_DAPID/data/LD_EUR.tsv ./data/indexSNP.txt 1 ./result/LD_EUR.txt

for f in ./result/LD_AFR.txt ./result/LD_AMR.txt ./result/LD_ASN.txt ./result/LD_EUR.txt; do sed -i "s/$/\t$(echo $f | cut -d'.' -f3|cut -d'_' -f2)/" $f; done

cat ./result/LD_AFR.txt ./result/LD_AMR.txt ./result/LD_ASN.txt ./result/LD_EUR.txt > ./result/LD_comb.txt

python ./script/process_LD_filter_3.py ./result/LD_comb.txt ./result/sig_LD_comb.txt

cd ./result

Rscript ../script/collapse.R

cut -f1 ./collapse-sig_LD_comb.txt | sort | uniq> ./LinkedSNP.txt

ExtractSnpLD /oak/stanford/groups/khavari/users/yzhao8/Jeet_DAPID/data/haploreg_v4.0_20151021_loc.vcf ./LinkedSNP.txt 3 ./LinkedSNP_location.vcf





