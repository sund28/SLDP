# this folder is to do SDLM, and change for 4 SNPs for correspond genes first, then if work do it on more SNPs. 
Use TSS(transcript start site) provide by gencode v12 as center of input sequence, and pass it to Enformer, use reference to compute output_ref and then change the SNP assigned to get output_alt, and compute output_diff = output_alt - output_ref, then sum each track of 896 bins to get summary value for each track, 5313 track in total. 

Huang NG23 and Gauvadis both use gencode v12, and huang NG use hg19 as reference seq. (but original Enformer training use gene code v32)
onek1k dataset on cluster is based on hg19 as well, like those vcf and bed/bam/fam SNP files. 

input seq length is 196608 bp, output is 114688 bp. so for all 8 SNP-gene provide, only 4 of SNP are within input region of Enformer, and are computed. and for some gene, output region of enformer may only cover part of gene region like 70%. 





## description for all files under folder: 

1. SLDP.py: main script used to compute difference matrix between ref and alt sequence, results saved in h5 file.  

2. format of results: each h5 file is named as gene_name_SNP_name.h5, and contain reference sequence result output_ref, alt sequence output_alt, difference sequence output_diff which is output_alt-output_ref, all of those three are (1,896,5313) matrix. sum896_diff, sum of all 896 bins, shape is (1, 5313) correspond to 5313 tracks. sum_gene_diff only sum bins correspond to gene region, shape is (1, 5313). 

3. check_result.py: input gene name, SNP name, and index of track of 5313 tracks want to check, it will print correspond track info and computed value by enformer. 

4. targets_human.txt: 5313 track correspondance. 

### 5313 tarck corerspond can found here 
https://github.com/calico/basenji/blob/master/manuscripts/cross2020/targets_human.txt

