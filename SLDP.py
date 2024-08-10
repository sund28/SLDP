'''
main script to do SLDP analysis, get ref sequence hg19 for target gene use genecode v12 as referecne to get TSS, and use SNP coordinate to get the ref and alt sequence. 
run Enformer for both seq, and save both result, and also try several ways to summary each track results to one value which is gene expression, may need to consider gene's exon region. 
'''

from pyfaidx import Fasta
import pyfaidx
import kipoiseq
from kipoiseq import Interval
import time 
import numpy as np
import pandas as pd 
import sys 
import os 
import argparse
import pickle 
import torch 
import h5py 
import pysam # to read vcf file. 

class FastaStringExtractor:
    
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()

def one_hot_encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

class Gene_SNP: # define a gene snps pair. 
    def __init__(self, gene_name, start, end, TSS, strand, chr, SNP_chr, SNP_pos, SNP_ref, SNP_alt, SNP_id):
        self.gene_name = gene_name
        self.start = start
        self.end = end
        self.TSS = TSS
        self.strand = strand
        self.chr = chr
        self.SNP_chr = SNP_chr
        self.SNP_pos = SNP_pos
        self.SNP_ref = SNP_ref
        self.SNP_alt = SNP_alt
        self.SNP_id = SNP_id 

gene_snp1 = Gene_SNP(  # rs2155219
    gene_name="LRRC32",
    start=76368568,
    end=76381791,
    TSS=76381791,
    strand="-",
    chr="11",
    SNP_chr="11",
    SNP_pos=76299194,
    SNP_ref="G", 
    SNP_alt="T",
    SNP_id="rs2155219" 
)

gene_snp2 = Gene_SNP(  # rs73127651
    gene_name="PPHLN1",
    start=42632249,
    end=42853517,
    TSS=42632249,
    strand="+",
    chr="12",
    SNP_chr="12",
    SNP_pos=42376118,
    SNP_ref="A",
    SNP_alt="G",
    SNP_id="rs73127651"
)

gene_snp3 = Gene_SNP(  # rs2075620
    gene_name="APOE",
    start=45408956,
    end=45412650,
    TSS=45408956,
    strand="+",
    chr="19",
    SNP_chr='19',
    SNP_pos=45480037,
    SNP_ref="A",
    SNP_alt="G", 
    SNP_id="rs2075620"   
)

gene_snp4 = Gene_SNP(  # rs2075620
    gene_name="MARK4",
    start=45754550,
    end=45808541,
    TSS=45754550,
    strand="+",
    chr="19",
    SNP_chr="19",
    SNP_pos=45480037,
    SNP_ref="A",
    SNP_alt="G", 
    SNP_id="rs2075620" 
)

gene_snp5 = Gene_SNP(  # rs2075620
    gene_name="CLASRP",
    start=45542298,
    end=45574214,
    TSS=45542298,
    strand="+",
    chr="19",
    SNP_chr="19",
    SNP_pos=45480037,
    SNP_ref="A",
    SNP_alt="G",
    SNP_id="rs2075620" 
)

gene_snp6 = Gene_SNP( # rs6690215
    gene_name="CD55",
    start=207494853,
    end=207534311,
    TSS=207494853,
    strand="+",
    chr="1",
    SNP_chr="1",
    SNP_pos=207656050,
    SNP_ref="C",
    SNP_alt="T",
    SNP_id="rs6690215" 
)

gene_snp7 = Gene_SNP(  # rs6690215
    gene_name="CR1",
    start=207669492,
    end=207813992,
    TSS=207669492,
    strand="+",
    chr="1",
    SNP_chr="1",
    SNP_pos=207656050,
    SNP_ref="C",
    SNP_alt="T",
    SNP_id="rs6690215" 
)

gene_snp8 = Gene_SNP(  # rs6690215
    gene_name="YOD1",
    start=207217194, # start of gene 
    end=207226325,
    TSS=207226325,
    strand="-",
    chr="1",
    SNP_chr="1",
    SNP_pos=207656050,
    SNP_ref="C",
    SNP_alt="T",
    SNP_id="rs6690215" 
)


# Setup
consensus_dir = '/projects/zhanglab/users/dongbo/personalized-expression-benchmark/data/ref_fasta' # TODO path to consensus sequences (reference seq)
genes_file = '/projects/zhanglab/users/dongbo/gene_3000.csv'

SEQUENCE_LENGTH = 196_608 # 
# SEQUENCE_LENGTH = 10 # test now, TSS-5, TSS+5, only one side is included, 5th [4] is TSS then. 



'''
input gene position and correpond SNP position here. 
'''
gene_snp = gene_snp8 # choose object to use, only 1,3,5,7 within enformer region now. 

gene_name = gene_snp.gene_name 
TSS = gene_snp.TSS 
ref = gene_snp.SNP_ref 
alt = gene_snp.SNP_alt 
chr = gene_snp.SNP_chr 
SNP_pos = gene_snp.SNP_pos 
SNP_id = gene_snp.SNP_id 
assert gene_snp.chr == gene_snp.SNP_chr, "gene and SNP should be in the same chromosome" 

print("TSS: ",TSS, " ref: ",ref, " alt: ",alt) 
start = int(TSS) - SEQUENCE_LENGTH // 2
end = int(TSS) + SEQUENCE_LENGTH // 2 
print("start of enformer input: ",start, " end of enformer input: ",end, 'chr: ',chr) 
ref_gene_file = os.path.join(consensus_dir, f'chr{chr}.fa')
fasta_extractor = FastaStringExtractor(ref_gene_file) 
target_interval = kipoiseq.Interval(f'chr{chr}', start, end)
# sequence_letter = fasta_extractor.extract(target_interval.resize(SEQUENCE_LENGTH)) # 10 now,TSS-5, TSS+4, only one side is included. 
sequence_letter = fasta_extractor.extract(target_interval) # 196608. mind rechange seq length, maybe make coordinate wrong. try make sure for more genes. 
TSS_letter_Index = SEQUENCE_LENGTH // 2 -1 # use to index TSS in extracted seq .
snp_inds = SNP_pos - TSS + TSS_letter_Index # index of SNP in extracted seq. make sure >=0 and < SEQUENCE_LENGTH. 
if snp_inds < 0 or snp_inds >= SEQUENCE_LENGTH:
    raise ValueError(f"Error: SNP index {snp_inds} out of bounds! It should be >= 0 and < {SEQUENCE_LENGTH}.")


print("gene length is: ",gene_snp.start - gene_snp.end) # print gene length. 
print("index of SNP reletive to TSS: ",SNP_pos - TSS) # print index of TSS reletive to TSS. 
print("distance of SNP to start and end of gene: ",SNP_pos - gene_snp.start, SNP_pos - gene_snp.end) # print distance of SNP to start and end of extracted seq. 
print("SNP index in extracted seq is: ",snp_inds) # print index of SNP in extracted seq. 
print("extracted SNP ref letter is:",sequence_letter[snp_inds]) # print SNP ref letter extracted from refseq, sholud be in line with it in SNP file. 
print("sequence_letter around SNP: ",sequence_letter[snp_inds-2:snp_inds+2]) # print around SNP. 


'''
load enformer and predict for ref and alt of SNP, and also save differene matrix and sum of each track, all in one h5 file for 3 keys. 
'''
refseq_letter = sequence_letter  # 196608 
altseq_letter = sequence_letter[:snp_inds] + alt + sequence_letter[snp_inds+1:] # replace SNP with ref letter. 
print("replaced SNP letter is:",altseq_letter[snp_inds]) # print SNP ref letter extracted from refseq, sholud be in line with it in SNP file. 
print("sequence_letter around SNP after replace: ",altseq_letter[snp_inds-2:snp_inds+2]) # print around SNP. 
print("length of refseq_letter: ",len(refseq_letter), " length of altseq_letter: ",len(altseq_letter)) # print length of ref and alt seq. 


import subprocess
import time
import torch
from enformer_pytorch import Enformer, seq_indices_to_one_hot
import subprocess 
from torch.cuda.amp import autocast as autocast # mix precision, after torch 1.6 just use it in cuda. 
from enformer_pytorch import from_pretrained

def print_gpu_usage():
    result = subprocess.run(['nvidia-smi', '-q', '-d', 'MEMORY'], stdout=subprocess.PIPE)
    print(result.stdout.decode())

print("loading enformer torch model...")
start_time = time.time()
enformer = from_pretrained('EleutherAI/enformer-official-rough', use_tf_gamma = False).cuda()
end_time = time.time()
print("time for loading model: ", end_time - start_time) 

print("finished loading model")

print("enformer: ",enformer)


# Check if CUDA is in use
if torch.cuda.is_available():
    print("CUDA is in use.")
else:
    print("CUDA is not in use.")

model = enformer # cuda 
model.eval() # eval mode. 

input_seq_ref = one_hot_encode(refseq_letter) # 196608, 4 
input_seq_alt = one_hot_encode(altseq_letter) # 196608, 4 
input_seq_ref = torch.tensor(input_seq_ref).unsqueeze(0).cuda() # 1, 196608, 4 (1 is batch size )
input_seq_alt = torch.tensor(input_seq_alt).unsqueeze(0).cuda() # 1, 196608, 4 

with torch.no_grad():
    output_ref = model(input_seq_ref)['human'] # (1, 896, 5313), human head 
    output_alt = model(input_seq_alt)['human'] # (1, 896, 5313)
# print("output_ref: ",output_ref.shape) # 

output_diff = output_alt - output_ref # (1, 896, 5313)

# Sum along the second dimension (axis=1)
sum896_diff = torch.sum(output_diff, dim=1) # sum all 896 bins 



'''
read in gene region annotation and only sum gene region. both gene and exon use gencode v12. 
'''
path_v12_geneanno = '/projects/zhanglab/users/dongbo/Enformer_retrain/gencode_v12/gene_annotation_v12.csv' # contain 50000 genes. 
# Read the CSV file
df = pd.read_csv(path_v12_geneanno)
# Filter out rows where the 'gene_name' column is equal to 'gene_name'
filtered_df = df[df['gene_name'] == gene_name]
row = filtered_df.iloc[0] # only one row. 
# print("filtered gene anno: ",row) # only one row. 
gene_start = row['start']    
gene_end = row['end'] 
gene_chr = row['seqname']  # like chr1 
assert gene_chr == 'chr'+str(chr), "gene and SNP should be in the same chromosome" 
bin_start = (gene_start - start) // 128 # 
bin_end = (gene_end - start) // 128 # index of bin start and end. 
length_bins = bin_end - bin_start + 1 # length of gene region in bins. use to compute after clip how many percent left. 
print("bin start and end of gene: ",bin_start, bin_end) # print bin start and end. 0 - 895 
if bin_start < 0: 
    bin_start = 0 
if bin_end > 895: 
    bin_end = 895 
percentage_cover = (bin_end - bin_start)/length_bins # percentage of gene region cover in 896 bins. 
print("percentage of gene region cover in 896 bins: ",percentage_cover) # print percentage of gene region cover in 896 bins. 


'''
TODO only sum exon region as well. use the finer annotation later after do base gene annotation version. 
'''

'''
save result to h5 file 
'''
output_path = f'/projects/zhanglab/users/dongbo/SLDP_project/result/{gene_name}_{SNP_id}.h5'

# Check if the file exists and delete it if it does
if os.path.exists(output_path):
    os.remove(output_path)
    print(f"Existing file '{output_path}' deleted.")

output_ref_np = output_ref.cpu().numpy()
output_alt_np = output_alt.cpu().numpy()
output_diff_np = output_diff.cpu().numpy()
sum896_diff_np = sum896_diff.cpu().numpy()

gene_region_diff_np = output_diff_np[:, bin_start:bin_end]
sum_gene_diff_np = gene_region_diff_np.sum(axis=1) #  only sum gene region.  (1, 127, 5313)
print("shape of sum region_diff_np: ",gene_region_diff_np.shape) # print shape of sumregion_diff_np. 
# The shape of vector_5313 will be (1, 5313)


print(sum896_diff.shape)  # Should output torch.Size([1, 5313])
print("sum of all bins:", sum896_diff)  # Should output a tensor of shape (1, 5313) # sum of all 896 bins 
print("sum of gene region:", sum_gene_diff_np)  # Should output a tensor of shape (1, 5313) # sum of gene region. 


with h5py.File(output_path, 'w') as h5f:
    h5f.create_dataset('output_ref', data=output_ref_np)
    h5f.create_dataset('output_alt', data=output_alt_np)
    h5f.create_dataset('output_diff', data=output_diff_np)
    h5f.create_dataset('sum896_diff', data=sum896_diff_np)
    h5f.create_dataset('sum_gene_diff', data=sum_gene_diff_np) 
    h5f.create_dataset('gene_coverage', data=percentage_cover) 

print(f"new Data saved to {output_path}")