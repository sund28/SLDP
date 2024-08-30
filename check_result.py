import argparse
import h5py 
import pandas as pd 
import numpy as np 

# Create an argument parser
parser = argparse.ArgumentParser(description="Process gene name and SNP name.")

# Add arguments for gene name and SNP name
parser.add_argument('--gene_name', type=str, default="LRRC32", help="The name of the gene.")
parser.add_argument('--snp_name', type=str, default="rs1149596", help="The name of the SNP.")
parser.add_argument('--index', type=int, default=1, help="The index of track number, 0 - 5312.")
parser.add_argument('--show_matrix', type=int, default=0, help="show the matrix (896, 5313) or not, 0 is not show 1 is show.")


# Parse the arguments
args = parser.parse_args()


'''
read in 5313 track annotation, and print out the annotation of track assigned in parser. 
'''
path_annotation = './targets_human.txt'
track_anno = pd.read_csv(path_annotation, sep='\t') 
print(track_anno)
description = track_anno.loc[args.index, 'description']



# path_result = f'./result/{args.gene_name}_{args.snp_name}.h5'
path_result = f'./result_reference/{args.gene_name}_{args.snp_name}.h5' 

print(f"checking the result for gene {args.gene_name} and SNP {args.snp_name} at {path_result}") 

# Open the HDF5 file and print all keys
with h5py.File(path_result, 'r') as h5f:
	print("Keys in the HDF5 file:")
	for key in h5f.keys():
		dataset = h5f[key]
		# print(f"Key: {key}, Shape: {dataset.shape}")

with h5py.File(path_result, 'r') as h5f:
	sum_896bins = h5f['sum896_diff'][0]
	# print(sum_896bins) # debug 
	print("gene value computed by summing all 896 bins: ",sum_896bins[args.index]) 
	sum_gene_region = h5f['sum_gene_diff'][0] 
	print("gene value computed by summing gene region: ",sum_gene_region[args.index]) 
	

	ref_matrix = h5f['output_ref'][()]	 # (1, 896, 5313) 
	alt_matrix = h5f['output_alt'][()]	 # (1, 896, 5313) 
	diff_matrix = h5f['output_diff'][()] 	# (1, 896, 5313) 	 
	# print("shape of diff matrix: ",diff_matrix.shape) 

	summed_diff = np.sum(diff_matrix[0, 447:450, :], axis=0)  # TSS3, 447, 448, 449, output is 5313 track. 
	# print("summed diff: ",summed_diff, "shape: ",summed_diff.shape) # shape is 5313, 
	tss3_diff = summed_diff[args.index] # choose track from 5313. 
	print("gene value computed by sum of middle 3 bins of chosen track of difference between ref and alt: ",tss3_diff) 

	summed_diff = np.sum(diff_matrix[0, 443:453, :], axis=0)  # TSS10 
	tss10_diff = summed_diff[args.index] # choose track from 5313.
	print("gene value computed by sum of middle 10 bins of chosen track of difference between ref and alt diff: ",tss10_diff) 	 

	gene_coverage = h5f['gene_coverage'][()]
	print("percentage of gene region cover in 896 bins: ",gene_coverage) 

	if args.show_matrix == 1:
		print("ref matrix: ",ref_matrix) 
		print("alt matrix: ",alt_matrix) 
		print("diff matrix: ",diff_matrix) 
		


print("description of track: ",description) 