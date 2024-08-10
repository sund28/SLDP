import argparse
import h5py 
import pandas as pd 

# Create an argument parser
parser = argparse.ArgumentParser(description="Process gene name and SNP name.")

# Add arguments for gene name and SNP name
parser.add_argument('--gene_name', type=str, default="LRRC32", help="The name of the gene.")
parser.add_argument('--snp_name', type=str, default="rs2155219", help="The name of the SNP.")
parser.add_argument('--index', type=int, default=5312, help="The index of track number, 0 - 5312.")


# Parse the arguments
args = parser.parse_args()


'''
read in 5313 track annotation, and print out the annotation of track assigned in parser. 
'''
path_annotation = './targets_human.txt'
track_anno = pd.read_csv(path_annotation, sep='\t') 
print(track_anno)
description = track_anno.loc[args.index, 'description']






path_result = f'./result/{args.gene_name}_{args.snp_name}.h5'
print(f"checking the result for gene {args.gene_name} and SNP {args.snp_name} at {path_result}") 

# Open the HDF5 file and print all keys
with h5py.File(path_result, 'r') as h5f:
	print("Keys in the HDF5 file:")
	for key in h5f.keys():
		dataset = h5f[key]
		print(f"Key: {key}, Shape: {dataset.shape}")

with h5py.File(path_result, 'r') as h5f:
	sum_896bins = h5f['sum896_diff'][0]
	# print(sum_896bins) # debug 
	print("value computed by summing all 896 bins: ",sum_896bins[args.index]) 
	sum_gene_region = h5f['sum_gene_diff'][0] 
	print("value computed by summing gene region: ",sum_gene_region[args.index]) 
	gene_coverage = h5f['gene_coverage_rate'][()]
	print("percentage of gene region cover in 896 bins: ",gene_coverage) 

print("description of track: ",description) 