import pandas as pd
import numpy as np
import pickle
import os
import sys
import subprocess
import argparse 

from io import StringIO


## INPUT parameters
parser = argparse.ArgumentParser()
parser.add_argument("--gtf", help="input gtf files")
parser.add_argument("--output_path", help="outpath paths")

args = parser.parse_args()

gtf_file = args.gtf

# output                                       
gene_obj_file = args.output_path + '/gene_obj.pickle'
gtf_df_file   = args.output_path + '/gtf_df.pickle'
gene_info_file= args.output_path + '/gene_info.txt'


class Gene:
    def __init__(self, chr_id, start, end, strand,  gene_id, gene_name, gene_type, transcript_list, transcript_dict, gtf_idx):
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.transcript_list = transcript_list
        self.transcript_dict = transcript_dict
        self.gtf_idx = gtf_idx # a 2-element list, first is start row idx, second is the last row idx

class Transcript:
    def __init__(self, start, end, gene_id, transcript_id, transcript_name, transcript_type, exon_list, exon_dict, gtf_idx):
        self.start = start
        self.end = end
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.transcript_type = transcript_type
        self.exon_list = exon_list
        self.exon_dict = exon_dict
        self.gtf_idx = gtf_idx

class Exon:
    def __init__(self, start, end, gene_id, transctipt_id, exon_id):
        self.start = start
        self.end = end
        self.gene_id = gene_id
        self.transctipt_id = transctipt_id
        self.exon_id = exon_id

def parse_info(gtf):
    gene_idx = gtf[gtf[2] == "gene"].index.tolist()
    gene_num = len(gene_idx)
    # add a psuedo index as the very last row index, for slicing purpose
    gene_idx.append(gtf.shape[0]) 

    gene_list = [None]*gene_num
    gene_dict = {}


    for i in range(gene_num):
        gene_gtf_start = gene_idx[i]
        gene_gtf_end = gene_idx[i+1]-1

        gene_gtf = gtf.loc[gene_gtf_start:gene_gtf_end] # different to  gtf[a:b]
        gene_info= gene_gtf[gene_gtf[2] == "gene"]

        # for gene record, do sth
        #print(gene_info.iloc[:, [0,3,4,6,8]].values)
        chr_id, gene_start, gene_end, strand, gene_notes = gene_info.iloc[:, [0,3,4,6,8]].values[0]
        # print(gene_info)
        gene_notes_list = gene_notes.split()
        gene_id = gene_notes_list[1].split("\"")[1]
        gene_type = gene_notes_list[3].split("\"")[1]
        gene_name = gene_notes_list[5].split("\"")[1]
        gene_gtf_idx = [gene_gtf_start, gene_gtf_end]


        transcript_idx = gene_gtf[gene_gtf[2] == "transcript"].index.tolist()
        transcript_num = len(transcript_idx)
        transcript_idx.append(gene_gtf.index[-1] + 1)

        transcript_dict = {}
        transcript_list = []
        transcript_exon_list = []

        for j in range(transcript_num):
            transcript_gtf_start = transcript_idx[j]
            transcript_gtf_end = transcript_idx[j+1] -1

            transcript_gtf = gene_gtf.loc[transcript_gtf_start:transcript_gtf_end]
            transcript_info = transcript_gtf[transcript_gtf[2] == "transcript"]

            # for gene record, do sth
            transcript_start, transcript_end, transript_strand, transcript_notes = transcript_info.iloc[:, [3,4,6,8]].values[0]
            # print(transcript_info)
            transcript_notes_list = transcript_notes.split()
            # print(transcript_notes_list)
            transcript_id = transcript_notes_list[3].split("\"")[1]
            transcript_type = transcript_notes_list[9].split("\"")[1]
            transcript_name = transcript_notes_list[11].split("\"")[1]
            transcript_gtf_idx = [transcript_gtf_start, transcript_gtf_end]


            exon_gtf = transcript_gtf[transcript_gtf[2] == "exon"]
            exon_idx = exon_gtf.index.tolist()

            exon_dict = {}
            exon_list = []
            exon_boundary_list = []

            for k in exon_idx:
                exon_info = transcript_gtf.loc[k]
                exon_start = exon_info[3]
                exon_end = exon_info[4]
                exon_notes = exon_info[8]
                exon_notes_list = exon_notes.split()
                exon_id = exon_notes_list[15].split("\"")[1]

                exon_dict[exon_id] = Exon(exon_start, exon_end, gene_id, transcript_id, exon_id)
                exon_list.append(exon_id)
                exon_boundary_list.append([exon_start, exon_end])

            transcript_list.append(transcript_id)
            transcript_dict[transcript_id] = Transcript(transcript_start, transcript_end, gene_id, transcript_id, transcript_name, transcript_type, exon_list, exon_dict, transcript_gtf_idx)
            transcript_exon_list.append(exon_boundary_list)


        gene_list.append(gene_id)
        gene_dict[gene_id] = Gene(chr_id, gene_start, gene_end, strand, gene_id, gene_name, gene_type, transcript_list, transcript_dict, gene_gtf_idx)

    return gene_dict


gtf_all = pd.read_csv(gtf_file, header = None, comment='#', sep = "\t")
# consider genes on autosome chromosome, chrX, chrY and chrM
chr_list = [chr_id for chr_id in list(set(gtf_all[0]))]
gtf = gtf_all.loc[gtf_all[0].isin(chr_list)] 
gene_dict = parse_info(gtf)


gtf_gene_info = gtf[gtf[2]=="gene"]
gtf_gene_attr = gtf_gene_info[8].str.split(";| ",expand=True)[[1,4,7]].apply(lambda s:s.str.replace('"', ""))
gtf_gene_info[[9,10,11]] = gtf_gene_attr            # gene_id, gene_type, gene_name
gtf_gene_info[12] = gtf_gene_info.index.tolist()    # index in the gtf dataframe

# you can filter some genes based on your needs here
gtf_keep_gene_info = gtf_gene_info

# save gene information files
gtf_keep_gene_info2 = gtf_keep_gene_info[[0,1,2,3,4,6,9,10,11,12]]
gtf_keep_gene_info2.columns=["chr","source","feature","start","end","strand","gene_id","gene_type","gene_name","gtf_index"]
gtf_keep_gene_info2.to_csv(gene_info_file, sep='\t', header=False, index=False)


with open(gene_obj_file, 'wb') as gtf_handle:
    pickle.dump(gene_dict, gtf_handle, protocol=pickle.HIGHEST_PROTOCOL) 
    
with open(gtf_df_file, 'wb') as gtf_handle2:
    pickle.dump(gtf, gtf_handle2, protocol=pickle.HIGHEST_PROTOCOL) 

