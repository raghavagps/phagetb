#################################################################
# This module contains the model class and the functions to run #
# The model to predict phage which will lyse the given host     #
# based on genomes sequence of the bacteria                     #
#################################################################

import warnings
warnings.filterwarnings('ignore')
import pickle as pkl
import pandas as pd
import numpy as np
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import Bio
from Bio import SeqIO
from readFasta import readFasta, write_fasta
import argparse
from model import *

temp_path = '/usr1/webserver/cgidocs/tmp/phagetb/'
base_path = '/usr1/webserver/cgibin/phagetb'
DB_PREFIX =  base_path + '/blastdb/refdb/reference_interactions'
blastn_path = '/usr1/webserver/cgidocs/raghava/humcfs/FragileSequences/ncbi-blast-2.2.29+/bin/blastn'

def blastCallSingle(query_file, db_prefix, numThread = 4, e_value = 0.01):
    # if not os.path.exists(temp_path + "/temp"):
    #     os.makedirs(temp_path + "/temp")
    output_file = temp_path + os.path.basename(query_file) + '.txt'
    blastcall = NcbiblastnCommandline(cmd = blastn_path, query=query_file, db=db_prefix, out=output_file, outfmt="6", word_size = 11, evalue = e_value, reward = 1, penalty = -2, gapopen = 0, gapextend = 0, perc_identity = 90, num_threads=numThread)
    blastcall()
    return output_file

def get_ref_bacteria(dna_path, num_of_results = 1):
    output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_PREFIX, numThread=4, e_value=0.01)
    try:
        df = pd.read_table(output_file, header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        df_grouped = df.loc[df.groupby('sseqid')['evalue'].idxmin()]
        df_grouped = df_grouped.reset_index()
    except:
        return None
    os.remove(output_file)
    return df_grouped.iloc[:num_of_results, :]

def align_genome_sequences(seq1, seq2):
    '''
    Align two sequences using blastn
    '''
    ran = np.random.randint(100, 100000)
    outfile = temp_path + "/temp_" + str(ran) + "/temp_" + str(ran) + ".csv"
    blastn_cline = NcbiblastnCommandline(cmd = blastn_path, query=seq1, subject=seq2, outfmt="6", out= outfile)
    stdout, stderr = blastn_cline()
    return outfile

if __name__ == '__main__':
    ## Read Arguments from command
    parser = argparse.ArgumentParser(description='Please provide following arguments') 
    parser.add_argument("-i", "--input", type=str, required=True, help="Input: genome sequence of the bacteria in FASTA format or single sequence per line in single letter code")
    parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
    parser.add_argument("-l", "--levels", type=int, nargs='+', default=[1,2,3,4], help="Levels: 1: Blast against phage reference DB, 2: Blast against host reference DB, 3: Integrated model, 4: CRISPR")
    parser.add_argument("-n", "--num_of_ref_hosts", type=int, default=1, help="Number of reference hosts to consider")
    parser.add_argument("-t", "--threshold", type=float, default=0.01,help="Threshold: evalue threshold for similarity score")
    args = parser.parse_args()
    if args.output:
        outfile = args.output
    else:
        ran = np.random.randint(100, 100000)
        outfile = temp_path + '/output/outfile' + str(ran) + '.csv'
    try:
        check_levels(args.levels)
        levels = args.levels
    except Exception as e:
        print(e)
        print("Please provide valid levels")
        exit()

    input_file = args.input
    ## if input file contains multiple sequences then convert to one sequence per file
    records = readFasta(input_file)
    file_list = []
    for i in range(len(records)):
        write_fasta(records[i], input_file + "-" + str(i) + '.fasta')
        file_list.append(input_file + "-" + str(i) + '.fasta')
    result_dfs = []
    for input_file in file_list:
        threshold = args.threshold
        num_of_ref_hosts = args.num_of_ref_hosts
        blasthit_results = get_ref_bacteria(input_file, num_of_ref_hosts)
        if blasthit_results is None:
            df = pd.DataFrame()
            df['Predicted Phage'] = ["No Reference Phage Found"]
            df['Reference Host of the phage'] = ["No Reference Host Found"]
            result_dfs.append(df)
        else:
            ref_interaction_ids = blasthit_results['sseqid'].tolist()
            ref_interaction_ids = [x.split('_') for x in ref_interaction_ids]
            ## Each entry in ref_interaction_ids is a 3 tuple - (virus_id, host_id, source)
            ## source is either 'ref' or 'approved/actual'
            ref_virus_ids = [x[0] for x in ref_interaction_ids]
            ref_host_ids = ["_".join(x[1:-1]) for x in ref_interaction_ids]
            ref_ids = zip(ref_virus_ids, ref_host_ids)
            ref_ids = list(set(ref_ids))
            saved_results = None
            for ref_vir, ref_host in ref_ids:
                df = pd.DataFrame()
                df['Predicted Phage'] = [ref_vir]
                df['Reference Host of the phage'] = [ref_host]
                result_dfs.append(df)

    ## Combine all the results and remove the intermediate files
    for file_temp in file_list:
        os.remove(file_temp)
    result_df = pd.concat(result_dfs)
    result_df.to_csv(outfile, index=False)
