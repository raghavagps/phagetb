#################################################################
# This module contains the model class and the functions to run #
# The model predict the host for a phage using the genome       #
# sequence of the phage only                                    #
#################################################################

import warnings
warnings.filterwarnings('ignore')
import pickle as pkl
import pandas as pd
import numpy as np
import os
from sklearn.linear_model import LogisticRegression
from Bio.Blast.Applications import NcbiblastnCommandline
import Bio
from Bio import SeqIO
from getKmerProfile import getKmerProfile
from readFasta import readFasta, write_fasta
import argparse

temp_path = '/usr1/webserver/cgidocs/tmp/phagetb'
base_path = '/usr1/webserver/cgibin/phagetb'
blastn_path = '/usr1/webserver/cgidocs/raghava/humcfs/FragileSequences/ncbi-blast-2.2.29+/bin/blastn'
## Blast DB paths ##
DB_PHAGE_PREFIX = base_path + '/blastdb/train-blastdb/db'
DB_HOST_PREFIX = base_path + '/blastdb/host-blastdb/db'
DB_CRISPR_PREFIX = base_path + '/blastdb/crispr-blastdb/allCRISPRs'

## Helper functions ##
def save_pickle_file(file_path, data):
    with open(file_path, 'wb') as f:
        pkl.dump(data, f)

def read_pickle_file(file_path):
    with open(file_path, 'rb') as f:
        return pkl.load(f)

def blastCallSingle(query_file, db_prefix, numThread = 4, e_value = 0.01):
    if not os.path.exists(temp_path + "/temp"):
        os.makedirs(temp_path + "/temp")
    output_file = temp_path + "/temp/" + os.path.basename(query_file) + '.hit'
    blastcall = NcbiblastnCommandline(cmd = blastn_path, query=query_file, db=db_prefix, out=output_file, outfmt="6", word_size = 11, evalue = e_value, reward = 1, penalty = -2, gapopen = 0, gapextend = 0, perc_identity = 90, num_threads=numThread)
    blastcall()
    return output_file

def crisprSingle(query_file, db_prefix, numThread = 4):
    if not os.path.exists(temp_path + "/temp"):
        os.makedirs(temp_path + "/temp")
    output_file = temp_path + "/temp/" + os.path.basename(query_file) + '.crispr'
    crispr_call = NcbiblastnCommandline(cmd = blastn_path, query=query_file, db=db_prefix, out=output_file, outfmt="6", evalue=1, gapopen=10, penalty=-1, gapextend=2, word_size=7, dust='no', task='blastn-short', perc_identity=90, num_threads=numThread)
    crispr_call()
    return output_file

class HierarchicalModel:
    def __init__(self, levels = [1, 2, 3, 4]):
        self.clf = read_pickle_file(base_path + '/pretrained/model.pkl')
        self.e1 = 1e-75
        self.e2 = 1e-16
        self.alpha = 0.9
        self.gamma = 0.5
        self.threshold = 0.6
        self.taxa = read_pickle_file(base_path + '/extras/hostTaxa.pkl')
        self.taxa.columns = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Name']
        self.taxa_list = self.taxa.columns[:-1]
        self.acc_to_id = read_pickle_file(base_path + '/extras/acc_to_id.pkl')
        self.hosts = read_pickle_file(base_path + '/extras/hosts_list.pkl')
        self.virus = read_pickle_file(base_path + '/extras/ref_virus_list.pkl')
        self.host = read_pickle_file(base_path + '/extras/ref_host_list.pkl')
        self.host2index = {self.hosts[i]: i for i in range(len(self.hosts))}
        self.virus2index = {self.virus[i]: i for i in range(len(self.virus))}
        self.host_labels = np.unique(self.host)
        self.levels = levels
                
    def calculate_feats_for_LR(self, dna):
        return getKmerProfile(dna=dna)

    def read_dna_from_file(self, file_path):
        seq_req = readFasta(fasta_file=file_path)
        return seq_req[0][2]

    def get_similarity_score(self, dna_path):
        output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_HOST_PREFIX, numThread=4, e_value=10)
        try:
            df = pd.read_table(output_file, header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df['sseqid'] = df['sseqid'].apply(lambda x: self.acc_to_id[x])
            df_grouped = df.loc[df.groupby('sseqid')['bitscore'].idxmax()]
            df_grouped.bitscore = (df_grouped.bitscore - df_grouped.bitscore.min())/(df_grouped.bitscore.max() - df_grouped.bitscore.min())
            df_grouped = df_grouped.reset_index()
            results = np.zeros(len(self.host_labels))
            for i in range(len(df_grouped)):
                target_host_index = self.host2index[df_grouped.sseqid[i]]
                results[target_host_index] = df_grouped.bitscore[i]
        except:
            results = np.zeros(len(self.host_labels))
        os.remove(output_file)
        return results

    def get_vir_vir_sim(self, dna_path):
        output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_PHAGE_PREFIX, numThread=4, e_value=10)
        try:
            df = pd.read_table(output_file, header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df_grouped = df.loc[df.groupby('sseqid')['bitscore'].idxmax()]
            df_grouped.bitscore = (df_grouped.bitscore - df_grouped.bitscore.min())/(df_grouped.bitscore.max() - df_grouped.bitscore.min())
            df_grouped = df_grouped.reset_index()
            results = np.zeros(len(self.virus))
            for i in range(len(df_grouped)):
                target_virus_index = self.virus2index[df_grouped.sseqid[i].split('.')[0]]
                results[target_virus_index] = df_grouped.bitscore[i]
        except:
            results = np.zeros(len(self.virus))
        os.remove(output_file)
        return results

    def get_host_host_sim(self, dna_path):
        output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_HOST_PREFIX, numThread=4, e_value=10)
        try:
            df = pd.read_table(output_file, header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df['sseqid'] = df['sseqid'].apply(lambda x: self.acc_to_id[x])
            df_grouped = df.loc[df.groupby('sseqid')['bitscore'].idxmax()]
            df_grouped.bitscore = (df_grouped.bitscore - df_grouped.bitscore.min())/(df_grouped.bitscore.max() - df_grouped.bitscore.min())
            df_grouped = df_grouped.reset_index()
            results = np.zeros(len(self.host_labels))
            for i in range(len(df_grouped)):
                target_host_index = self.host2index[df_grouped.sseqid[i]]
                results[target_host_index] = df_grouped.bitscore[i]
        except:
            results = np.zeros(len(self.host_labels))
        os.remove(output_file)
        return results

    def getBlastH(self, dna_path):
        output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_HOST_PREFIX, numThread=4, e_value=0.01)
        try:
            df = pd.read_table(output_file, header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df['sseqid'] = df['sseqid'].apply(lambda x: self.acc_to_id[x])
            df_grouped = df.loc[df.groupby('sseqid')['evalue'].idxmin()]
            df_grouped = df_grouped.reset_index()
            results = np.ones(len(self.host_labels))*1e9
            for i in range(len(df_grouped)):
                target_host_index = self.host2index[df_grouped.sseqid[i]]
                results[target_host_index] = df_grouped.evalue[i]
        except:
            results = np.ones(len(self.host_labels))*1e9
        os.remove(output_file)
        return results

    def getBlastP(self, dna_path):
        output_file = blastCallSingle(query_file=dna_path, db_prefix=DB_PHAGE_PREFIX, numThread=4, e_value=0.01)
        try:
            df = pd.read_table(output_file, header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            df_grouped = df.loc[df.groupby('sseqid')['evalue'].idxmin()]
            df_grouped = df_grouped.reset_index()
            results = np.ones(len(self.virus))*1e9
            for i in range(len(df_grouped)):
                target_virus_index = self.virus2index[df_grouped.sseqid[i].split('.')[0]]
                results[target_virus_index] = df_grouped.evalue[i]
        except:
            results = np.ones(len(self.virus))*1e9
        os.remove(output_file)
        return results

    def getCRISPR(self, dna_path):
        output_file = crisprSingle(query_file=dna_path, db_prefix=DB_CRISPR_PREFIX, numThread=4)
        try:
            df = pd.read_table(output_file, header = None)
            hosts =  df[1].apply(lambda x: x.split('|')[1]).tolist()
        except:
            hosts = []
        os.remove(output_file)
        return hosts

    def predict(self, file_path):
        predicted_host = None
        if 1 in self.levels:
            blastP_results = self.getBlastP(file_path)
        else:
            blastP_results = []
        if 2 in self.levels:
            blastH_results = self.getBlastH(file_path)
        else:
            blastH_results = []
        if 4 in self.levels:
            CRISPR_results = self.getCRISPR(file_path)
        else:
            CRISPR_results = []
        dna = self.read_dna_from_file(file_path)
        feats = self.calculate_feats_for_LR(dna)
        y_pred_probs = self.clf.predict_proba(np.reshape(feats, (1, -1)))
        similarity_scores = self.get_similarity_score(file_path)
        vir_vir_sim_most_sim = self.get_vir_vir_sim(file_path)
        host_host_sim = self.get_host_host_sim(file_path)
        predicted_host = None
        ## Predictions using Blast against phage reference DB
        if 1 in self.levels and blastP_results[np.argmin(blastP_results)] < self.e2:
            predicted_host = self.host[np.argmin(blastP_results)]
        else:
            ## Predictions using Blast against host reference DB
            if 2 in self.levels and blastH_results[np.argmin(blastH_results)] < self.e1:
                predicted_host = np.argmin(blastH_results)
            else:
                ## Predictions using integrated model
                prob_of_pred = y_pred_probs
                x = similarity_scores
                most_sim_virus = np.argmax(vir_vir_sim_most_sim)
                host_of_most_sim_virus = self.host[most_sim_virus]
                sim_with_ref_host = host_host_sim[host_of_most_sim_virus]
                x = prob_of_pred * (1 - self.gamma) + (x * (1-self.alpha) + sim_with_ref_host * self.alpha) * self.gamma
                ## Predictions from integrated model if confidence is high enough else predictions from Blast against CRISPR reference DB
                if len(CRISPR_results) > 0:
                    if 3 in self.levels and np.max(x) > self.threshold:
                        predicted_host = np.argmax(x)
                    elif 4 in self.levels:
                        predicted_host = self.host2index[CRISPR_results[0]]
                elif 3 in self.levels:
                    predicted_host = np.argmax(x)
        ## Host (ID) wrt to reference host database
        if predicted_host is not None:
            results = {'Host': predicted_host}
            for t in self.taxa_list:
                results[t] = self.taxa.loc[self.hosts[predicted_host]][t]
            return results
        else:
            return None

def check_levels(levels):
    if not isinstance(levels, list):
        raise TypeError('levels must be a list of integers')
    for i in levels:
        if not isinstance(i, int):
            raise TypeError('levels must be a list of integers')
        if i < 1 or i > 4:
            raise ValueError('levels must be integers between 1 and 4')

if __name__ == '__main__':
    ## Read Arguments from command
    parser = argparse.ArgumentParser(description='Please provide following arguments') 
    parser.add_argument("-i", "--input", type=str, required=True, help="Input: genome sequence of the phage in FASTA format or single sequence per line in single letter code")
    parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
    parser.add_argument("-l", "--levels", type=int, nargs='+', default=[1,2,3,4], help="Levels: 1: Blast against phage reference DB, 2: Blast against host reference DB, 3: Integrated model, 4: CRISPR")
    args = parser.parse_args()

    if args.input:
        input_file = args.input
    else:
        print("Please provide input file")
        exit()

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

    model = HierarchicalModel(args.levels)
    input_file = args.input
    ## if input file contains multiple sequences then convert to one sequence per file
    records = readFasta(input_file)
    file_list = []
    for i in range(len(records)):
        write_fasta(records[i], input_file + "-" + str(i) + '.fasta')
        file_list.append(input_file + "-" + str(i) + '.fasta')
    result_dfs = []

    for input_file in file_list:
        res = model.predict(input_file)
        if res is not None:
            df = pd.DataFrame(res, index=[0])
            result_dfs.append(df)
            # df.to_csv(outfile, index=False)
        else:
            df = pd.DataFrame()
            df['Host'] = ["No Prediction Found"]
            df['Domain'] = ['No Prediction Found']
            df['Phylum'] = ['No Prediction Found']
            df['Class'] = ['No Prediction Found']
            df['Order'] = ['No Prediction Found']
            df['Family'] = ['No Prediction Found']
            df['Genus'] = ['No Prediction Found']
            df['Species'] = ['No Prediction Found']
            result_dfs.append(df)
            # df.to_csv(outfile, index=False)
            print("No prediction found")

    ## Combine all the results and remove the intermediate files
    for file_temp in file_list:
        os.remove(file_temp)
    result_df = pd.concat(result_dfs)
    result_df.to_csv(outfile, index=False)


