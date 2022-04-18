# PhageTB

![Pics](./extras/PhageTB.png)

PhageTB is a multilevel prediction method for predicting interactions between bacteriophages and pathogenic bacterial hosts. This study develops a novel host prediction method for predicting hosts of query phages by their genome sequences utilizing alignment-based and alignment-free features.

## Installation
```
git clone https://github.com/raghavagps/phagetb.git
cd phagetb/
pip install -r requirements.txt
python model.py -i input.fasta -o output.csv -l 1 2 3 4
```

## Getting Started

This project is hosted on [PhageTB](https://webs.iiitd.edu.in/raghava/phagetb/) and can be accessed by clicking on the link above or can be used as a standalone application by downloading the source code from this GitHub repository.

There are 3 prediction methods available in this model. 
1. Predict The Bacterial Host For A Query Phage (model.py)
This Module Allows Users To Predict The Bacterial Hosts Corresponding To The Query Phages using the genome sequence of the phage. 
1. Predict Interaction Of Query Phage-Bacteria Pair (model_phage_host_pair.py)
This Module Allows Users To Predict Whether Given Phage And Bacterial Hosts Are Likely To Interact With One Another. The prediction from this module for the phage is used as a query for the BLAST task (blastn) against the query bacterial host. The BLAST task is performed using the NCBI BLAST+ tool. The BLAST output is parsed and if the predicted host and query host have a similarity higher than the threshold, then the phage-host pair is predicted to interact.
1. Predict The Lytic Phage For Query Bacteria (model_bacteria.py)
This Module Allow Users To Predict The Target Phage Likely To Infect Query Bacteria.

## Usage
Following is the complete list of all options, (with default values) that can be used to run the model. you may get these options by "python model.py -h" (and similarly for other modules).

You need to change the following paths in the respective files:
```
temp_path: where the temporary files are stored
base_path: where the source files are stored
blastn_path: where the blastn executable is stored
```

#### model.py
```
usage: model.py [-h] -i INPUT [-o OUTPUT] [-l LEVELS [LEVELS ...]]

Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: genome sequence of the phage in FASTA format or
                        single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -l LEVELS [LEVELS ...], --levels LEVELS [LEVELS ...]
                        Levels: 1: Blast against phage reference DB, 2: Blast
                        against host reference DB, 3: Integrated model, 4:
                        CRISPR
```

#### model_phage_host_pair.py
```
usage: model_phage_host_pair.py [-h] -v INPUT_PHAGE -b INTPUT_BACTERIA
                                [-o OUTPUT] [-l LEVELS [LEVELS ...]]
                                [-t THRESHOLD]

Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -v INPUT_PHAGE, --input_phage INPUT_PHAGE
                        Input: genome sequence of the phage in FASTA format or
                        single sequence per line in single letter code
  -b INTPUT_BACTERIA, --intput_bacteria INTPUT_BACTERIA
                        Input: genome sequence of the bacteria in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -l LEVELS [LEVELS ...], --levels LEVELS [LEVELS ...]
                        Levels: 1: Blast against phage reference DB, 2: Blast
                        against host reference DB, 3: Integrated model, 4:
                        CRISPR
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: evalue threshold for similarity score
```

#### model_bacteria.py
```
usage: model_bacteria.py [-h] -i INPUT [-o OUTPUT] [-l LEVELS [LEVELS ...]]
                         [-n NUM_OF_REF_HOSTS] [-t THRESHOLD] [--only_blast]

Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: genome sequence of the bacteria in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -l LEVELS [LEVELS ...], --levels LEVELS [LEVELS ...]
                        Levels: 1: Blast against phage reference DB, 2: Blast
                        against host reference DB, 3: Integrated model, 4:
                        CRISPR
  -n NUM_OF_REF_HOSTS, --num_of_ref_hosts NUM_OF_REF_HOSTS
                        Number of reference hosts to consider
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: evalue threshold for similarity score
```

## File descriptions
```
blastdb: The database of the hosts, phages and CRISPR sequences.
extras: The directory containing the extra files required for predictions.
pretrained: The directory containing the pretrained models.
genome_data: The directory containing the genome data of reference hosts
input.fasta: The input file containing the query sequences.
getKmerProfiles.py: The python script for generating kmer profiles.
readFasta.py: The python script for reading fasta files.
model.py: The python script for predicting hosts.
model_bacteria.py: The python script for predicting target phages for a bacteria.
model_phage_host_pair.py: The python script for predicting interaction for a phage and host pair.
``` 

## Address for contact
In case of any query, feel free to reach out to us at 
```
Prof. G. P. S. Raghava, Head Department of Computational Biology,            
Indraprastha Institute of Information Technology (IIIT), 
Okhla Phase III, New Delhi 110020 ; Phone:+91-11-26907444; 
Email: raghava@iiitd.ac.in  Web: http://webs.iiitd.edu.in/raghava/
```
