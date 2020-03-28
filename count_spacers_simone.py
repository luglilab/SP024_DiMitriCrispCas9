from Bio import SeqIO
import gzip
import csv
import pandas as pd
import re
import argparse
import numpy as np
## input
path = "/Users/simone/Desktop/Dimitri"
inputfile = "a.csv"

parser = argparse.ArgumentParser(
    description='Analyze sequencing data for sgRNA library distribution')
parser.add_argument('-f', '--fastq', type=str, dest='fastq_file',
                    help='fastq file name', default='NGS.fastq')
parser.add_argument('-o', '--output', type=str, dest='output_file',
                    help='output file name', default='library_count.csv')
parser.add_argument('-i', '--input', type=str, dest='input_file',
                    help='input file name', default='library_sequences.csv')
args = parser.parse_args()

# occurrence variable store the number of reads that contain the key substring
occurrence = 0
# notvalidsequence variable store the number of reads that does not contain the key substring
notvalidsequence = 0
# key used to exctract the guide (Attenction the guanine is present and is specific for Gecko Add gene kit)
key = "CGAAACACCG"
# variable that store the number of guide that match with library_sequences
perfect_matches = 0
# variable that store the number of guide that partially match with library_sequences
non_perfect_matches = 0
# numb of total reads
num_read = 0
# import file with library_sequences
df = pd.read_csv(args.input_file, header=0)
# create dataframe
df.insert(3, 'Count', 0)
# decompress fastq.gz
with gzip.open(args.fastq_file,"rt") as handle:
    # extract sequence id and sequence nucleotides
    for record in SeqIO.parse(handle, "fastq"):
        num_read +=1
        print(num_read)
        # add to seqid and seqnt list the sequence id and nucleotide sequence
        if key in record.seq:
            # add 1 to occurrence count because the reads is valid and could be analyzed
            occurrence += 1
            # store the position the the key has the first match on the reads
            if (df['seq'] == str(record.seq)[int(str(record.seq).find(key))+len(key):int(str(record.seq).find(key))+len(key)+20]).any():
                # update total counting
                perfect_matches += 1
                # update the complete-match guide count
                df.loc[df['seq'].str.contains(str(record.seq)[int(str(record.seq).find(key))+len(key):int(str(record.seq).find(key))+len(key)+20]), 'Count'] = df['Count'] + 1
            else:
                # update the partially-match guide count
                non_perfect_matches += 1
        else:
            # count the number of reads discarded
            notvalidsequence += 1

# sort using Count column
df.sort_values(by=['Count'], ascending=False,inplace=True)
#export the dataframe with count
df.to_excel(".".join([args.output_file,"xlsx"]), index = False,header=True)

# percentage of guides that matched perfectly
percent_matched = round(perfect_matches / float(perfect_matches + non_perfect_matches) * 100, 1)
# percentage of undetected guides with no read counts
guides_with_reads = np.count_nonzero(df['Count'])
guides_no_reads = len(df['Count']) - guides_with_reads
percent_no_reads = round(guides_no_reads / float(len(df['Count'])) * 100, 1)
# skew ratio of top 10% to bottom 10% of guide counts
top_10 = np.percentile(df['Count'], 90)
bottom_10 = np.percentile(df['Count'], 10)
if top_10 != 0 and bottom_10 != 0:
    skew_ratio = top_10 / bottom_10
else:
    skew_ratio = 'Not enough perfect matches to determine skew ratio'

# write analysis statistics to statistics.txt
with open("".join([args.output_file,'statistics.txt']), 'w') as infile:
    infile.write('Number of perfect guide matches: ' + str(perfect_matches) + '\n')
    infile.write('Number of nonperfect guide matches: ' + str(non_perfect_matches) + '\n')
    infile.write('Number of reads where key was not found: ' + str(notvalidsequence) + '\n')
    infile.write('Number of reads processed: ' + str(num_read) + '\n')
    infile.write('Percentage of guides that matched perfectly: ' + str((perfect_matches*100)/occurrence) + '\n')
    infile.write('Percentage of undetected guides: ' + str((non_perfect_matches*100)/occurrence) + '\n')
    infile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio))
    infile.close()



