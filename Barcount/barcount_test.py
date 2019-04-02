import numpy as np
import pandas as pd
from StringIO import StringIO
from datetime import datetime
import os
from subprocess import check_call

print 'Welcome to the test script for barcount. This script will run barcount on a simulated fastq files using a range of parameters and check that the output matches the expected. Please run this after installation to ensure bacrount is running correctly. Before running, please ensure that barcount is added to your PATH.'

directory = datetime.now().strftime("%m%d%Y_%H-%M-%S_barcountTest")
os.mkdir(directory)
os.chdir(directory)
print 'Created folder for testing'

bases = np.array(['A','G','T', 'C'])

#Simulated reads are inspired by Catalina's read structure
def generate_random_seq(length):
    return ''.join(np.random.choice(bases, length, replace=True))

def generate_random_sequences(length, n):
    return [generate_random_seq(length) for i in range(n)]

# We will use 3 barcodes
np.random.seed(17)
barcodes = generate_random_sequences(20,3)
#bc0: 'CGCTTGAGTCACAACGTGGA'
#bc1: 'TTCCTCTCGGTGTGGGAGGT'
#bc2: 'TTATGGACTCTGAGACCTCC'


#Write those into a database file
with open('test_barcode_database.csv', 'w') as bdb:
    for i, b in enumerate(barcodes):
        bdb.write('gene%i,%s\n'%(i,b))

#Simulate a fastq file
test_sequences = pd.DataFrame(columns=['barcode_index', 'barcode', 'barcode_mismatches', 'umiA', 'umiB', 'left_flank','left_flank_mismatches', 'right_flank', 'right_flank_mismatches'])
index_counter = 0

#Add 5 correct reads for each barcode
for b in barcodes:
    for c in range(5):
        test_sequences.loc[index_counter] = [index_counter, b, 0, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
        index_counter += 1

#For bc0, add 1 read with 1 mismatch in the barcode. Use string replace with the count option
test_sequences.loc[index_counter] = [index_counter, barcodes[0].replace('G', 'A', 1), 1, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1
#For bc0, add 1 read with 2 mismatch in the barcode
test_sequences.loc[index_counter] = [index_counter, barcodes[0].replace('G', 'A', 2), 2, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1
#For bc0, add 1 read with 3 mismatch in the barcode
test_sequences.loc[index_counter] = [index_counter, barcodes[0].replace('G', 'A', 3), 3, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1
#For bc0, add 1 read with 4 mismatch in the barcode
test_sequences.loc[index_counter] = [index_counter, barcodes[0].replace('G', 'A', 4), 4, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1

#For bc0, add one read with a short barcode of length 19. Should still be picked up except when min read length filter is not set
test_sequences.loc[index_counter] = [index_counter, barcodes[0].replace('G', '', 1), 1, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1

#For bc1, add 1 read with 1 mismatch in left flank and 1 in right flank
test_sequences.loc[index_counter] = [index_counter, barcodes[1], 0, generate_random_seq(4), generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATGTC', 1, 'TTTAAATGCGAAGTAAGGTGGGAGCG', 1]
index_counter += 1

#For bc2, add 1 read with same umiA
test_sequences.loc[index_counter] = [index_counter, barcodes[2], 0, 'GCGC', generate_random_seq(4), 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1

#For bc2, add 1 read with same umiA and same umiB
test_sequences.loc[index_counter] = [index_counter, barcodes[2], 0, 'GCGC', 'CGTA', 'GGGGACGAGGCAAGCTAAGATATC', 0, 'TTTAAATGCGAAGTAAGGCGGGAGCG', 0]
index_counter += 1

test_sequences['simulated_read'] = 'AGTA' + test_sequences['umiA'] + test_sequences['left_flank'] + test_sequences['barcode'] + test_sequences['right_flank'] + test_sequences['umiB'] + 'TGAC'
test_sequences.to_csv('test_simulated_read_data.csv')

def write_fastq_entry(idcounter, read, handle):
    handle.write('@testseq-%i\n'%idcounter)
    handle.write(read + '\n')
    handle.write('+\n')
    handle.write(''.join(['I' for i in range(len(read))]) + '\n')
    
fo = open('test_sequences.fastq', 'w')
for i,s in test_sequences['simulated_read'].iteritems():
    write_fastq_entry(i, s, fo)
fo.close()

print 'Simulated read file and barcount database'


#Run barcount with the following parameters
parameters = {'strict':{'min_read_length':86, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':0, 'max_distance_barcode':0, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},
             '1bc_mismatch':{'min_read_length':86, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':0, 'max_distance_barcode':1, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},
             '2bc_mismatch':{'min_read_length':86, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':0, 'max_distance_barcode':2, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},
             '3bc_mismatch':{'min_read_length':86, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':0, 'max_distance_barcode':3, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},
             '85bpallowed':{'min_read_length':85, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':0, 'max_distance_barcode':0, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},
             '1flankmismatch':{'min_read_length':86, 'max_read_length':86, 'flanking_left': 'GGGGACGAGGCAAGCTAAGATATC', 'flanking_right' : 'TTTAAATGCGAAGTAAGGCGGGAGCG', 'max_distance_flanks':1, 'max_distance_barcode':0, 'umiA_position':'" 4:8"', 'umiB_position':'" -8:-4"'},

             }
for pi,pval in parameters.items():
    check_call('barcount --barcode_table test_barcode_database.csv --debug --verbose --save_extracted_barcodes --fastq test_sequences.fastq --out test_parameters_%s %s'%(pi, ' '.join(['--%s %s'%(k,v) for k,v in pval.items()])), shell=True)
#Ran barcount using 6 different parameter sets


#Read in the results count file
results = pd.Series()
for i in parameters:
    results[i] = pd.read_csv('test_parameters_%s.csv'%i, index_col=0)
    results[i] = results[i].drop('gene', axis=1)
results = pd.concat(results.to_dict(), axis=1, sort=True)
results.to_csv('test_result.csv')
print 'Read in barcount result files'

result_string = ''',1bc_mismatch,1bc_mismatch,1bc_mismatch,1bc_mismatch,1flankmismatch,1flankmismatch,1flankmismatch,1flankmismatch,2bc_mismatch,2bc_mismatch,2bc_mismatch,2bc_mismatch,3bc_mismatch,3bc_mismatch,3bc_mismatch,3bc_mismatch,85bpallowed,85bpallowed,85bpallowed,85bpallowed,strict,strict,strict,strict
,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths,count_postFilters,failed_sameUMI,mean_matching_distance,mean_barcode_lengths
CGCTTGAGTCACAACGTGGA,6,0,0.16666666666666666,20.0,5,0,0.0,20.0,7,0,0.42857142857142855,20.0,8,0,0.75,20.0,5,0,0.0,20.0,5,0,0.0,20.0
TTATGGACTCTGAGACCTCC,6,1,0.0,20.0,6,1,0.0,20.0,6,1,0.0,20.0,6,1,0.0,20.0,6,1,0.0,20.0,6,1,0.0,20.0
TTCCTCTCGGTGTGGGAGGT,5,0,0.0,20.0,6,0,0.0,20.0,5,0,0.0,20.0,5,0,0.0,20.0,5,0,0.0,20.0,5,0,0.0,20.0
unassigned,0,0,,,0,0,,,0,0,,,0,0,,,0,0,,,0,0,,'''
expected = pd.read_csv(StringIO(result_string), header=[0,1], index_col=0)
assert np.isclose(results, expected, equal_nan=True).all(), 'Test failed. Please raise an issue on github.'

print 'Results match expected. Congratulations, barcount appears to be running as expected (on test simulated data at least)'


