import pandas as pd
import re
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import gc

#Global Variables
Lib_Name = 'TL4S1_czb_NEW' #UPDATE with library name or base name you want for your files
Fig_Format = 'jpeg' #UPDATE default graph file format
Threshold = 100 #UPDATE with the min number of reads you want for tiles not in designed file  

Output_Directory = f'{Lib_Name}_Maps_and_Graphs'
os.makedirs(Output_Directory, exist_ok=True)

#create summary tabble
summary_dict = {'Category': [], 'Read Count': []}

#functions to process the input files
def find_designed(des):
    """Creates a lookup dictionary of all designed tiles from a file."""
    dt = []
    with open(des, 'r') as f_des:
        for line in f_des:
            # Remove the left primer
            left_trimmed = line.replace("CCCAGCTTAAGCCACCATG", "") #UPDATE with the left primer seq that is the same for all seq in the design file you want to remove
            
            # Remove everything after (and including) the right sequence
            right_trimmed = left_trimmed.split("GGATCCGAGCTCG")[0] #UPDATE with the right seq you want to remove and everything after it 
            
            dt.append(right_trimmed.strip())
    return {tile: 1 for tile in dt}

def getmid(seq, pre, post):
    """Extracts the sequence between pre and post substrings."""
    match = re.search(f"{pre}(.*){post}", seq)
    return match.group(1) if match else "X"

def tilebc_mapper(readfile, dtd, t_len=6, bc1_len=9, rtbc_len=16, designed_len=162, #UPDATE with correct barcode lengths and flanking sequences
                  tile_pre="CTCGAGATAACTTCGTATAATGTATGCTAT", tile_post="GGCCGGCCATAGGGCCCC",
                  bc1_pre="GAGCTCGCTAGC", bc1_post="CTCGAGATAA",
                  rtbc_pre="GGCCGGCCATAGGGCCCC", rtbc_post="GCGGTCCA",
                  designed_pre="CACCATG", designed_post="GGATCCG"): #CHANGED POST TO JUST BE AD LENGTH
    """Processes input sequences to map tiles, HawkBCs, RTBCs, and Designed sequences."""

    # Lists to store extracted data
    tile_list, tile_lengths, tq_list = [], [], []
    bc1_list, bc1_lengths, bc1q_list = [], [], []
    rtbc_list, rtbc_lengths, rtbcq_list = [], [], []
    designed_list, designed_lengths, designedq_list, in_designlist = [], [], [], []
    sequences = []
    total_sequences = 0

    with open(readfile, 'r') as fin:
        for line in fin:
            if line.startswith('@'):
                seq = next(fin).strip()
                sequences.append(seq)
                total_sequences += 1

                # Extract Tile
                tile = getmid(seq, tile_pre, tile_post)
                tile_len = len(tile)
                tile_quality = 1 if tile_len == t_len else 0

                # Extract HawkBC
                adBC = getmid(seq, bc1_pre, bc1_post)
                adBC_len = len(adBC)
                adBC_quality = 1 if adBC_len == bc1_len else 0

                # Extract RTBC
                rtbc = getmid(seq, rtbc_pre, rtbc_post)
                rtbc_len_actual = len(rtbc)
                rtbc_quality = 1 if rtbc_len_actual == rtbc_len else 0

                # Extract Designed
                designed = getmid(seq, designed_pre, designed_post)
                designed_len_actual = len(designed)
                designed_quality = 1 if designed_len_actual == designed_len else 0
                in_design = 1 if designed in dtd else 0

                # Store all values
                tile_list.append(tile)
                tile_lengths.append(tile_len)
                tq_list.append(tile_quality)

                bc1_list.append(adBC)
                bc1_lengths.append(adBC_len)
                bc1q_list.append(adBC_quality)

                rtbc_list.append(rtbc)
                rtbc_lengths.append(rtbc_len_actual)
                rtbcq_list.append(rtbc_quality)

                designed_list.append(designed)
                designed_lengths.append(designed_len_actual)
                designedq_list.append(designed_quality)
                in_designlist.append(in_design)

    # Create DataFrame
    tileBC_df = pd.DataFrame({
        "Reads": sequences,
        "ADBC2": tile_list,
        "ADBC2 Len": tile_lengths,
        "ADBC2 Qual": tq_list,
        "HawkBCs": bc1_list,
        "HawkBC Len": bc1_lengths,
        "HawkBC Qual": bc1q_list,
        "RTBC": rtbc_list,
        "RTBC Len": rtbc_lengths,
        "RTBC Qual": rtbcq_list,
        "Designed": designed_list,
        "Designed Len": designed_lengths,
        "Designed Qual": designedq_list,
        "In Des File" : in_designlist
    })

    return tileBC_df

def process_maps(input_file, design_file):
    designed_tile_dict = find_designed(design_file)
    map1 = tilebc_mapper(input_file, designed_tile_dict)
    return map1, designed_tile_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine and optimize scripts.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input paired seq file path')
    parser.add_argument('-d', '--design', type=str, required=True, help ='Input design file')
    args = parser.parse_args()
    map1, designed_tile_dict = process_maps(args.input, args.design)

#export the Map1 LUT
map1.to_csv(os.path.join(Output_Directory, f'{Lib_Name}_map1.csv'), index=False)

###Adding different quality counts to summary file
count_rows_a = ((map1['HawkBC Qual'] == 1)).sum()
summary_dict['Category'].append(f'Reads with Correct HawkBC Length')
summary_dict['Read Count'].append(count_rows_a)

count_rows_t = ((map1['ADBC2 Qual'] == 1)).sum()
summary_dict['Category'].append(f'Reads with Correct ADBC2 Length')
summary_dict['Read Count'].append(count_rows_t)

count_rows_r = ((map1['RTBC Qual'] == 1)).sum()
summary_dict['Category'].append(f'Reads with Correct RTBC Length')
summary_dict['Read Count'].append(count_rows_r)

count_rows_att = ((map1['Designed Qual'] == 1)).sum()
summary_dict['Category'].append(f'Reads with Correct Tile Length')
summary_dict['Read Count'].append(count_rows_att)

count_rows_attf = ((map1['In Des File'] == 1)).sum()
summary_dict['Category'].append(f'Reads in Design File')
summary_dict['Read Count'].append(count_rows_attf)

count_rows_one = len(map1[(map1['ADBC2 Qual'] == 1) & (map1['HawkBC Qual'] == 1)])
summary_dict['Category'].append(f'Rows with correct ADBC2 length and HawkBC length')
summary_dict['Read Count'].append(count_rows_one)

count_rows_2 = len(map1[(map1['ADBC2 Qual'] == 1) & (map1['RTBC Qual'] == 1)])
summary_dict['Category'].append(f'Rows with correct ADBC2 length and RTBC length')
summary_dict['Read Count'].append(count_rows_2)

count_rows_3 = len(map1[(map1['HawkBC Qual'] == 1) & (map1['RTBC Qual'] == 1)])
summary_dict['Category'].append(f'Rows with correct HwkBC length and RTBC length')
summary_dict['Read Count'].append(count_rows_3)

count_rows_4 = len(map1[(map1['HawkBC Qual'] == 1) & (map1['RTBC Qual'] == 1) & (map1['ADBC2 Qual'] == 1)])
summary_dict['Category'].append(f'Rows with correct HwkBC length, RTBC length, and AD BC2 length')
summary_dict['Read Count'].append(count_rows_4)

count_rows_5 = len(map1[(map1['HawkBC Qual'] == 1) & (map1['RTBC Qual'] == 1) & (map1['ADBC2 Qual'] == 1)& (map1['Designed Qual'] == 1)])
summary_dict['Category'].append('Rows with correct HwkBC length, RTBC length, AD BC2 length and Tile Length')
summary_dict['Read Count'].append(count_rows_5)

count_rows_59 = len(map1[(map1['HawkBC Qual'] == 1) & (map1['RTBC Qual'] == 1) & (map1['ADBC2 Qual'] == 1) & (map1['Designed Qual'] == 1) & (map1['In Des File'] == 1)])
summary_dict['Category'].append('Rows with correct HwkBC length, RTBC length, AD BC2 length, Tile length, AND in designed')
summary_dict['Read Count'].append(count_rows_59)

###Create Histograms of quality and length of barcodes and Tiles#plot tile length histogram
def map1_graphs(map1, Lib_Name, Fig_Format):
    plt.hist(map1['HawkBC Len'], bins=100)
    plt.xlim([0, 12]) #UPDATE if your bc1 length won't fit in this range
    plt.title(f'{Lib_Name } HawkBC Length Frequency')
    plt.xlabel('HawkBC Length')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_HawkBC_length.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['HawkBC Qual'])
    plt.title(f'{Lib_Name } HawkBC Quality Frequency')
    plt.xlabel('HawkBC Qual')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_HawkBC_Quality.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['ADBC2 Len'])
    plt.xlim([0, 10]) #UPDATE if your tile length won't fit in this range
    plt.title(f'{Lib_Name } ADBC2 Length Frequency')
    plt.xlabel('ADBC2 Length')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_ADBC2_Length.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['ADBC2 Qual'])
    plt.title(f'{Lib_Name } ADBC2 Quality Frequency')
    plt.xlabel('ADBC2 Qual')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_ADBC2_Quality.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['RTBC Len'], bins=100)
    plt.xlim([0, 20]) #UPDATE if your bc1 length won't fit in this range
    plt.title(f'{Lib_Name } RTBC Length Frequency')
    plt.xlabel('RTBC Length')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_RTBC_length.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['RTBC Qual'])
    plt.title(f'{Lib_Name } RTBC Quality Frequency')
    plt.xlabel('RTBC Qual')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_RTBC_Quality.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['Designed Len'], bins=100)
    plt.xlim([0, 250]) #UPDATE if your bc1 length won't fit in this range
    plt.title(f'{Lib_Name } Tile Length Frequency')
    plt.xlabel('Tile Length')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_Tile_length.{Fig_Format}'))
    plt.clf()

    plt.hist(map1['Designed Qual'])
    plt.title(f'{Lib_Name } Tile length Quality Frequency')
    plt.xlabel('Tile length Qual')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_Tile_length_Quality.{Fig_Format}'))
    plt.clf()

map1_graphs(map1, Lib_Name, Fig_Format) #UPDATE run if you want to make qulaity and length graphs for map1
def not_in_des_file_extract(map1,thresh):
    not_in_design = map1[(map1['Designed Qual'] == 1) & (map1['HawkBC Qual'] == 1) & (map1['RTBC Qual'] == 1) & (map1['ADBC2 Qual'] == 1) & (map1['In Des File'] == 0)]

    nd_labels = ['Reads', 'ADBC2 Len','ADBC2 Qual', 'HawkBC Len','HawkBC Qual', 'RTBC Qual', 'RTBC Len', 'Designed Len', 'Designed Qual', 'In Des File']
    not_in_design = not_in_design.drop(nd_labels, axis = 1)
    print(f'reads in not in design but all other quality 1 : {not_in_design.shape[0]}')

    lut_not_in_design = not_in_design.copy()
    lut_not_in_design['Cat'] = lut_not_in_design['HawkBCs'].str.cat([lut_not_in_design['ADBC2'], lut_not_in_design['RTBC'], lut_not_in_design['Designed']], sep='-')
    lut_not_in_design['HA'] = lut_not_in_design['HawkBCs'].str.cat([lut_not_in_design['ADBC2']], sep='-')
    lut_not_in_design['HAR'] = lut_not_in_design['HawkBCs'].str.cat([lut_not_in_design['ADBC2'], lut_not_in_design['RTBC']], sep='-')

    # filter to only keep reads that appear >= Threshold times 
    lut_not_in_designcat_counts = lut_not_in_design['Cat'].value_counts()
    valid_catz = lut_not_in_designcat_counts[lut_not_in_designcat_counts >= f'{thresh}'].index
    nd_g5 = lut_not_in_design[lut_not_in_design['Cat'].isin(valid_catz)].copy()
    print(f'reads in not in design but all other quality 1 >= {thresh} reads: {nd_g5.shape[0]}')

    #create table of only unique that still keeps track of the origional read count 
    cat_countz = nd_g5['Cat'].value_counts()
    nd_final = nd_g5.drop_duplicates(subset='Cat').copy()
    nd_final['Cat_Counts'] = nd_final['Cat'].map(cat_countz)
    print(f'Unique Cat in not in designed >= {thresh} reads: {nd_final.shape[0]}')

    nd_final.to_csv(os.path.join(Output_Directory,f'{Lib_Name}_not_in_design_unique_Cat_min_ge{thresh}_reads.csv'), index=False)

###Run this part below if you want to look at reads where all BC quality = 1 and Tile length is correct but the Tile is not in the design file 
not_in_des_file_extract(map1,f'{Threshold}') #UPDATE uncomment if you want to make a seperate file that contains reads with all quality = 1 but not in the designed file


###removing reads that do not have quality = 1 
#Replace all 0s in map1 with NaN to filter out any Qual=0 reads
map1_nans = map1.replace(0, np.nan)
map2 = map1_nans.dropna().reset_index()

#get rid of some now useless columns
clabels = ['index','Reads', 'ADBC2 Len','ADBC2 Qual', 'HawkBC Len','HawkBC Qual', 'RTBC Qual', 'RTBC Len', 'Designed Len', 'Designed Qual', 'In Des File']
map2 = map2.drop(clabels, axis = 1)



summary_dict['Category'].append(f'New Section')
summary_dict['Read Count'].append('Quality = 0 rows removed')

map3 = map2.copy()
del map1, map2 #doing to save memory 
gc.collect()

###Add column that connects barcodes to the Tile it is paired with (Cat for concatenation) 
# Create the Cat column by concatenating HawkBCs, ADBC2, and RTBC
# Create the Cat column by concatenating HawkBCs, ADBC2, RTBC, and Designed
map3['Cat'] = map3['HawkBCs'].str.cat([map3['ADBC2'], map3['RTBC'], map3['Designed']], sep='-') #UPDATE if you want different combinations or do not like the concatenation column name
map3['HA'] = map3['HawkBCs'].str.cat([map3['ADBC2']], sep='-')
map3['HAR'] = map3['HawkBCs'].str.cat([map3['ADBC2'], map3['RTBC']], sep='-')

summary_dict['Category'].append('Map3 Shape')
summary_dict['Read Count'].append(map3.shape[0])

### export Map3 to csv 
#make csv of map3
map3.to_csv(os.path.join(Output_Directory, f'{Lib_Name}_map3_all_good_qual_reads.csv'), index=False)

###Add counts to  summary table 
#Hawkings BC
abcov = map3['HawkBCs'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique HawkBC')
summary_dict['Read Count'].append(abcov.shape[0])

#ADBC2
tcov = map3['ADBC2'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique ADBC2 coverage')
summary_dict['Read Count'].append(tcov.shape[0])

#RTBC
Rbcov = map3['RTBC'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique RTBC')
summary_dict['Read Count'].append(Rbcov.shape[0])

#Tile
tcovArd3t = map3['Designed'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Map3 Unique Tiles ')
summary_dict['Read Count'].append(tcovArd3t.shape[0])

#HA
tcovA = map3['HA'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique Hawk+ADBC (HA) combos coverage')
summary_dict['Read Count'].append(tcovA.shape[0])

#HAR
tcovAr = map3['HAR'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique Hawk+ADBC+ RTBC (HAR) combos coverage')
summary_dict['Read Count'].append(tcovAr.shape[0])

#CAT
tbcov = map3['Cat'].value_counts().to_frame().reset_index()
summary_dict['Category'].append(f'Unique Cat')
summary_dict['Read Count'].append(tbcov.shape[0])

###Create a smaller file that is  the look up table of unique Cat that keeps track of the number of times each unique Cat occured in the origional file
cat_countss = map3['Cat'].value_counts()
map4 = map3.drop_duplicates(subset='Cat').copy()
map4['Cat_Count'] = map4['Cat'].map(cat_countss)
map4 = map4.sort_values(by='Cat_Count', ascending=False)
map4.to_csv(os.path.join(Output_Directory, f'{Lib_Name}_map4_unique_cat.csv'), index=False)

#export summary table
summary_dict_df = pd.DataFrame.from_dict(summary_dict)
summary_dict_df.to_csv(os.path.join(Output_Directory, f'{Lib_Name}_Initial_Summary.csv'), index=False)

def map3_graphs(abcov,tcov,Rbcov,tbcov,):
    plt.hist(abcov['HawkBCs'], bins=100)
    plt.title(f' Unique HawkBC Read Coverage Frequency')
    plt.xlabel('Coverage')
    plt.ylabel('Counts')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_HawkBC_cov.{Fig_Format}'))
    plt.clf()

    plt.hist(tcov['ADBC2'], bins=100)
    plt.title(f'Unique ADBC2 Coverage Frequency')
    plt.xlabel('Coverage')
    plt.ylabel('Counts')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_ADBC2_cov.{Fig_Format}'))
    plt.clf() 

    plt.hist(tcov['ADBC2'], bins=100)
    plt.title(f' Unique ADBC2 Coverage Frequency LOG')
    plt.xlabel('Coverage')
    plt.ylabel('Counts')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_ADBC2_log_cov.{Fig_Format}'))
    plt.clf()

    plt.hist(Rbcov['RTBC'], bins=100)
    plt.title(f' Unique RTBC Read Coverage Frequency')
    plt.xlabel('Coverage')
    plt.ylabel('Counts')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_RTBC_cov.{Fig_Format}'))
    plt.clf()

    plt.hist(Rbcov['RTBC'], bins=100)
    plt.title(f' Unique RTBC Read Coverage Frequency LOG')
    plt.xlabel('Coverage')
    plt.ylabel('Counts')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(os.path.join(Output_Directory, f'{Lib_Name}_RTBC_cov_log.{Fig_Format}'))
    plt.clf()

map3_graphs(abcov,tcov,Rbcov,tbcov,) #UPDATE if you want to make map3 graphs of read coverage for unique barcodes