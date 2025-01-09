#! /usr/bin/env python3

import os
import sys
import re
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

### Filter seqeunces containing intact recorder seq (non-edited and edited and with and without 8N insertion) and create lists of these sequences) 

total = []
reads_rec = []
reads_bad = []
with open ("WH85_ZF1412_3_CKDL240013352-1A_H5NVGDSXC_L1.txt", "r") as recorder_seq:
    for line in recorder_seq:
        line = line.rstrip()
        total.append(line)
        found = re.search(r"TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TTAAATT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G](?=C)", line)
        if found:
            reads_rec.append(line)
        else:
            reads_bad.append(line)
all_reads = len(total)
good = len(reads_rec)
bad = len(reads_bad)
#print(reads_empty)


# Filter all reads having the same 5' ends and 8mer insertion



seqs_all = []
motifs = []
recorders = []
total_edits = []

for sequence in reads_rec:
    if 'N' in sequence:
        continue  # Skip this iteration if "N" is present

    found = re.search(r"^GCCAACGCTATTCTGGCTGA(.{7})ATCTATCCAACGCAATTCTGGCACG", sequence)
    if found:
        seqs_all.append(sequence)
        motif = found.group(1)  # Extract the matched motif
        motifs.append(motif)
    else:
        foundE = re.search(r"^GCCAACGCTATTCTGGCTGAATTCTATCCAACGCAATTCTGGCACG", sequence)
        if foundE:
            seqs_all.append(sequence)
            motifE = "EMPTY"  # Assign "EMPTY" as the motif
            motifs.append(motifE)






for sequence in seqs_all:
    found1 = re.search("TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G](?=C)", sequence)
    if found1:
        recorder = found1.group()
        recorders.append(recorder)

        start = found1.start()
        end = found1.end()
        region = sequence[start:end]

        # Count occurrences of "TTGGA" in the extracted region
        count1 = region.count("TTGGA")
        print(count1)
        total_edits.append(count1)
#print(total_edits)

### Count all edited reads. Create a list with either 1 or 0 values. 1 is assigned if more then 1 edit is present in a read

molecules = []
for x in total_edits:
    if x >= 1:
        molecules.append(1)
    else:
        molecules.append(0)

### Count abundance of each motif and make a dictinary 
#print(len(seqs_all))
#print(len(motifs))
#print(len(total_edits))
#print(len(molecules))




occurrence = {}
for item in motifs:
    if item in occurrence:
        occurrence[item] +=1
    else:
        occurrence[item] =1
#occurrence_sorted = sorted(occurrence) 
#print(occurence)



### Create a dictonary of all motifs and correspond reads that have at least one edit in one molecule 

dict_0 = {}
s = iter(motifs)
p = iter(molecules)
m = list(zip(s,p))
for (s,p) in m:
    if s in dict_0:
        dict_0[s] = dict_0[s]+p
    else:
        dict_0[s] = p
#print(dict_0)

		
### Build a dictonary of all motifs and total number of editing events for each ot them 

dict_1 = {}
i = iter(motifs)
j = iter(total_edits)
k = list(zip(i,j))
for (x,y) in k:
    if x in dict_1:
        dict_1[x] = dict_1[x]+y
    else:
        dict_1[x] = y
#print(dict_1)


### Building combined dictonary:

ds = [occurrence, dict_0, dict_1]
d = {}
for k in occurrence.keys():
    d[k] = tuple(d[k] for d in ds)
#print(d)




### Convert dictonary in pands dataframe

data = pd.DataFrame.from_dict(d, orient='index')
data.columns = ['occurrence','reads_edited','total_edits']
number_of_good_reads = len(seqs_all)
data['all_reads'] = number_of_good_reads

data['expression_percent']=(data['occurrence']/data['all_reads'])*100

data['efficiency_sq'] =(data['reads_edited']/data['occurrence'])*(data['total_edits']/data['occurrence'])
#data['efficiency1_new_norm_all'] =((data['reads_edited']/data['occurrence'])*(data['total_edits']/data['all_reads']))*1000000
data['percent_edited']=(data['total_edits']/(data['occurrence']*6))*100


data_filtered = data[data['occurrence'] >= 10]



# Sort the DataFrame by Z-score in descending order
final = data_filtered.sort_values(by=['percent_edited'], ascending=False)
final = final.rename_axis("Motif").reset_index()
final = final[(final['Motif'].str.len() >2) | (final['Motif'] == "EMPTY")]
print(final.to_string())

final.to_csv("ZF_7N_1412_3_6A.csv")

# Calculate mean efficiency and standard deviation based on filtered data
#mean_efficiency_filtered = np.mean(data_filtered['efficiency'])
#std_dev_efficiency_filtered = np.std(data_filtered['efficiency'])
#print(mean_efficiency_filtered)
#print(std_dev_efficiency_filtered)
#data_filtered = data_filtered.copy()
#data_filtered.loc[:, 'z_score'] = (data_filtered['efficiency'] - mean_efficiency_filtered) / std_dev_efficiency_filtered

# Perform Min_Max_normalization
#min_val = min(data_filtered['efficiency'])
#max_val = max(data_filtered['efficiency'])
#print(min_val)
#print(max_val)
#data_filtered['min_max'] = (data_filtered['efficiency'] - min_val) / (max_val - min_val)




#final = final.dropna(subset=['Motif', 'efficiency'])
final = final.sort_values(by='occurrence', ascending=False)



#sns.set(style="whitegrid")
#plt.figure(figsize=(10, 6))
filtered_final = final[final['Motif'] != 'EMPTY']
# Create a density plot
sns.kdeplot(data=filtered_final['occurrence'], fill=True)

# Add labels and title
plt.title('Density Plot of Reads per motif')
plt.xlabel('Motifs')
plt.ylabel('Reads')


#plt.savefig('ZF_1534_1.png', format='png')
#plt.show()




#sns.heatmap(final)
#plt.show()
#ax = final.plot.bar(x="Motif", y = "z_score", rot =0,color='blue')
#ax.set_xticklabels([])
#plt.savefig('1534_7N_1331_1_z-score.png', format='png')
#plt.show()
# Create the bar plot with efficiency values (blue)
#plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
#sns.barplot(x="Motif", y="occurrence", data=final, color='blue', alpha=0.5)
#sns.barplot(x="Motif", y="efficiency", data=final, color='red')


#plt.xticks(rotation=0)  # Display x-axis labels
#plt.tight_layout()  # Adjust layout for better display

#plt.savefig('ZF_1534_1000pum_1338-E_3.png', format='png')
#plt.show()
