
# coding: utf-8

# In[1]:

import pandas as pd
import os
import numpy as np


# In[2]:

def impact(x): 
    if x == 'HIGH':
        return 3
    elif x == 'MODERATE':
        return 2
    else:
        return 1

def resolve_duplicates(dataframe, key):
    dupes = dataframe[dataframe.duplicated([key],keep=False)]
    dupes = dupes.sort_values(by=key)
    dupes.to_csv('duplicated_variants.txt', sep="\t", mode='a',index=False)
    dataframe.drop_duplicates(subset=key, keep=False, inplace=True)
    i = 0
    while i < len(dupes.index):
        dupeframe = dupes[dupes[key] == dupes.iloc[i][key]]
        if key == "Hugo_Symbol":
            dupeframe['IMPACT'] = pd.Categorical(dupeframe['IMPACT'], ["HIGH", "MODERATE","LOW","MODIFIER"])
            dupeframe.sort_values("IMPACT", inplace=True)
        else:
            dupeframe.sort_values("ploidy_adj_cn", inplace=True, ascending=False)
        #print(dupeframe)
        #print(dupeframe.iloc[0])
        dataframe = dataframe.append(dupeframe.iloc[0])
        i = i + len(dupeframe.index)
    return dataframe


# In[3]:

maf_dir = "MAF/"
maf_df = pd.DataFrame()
print("Processsing MAF files:")
for i in os.listdir(maf_dir):
    if i.endswith('.maf') or i.endswith('.MAF'):
        print("Reading in file: " + i)
        input_df = pd.read_csv(maf_dir + "/" + i,sep="\t")
        input_df = resolve_duplicates(input_df, "Hugo_Symbol")
        maf_df = maf_df.append(input_df)
        
#filter for HIGH or MODERATE impact rows
maf_df = maf_df[(maf_df["IMPACT"]=="HIGH") | (maf_df["IMPACT"]=="MODERATE")]


# In[61]:
ploidy_exists = False
if os.path.isfile('metadata.txt') :
    metadata = pd.read_csv("metadata.txt",sep="\t")
    if 'Ploidy' in metadata.columns and 'Gender' in metadata.columns:
        print("metadata.txt with 'Gender' exists")
        ploidy_exists = True
    
bed_dir = "BED/"
bed_df = pd.DataFrame()
print("Processing BED files")
for i in os.listdir(bed_dir):
    if i.endswith('.BED') or i.endswith('.bed'):
        print("Reading in file: " + i)
        input_df = pd.read_csv(bed_dir + "/" + i,sep="\t")
        #adjust for ploidy
        if ploidy_exists:
            samplename = str(input_df['tumor'][0]) + "--" + str(input_df['normal'][0])
            sample_ploidy = metadata[metadata['Sample'] == samplename]['Ploidy'].iloc[0]
            sample_ploidy = round(sample_ploidy * 2) / 2
            input_df['ploidy_adj_cn'] = abs(input_df['cnv_total_cn'] - sample_ploidy)
            if(metadata[metadata['Sample'] == samplename]['Gender'].iloc[0] == 'male'):
                print('male sample')
                input_df.loc[input_df['cnv_chr'] == 'chrX', 'ploidy_adj_cn'] = abs(input_df.loc[input_df['cnv_chr'] == 'chrX', 'cnv_total_cn'] - (sample_ploidy/2))
                input_df.loc[input_df['cnv_chr'] == 'chrY', 'ploidy_adj_cn'] = abs(input_df.loc[input_df['cnv_chr'] == 'chrY', 'cnv_total_cn'] - (sample_ploidy/2))
        else:
            input_df['ploidy_adj_cn'] = abs(input_df['cnv_total_cn'] - 2) #adjust for ploidy 2 for now
        input_df = resolve_duplicates(input_df, "gene_name")
        input_df['cnv_scale'] = 'largescale'
        input_df.loc[input_df['cnv_end'] - input_df['cnv_start'] <= 3000000, 'cnv_scale'] = 'focal'
        bed_df = bed_df.append(input_df)
        

# In[184]:

#sample	gene	impact	cnv	ploidy(?)	metadata(comma separated?)

df = pd.merge(maf_df, bed_df, how='outer', left_on=['Hugo_Symbol','Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode'],right_on=['gene_name','tumor','normal'], validate="one_to_one")


# In[185]:

def merge_keys(row):
    if pd.isnull(row['Hugo_Symbol']):
        gene = row['gene_name']
    else:
        gene = row['Hugo_Symbol']
    if pd.isnull(row['Tumor_Sample_Barcode']):
        sample = str(row['tumor']) + '--' + str(row['normal'])
    else:
        sample = str(row['Tumor_Sample_Barcode']) + '--' + str(row['Matched_Norm_Sample_Barcode'])
    
    return gene, sample

newcols = df.apply(merge_keys, axis=1, result_type='expand')
df['Gene'] = newcols[0]
df['Sample'] = newcols[1]


#df[['Gene','Sample']] = test


# In[186]:

geneList = df['Gene'].unique()


# In[190]:

#python implementation of the app side ranking? not sure to include it here or in app

#ploidy = 2 #harcode for now
#rankList = pd.DataFrame(columns = ["gene","score"])
#rankList
#for i in geneList:
#    geneData = df[df['Gene'] == i]
#    high_count = len(geneData[geneData['IMPACT'] == "HIGH"])
#    moderate_count = len(geneData[geneData['IMPACT'] == "MODERATE"])
#    high_amp_count = len(geneData[geneData['cnv_total_cn'] - ploidy >= 2])
#    amp_count = len(geneData.query('cnv_total_cn - @ploidy < 2 & cnv_total_cn - @ploidy >= 1'))
#    deep_del_count = len(geneData[geneData['cnv_total_cn'] - ploidy <= -2])
#    del_count = len(geneData.query('cnv_total_cn - @ploidy > -2 & cnv_total_cn - @ploidy <= -1'))
#    score = high_count * 20 + moderate_count * 10 + deep_del_count * 4 + high_amp_count * 3 + del_count * 2 + amp_count 
#    rankList = rankList.append({'gene':i, 'score': score}, ignore_index = True)
#    
#rankList.sort_values('score', ascending=False, inplace=True)
#print(rankList)


# In[172]:

print("Making oncomatrix...")
output = df[["Sample","Gene","IMPACT","cnv_total_cn",'cnv_scale']]
output.columns = ["Sample","Gene","Impact","CNV","CNV_Scale"]
output.to_csv("./oncomatrix.txt",sep="\t",index=False)

print("Done.")


# In[ ]:



