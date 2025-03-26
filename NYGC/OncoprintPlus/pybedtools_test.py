import pandas as pd
import os
import numpy as np
import math
import pybedtools

print("Reading in file: ")
input_df = pd.read_csv("BED/MILD-B58G-TTP1-C-1-1-D-A865-36.segtab.txt",sep="\t")
input_df[input_df['End.bp'] - input_df['Start.bp'] <= 3000000]
#input_df.columns = ['gene_chr','gene_start','gene_end','gene_name','cnv_chr','cnv_start','cnv_end','cnv_total_cn','overlap_bp','tumor','normal']
#adjust for ploidy
#- If the fractional value of the ploidy is below 0.4, round down (e.g. 3.2 -> 3)
#- If the fractional value of the ploidy is above 0.6, round up (e.g. 3.6 -> 4)
#- Anything above neutral CN is a gain
#- Anything below neural CN is a deletion
#CN=0 complete loss, and 2X ploidy amplification
#- If the fractional value of the ploidy is bewteen 0.4 and 0.6, round up and down (neutral CN for 3.5 is 3 and 4)
#For example, if the ploidy is 3.5, neutral copy numbers are 3 and 4. Anything greater than 4 is a gain, anything less than 3 is a deletion
samplename = str(input_df['sample'][0])
if ploidy_exists:
    sample_ploidy = metadata[metadata['tumor_submitter_id'] == samplename]['Consensus Ploidy Value'].iloc[0]
else:
    sample_ploidy = 2.0
ploidy_fraction = sample_ploidy % 1
if ploidy_fraction > 0.6:
    rounded_ploidy = math.ceil(sample_ploidy)
    gain_thresh = rounded_ploidy
    del_thresh = rounded_ploidy
elif ploidy_fraction < 0.4:
    rounded_ploidy = math.floor(sample_ploidy)
    gain_thresh = rounded_ploidy
    del_thresh = rounded_ploidy
else:
    rounded_ploidy = math.floor(sample_ploidy)
    gain_thresh = rounded_ploidy+1
    del_thresh = rounded_ploidy

input_df.loc[input_df['cnv_total_cn'] < del_thresh, 'cn_class'] = 'deletion'
input_df.loc[input_df['cnv_total_cn'] > gain_thresh, 'cn_class'] = 'gain'
input_df.loc[input_df['cnv_total_cn'] == 0, 'cn_class'] = 'loss'
input_df.loc[input_df['cnv_total_cn'] >= (del_thresh * 2), 'cn_class'] = 'amplification'

print(input_df)

#bedtools intersect -wo -a  -b stdin
#ensembl_genes.intersect()