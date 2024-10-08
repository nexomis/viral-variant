#!/usr/bin/env python3


import pandas as pd
import os

def filter_variants(tsv_files, csv_positions, prefix):
    # import interest positions
    positions_df = pd.read_csv(csv_positions, sep='\\t', header=None, names=['chr', 'pos'])
    positions_set = set(zip(positions_df['chr'], positions_df['pos']))
    
    output_dir = prefix + "_batchFiltered"
    os.makedirs(output_dir)
    common_columns = ['REGION', 'POS', 'REF', 'ALT', 'GFF_FEATURE', 'REF_CODON', 'REF_AA', 'ALT_CODON', 'ALT_AA', 'POS_AA']
    all_samples_filtered = []
    for tsv_file in tsv_files.split():
        sample_name = os.path.basename(tsv_file).split('.')[0]

        # import, filter and save in specific file tsv of sample
        tsv_df = pd.read_csv(tsv_file, sep='\\t', low_memory=False)     # 'low_memory=False' to avoid the following warning (which seems unfounded): '<stdin>:1: DtypeWarning: Columns (3,13,14,15,16,17,18) have mixed types. Specify dtype option on import or set low_memory=False.'. But it may have a significant impact on memory in the case of large genomes.
        filtered_df = tsv_df[tsv_df.apply(lambda row: (row['REGION'], row['POS']) in positions_set, axis=1)]
        filtered_output_path = os.path.join(output_dir, f"{sample_name}_batchFiltered.tsv")
        filtered_df.to_csv(filtered_output_path, sep='\\t', index=False)
        
        # reordone and rename specific columns using sample_name and append to list of df to merge
        variable_columns = [col for col in tsv_df.columns if col not in common_columns]
        filtered_df = filtered_df[common_columns + variable_columns]
        renamed_columns = {col: f"{col}({sample_name})" for col in variable_columns}
        filtered_df = filtered_df.rename(columns=renamed_columns)
        all_samples_filtered.append(filtered_df)

    # merge filtered df of all sample on one unique file
    global_output_file = prefix + "_summary_all_iSNVs.tsv"
    global_df = all_samples_filtered[0]
    for df in all_samples_filtered[1:]:
        global_df = pd.merge(global_df, df, on=common_columns, how='outer')
    global_df.to_csv(global_output_file, sep='\\t', index=False)

    # TODO: for each pos in global_df, for each smpl, increade value of REF_DP, ALT_DP, ALT_FREQ and TOAL_DP using dictionay constructed about mpileup result.



def main():
    tsv_files="${ivar_tsv_files}"
    csv_positions="${pos_file}"
    prefix="${meta.id}"

    filter_variants(tsv_files, csv_positions, prefix)


if __name__ == "__main__":
    main()
