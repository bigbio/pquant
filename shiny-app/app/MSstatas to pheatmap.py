import pandas as pd

def MS_change_pht(MS):
    MS = MS[~MS['log2FC'].isin(['Inf'])]  # By ~ inverse;Without the outermost log2FC[], the output would be True/False
    MS = MS[~MS['log2FC'].isin(['-Inf'])]

    Protein = MS['Protein'].drop_duplicates().tolist()
    Label = MS['Label'].drop_duplicates().tolist()

    pHT_Protein = []
    pHT_log2FC = []

    for i in range(len(Protein)):
        tmp = MS[MS['Protein'].isin([Protein[i]])]
        len_tmp = len(tmp['log2FC'])
        if len_tmp == len(Label):
            pHT_Protein.append(Protein[i])
            pHT_log2FC.append(tmp['log2FC'].tolist())

    pHT = pd.DataFrame(pHT_Protein, columns=['Protein'])
    for i in range(len(Label)):
        pHT[Label[i]] = [j[i] for j in pHT_log2FC]    # Get every i-th element in a two-dimensional list
    
    pHT.drop(index=[15], inplace=True) # 'sp|P40926|MDHM_HUMAN' This row of data is an outlier. At present, we manually delete the data.
    return pHT


if __name__ == '__main__':

    # TODO This code implements data processing
    now_dir = r'D:\dataset\R downstream analysis\shiny\shiny-app'   # Your own work path
    MS_ouput = now_dir + '\\' + 'MSstats_output.csv'

    df_MS = pd.read_csv(MS_ouput)
    
    # pHT = pheatmap
    # We only need 3 columns of Protein, Label, and log2FC
    # Output file format: the protein column remains unchanged (the ordinate), log2FC becomes the abscissa, find the log2FC value corresponding to each protein to fill the matrix

    pHT = MS_change_pht(df_MS)   

    pHT.to_csv(now_dir + '\\' + 'pheatmap_input.csv', index=False)