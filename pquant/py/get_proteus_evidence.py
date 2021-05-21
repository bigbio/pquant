import re
import pandas as pd
import numpy as np
from collections import Counter
import os


def function(sequence, reference, evi_protein_group):
    evi_sequence = []  # Total cleaning data
    evi_modified_sequence = []
    out_C = []
    evi_modifications = []
    evi_experiment = []

    for i in range(0, len(csv_PeptideSequence)):
        top = csv_PeptideSequence[i]
        while "(" in top:
            a = top.find('(')
            b = top.find(')')
            top = top.replace(top[a:b + 1], '')  # Replace the data in brackets with empty
        result = re.sub('[\W_\d]+', '', top)  # Filter character
        evi_sequence.append(result)
        evi_experiment.append(reference[i][:-5])
    #print(len(evi_sequence))
    #print(len(evi_experiment))

    for i in range(0, len(csv_PeptideSequence)):
        top = csv_PeptideSequence[i]
        k = csv_PeptideSequence[i]
        date = "_"
        #row = []
        #col = []
        next = 0
        while "(" in top:
            a = top.find('(')
            b = top.find(')')
            c = a + next
            d = b + next
            temp = k[c + 1:c + 3]  # Lowercase the first two strings
            temp = temp.lower()
            k = k.replace(k[c + 1:d], temp)
            top = top.replace(top[a:b + 1], '')
            next = 4
        k = k.replace(".", '')
        evi_modified_sequence.append(date + k + date)

    
    for i in range(0, len(csv_PeptideSequence)):
        top = csv_PeptideSequence[i]
        flag = True   # True means this line is unmodified
        out_C = {}
        while "(" in top:
            flag = False
            a = top.find('(')
            b = top.find(')')
            temp = top[a + 1: b]  # Storage modification
            lent = len(top.split(temp)) - 1  # Modification occurrences
            top = top.replace(top[a:b + 1], '')
            out_C[temp] = lent
        
        # Convert dictionary to array
        tem = []
        for key, value in out_C.items():
            if value == 1:  # If it only appears once, no need to display the number of occurrences
                tem.append(key)
            else:
                tem.append(str(value) + ' ' + key)
        tem = ','.join(tem)

        if flag:
            evi_modifications.append("Unmodified")
        else:
            evi_modifications.append(tem)

    evi_protein = []
    for i in range(0, len(evi_protein_group)):
        protein_group = evi_protein_group[i]
        if ";" in protein_group:
            semicolon = protein_group.index(';')
            evi_protein.append(protein_group[:semicolon])
        else:
            evi_protein.append(protein_group)

    return evi_sequence, evi_modified_sequence, evi_modifications, evi_experiment, evi_protein


def get_mztab(pri_mztab, data_dir):
    df_pri = pd.read_csv(pri_mztab)

    tmp_peh = (df_pri[df_pri['MTD'].isin(['PEH'])].values)[0]
    tmp_pep = df_pri[df_pri['MTD'].isin(['PEP'])].values

    tmp = pd.DataFrame(
        data = tmp_pep,
        columns = tmp_peh,
        index = range(len(tmp_pep))
    )

    tmp.to_csv(data_dir + '\\pep.csv', index=False)

    return ''





if __name__ == "__main__":
    
    # TODO These code implements data processing
    now_dir = os.getcwd()
    data_dir = now_dir

    csv = data_dir  + "\\out_msstats.csv"
    pri_mztab = data_dir  + "\\out_mzTab.csv"

    get_mztab(pri_mztab, data_dir)

    mztab = data_dir + "\\pep.csv"
    
    df_csv = pd.read_csv(csv)
    df_mztab = pd.read_csv(mztab)

    csv_PeptideSequence = df_csv['PeptideSequence']
    csv_ProteinName = df_csv['ProteinName']
    csv_experiment = df_csv['Reference']
    csv_PrecursorCharge = df_csv['PrecursorCharge']
    csv_intensity = df_csv['Intensity']

    evi_protein_group = csv_ProteinName
    evi_charge = csv_PrecursorCharge
    evi_intensity = csv_intensity


    evi_sequence, evi_modified_sequence, evi_modifications, evi_experiment, evi_protein = function(csv_PeptideSequence, csv_experiment, evi_protein_group)
    evidence = pd.DataFrame({
        "PeptideSequence": evi_sequence,
        "modified_sequence": evi_modified_sequence,
        "modifications": evi_modifications,
        "protein_group": evi_protein_group,
        "protein": evi_protein,
        "experiment": evi_experiment,
        "charge": evi_charge,
        #"reverse": [],
        #"contaminant": [],
        "intensity": evi_intensity})
    evidence.to_csv(data_dir + "\\result_1.csv")
    
    # TODO：The following code implements the VLOOKUP function operation
    data_text = data_dir + "\\result_1.csv"
    pep_text = data_dir + "\\pep.csv"
    pep = pd.read_csv(pep_text)
    df = pd.read_csv(data_text)
    data = df["PeptideSequence"]
    Pep = pep[["sequence", "accession"]]
    Pep = Pep.drop_duplicates(subset="sequence")  # De-duplicate the second table
    
    df_merge = pd.merge(left=df, right=Pep, left_on="PeptideSequence", right_on="sequence", how='left', )


    if df_merge.columns[0] != 'PeptideSequence':
        tmp = df_merge.iloc[:,1:len(df_merge.columns)]  # Deletes additional data in the first column that is not known when it was generated
        df_merge = tmp
        #print(df_merge)

    df_merge.to_csv(data_dir + "\\out_proteus.csv", index=False)
