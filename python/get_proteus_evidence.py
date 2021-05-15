import re
import pandas as pd
import numpy as np
from collections import Counter


def function(sequence, reference):
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
            top = top.replace(top[a:b + 1], '')  # Replace the data in parentheses with null
        result = re.sub('[\W_\d]+', '', top)  # It means a filter character
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
            temp = k[c + 1:c + 3]  # The first two strings are lowercase
            temp = temp.lower()
            k = k.replace(k[c + 1:d], temp)
            top = top.replace(top[a:b + 1], '')
            next = 4
        k = k.replace(".", '')
        evi_modified_sequence.append(date + k + date)

    
    for i in range(0, len(csv_PeptideSequence)):
        top = csv_PeptideSequence[i]
        flag = True   # True means this sequence is unmodification 
        out_C = {}
        while "(" in top:
            flag = False
            a = top.find('(')
            b = top.find(')')
            temp = top[a + 1: b]  # Store modification
            lent = len(top.split(temp)) - 1  # Number of modifications' occurrences
            top = top.replace(top[a:b + 1], '')
            out_C[temp] = lent
        
        # Convert a dictionary into an array
        tem = []
        for key, value in out_C.items():
            if value == 1:  # If the modification occurs only once, the number need not be displayed
                tem.append(key)
            else:
                tem.append(str(value) + ' ' + key)
        tem = ','.join(tem)

        if flag:
            evi_modifications.append("Unmodified")
        else:
            evi_modifications.append(tem)
    return evi_sequence, evi_modified_sequence, evi_modifications, evi_experiment


if __name__ == "__main__":
    
    # TODO data processing
    now_dir = r"D:\dataset\R downstream analysis\proteus\code"
    csv = now_dir + '\\' + "out.csv"
    mztab = now_dir + '\\' + "pep.csv"
    
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

    evi_sequence, evi_modified_sequence, evi_modifications, evi_experiment = function(csv_PeptideSequence, csv_experiment)
    evidence = pd.DataFrame({
        "PeptideSequence": evi_sequence,
        "modified_sequence": evi_modified_sequence,
        "modifications": evi_modifications,
        "protein_group": evi_protein_group,
        #"protein": [],
        "experiment": evi_experiment,
        "charge": evi_charge,
        #"reverse": [],
        #"contaminant": [],
        "intensity": evi_intensity})
    evidence.to_csv(r"result_1.csv")
    
    # TODO：Implements the VLOOKUP function
    data_text = "result_1.csv"
    pep_text = "pep.csv"
    pep = pd.read_csv(pep_text)
    df = pd.read_csv(data_text)
    data = df["PeptideSequence"]
    Pep = pep[["sequence", "accession"]]
    Pep = Pep.drop_duplicates(subset="sequence")  # Delete the duplicate values in the second table
    #print(data.shape)
    df_merge = pd.merge(left=df, right=Pep, left_on="PeptideSequence", right_on="sequence", how='left', )
    #print(df_merge.shape)
    df_merge.to_csv("result.csv", index=False)
