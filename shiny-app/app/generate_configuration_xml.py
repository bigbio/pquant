import os
import pandas as pd
import xml.dom.minidom

def get_data(disease, information, data_dir):

    len_d = len(disease)
    len_i = len(information)

    assay_label = []
    assay_id = []
    assay = []

    contrast_name = []
    contrast_id = []

    g = 1
    g_name = []
    g_id = []
    if len_d >= len_i:
        for i in range(len_i):
            for j in range(len_d):

                atmp = df_csv.loc[(df_csv.iloc[:,-1] == information[j]) &
                            (df_csv['characteristics[disease]'] == disease[i])]
                if atmp.empty:
                    continue

                label = information[i] + '; ' + disease[j]
                assay_label.append(label)
                
                aid = 'g' + str(g)
                assay_id.append(aid)
                g += 1

                assay.append(list(atmp['Source Name']))

                if g % 2 == 0:
                    cid = 'g' + str(g) + '_g' + str(g-1)
                    contrast_id.append(cid)

                    cname = "\'" + disease[j] + "\' vs \'" + disease[j-1] + "\' in \'" + information[i] + "\'"
                    contrast_name.append(cname)

                    gname = disease[j] + '-' + disease[j-1]
                    g_name.append(gname)
                    g_id.append(cid)
                    gname = disease[j-1] + '-' + disease[j]
                    g_name.append(gname)
                    g_id.append(cid)

    else:
        for i in range(len_d):
            for j in range(len_i):

                atmp = df_csv.loc[(df_csv.iloc[:,-1] == information[j]) &
                            (df_csv['characteristics[disease]'] == disease[i])]
                if atmp.empty:
                    continue

                label = disease[i] + '; ' + information[j]
                assay_label.append(label)
                
                aid = 'g' + str(g)
                assay_id.append(aid)
                g += 1
                
                assay.append(list(atmp['Source Name']))

                if g % 2 == 0:
                    cid = 'g' + str(g) + '_g' + str(g-1)
                    contrast_id.append(cid)

                    cname = "\'" + information[j] + "\' vs \'" + information[j-1] + "\' in \'" + disease[i] + "\'"
                    contrast_name.append(cname)

                    gname = information[j] + '-' + information[j-1]
                    g_name.append(gname)
                    g_id.append(cid)
                    gname = information[j-1] + '-' + information[j]
                    g_name.append(gname)
                    g_id.append(cid)

    gg = g_name + g_id
    gg = ','.join(gg)

    f = open(data_dir + 'g_g_name.txt', 'w')
    f.write(gg)
    f.close()

    return assay_label, assay_id, assay, contrast_name, contrast_id




def get_xml(assay_label, assay_id, assay, contrast_name, contrast_id, data_dir):
    doc = xml.dom.minidom.Document()

    root = doc.createElement('Configuration')
    root.setAttribute('experimentType', 'rnaseq_mrna_differential')
    root.setAttribute('r_data', '0')
    doc.appendChild(root)

    node_analytics = doc.createElement('analytics')
    root.appendChild(node_analytics)

    node_assayGs = doc.createElement('assay_groups')
    node_analytics.appendChild(node_assayGs)

    node_contrasts = doc.createElement('contrasts')
    node_analytics.appendChild(node_contrasts)

    for i in range(len(assay_id)):
        node_assayG = doc.createElement('assay_group')
        node_assayG.setAttribute('id', assay_id[i])
        node_assayG.setAttribute('label', assay_label[i])
        node_assayGs.appendChild(node_assayG)
        
        for j in range(len(assay[i])):
            node_assay = doc.createElement('assay')
            node_assay.appendChild(doc.createTextNode(assay[i][j]))
            node_assayG.appendChild(node_assay)
        


    for i in range(len(contrast_id)):
        node_contrast = doc.createElement('contrast')
        node_contrast.setAttribute('id', contrast_id[i])
        node_contrast.setAttribute('cttv_primary', '1')
        node_contrasts.appendChild(node_contrast)

        node_c_name = doc.createElement('name')
        node_c_name.appendChild(doc.createTextNode(contrast_name[i]))

        node_c_ref = doc.createElement('reference_assay_group')
        node_c_ref.appendChild(doc.createTextNode(contrast_id[i][:2]))

        node_c_test = doc.createElement('test_assay_group')
        node_c_test.appendChild(doc.createTextNode(contrast_id[i][-2:]))

        node_contrast.appendChild(node_c_name)
        node_contrast.appendChild(node_c_ref)
        node_contrast.appendChild(node_c_test)


    fp = open(data_dir + 'NAME_configuration.xml', 'w')
    doc.writexml(fp, indent='\t', addindent='\t', newl='\n', encoding="utf-8")
    fp.close()






if __name__ == "__main__":
    
    data_dir = 'D:\\dataset\\R downstream analysis\\shiny\\data\\' # change to the path where you put data
    #data_dir = os.getcwd()
    csv = data_dir  + "PXD015270-cell-lines.sdrf.csv"
    df_csv = pd.read_csv(csv)

    disease = df_csv['characteristics[disease]']
    information = df_csv.iloc[:,-1]  # you can change the column names that need to be compared

    disease = sorted(disease.unique())
    information = sorted(information.unique())

    

    assay_label, assay_id, assay, contrast_name, contrast_id = get_data(disease, information, data_dir)

    get_xml(assay_label, assay_id, assay, contrast_name, contrast_id, data_dir)   
