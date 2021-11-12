import os, sys
import pandas as pd
import xml.dom.minidom


def get_data1(disease, information, data_dir):
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

        atmp = df_csv.loc[(df_csv.iloc[:, -1] == information[j]) &
                          (df_csv['characteristics[disease]'] == disease[i])]
        if atmp.empty:
          continue

        label = information[i] + '; ' + disease[j]
        assay_label.append(label)

        aid = 'g' + str(g)
        assay_id.append(aid)
        g += 1

        assay.append(list(atmp['source name']))

        if g % 2 == 0:
          cid = 'g' + str(g) + '_g' + str(g - 1)
          contrast_id.append(cid)

          cname = "\'" + disease[j] + "\' vs \'" + disease[j - 1] + "\' in \'" + information[i] + "\'"
          contrast_name.append(cname)

          gname = disease[j] + '-' + disease[j - 1]
          g_name.append(gname)
          g_id.append(cid)
          gname = disease[j - 1] + '-' + disease[j]
          g_name.append(gname)
          g_id.append(cid)

  else:
    for i in range(len_d):
      for j in range(len_i):

        atmp = df_csv.loc[(df_csv.iloc[:, -1] == information[j]) &
                          (df_csv['characteristics[disease]'] == disease[i])]
        if atmp.empty:
          continue

        label = disease[i] + '; ' + information[j]
        assay_label.append(label)

        aid = 'g' + str(g)
        assay_id.append(aid)
        g += 1

        assay.append(list(atmp['source name']))

        if g % 2 == 0:
          cid = 'g' + str(g) + '_g' + str(g - 1)
          contrast_id.append(cid)

          cname = "\'" + information[j] + "\' vs \'" + information[j - 1] + "\' in \'" + disease[i] + "\'"
          contrast_name.append(cname)

          gname = information[j] + '-' + information[j - 1]
          g_name.append(gname)
          g_id.append(cid)
          gname = information[j - 1] + '-' + information[j]
          g_name.append(gname)
          g_id.append(cid)

  gg = g_name + g_id
  gg = ','.join(gg)
  print(gg)

  f = open(data_dir + '/g_g_name.txt', 'w')
  f.write(gg)
  f.close()

  return assay_label, assay_id, assay, contrast_name, contrast_id


def get_data2(disease, information, data_dir, tmp_info0, tmp_info1, judge_information):
  len_d = len(disease)
  len_i = len(information)
  len_i0 = len(tmp_info0)
  len_i1 = len(tmp_info1)

  assay_label = []
  assay_id = []
  assay = []

  contrast_name = []
  contrast_id = []

  g = 1
  g_name = []
  g_id = []
  if len_d >= len_i:
    for i in range(len_i0):
      for x in range(len_i1):
        for j in range(len_d):

          atmp = df_csv.loc[(df_csv[judge_information[0]] == tmp_info0[i]) &
                            (df_csv[judge_information[1]] == tmp_info1[x]) &
                            (df_csv['characteristics[disease]'] == disease[j])]
          if atmp.empty:
            continue

          label = disease[j] + '; ' + tmp_info0[i] + '|' + tmp_info1[x]
          assay_label.append(label)

          aid = 'g' + str(g)
          assay_id.append(aid)
          g += 1

          assay.append(list(atmp['source name']))

          if g % 2 == 0:
            cid = 'g' + str(g) + '_g' + str(g - 1)
            contrast_id.append(cid)

            tmp_info_g = tmp_info0[i] + '|' + tmp_info1[x]
            tmp_info_g_1 = tmp_info0[i] + '|' + tmp_info1[x - 1]

            cname = "\'" + tmp_info_g + "\' vs \'" + tmp_info_g_1 + "\' in \'" + disease[i] + "\'"
            contrast_name.append(cname)

            gname = tmp_info_g + '-' + tmp_info_g_1
            g_name.append(gname)
            g_id.append(cid)
            gname = tmp_info_g_1 + '-' + tmp_info_g
            g_name.append(gname)
            g_id.append(cid)

            g += 2

            cid = 'g' + str(g) + '_g' + str(g - 1)
            contrast_id.append(cid)

            tmp_info_g = tmp_info0[j] + '|' + tmp_info1[x]
            tmp_info_g_1 = tmp_info0[j - 1] + '|' + tmp_info1[x - 1]

            cname = "\'" + tmp_info_g + "\' vs \'" + tmp_info_g_1 + "\' in \'" + disease[i] + "\'"
            contrast_name.append(cname)

            gname = tmp_info_g + '-' + tmp_info_g_1
            g_name.append(gname)
            g_id.append(cid)
            gname = tmp_info_g_1 + '-' + tmp_info_g
            g_name.append(gname)
            g_id.append(cid)


  else:
    for i in range(len_d):
      for j in range(len_i0):
        for x in range(len_i1):

          atmp = df_csv.loc[(df_csv[judge_information[0]] == tmp_info0[j]) &
                            (df_csv[judge_information[1]] == tmp_info1[x]) &
                            (df_csv['characteristics[disease]'] == disease[i])]
          if atmp.empty:
            continue

          label = disease[i] + '; ' + tmp_info0[j] + '|' + tmp_info1[x]
          assay_label.append(label)

          aid = 'g' + str(g)
          assay_id.append(aid)
          g += 1

          assay.append(list(atmp['source name']))

          if g % 2 == 0:
            cid = 'g' + str(g) + '_g' + str(g - 1)
            contrast_id.append(cid)

            tmp_info_g = tmp_info0[j] + '|' + tmp_info1[x]
            tmp_info_g_1 = tmp_info0[j] + '|' + tmp_info1[x - 1]

            cname = "\'" + tmp_info_g + "\' vs \'" + tmp_info_g_1 + "\' in \'" + disease[i] + "\'"
            contrast_name.append(cname)

            gname = tmp_info_g + '-' + tmp_info_g_1
            g_name.append(gname)
            g_id.append(cid)
            gname = tmp_info_g_1 + '-' + tmp_info_g
            g_name.append(gname)
            g_id.append(cid)

            g += 2

            cid = 'g' + str(g) + '_g' + str(g - 1)
            contrast_id.append(cid)

            tmp_info_g = tmp_info0[j] + '|' + tmp_info1[x]
            tmp_info_g_1 = tmp_info0[j - 1] + '|' + tmp_info1[x - 1]

            cname = "\'" + tmp_info_g + "\' vs \'" + tmp_info_g_1 + "\' in \'" + disease[i] + "\'"
            contrast_name.append(cname)

            gname = tmp_info_g + '-' + tmp_info_g_1
            g_name.append(gname)
            g_id.append(cid)
            gname = tmp_info_g_1 + '-' + tmp_info_g
            g_name.append(gname)
            g_id.append(cid)

  gg = g_name + g_id
  gg = ','.join(gg)

  f = open(data_dir + '/g_g_name.txt', 'w')
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

  fp = open(data_dir + '/NAME_configuration.xml', 'w')
  doc.writexml(fp, indent='\t', addindent='\t', newl='\n', encoding="utf-8")
  fp.close()


# def get_convert_exp_values(data_dir):
#   name_data = pd.read_csv(data_dir + 'g_g_name.txt', sep=',')
#
#   for (i in 1:(length(name_data) / 2)){
#   data = read.csv("./MSstats_output.csv")
#
#   temp.name = name_data[, i]
#   data < - data[which(data$Label % in % temp.name), ]
#
#   if (any( is.na(data))) {
#   matrix_final = data.frame((Protein = data$Protein))
#   colnames(matrix_final)[1] < - 'Protein'
#   break
#
# }
# }
#
#
#
# for (i in 1:(length(name_data) / 2)){
# data < - read.csv("./MSstats_output.csv")
#
# temp.name = name_data[, i]
# data < - data[which(data$Label % in % temp.name), ]
#
# if (any( is.na(data))) {
# num = 1
#
# name_p = paste(name_data[, i + length(name_data) / 2], '.p-value', sep = '')
# name_log2fc = paste(name_data[, i + length(name_data) / 2], '.log2foldchange', sep = '')
#
# matrix_tmp = data.frame(pp=data$pvalue, llog2 = data$log2FC)
# colnames(matrix_tmp)[1] < - name_p
# colnames(matrix_tmp)[2] < - name_log2fc
#
# matrix_final = cbind(matrix_final, matrix_tmp)
# }
# else {
# next
# }
# }
#
#
# ### Delete lines containing 'NA'
# matrix_final = na.omit(matrix_final)
#
# write.csv(matrix_final, file='./NAME_analytics.csv', quote=F, row.names = F)
#
# }

data_dir = '../data/'

sdrfName = 'PXD015270-cell-lines.sdrf.tsv'
csv = data_dir + "/" + sdrfName
df_csv = pd.read_csv(csv, sep='\t')

rownames = list(df_csv.columns)
low_rownames = [i.lower() for i in rownames]
df_csv.columns = low_rownames

disease = df_csv['characteristics[disease]']
disease = sorted(disease.unique())

all_information = list(df_csv)
judge_information = list(filter(lambda x: len(x) != len(x.replace('factor', '')), all_information))

if len(judge_information) == 1:
  information = sorted(df_csv[judge_information[0]].unique())

  assay_label, assay_id, assay, contrast_name, contrast_id = get_data1(disease, information, data_dir)

elif len(judge_information) == 2:
  for i in range(2):
    exec("tmp_info%s = sorted(df_csv[judge_information[%s]].unique())" % (i, i))

  information = []
  for i in range(len(tmp_info0)):
    for j in range(len(tmp_info1)):
      information.append(tmp_info0[i] + '|' + tmp_info1[j])

  assay_label, assay_id, assay, contrast_name, contrast_id = get_data2(disease, information, data_dir, tmp_info0,
                                                                       tmp_info1, judge_information)

else:
  sys.exit(1)

get_xml(assay_label, assay_id, assay, contrast_name, contrast_id, data_dir)
