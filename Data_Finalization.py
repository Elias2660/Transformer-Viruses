#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#mount google drive (in colab) to get the files

import os

from google.colab import drive
drive.mount('/content/drive')

os.listdir("/content/drive/MyDrive/Transformer-Viruses")

# In[ ]:


#making the dataframe of transmission levels

import pandas as pd
virus_data = pd.read_excel("/content/drive/MyDrive/Transformer-Viruses/data/Woolhouse_and_Brierley_RNA_virus_database.xlsx")

# In[ ]:


virus_data.head()

# In[ ]:


genome_data = pd.DataFrame(columns=["Locus", "Position/Length Indicator?", "Virus Name", "Genome", "Estimated Transmission Level", "temp ratio"])

# In[ ]:


#making the dataframe of genomes

index = 0
stor = []
with open("/content/drive/MyDrive/Transformer-Viruses/data/humanVirusSeq.fa") as d:
    lines = d.readlines()
    print(len(lines))
    for line in lines:
        if index % 2 == 0:
            stor = line.strip().replace("\n","").replace(">","").split("|")
            index += 1
        else:
            genome_data.loc[len(genome_data.index)] = [stor[0], stor[1], stor[2], line.strip(), None, 0]
            index += 1

# In[ ]:


genome_data

# In[ ]:


len(genome_data.iloc[0]["Genome"].strip())

# In[ ]:


#test of the normalizer used to distinguish genomes with and without corresponding transmission levels

i = 0
j = 0

while True:
    compValues.append(df.SequenceMatcher(None, virus_data.iloc[i,0], genome_data.iloc[j, 2]).ratio())
    i += 1
    if i == 213:
      i = 0
      compValues.append(1.0)
      norm_compValues = [(float(i)-min(compValues))/(max(compValues)-min(compValues)) for i in compValues]
      #print(norm_compValues)
      for k in range(213):
        if norm_compValues[k] > 0.75:
          #print(k, virus_data.iloc[k,20])
          break
      break

# In[ ]:


# adding transmission levels to all genomes that have corresponding transmission levels

import difflib as df
i = 0
j = 0
compValues = []

for j in range(1620):
  while True:
    compValues.append(df.SequenceMatcher(None, virus_data.iloc[i,0], genome_data.iloc[j, 2]).ratio())
    i += 1
    if i == 213:
      i = 0
      compValues.append(1.0)
      norm_compValues = [(float(i)-min(compValues))/(max(compValues)-min(compValues)) for i in compValues]
      for k in range(213):
        if norm_compValues[k] > 0.75:
          genome_data.loc[j, "Estimated Transmission Level"] = virus_data.iloc[k,20]
          break
      break
  compValues = []
  j += 1
  print(j)

genome_data

# In[ ]:


genome_data.loc[52]

# In[ ]:


#delete rows that dont have any transmission level recorded

deleteList = []

for i in range(len(genome_data)):
  if genome_data.loc[i, "Estimated Transmission Level"] == None:
    deleteList.append(i)

print(deleteList)

genome_data = genome_data.drop(deleteList)
genome_data = genome_data.reset_index(drop=True)

# In[ ]:


genome_data = genome_data.drop(columns = "temp ratio")
genome_data

# In[ ]:


#download the resulting spreadsheet

genome_data.to_excel("genome_data_with_transmission_levels.xlsx")
!cp genome_data_with_transmission_levels.xlsx "drive/My Drive/Transformer-Viruses"
