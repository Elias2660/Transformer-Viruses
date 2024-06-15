#!/usr/bin/env python
# coding: utf-8

# ## Mount Google Drive so that we can access the genomic sequence data

# In[ ]:


import os

# from google.colab import drive # for the google colab
# drive.mount('/content/drive')

# now with os.listdir("/content/drive/...") we can get the list of files in any directory in our Google Drive

# ## Load our dataset into memory

# In[ ]:


import pandas as pd
virus_data = pd.read_excel("data/Woolhouse_and_Brierley_RNA_virus_database.xlsx")

# In[ ]:


virus_data.head()

# In[ ]:


genome_data = pd.DataFrame(columns=["Locus", "Position/Length Indicator?", "Virus Name", "Genome"])

# In[ ]:


index = 0
stor = []
with open("data/humanVirusSeq.fa") as d:
    lines = d.readlines()
    print(len(lines))
    for line in lines:
        if index % 2 == 0:
            stor = line.strip().replace("\n","").replace(">","").split("|")
            index += 1
        else:
            genome_data.loc[len(genome_data.index)] = [stor[0], stor[1], stor[2], line.strip()]
            index += 1

# In[ ]:


genome_data

# In[ ]:


len(genome_data.iloc[0]["Genome"].strip())

# In[ ]:


# assuming that we have a CSV with rows like this:
# genomic_data_filepath.txt, R0_value
import csv

list_of_genomic_filepaths = []
list_of_genomic_sequences = []
list_of_r0_values = []

with open('{CSV FILE HERE}', 'r') as csv_file:
  csvreader = csv.reader(csv_file)
  for row in csvreader:
    # add the filepaths to list_of_genomic_filepaths
    # add the strings of our genomic file sequences to the list_of_genomic_sequences
    # add the r0 values to list_of_r0_values

# ## Tokenize our sequences

# In[ ]:


# define a vocabulary, this should be a mapping from all of the possible characters
# that can appear in our sequence
# we can either do this manually by defining a dictionary and writing {'A': 0, 'G': 1, 'T': 2, ...}
# or you can use torchtext.vocab.build_vocab_from_iterator (this is probably easier)

sample_sequence = list_of_genomic_file_contents[0]
vocab = None # write code here
print("Vocab is: {}".format(vocab)) # sanity check that all of the possible values are reflected

# now we builder our tokenizer
tokenizer = torch.classes.torchtext.Tokenizer(vocab=vocab, split='character')
sample_tokens = tokenizer(sample_sequence) # sanity check that it looks good

# ## Install some necessary libraries

# In[ ]:


!pip install torch

# ## Create our dataset

# In[ ]:


import torch
from torch.utils.data import Dataset

# Pytorch defines a nice dataset class that only requires we implement two functions:
# 1. __len__
# 2. __getitem__
# Our dataset is just the set of all samples that we want to train our model on
# __len__ should get us the total number of samples
# __getitem__ takes in an integer and should give us the corresponding element in the dataset

# You might be wondering why they have those weird underscores? That enables us to
# call those functions on the class instances directly
#
# For instance if we define an instance of our Dataset class, then we can get the length of it
# by invoking `len` on it directly
# new_dataset = GenomicR0ValueDataset(some_sample_sequences, some_sample_r0_values, new_tokenizer)
# len(new_dataset) # this will call GenomicR0ValueDataset.__len__

class GenomicR0ValueDataset(Dataset):
  def __init__(self, sequences, r0_values, tokenizer):
    # here we should initialize some class variables using `self.` so that
    # we can access them further along
    # make sure all of our sequences are tokenized so that when we return them in __getitem__
    # we don't have to do any post-processing on them during training

  def __len__(self):
    # this needs to return the number of elements in our dataset
    # so this should be the total number of sequences
    pass

  def __getitem__(self, index):
    # get item should take in an index and return the corresponding genomic sequence
    # AND its r0 value. we need to return both because every time we pass a sequence
    # through our model, we need to compare its predicted r0 value to the correct r0
    # value. so the format of this output should be returning two things like this:
    # return _, _
    pass

# ## Create our dataloader

# In[ ]:


import torch
from torch.utils.data import random_split, DataLoader
from torchvision import datasets, transforms
from torch import nn, optim

# pytorch dataloaders: https://pytorch.org/tutorials/beginner/basics/data_tutorial.html
# dataloaders take in a pytorch.Dataset as an argument (like the one we defined above!)

dataset = None # lets instantiate that GenomicR0ValueDataset we defined above

# define some train test split
train_test_ratio = 0.9
train_size = int(train_test_ratio * len(dataset)) # how convenient that we can called `len` on our dataset!!
test_size = len(dataset) - train_size

# now use random_split from torch.utils.data to define our two datasets
train_dataset, test_dataset = None, None

# define our batch_size
batch_size = 4 # batch size defines how many sequences our model will be processing at once
# higher batch sizes mean training will be faster, but will make updates slightly less precisely (will talk about this in person)
# another thing to keep in mind is the amount of MEMORY that we have! our batches can't get too big!

train_dataloader = DataLoader(train_dataset, batch_size=batch_size)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size)


# ## Define our model architecture

# ### Positional Encoding

# In[ ]:


# positional encoding explanation
# https://machinelearningmastery.com/a-gentle-introduction-to-positional-encoding-in-transformer-models-part-1/

import torch
import torch.nn as nn

class PositionalEncoding(nn.Module):

    def __init__(self, d_model: int, max_len: int = 5000):
        super().__init__()
        self.dropout = nn.Dropout(p=0.1)

        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
        pe = torch.zeros(max_len, 1, d_model)
        pe[:, 0, 0::2] = torch.sin(position * div_term)
        pe[:, 0, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        """
        Arguments:
            x: Tensor, shape ``[seq_len, batch_size, embedding_dim]``
        """
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)


# positional encoding is critical for the transformer model
# d_model is the embedding dimension of the model
# max_len here is our context window size

# In[ ]:


import math

import torch
from torch import nn, Tensor
from torch.nn import TransformerEncoder, TransformerEncoderLayer


class GenomeR0ValueModel(nn.Module):
    def __init__(self, ntoken: int, d_model: int, nhead: int, d_hid: int,
                 nlayers: int):
        super().__init__()
        self.model_type = 'Transformer'
        self.pos_encoder = PositionalEncoding(d_model)
        encoder_layers = TransformerEncoderLayer(d_model, nhead, d_hid)
        self.transformer_encoder = TransformerEncoder(encoder_layers, nlayers)
        self.embedding = nn.Embedding(ntoken, d_model)
        self.d_model = d_model
        self.linear = nn.Linear(d_model, ntoken)

        self.init_weights()

    def init_weights(self) -> None:
        initrange = 0.1
        self.embedding.weight.data.uniform_(-initrange, initrange)
        self.linear.bias.data.zero_()
        self.linear.weight.data.uniform_(-initrange, initrange)

    def forward(self, src: Tensor, src_mask: Tensor = None) -> Tensor:
        """
        Arguments:
            src: Tensor, shape ``[seq_len, batch_size]``
            src_mask: Tensor, shape ``[seq_len, seq_len]``

        Returns:
            output Tensor of shape ``[seq_len, batch_size, ntoken]``
        """
        src = self.embedding(src) * math.sqrt(self.d_model)
        src = self.pos_encoder(src)
        if src_mask is None:
            """Generate a square causal mask for the sequence. The masked positions are filled with float('-inf').
            Unmasked positions are filled with float(0.0).
            """
            src_mask = nn.Transformer.generate_square_subsequent_mask(len(src)).to(device)
        output = self.transformer_encoder(src, src_mask)
        output = self.linear(output)
        return output


# here, we define the model class, it is based on a transformer architecture
# and it takes in a ntoken, which is the size of the vocabulary (number of unique characters in our input)
# d_model, which is the embedding size of the tokens
# nhead is the number of heads in our self-attention setup
# d_hid is the dimension of the hidden layer
# n_layers is the number of hidden layers

# ## Create our train and test loop and begin trainin

# In[ ]:


# initialize our model
import torch
import torch.nn as nn

model = None

import torch
from torch.utils.data import DataLoader
from tqdm import tqdm  # Progress bar (optional, but very helpful)

def train_model(model, train_loader, test_loader, criterion, optimizer, num_epochs, device="cuda" if torch.cuda.is_available() else "cpu"):
  # move the
  model.to(device)
  best_val_loss = float("inf")
  for epoch in range(num_epochs):
    # Training Phase
    model.train()
    train_loss = 0.0
    pbar = tqdm(train_loader, desc=f"Epoch {epoch + 1}/{num_epochs}")
    for inputs, targets in pbar:
      inputs, targets = inputs.to(device), targets.to(device)

        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs, targets.unsqueeze(1))  # Ensure target shape matches output
        loss.backward()
        optimizer.step()

        train_loss += loss.item()
        pbar.set_postfix({"Train Loss": loss.item()})

    train_loss /= len(train_loader)

    # Test Phase
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
      for inputs, targets in val_loader:
        inputs, targets = inputs.to(device), targets.to(device)
        outputs = model(inputs)
        loss = criterion(outputs, targets.unsqueeze(1))
        val_loss += loss.item()
      val_loss /= len(val_loader)

        print(f"Epoch {epoch + 1}/{num_epochs}, Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}")

        # Save the best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), "best_model.pth")


# Example usage
train_loader = DataLoader(...)  # Your training dataloader (we defined these above)
val_loader = DataLoader(...)   # Your validation dataloader

criterion = torch.nn.MSELoss()
optimizer = torch.optim.Adam(learning_rate=0.1)

# Train the model
train_model(model, train_loader, val_loader, criterion, optimizer, num_epochs=10)





# ## Begin training

# In[ ]:




