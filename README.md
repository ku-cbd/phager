# Phager

This is Phager - a novel machine learning based tool that recognizes phage contigs. 

# General information

Phager was trained on biological features of phage and non-phage genes. The script phager.py takes a raw FASTA file as input, performs gene calling, and generates biological features for each gene. These features are then combined into feature triplets, corresponding to gene triplets along a shifting reading frame. The gene feature triplets are then used by a LightGBM model to determine whether the input genome is a phage.
![Phager_scheme  001](https://github.com/user-attachments/assets/85381369-d621-45f6-b9d4-678c86db1f62)


# Getting started 
## Requirements: 

- Python version 3.8 or later
- Biopython version 1.79 or later

## Installation 
### From GitHub

```
conda create -y -n Phager-env python=3.8 
conda activate Phager-env 
git clone git@github.com:ku-cbd/Phager.git 
```

## Example

```
phager.py -h
phager.py -c 10000 -a contigs.fa -d myresults -v
```

## Citation
