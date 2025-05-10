# Phager

This is Phager - a novel machine learning based tool that recognizes phage contigs. 

# General information

Phager was trained on biological features of phage and non-phage genes. The script phager.py takes a raw FASTA file as input, performs gene calling, and generates biological features for each gene. These features are then combined into feature triplets, corresponding to gene triplets along a shifting reading frame. The gene feature triplets are then used by a LightGBM model to determine whether the input genome is a phage.
![Phager_img 001](https://github.com/user-attachments/assets/69961a06-c552-4315-a08c-b563a1dad561)



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
## References

Pyrodigal (Larralde, 2022), a Python library binding to Prodigal (Hyatt et al., 2010). (https://github.com/althonos/pyrodigal)

Sirén,K., Millard,A., Petersen,B., Gilbert,M.T.P., Clokie,M.R.J. and Sicheritz-Pontén,T. (2020) Rapid discovery of novel prophages using biological feature engineering and machine learning. 10.1101/2020.08.09.243022. https://www.biorxiv.org/content/10.1101/2020.08.09.243022v1.abstract

## Citation

