# Phager

Phager - a rapid phage contig predictor using biological feature based machine learning. 

### Authors
Anastasiya Gæde and Thomas Sicheritz-Pontén

# General information
Phager was trained on biological features of phage and non-phage genes. The script phager.py takes a raw or gzipped FASTA file as input, performs gene calling, and generates biological features for each gene. These features are then combined into feature triplets, corresponding to gene triplets along a shifting reading frame. The gene feature triplets are then used by a LightGBM model to determine whether the input genome is a phage.
<img width="6000" height="4200" alt="phager_overview" src="https://github.com/user-attachments/assets/1ade0697-97e8-49bf-8c09-d7d3beda2bfe" />




# Getting started 
## Requirements: 

- Python version 3.8 or later
- Biopython version 1.79 or later
- Pandas
- LighGBM


## Installation 
### From GitHub

```
conda create -n phager pandas
conda activate phager
pip install tqdm icecream colorama pyrodigal biopython more_itertools lightgbm scikit-learn pyarrow
git clone https://github.com/ku-cbd/phager.git
cd phager


```

## Example

```
python phager.py -h
python phager.py -a example/example_contigs.fasta -o myresults -v
```
## References

Pyrodigal (Larralde, 2022), a Python library binding to Prodigal (Hyatt et al., 2010). (https://github.com/althonos/pyrodigal)

Sirén,K., Millard,A., Petersen,B., Gilbert,M.T.P., Clokie,M.R.J. and Sicheritz-Pontén,T. (2020) Rapid discovery of novel prophages using biological feature engineering and machine learning. 10.1101/2020.08.09.243022. https://www.biorxiv.org/content/10.1101/2020.08.09.243022v1.abstract

## Citation
**Large-scale analysis of bacterial genomes reveals thousands of lytic phages**
***Alexander Perfilyev, Anastasiya Gæde, Steve Hooton, Sara A. Zahran, Panos G. Kalatzis, Caroline Sophie Winther-Have, Rodrigo Ibarra Chavez, Rachael C. Wilkinson, Anisha M. Thanki, Zhengjie Liu, Qing Zhang, Qianghua Lv, Yuqing Liu, Adriano M. Gigante, Robert J. Atterbury, Bent Petersen, Andrew D. Millard, Martha R. J. Clokie & Thomas Sicheritz-Pontén***
Nature Microbiology volume 11, pages42–52 (2026); doi: https://doi.org/10.1038/s41564-025-02203-4
