#!/usr/bin/env python3
# Created: Wed Jul 21 17:55:13 2021
# Last changed: Time-stamp: <Last changed 2025-05-11 12:40:22 by Thomas Sicheritz-Pontén, thomas>

import string, re
import os, sys

ENVIRONMENT = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
               "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS")

for e in ENVIRONMENT:
    os.environ[e] = '1'

sys.path.insert(0, '.')
import glob, gzip
import pandas as pd
import numpy as np
import time
from icecream import ic
from tqdm import tqdm
import pickle
import shutil
import logging
from colorama import Fore, Back, Style
PHAGER_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, PHAGER_DIR)

from phager_features import create_triplets_from_fna
from Bio.SeqIO.FastaIO import SimpleFastaParser

def secure_filename(name):
    import unicodedata
    name = unicodedata.normalize('NFKD', str(name)).encode('ascii', 'ignore').decode('ascii')
    name = re.sub('[^a-zA-Z0-9_.-]','_', name)
    return name

def exists(file): return os.path.exists(file)

            
class Phager:
    def __init__(self, outdir, model_file, features_dir='.', min_length=10, verbose = True, plot_predictions=False, logger=None, prefix='', threads=1):
        self.__version__ = '2.1'
        self.verbose = verbose
        if not logger:
            loggfile = os.path.join(outdir, 'phager.log')
            logging.basicConfig(filename=loggfile, level=logging.INFO, filemode='w', format='%(asctime)s %(message)s', datefmt='%m/%d/%Y%I:%M:%S %p')

            logger=logging.getLogger() 
        self.logger = logger

        self.threads = threads
        self.outdir = outdir
        self.features_dir = features_dir
        self.plot_predictions = plot_predictions
        self.__files_to_delete = []
        self.__files_to_compress = []
        self.__dirs_to_delete = []
        self.model_file = model_file
        self.last_status = ''
        self.prefix = prefix
        self.min_length = min_length
        #self.frac_th = frac_th
        self.load_model()
        
    def exit(self):
        ic(f"Cleaning up {len(self.__files_to_delete)} files and {len(self.__dirs_to_delete)} directories")
        for file in self.__files_to_delete:
            if os.path.exists(file):
                os.remove(file)
        for dir in self.__dirs_to_delete:
            if os.path.exists(dir):
                os.rmdir(dir)
        ic(f"Cleaning up compressing {len(self.__files_to_delete)}")
        for file in self.__files_to_compress:
            if os.path.exists(file):
                os.remove(f"gzip {file}")

    def load_model(self):
        ic("Loading model from", self.model_file)
        self.model = pickle.load(open(self.model_file, 'rb'))
        self.model.set_params(n_jobs=1)
        self.threshold = self.model.optimal_threshold
        self.frac_th = self.model.optimal_frac
        ic(f"Got {self.threshold=} and {self.frac_th=}")
        self.features = list(self.model.feature_name_)
        
    def runDNA(self, file):
        self.fasta_file = file
        self.genome_length = 0
        self.base = os.path.basename(self.fasta_file).split('.fasta')[0].split('.fsa')[0].split('.fna')[0]
        ic("Creating features for", file)
        try:
            df = create_triplets_from_fna(file, outdir=self.features_dir, return_df=True)
        except ValueError:
            return (self.base, 0, -1, -1, -1)
        if df is None: return (self.base, 0, -1, -1, -1)
        self.genome_length = df['genome_length-MID'].unique()[0]
        df = self.engineer_features(df)
        return self.predict(df)

    def runAA(self, file, entries=None, name = None, spades=True):
        if file:
            self.base = os.path.basename(file).split('.faa')[0].split('.fasta')[0].split('.fsa')[0].split('.fna')[0]
            handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
            entries = list(SimpleFastaParser(handle))
        else:
            assert len(entries) >0
            self.base = name
        if len(entries) < 3:
            self.last_status = f"Too few genes ({len(entries)}) in {file}"
            return "failed"
            
        self.faa_to_features(entries=entries, spades=spades)
        self.make_triplets()
        self.engineer_features()
        
    def engineer_features(self, df):
        ic("Engineering more features")
        df = df.drop([x for x in df.columns if x.startswith('HTS')], axis=1)
        df = df.drop([x for x in df.columns if x.find('Nratio') >-1], axis=1)
        df = df.drop([x for x in df.columns if x.find('genome_length') >-1], axis=1)
        df = df.drop(['aa-LEFT', 'aa-MID', 'aa-RIGHT', 'dna-LEFT', 'dna-MID', 'dna-RIGHT', 'gene-LEFT', 'gene-MID', 'gene-RIGHT', 'tmp_acc-LEFT', 'tmp_acc-MID', 'tmp_acc-RIGHT'], axis=1, errors='ignore')

        df['inter12'] = df['start-MID'] - df['end-LEFT']
        df['inter23'] = df['start-RIGHT'] - df['end-MID']

        df['length_inter'] = df['inter12'] + df['inter23']
        df['mean_inter'] = df['length_inter']/2
        df['length_triplet'] = df['end-RIGHT'] - df['start-LEFT']
        df['ratio_genic_intergenic'] = df['length_triplet'] / (df['inter12'] + df['inter23'] + 0.0001)

        df = df.drop([x for x in df.columns if x.find('end') >-1], axis=1)
        df = df.drop([x for x in df.columns if x.find('start') >-1], axis=1)
        to_drop = ['genome_GC-LEFT','genome_GC-RIGHT', 'GC-LEFT',
                   'GC-MID', 'GC-RIGHT',
                   'genome_coding_capacity-LEFT','genome_coding_capacity-RIGHT',
                   'coding_capacity-LEFT',
                   'coding_capacity-MID','coding_capacity-RIGHT']

        df = df.drop(to_drop, axis=1)

        features = [x for x in df.columns if x not in 'name cloud label source partition'.split()]
        base_features = [x.replace('-LEFT','') for x in features if x.endswith('-LEFT')]

        ic("Features log/sqrt/exp")
        for feature in features:
            s = {}
            s[f"{feature}_log"] = np.log(df[feature].mask(df[feature] <=0)).fillna(0).astype(np.float32)
            s[f"{feature}_sqrt"] = np.sqrt(df[feature].mask(df[feature] <=0)).fillna(0).astype(np.float32)

            if feature.find('molar_extinction_coefficient') >-1: continue
            if feature.find('weight') >-1: continue
            if feature.find('length') >-1: continue
            if feature.find('inter') >-1: continue

            s[f"{feature}_exp"] = np.exp(df[feature]).astype(np.float32)

        ic("concatenating ...")
        df = pd.concat([df, pd.DataFrame(s)], axis=1)
        del s

        s = {}
        ic("features mean/min/max")
        for base_feature in base_features:
            s[f"mean_{base_feature}"] = ((df[f"{base_feature}-LEFT"] + df[f"{base_feature}-MID"] + df[f"{base_feature}-RIGHT"])/3).astype(np.float32)
            s[f"max_{base_feature}"] = (df[[f"{base_feature}-LEFT", f"{base_feature}-MID", f"{base_feature}-RIGHT"]].max(axis=1)).astype(np.float32)
            s[f"min_{base_feature}"] = (df[[f"{base_feature}-LEFT", f"{base_feature}-MID", f"{base_feature}-RIGHT"]].min(axis=1)).astype(np.float32)
            s[f"sum13_{base_feature}"] = (df[f"{base_feature}-LEFT"] + df[f"{base_feature}-RIGHT"]).astype(np.float32)
            s[f"diff13_{base_feature}"] = (df[f"{base_feature}-LEFT"] - df[f"{base_feature}-RIGHT"]).abs().astype(np.float32)        
            s[f"multiply_{base_feature}"] = (df[f"{base_feature}-LEFT"] * df[f"{base_feature}-MID"] * df[f"{base_feature}-RIGHT"]).astype(np.float32)

        ic("concatenating ...")
        df = pd.concat([df, pd.DataFrame(s)], axis=1).reset_index(drop=True)
        del s

        ic("downcasting ...")
        df = df.astype({x:'float32' for x in df.select_dtypes(include=['float64', 'float16']).columns})
        if 'label' in df.columns: df.label = df.label.astype('int8')
        return df

    def predict_genes(self, df):
        ic("Predicting genes")
        y_pred = self.model.predict_proba(df[self.features])
        self.y_pred = y_pred
        self.predictions = [int(x[1]>=self.threshold) for x in y_pred]
        
    def ascii_plot(self, phage_frac):
        print('PPP PLOT', round(phage_frac,2), self.base, ''.join(['.' if x == 0 else 'x' for x in self.predictions]), flush=True)

    def plot(self):
        import plotly.express as px
        a = [x[1] for x in self.y_pred]
        fig = px.line(y=a, x=range(len(a)), range_y=(0,1), title=f'Phager prediction for {self.base}')
        outfile = f"{self.outdir}/{self.base}_phager.html"
        fig.write_html(outfile)
        print(outfile)

    def predict(self, df):
        self.predict_genes(df)
        predictions = self.predictions
        L = len(predictions)
        #if L < self.min_length: return (self.base, L, -1, 0, self.genome_length)
        phage_frac = sum(self.predictions)/L
        is_phage = self.decide_if_phage()
        return (self.base, L, is_phage, phage_frac, self.genome_length)

    def decide_if_phage(self, min_value=0.1, min_value_th = 0.035):
        min_length = self.min_length

        L = len(self.predictions)
        phage_frac = sum(self.predictions)/L
        is_phage = int(phage_frac >= self.frac_th)

        #print(f"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  {is_phage=} {phage_frac=} {self.frac_th=} {phage_frac >= self.frac_th}")

        color = Fore.GREEN if is_phage else Fore.RED

        what = {1:'Phage', 0:'Bacteria', -1:'TooShort'}
        if len(self.predictions) < min_length:
            if color == Fore.GREEN: color = Fore.BLUE
            is_phage = -1
            
        min_n_value = len([x[1] for x in self.y_pred if x[1] <= min_value])

        ##############################################################################
        # the fraction of really bad predictions - will override lower frac_th values
        #
        if min_n_value/L > min_value_th:
            color = Fore.RED
            is_phage = 0

        if self.plot_predictions:
            blocks = " ▁▂▃▄▅▆▇████"
            plot1 = ''.join(['.' if x == 0 else 'x' for x in self.predictions])
            plot2 = ''.join([blocks[int(x[1]*10)] for x in self.y_pred])
            txt =  "\n"
            #txt += f"{color}{self.base} {round(phage_frac,4)} {what[is_phage]} ({len(self.y_pred)} genes){Style.RESET_ALL} {self.frac_th=} {self.threshold=}\n"
            txt += f"{color}{self.base} {round(phage_frac,4)} {what[is_phage]} ({len(self.y_pred)} genes){Style.RESET_ALL}\n"
            #txt += f"{color}{self.base}{Style.RESET_ALL} {phage_frac} {min_n_value}\n"
            txt += f"{Fore.LIGHTBLACK_EX}{plot1}{Style.RESET_ALL}\n"
            txt += f"{color}{plot2}{Style.RESET_ALL}\n"
            #txt += f"{self.base} {what[is_phage]} {round(phage_frac,2)}\n"
            txt += "="*80

            self.logger.info(txt)
            if self.verbose: print(txt, flush=True)

        return is_phage

    def cleanup(self, to_compress=None):
        self.__files_to_delete.append(self.eng_feature_file)
        self.__files_to_delete.append(self.feature_file)
        self.__files_to_delete.append(self.triplets_file)
        if to_compress:
            self.__files_to_compress.extend(to_compress)
        self.exit()

def Run(_file):
    file, outdir, model_file, features_dir, plot_predictions, verbose, logger, n = _file
    phager = Phager(outdir, model_file, features_dir=features_dir, verbose = verbose, logger=logger, plot_predictions=plot_predictions, prefix=str(n))
    return phager.runDNA(file)
    

def predict_on_spades_assembly(file, outfile, outdir, model_file, logger, threads=1,
                               min_length=5000, verbose=False, save_intermediate_files=False,
                               extract_phage_genomes=False, plot_predictions=False):
 
    timer_start = time.time()
    base = os.path.basename(file).split('.fa')[0]

    logger.info(f"Started Phager with: {' '.join(sys.argv)}")
    fasta_dir = f"{outdir}/fastas"
    os.makedirs(fasta_dir, exist_ok=True)

    fasta_files = glob.glob(f'{outdir}/fastas/*.fasta.gz')

    if not fasta_files:
        handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
        entries = list(SimpleFastaParser(handle))

        fasta_files = []
        for name, seq in tqdm(entries):
            if len(seq) < min_length: continue
            name = name.split()[0]
            _file = f"{fasta_dir}/{secure_filename(name)}.fasta.gz"
            if _file in fasta_files:
                print(f"Duplicated sequence name - dropping {name}")
                continue
            fasta_files.append(_file)
            if not os.path.exists(_file):
                with gzip.open(_file, 'wt') as fid:
                    print(f">{name}\n{seq}", file=fid)

    logger.info(f"Found {len(entries)} contigs in {file}, using {len(fasta_files)} contigs with min size >= {min_length:_}bp")
    features_dir = f"{outdir}/features"
    os.makedirs(features_dir, exist_ok=True)
            
    if extract_phage_genomes:
        outdir = f"{outdir}/{base}.predicted_phages"
        if not os.path.exists(outdir): os.mkdir(outdir)

    # ugly multiprocessing hack
    if threads >1:
        from multiprocessing import Pool, cpu_count

        files = [(file, outdir, model_file, features_dir, plot_predictions, verbose, logger, n) for n,file in enumerate(fasta_files)]
        num_cores = threads
        pool = Pool(processes=num_cores)
        results = [_ for _ in tqdm(pool.imap_unordered(Run, files), total=len(files))]
    else:
        phager = Phager(outdir, model_file, features_dir=features_dir, logger=logger, verbose = verbose, plot_predictions=plot_predictions)
        results = [phager.runDNA(file) for file in tqdm(fasta_files)]

        # dont use this yet, need to move extract function before the cleanup
        #if not save_intermediate_files: phager.cleanup()
            
    run_time = time.time() - timer_start
    
    dfRes = pd.DataFrame(results, columns="name genes is_phage score bp".split())
    dfRes.score = dfRes.score.round(2)
    dfRes.index.name = f"contig (runtime:{round(run_time):_} seconds)"
    dfRes.to_csv(outfile, sep='\t')

    return base, outfile


def extract_phage_genomes(phage_results, fasta_files, outdir, base, min_genes=20):
    global df, entries
    outdir = f"{outdir}/{base}.predicted_phages"
    if not os.path.exists(outdir): os.mkdir(outdir)

    entries = {}
    for file in fasta_files:
        handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
        entries.update({name.split()[0]:seq for name,seq in SimpleFastaParser(handle)})

    df = pd.read_table(phage_results, index_col=0)
    df = df[(df.is_phage==1) & (df.genes >= min_genes)]
    
    for contig, genes, is_phage, score, bp in df.values:
        name = re.sub('[^a-zA-Z0-9_.]','_', contig.split()[0])
        outfile = f"{outdir}/{name}_score_{score}_ngenes_{genes}.fasta.gz"
        ic("Writing", outfile)
        if os.path.exists(outfile): continue
        fid = gzip.open(outfile, 'wt')
        seq = entries.get(contig)
        if not seq:
            seq = entries.get(contig + '.fasta')
        if not seq:
            seq = entries.get(secure_filename(contig))
        print(f">{contig}\tscore={score}\tgenes={genes}\n{seq}", file=fid)
        fid.close()

def sanity_check_results_numbers(fasta_files, contig_min_size, outdir, result_file):
    handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
    n_entries = len([x for x in SimpleFastaParser(handle) if len(x[1]) >= contig_min_size])
    n_generated_fasta_files = len(glob.glob(f"{outdir}/fastas/*.fasta.gz"))
    n_generated_features = len(glob.glob(f"{outdir}/features/*dna_aa.pa"))
    n_results = len(pd.read_table(result_file))
    correct = len(set([n_entries, n_generated_fasta_files, n_generated_features, n_results])) == 1 and n_results >0
    color = Fore.GREEN if correct else Fore.RED
    txt = f"{color}CHECK {'OK' if correct else 'NOTOK'} {[n_entries, n_generated_fasta_files, n_generated_features, n_results]}{Style.RESET_ALL}"
    logger.info(txt)
    if verbose: print(txt, flush=True)

    
def create_logger(logfile, web=False):
    if not web:
        logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w', format='%(asctime)s %(message)s', datefmt='%m/%d/%Y%I:%M:%S %p')
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s<br>', datefmt='%m/%d/%Y%I:%M:%S %p')
    logger=logging.getLogger()
    return logger


if __name__ == '__main__':
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    usage = "%prog [options] file (or - for stdin)"
    parser = ArgumentParser(usage, formatter_class=RawDescriptionHelpFormatter, description="""Phager predicts if a contig is phage like.
The most frequent way to run phager is with an assembly file.\n\n
example:
    phager2.py -c 10000 -a contigs.fa -d myresults -v

Author: Thomas Sicheritz-Pontén, thomassp@sund.ku.dk
    """)
    
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-a", "--assembly", action="store_true", default=False)
    parser.add_argument("-d", "-o", "--outdir", action="store", default=None)
    parser.add_argument("-l", "--listfile", action="store", default=None)
    parser.add_argument("-t", "--threads", action="store", type=int, default=1)
    parser.add_argument("-c", "--contig_min_size", action="store", type=int, default=5_000)
    parser.add_argument("-V", "--view", action="store_true", default=False)
    parser.add_argument("-x", "--extract", action="store_true", default=False)
    parser.add_argument("-m", "--extract_min_genes", action="store", default=10)
    parser.add_argument("-s", "--save_intermediate_files", action="store_true", default=False)
    parser.add_argument("-P", "--plot", action="store_true", default=True)
    parser.add_argument("-M", "--model", action="store", type=str, default=None)
    parser.add_argument("-w", "--web", action="store_true", default=False)
    parser.add_argument("-C", "--check", action="store_true", default=False)

    args, nargs = parser.parse_known_args()

    ############################################################
    # Current model
    model_file = os.path.join(PHAGER_DIR, 'current_phager_model.pkl') if not args.model else args.model

    web = args.web
    verbose = 0 if web else args.verbose
    if not verbose: ic.disable()
    threads = args.threads
    run_on_assembly = args.assembly
    contig_min_size = args.contig_min_size
    outdir = args.outdir
    view_only = args.view
    save_intermediate_files = args.save_intermediate_files
    plot_predictions = args.plot
    check_number_of_results = args.check
    
    if args.listfile:
        fasta_files = [x.strip() for x in open(args.listfile).readlines()]
    else:
        fasta_files = nargs

    if not fasta_files:
        parser.print_help()
        sys.exit(0)
        
    if not outdir:
        outdir = os.path.basename(fasta_files[0]).split('.faa')[0].split('.spades_contigs.fasta')[0].split('.contigs.fasta')[0].split('.fasta')[0]
    outdir = secure_filename(outdir)
    if not exists(outdir): os.mkdir(outdir)

    if run_on_assembly:
        ic("Running on assembly")
        for file in fasta_files:
            base = os.path.basename(file).split('.fa')[0].split('.fasta')[0].split('.fna')[0]
            outfile = f'{outdir}/{base}.phager_results.csv.gz' if not web else None
            logfile = f'{outdir}/{base}.phager.log' if not web else None
            logger = create_logger(logfile, web)
            base, result_file = predict_on_spades_assembly(file, outfile, outdir, model_file, logger, min_length=contig_min_size, verbose=verbose, plot_predictions=plot_predictions, threads=threads)
            if args.extract:
                extract_phage_genomes(result_file, fasta_files, outdir, base, args.extract_min_genes)
        if check_number_of_results: sanity_check_results_numbers(fasta_files, contig_min_size, outdir, result_file)

        sys.exit(0)
        
    not_finished = [file for file in fasta_files if not os.path.exists(f"{outdir}/{os.path.basename(file).split('.fa')[0]}.phager_results.csv.gz")]
    fasta_files, fasta_files_orig = not_finished, fasta_files

    phager = Phager(outdir, model_file, verbose=verbose, plot_predictions=plot_predictions)
    results = []
    for fasta_file in fasta_files:
        results.append(phager.runAA(fasta_file)) if fasta_file.endswith('.faa.gz') else results.append(phager.runDNA(fasta_file))
    dfR = pd.DataFrame(results, columns = 'name size is_phage score bp'.split())
    dfR.to_csv(f'{outdir}/phager_results.csv.gz', sep='\t', index=False)
        
