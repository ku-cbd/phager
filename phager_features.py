#!/usr/bin/env python3
# Created: Mon Jul  9 20:28:32 2018
# Last changed: Time-stamp: <Last changed 2025-05-12 17:23:26 by Thomas Sicheritz-Pontén, thomas>

import string, re
import os, sys, subprocess

import math
import gzip
import time
from datetime import timedelta
import itertools
import functools
import hashlib
import io
aa_fast_cache_max_size = os.getenv('AA_FAST_CACHE_MAX_SIZE')
if aa_fast_cache_max_size is not None:
    AA_FAST_CACHE_MAX_SIZE = int(aa_fast_cache_max_size)
else:
    AA_FAST_CACHE_MAX_SIZE = None


from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
sys.path.insert(0, os.path.dirname(__file__))
from CodonUsage import CodonAdaptationIndex, CodonsDict
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils import GC_skew, GC123
from Bio.Seq import Seq
def GC(seq):
    gc = sum(seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])
    try:
        return gc * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0

from zlib import compress
import pandas as pd
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

import numpy as np
from math import log
from more_itertools import chunked
from collections import namedtuple
from icecream import ic


def msg(*txts):
    global verbose
    if verbose:
        print(*txts, file=sys.stderr)

def guess_file_format(file):
    if any([file.endswith(suf) for suf in '.gbk.gz .gbk .gb.gz .gb'.split()]): return 'gbk'
    if any([file.endswith(suf) for suf in '.fasta.gz .fasta .fna.gz .fna .fas.gz .fas .fa.gz .fa'.split()]): return 'fasta'
    return None

def split_protein(seq, Nterm = 30, Cterm=30, Mminsize=50):
    Nseq = seq[:Nterm]
    Cseq = seq[-Cterm:]
    Mseq = seq[Nterm:-Cterm]
    if len(Mseq) < Mminsize:
        x = (len(seq)-Mminsize)//2
        if x <=0:
            Mseq = seq
        else:
            Mseq = seq[x:-x]
    return Nseq, Mseq, Cseq


def extract_location(entry):
    header, seq = entry
    rx = re.compile(r'([0-9]+)\.\.([0-9]+)')
    try:
        location = rx.search(header).group()
    except:
        location = None

    if location:
        fields = location.split('..')
        start, end = map(int, fields)
        strand = -1 if header.find('complement') > -1 else 1
        return (start, end, strand)
    try:
        start = int(re.search('start:([0-9]+)').group(1))
        m = re.search('end:([0-9]+)').group(0)
        end = int(m.group(1)) if m else start + len(seq)
        m = re.search('strand:([0-9]+)').group(0)
        strand = int(m.group(1)) if m else 1
    except:
        return (0,len(seq),1)

    
class myCodonAdaptationIndex(CodonAdaptationIndex):
    def __init__(self, verbose = False):
        CodonAdaptationIndex.__init__(self)
        self.verbose = verbose
        
    def _count_codons(self, sequences): 
        # make the codon dictionary local 
        self.codon_count = CodonsDict.copy() 
   
        # iterate over sequence and count all the codons in the FastaFile. 
        for seq in sequences:
            # make sure the sequence is lower case
            dna_sequence = str(seq).upper() if str(seq).islower() else str(seq)
            
            for i in range(0, len(dna_sequence), 3): 
                codon = dna_sequence[i:i + 3]
                if len(codon) != 3: continue
                if codon in self.codon_count: 
                    self.codon_count[codon] += 1 
                else:
                    if self.verbose:
                        print("illegal codon %s in gene" % (codon), dna_sequence, Seq(dna_sequence).translate(), file=sys.stderr)

        # fix for ZeroDivisionError
        self.codon_count = {c:n if n else 0.5 for c,n in self.codon_count.items()}
        
    def cai_for_gene(self, dna_sequence):
        # adapted from Biopython not raise an error on illegal codons
        cai_value, cai_length = 0, 0
        if self.index == {}: self.set_cai_index(SharpEcoliIndex)
        if dna_sequence.islower(): dna_sequence = dna_sequence.upper()

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i : i + 3]
            if codon in self.index:
                # these two codons are always one, exclude them:
                if codon not in ["ATG", "TGG"]:
                    cai_value += math.log(self.index[codon])
                    cai_length += 1
            # # some indices may not include stop codons:
            # elif codon not in ["TGA", "TAA", "TAG"]:
            #     raise TypeError(
            #         "illegal codon in sequence: %s.\n%s" % (codon, self.index)
            #     )

        try:
            return math.exp(cai_value / (cai_length - 1.0))
        except ZeroDivisionError:
            return 0
        
def calculate_CAI(entries):
    name, sequences = list(zip(*entries))
    CAI = myCodonAdaptationIndex()
    CAI.generate_index(sequences)    
    cais = []
    for seq in sequences:
        try:
            cai = CAI.cai_for_gene(str(seq))
        except TypeError:
            cai = 1   # add temporarily 1 for non-calculatable sequences (and change to median later)
            #            print cai, self.name
        cais.append(cai)

    # change the non-calculatable sequences from 1 to median score
    _cais = [x for x in cais if x < 1]
    if not len(_cais):
        return

    medium_score = sum(_cais)/len(_cais)
    return pd.Series([medium_score if x == 1 else x for x in cais], name='CAI')


def calculate_GC(entries):
    gc, gc1, gc2, gc3 = [], [], [], []
    name, sequences = list(zip(*entries))

    for seq in sequences:
        _gc, _gc1, _gc2, _gc3 = GC123(seq)
        gc.append(_gc)
        gc1.append(_gc1)
        gc2.append(_gc2)
        gc3.append(_gc3)

    df = pd.DataFrame([[x/100.0 for x in gc], [x/100.0 for x in gc1], [x/100.0 for x in gc2], [x/100.0 for x in gc3]]).T
    df.columns = "gc_content gc1_content gc2_content gc3_content".split()
    return df

                      
prosite = namedtuple("prosite", "AC ID hits phits pattern taxonomy")
prosite_patterns = [prosite(AC='PS00010', ID='ASX_HYDROXYL', hits=1775, phits=430, pattern=re.compile('C.[DN].{4}[FY].C.C'), taxonomy='eukaryota;eukaryota only'), prosite(AC='PS00012', ID='PHOSPHOPANTETHEINE', hits=1245, phits=1113, pattern=re.compile('[DEQGSTALMKRH][LIVMFYSTAC][GNQ][LIVMFYAG][DNEKHS]S[LIVMST][^PCFY][STAGCPQLIVMF][LIVMATN][DENQGTAKRHLM][LIVMWSTA][LIVGSTACR][^LPIY][^VY][LIVMFA]'), taxonomy='bacteria;eukaryota'), prosite(AC='PS00018', ID='EF_HAND_1', hits=3890, phits=2050, pattern=re.compile('D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]'), taxonomy='bacteria;bacteriophage;eukaryota;viruses'), prosite(AC='PS00022', ID='EGF_1', hits=3326, phits=962, pattern=re.compile('C.C.{2}[^V].{2}G[^C].C'), taxonomy='eukaryota;viruses'), prosite(AC='PS00027', ID='HOMEOBOX_1', hits=1340, phits=1334, pattern=re.compile('[LIVMFYG][ASLVR].{2}[LIVMSTACN].[LIVM][^Y].{2}[^L][LIV][RKNQESTAIY][LIVFSTNKH]W[FYVC].[NDQTAH].{5}[RKNAIMW]'), taxonomy='eukaryota;viruses'), prosite(AC='PS00028', ID='ZINC_FINGER_C2H2_1', hits=13396, phits=2324, pattern=re.compile('C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H'), taxonomy='archaea;eukaryota;viruses'), prosite(AC='PS00039', ID='DEAD_ATP_HELICASE', hits=1001, phits=1001, pattern=re.compile('[LIVMF]{2}DEAD[RKEN].[LIVMFYGSTN]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00055', ID='RIBOSOMAL_S12', hits=1018, phits=1018, pattern=re.compile('[RK].PNS[AR].R'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00086', ID='CYTOCHROME_P450', hits=1153, phits=1151, pattern=re.compile('[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00097', ID='CARBAMOYLTRANSFERASE', hits=1091, phits=1091, pattern=re.compile('F.[EK].S[GT]RT'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00101', ID='HEXAPEP_TRANSFERASES', hits=1579, phits=1224, pattern=re.compile('[LIV][GAED].{2}[STAV].[LIV].{3}[LIVAC].[LIV][GAED].{2}[STAVR].[LIV][GAED].{2}[STAV].[LIV].{3}[LIV]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00103', ID='PUR_PYR_PR_TRANSFER', hits=1303, phits=1303, pattern=re.compile('[LIVMFYWCTA][LIVM][LIVMA][LIVMFC][DE]D[LIVMS][LIVM][STAVD][STAR][GAC].[STAR]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00107', ID='PROTEIN_KINASE_ATP', hits=3356, phits=3332, pattern=re.compile('[LIV]G[^P]G[^P][FYWMGSTNH][SGA][^PW][LIVCAT][^PD].[GSTACLIVMFY].{5,18}[LIVMFYWCSTAR][AIVP][LIVMFAGCKR]K'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00108', ID='PROTEIN_KINASE_ST', hits=3099, phits=3068, pattern=re.compile('[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00142', ID='ZINC_PROTEASE', hits=1400, phits=1391, pattern=re.compile('[GSTALIVN][^PCHR][^KND]HE[LIVMFYW][^DEHRKP]H[^EKPC][LIVMFYWGSPQ]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00152', ID='ATPASE_ALPHA_BETA', hits=2288, phits=2288, pattern=re.compile('P[SAP][LIV][DNH][^LKGN][^F][^S]S[^DCPH]S'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00178', ID='AA_TRNA_LIGASE_I', hits=4385, phits=4382, pattern=re.compile('P.{0,2}[GSTAN][DENQGAPK].[LIVMFP][HT][LIVMYAC]G[HNTG][LIVMFYSTAGPC]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00198', ID='4FE4S_FER_1', hits=2570, phits=1422, pattern=re.compile('C.[^P]C[^C].C[^CP].[^C]C[PEG]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00211', ID='ABC_TRANSPORTER_1', hits=4196, phits=3887, pattern=re.compile('[LIVMFYC][SA][SAPGLVFYKQH]G[DENQMW][KRQASPCLIMFW][KRNQSTAVM][KRACLVM][LIVMFYPAN][^PHY][LIVMFW][SAGCLIVP][^FYWHP][^KRHP][LIVMFYWSTA]'), taxonomy='archaea;bacteria;bacteriophage;eukaryota;viruses'), prosite(AC='PS00232', ID='CADHERIN_1', hits=1502, phits=309, pattern=re.compile('[LIV].[LIV].D.ND[NH].P'), taxonomy='archaea;eukaryota'), prosite(AC='PS00237', ID='G_PROTEIN_RECEP_F1_1', hits=2145, phits=2141, pattern=re.compile('[GSTALIVMFYWC][GSTANCPDE][^EDPKRH].[^PQ][LIVMNQGA][^RK][^RK][LIVMFT][GSTANC][LIVMFYWSTAC][DENH]R[FYWCSH][^PE].[LIVM]'), taxonomy='eukaryota;viruses'), prosite(AC='PS00296', ID='CHAPERONINS_CPN60', hits=1047, phits=1047, pattern=re.compile('A[AS][^L][DEQ]E[^A][^Q][^R].GG[GA]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00297', ID='HSP70_1', hits=1198, phits=1198, pattern=re.compile('[IV]DLGT[ST].[SC]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00301', ID='G_TR_1', hits=3119, phits=3119, pattern=re.compile('D[KRSTGANQFYW].{3}E[KRAQ].[RKQD][GC][IVMK][ST][IV].{2}[GSTACKRNQ]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00323', ID='RIBOSOMAL_S19', hits=1011, phits=1011, pattern=re.compile('[STDNQ]G[KRNQMHSI].{6}[LIVM].{4}[LIVMC][GSD].{2}[LFI][GAS][DE][FYM].{2}[ST]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00329', ID='HSP70_2', hits=1232, phits=1232, pattern=re.compile('[LIVMF][LIVMFY][DN][LIVMFS]G[GSH][GS][AST].{3}[ST][LIVM][LIVMFC]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS00600', ID='AA_TRANSFER_CLASS_3', hits=1035, phits=1035, pattern=re.compile('[LIVMFYWCS][LIVMFYWCAH].D[ED][IVA].{2,3}[GAT][LIVMFAGCYN].{0,1}[RSACLIH].[GSADEHRM].{10,16}[DH][LIVMFCAG][LIVMFYSTAR].{2}[GSA]K.{2,3}[GSTADNV][GSAC]'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00632', ID='RIBOSOMAL_S4', hits=1044, phits=1044, pattern=re.compile('[LIVM][DERA].R[LI].{3}[LIVMC][VMFYHQL][KRTS].{3}[STAGCVF].[ST].{3}[SAI][KRQ].[LIVMF]{2}'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS00678', ID='WD_REPEATS_1', hits=3357, phits=1796, pattern=re.compile('[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]'), taxonomy='bacteria;eukaryota;viruses'), prosite(AC='PS00962', ID='RIBOSOMAL_S2_1', hits=1026, phits=1026, pattern=re.compile('[LIVMFA].[^GPRV][LIVMFYC]{2}[^LPC][STAC][GSTANQEKR][STALV][HY][LIVMF]G'), taxonomy='archaea;bacteria;eukaryota'), prosite(AC='PS01036', ID='HSP70_3', hits=1193, phits=1186, pattern=re.compile('[LIVMY].[LIVMF].GG.[ST][^LS][LIVM]P.[LIVM].[DEQKRSTA]'), taxonomy='archaea;bacteria;eukaryota;viruses'), prosite(AC='PS01186', ID='EGF_2', hits=3513, phits=1057, pattern=re.compile('C.C.{2}[GP][FYW].{4,8}C'), taxonomy='eukaryota;viruses'), prosite(AC='PS01187', ID='EGF_CA', hits=1697, phits=415, pattern=re.compile('[DEQN].[DEQN]{2}C.{3,14}C.{3,7}C.[DN].{4}[FY].C'), taxonomy='eukaryota;eukaryota only'), prosite(AC='PS01278', ID='MTTASE_RADICAL', hits=1059, phits=1059, pattern=re.compile('[LIVM].[LIVMT].{2}GC.{3}C[STAN][FY]C.[LIVMT].{4}G'), taxonomy='archaea;bacteria;eukaryota')]
murphy_10_tab = {'A': 'A', 'C': 'C', 'E': 'E', 'D': 'E', 'G': 'G', 'F': 'F', 'I': 'L', 'H': 'H', 'K': 'K', 'M': 'L', 'L': 'L', 'N': 'E', 'Q': 'E', 'P': 'P', 'S': 'S', 'R': 'K', 'T': 'S', 'W': 'F', 'V': 'L', 'Y': 'F'}
dayhoff_freq = {'A': 8.6, 'C': 2.9, 'E': 6.0, 'D': 5.5, 'G': 8.4, 'F': 3.6, 'I': 4.5, 'H': 2.0, 'K': 6.6, 'M': 1.7, 'L': 7.4, 'N': 4.3, 'Q': 3.9, 'P': 5.2, 'S': 7.0, 'R': 4.9, 'U': 0.1, 'T': 6.1, 'W': 1.3, 'V': 6.6, 'Y': 3.4}
property_residues = {'Polar': ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'], 'Aliphatic': ['A', 'I', 'L', 'V'], 'Aromatic': ['F', 'H', 'W', 'Y'], 'Basic': ['H', 'K', 'R'], 'Small': ['A', 'B', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'], 'Acidic': ['B', 'D', 'E', 'Z'], 'Charged': ['B', 'D', 'E', 'H', 'K', 'R', 'Z'], 'Tiny': ['A', 'C', 'G', 'S', 'T'], 'Non-polar': ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']}

"""50 years of amino acid hydrophobicity scales: revisiting the capacity for peptide classification"""
AAH_scale_patterns = {'dc-helix': ['ALAA', 'AALA', 'ALLE', 'ALLDA', 'AALAA', 'ALDAA', 'AAALA'], 'dc-random': ['GSSG', 'SSGS', 'SGSS', 'HHHH', 'EEEE', 'GSSGS', 'SGSSG', 'SSGSS', 'HHHHH'], 'dc-sheet': ['VLLV', 'VLVNA', 'SDTVV', 'KGTVT', 'YLVNM'], 'dd-helix': ['EELL', 'LEEL', 'LTEEE', 'LTLEE', 'ELLAD'], 'dd-random': ['GSSG', 'SSGS', 'SGSS', 'SSGL', 'GSSGS', 'SGSSG', 'SSGSS', 'GSSGL', 'TILPL'], 'dd-sheet': ['GEVV', 'PDGT', 'DGSV', 'TLDGG', 'SVIDT', 'LTVTG'], 'no-helix': ['SGSS', 'PDGS', 'GDSGG'], 'no-random': ['VVGI', 'QELD', 'VGIVT', 'TGHSL'], 'no-sheet': ['LEAL', 'SSGSS', 'SGSSG'], 'random': ['GSSG', 'GPSS', 'SGPS', 'SSGS', 'SGSS', 'SGPSS', 'GSSGS', 'SGSSG', 'SSGSS', 'HHHHH'], 's-sheet': ['VKVI', 'CGGSL', 'GIVSW', 'YGGVT'], 'krtm-helix': ['LGLL', 'VLLV', 'GIAL', 'LLVGI', 'LAAVA', 'FLAVL', 'YVFFG', 'YPIVW'], 'tm-helix': ['LLLL', 'LILL', 'LLLV', 'ILLL', 'LILLL', 'LLLLV'], 'krtm-sheet': ['SIGA', 'TGTLE'], 'tm-sheet': ['SGPL', 'SLNL', 'LYGG', 'PTLDL', 'LYGKV', 'SASAG', 'RQFNV'], 'all': ['VIGGG', 'IIGGG', 'LADAG', 'IVGAG', 'GVDVV'], 's-helix': ['EELKK']}


@functools.lru_cache(maxsize=AA_FAST_CACHE_MAX_SIZE)
def entropy(seq, is_protein = False):
    """
    Calculate Shannon entropy of sequence.
    Args:
        seq (str): Nucleotide sequence
    Examples:
        >>> sequtils.entropy('AGGATAAG')
        1.40
v        >>> sequtils.entropy('AAAACCGT')
        1.75
    """
    alphabet = list(dayhoff_freq.keys()) if is_protein else "ACGT" 
    cnt = [seq.count(i) for i in alphabet]
    d = sum(cnt)
    ent = []
    for i in [float(i)/d for i in cnt]:
        # round corner case that would cause math domain error
        if i == 0:
            i = 1
        ent.append(i * log(i, 2))
    return -1 * sum(ent)

def fix_name(name):
    name = name.split()[0]
    if name.startswith('sp|') and name.count('|') >= 2:
        name = name.split('|')[2]
    return name

def create_validate_amino_acids_seq(seq):
    return seq.replace('B', 'D').replace('Z','E').replace('J','I').replace('X','L').replace('U','C').replace('*','').replace('O','K')

def validate_amino_acids(seq):
    for aa in seq:
        if aa not in dayhoff_freq:
                return False
    return True

def prosite_matches(entries):
    m = []
    for name, seq in entries:
        m.append([int(len(pp.pattern.findall(seq))>0) for pp in prosite_patterns])
    
    df = pd.DataFrame(m)
    df.columns = ['PROSITE:%s' % x.ID for x in prosite_patterns]
    return df

@functools.lru_cache(maxsize=AA_FAST_CACHE_MAX_SIZE)
def AAH_scale_patterns_score_pep(seq, name, splitted=False):
    for pep in AAH_scale_patterns[name]:
        if seq.find(pep) >-1: return 1
    return 0

    
def AAH_scale_patterns_score_all_peps(entries, splitted=False):
    L = []
    features = []
    
    for name, peps in AAH_scale_patterns.items():
        L.append([AAH_scale_patterns_score_pep(x[1], name) for x in entries])
        features.append('AAH_%s' % name)
    df = pd.DataFrame(L, index=features).T
    #df.index = df[0]
    #df = df.drop(0, axis=1).T
    df.columns = features
    ADF = df
    return df
    mean = df.T.mean()
    mean.name = 'AAH_scale_peps'
    return mean

    
def biopython_proteinanalysis_seq_splitted(seq, scaling = True):
    Nseq, Mseq, Cseq = split_protein(seq)
    d = {}

    dW = biopython_proteinanalysis_seq(seq, scaling=scaling)
    d.update(dW)
    
    dN = biopython_proteinanalysis_seq(Nseq, scaling=scaling)
    dN = dict(zip(['%s_Nterm' % x for x in dN.keys()], list(dN.values())))
    d.update(dN)
    
    dM = biopython_proteinanalysis_seq(Mseq, scaling=scaling)
    dM = dict(zip(['%s_Mregion' % x for x in dM.keys()], list(dM.values())))    
    d.update(dM)
    
    dC = biopython_proteinanalysis_seq(Cseq, scaling=scaling)
    dC = dict(zip(['%s_Cterm' % x for x in dC.keys()], list(dC.values())))    
    d.update(dC)

    return d

@functools.lru_cache(maxsize=AA_FAST_CACHE_MAX_SIZE)
def biopython_proteinanalysis_seq(seq, scaling=False):
    res = ProteinAnalysis(seq)
    d = {}
    d['length'] = res.length
    #d['monoisotopic'] = int(res.monoisotopic)

    d['molecular_weight'] = res.molecular_weight()
    #d['aromaticity'] = res.aromaticity() # skip, as we will have it in the property_residues
    d['instability_index'] = res.instability_index()
    d['isoelectric_point'] = res.isoelectric_point()
    r, c = res.molar_extinction_coefficient()
    d['molar_extinction_coefficient_reduced'], d['molar_extinction_coefficient_cysteines'] = r, c

    if scaling:
        d['molecular_weight'] = d['molecular_weight'] / 100000
        d['instability_index'] = d['instability_index'] / 100
        d['isoelectric_point'] = d['isoelectric_point'] / 10
        d['molar_extinction_coefficient_reduced'], d['molar_extinction_coefficient_cysteines'] = d['molar_extinction_coefficient_reduced']/ 100000.0, d['molar_extinction_coefficient_cysteines'] / 100000.0


    d['percent_helix_naive'],d['percent_turn_naive'],d['percent_strand_naive']   = res.secondary_structure_fraction()

    aap = res.get_amino_acids_percent()
    aas = sorted(aap.keys())
    d.update({'percent:%s' % aa:aap[aa] for aa in aas})

    flex = np.array(res.flexibility())
    d['flex:min'], d['flex:max'], d['flex:std'] = flex.min(), flex.max(), flex.std()
    d['gravy'] = res.gravy()

    d.update({'prop_res_%s' % key:sum([aap.get(x, 0) for x in value]) for key, value in list(property_residues.items())})

    return d

def biopython_proteinanalysis(entries, scaling=False, splitted=False):
    m = []
        
    for name, seq in entries:
        if splitted:
            d = biopython_proteinanalysis_seq_splitted(seq, scaling=scaling)
        else:
            d = biopython_proteinanalysis_seq(seq, scaling=scaling)
        m.append([d[x] for x in list(d.keys())])

    df = pd.DataFrame(m)
    df.columns = ["BP:" + x for x in list(d.keys())]
    return df

@functools.lru_cache(maxsize=AA_FAST_CACHE_MAX_SIZE)
def oligopetide_frequencies(seq, n=2, alphabet="ACEDGFIHKMLNQPSRTWVY"):
    L = float(len(seq))
    nmers = sorted([''.join(x) for x in list(itertools.product(list(alphabet), repeat=n))])
    return nmers, [seq.count(x)/L for x in nmers]


def dipeptide_frequencies(entries, prefix="DIPEP", alphabet="ACEDGFIHKMLNQPSRTWVY"):
    m = []
    for n, (name, seq) in enumerate(entries):
        nmers, freqs = oligopetide_frequencies(seq, n=2, alphabet=alphabet)
        values = freqs
        m.append(values)

    df = pd.DataFrame(m)
    df.columns = ['%s:%s' % (prefix, x) for x in nmers]
    return df

def tripeptide_frequencies(entries, prefix = "TRIPEP", alphabet="ACEDGFIHKMLNQPSRTWVY"):
    m = []
    for n, (name, seq) in enumerate(entries):
        nmers, freqs = oligopetide_frequencies(seq, n=3, alphabet=alphabet)
        values = freqs
        m.append(values)

    df = pd.DataFrame(m)
    df.columns = ['%s:%s' % (prefix, x) for x in nmers]        
    return df

HTS = {'V': 'H', 'I': 'H', 'Y': 'H', 'F': 'H', 'W': 'H', 'L': 'S', 'N': 'T', 'P': 'T', 'G': 'T', 'S': 'T', 'E': 'S', 'M': 'S', 'A': 'S'}
def count_HTS(s, char, min_n = 4, max_n=31):
    d = {}
    for i in range(min_n,max_n): d[i] = 0
    for i in range(min_n,max_n):
        match = '%s{%d}[%s]+' % (char, (i-1), char)
        found = len(re.findall(match, s))
        if found == 0:
            break
        d[i] = found
    return d

@functools.lru_cache(maxsize=AA_FAST_CACHE_MAX_SIZE)
def HTS_score(seq):
    s = ''.join([HTS.get(x,'.') for x in seq])
    dH = count_HTS(s, 'H', min_n = 4, max_n = 31)
    dT = count_HTS(s, 'T', min_n = 2, max_n = 65)
    dS = count_HTS(s, 'S', min_n = 4, max_n = 31)
    values = [*[dH[x] for x in sorted(dH.keys())], *[dT[x] for x in sorted(dT.keys())], *[dS[x] for x in sorted(dS.keys())]]
    header = [*['HTS-H:%d' % x for x in sorted(dH.keys())], *['HTS-T:%d' % x for x in sorted(dT.keys())], *['HTS-S:%d' % x for x in sorted(dS.keys())]]
    
    return (values, header)

def calculate_HTS(entries, splitted=False):
    name, sequences = list(zip(*entries))
    L = []
    for seq in sequences:
        values, header = HTS_score(seq)
        L.append(values)
    df = pd.DataFrame(L, columns=header)
    return df


def calculate_GCskew(entries):
    name, sequences = list(zip(*entries))
    total_seq = ''.join(sequences)
    L = len(total_seq)
    window = L//len(sequences)
    gc_skew_values = GC_skew(total_seq, window=window)
    #print(L, len(sequences), len(gc_skew_values), len(total_seq), window)
    return pd.Series(gc_skew_values[:len(sequences)], name='gc_skew')

def calculcate_melting_temperature(entries):
    name, sequences = list(zip(*entries))    
    return pd.Series([Tm_NN(x, check=True)/100.0 for x in sequences], name='melting_temp')

def calculcate_gene_lengths(entries, locations):
    # OBS log transformed
    name, sequences = list(zip(*entries))
    s = [locations[x+1][0]-locations[x][1] for x in range(len(sequences)-1)]
    s.append(np.mean(s))
    s = pd.Series(s)
    diffs = np.sign(s)*np.log(abs(s)+0.001)/10

    df = pd.DataFrame([[log(len(x))/14.0 for x in sequences], [x[2] for x in locations], list(diffs)]).T
    df.columns = "gene_lengths gene_strands gene_diffs".split()
    return  df

def calcultate_compressebility(entries):
    name, sequences = list(zip(*entries))        
    # OBS normalised to 50
    return pd.Series([len(compress(x[:50].encode('ASCII')))/50.0 for x in sequences], name='compressebility_minsize_50')

def calcultate_compressebility_frac(entries):
    name, sequences = list(zip(*entries))        
    return pd.Series([len(compress(x.encode('ASCII')))/len(x) for x in sequences], name='compressebility_frac')

def calculate_entropy(entries, is_protein=False):
    name, sequences = list(zip(*entries))        
    return pd.Series([entropy(x, is_protein=is_protein)/2.0 for x in sequences], name='entropyAA' if is_protein else 'entropy')

def RunAAandDNA(DNA_entries, AA_entries, locations, skip_pep_freq=False, verbose = False):
    #global DF, DF_AA, DF_DNA
    DF_DNA = RunDNA(DNA_entries, locations, verbose = verbose)
    DF_AA = Run(AA_entries, verbose = verbose, skip_pep_freq=skip_pep_freq)
    DF_AA.drop('length', axis=1, inplace=True, errors='ignore')
    DF = pd.concat([DF_AA, DF_DNA], axis=1)
    return DF, DF_AA, DF_DNA

def RunDNA(entries, locations=None, verbose = False):
    names, sequences = list(zip(*entries))
    df_CAI = calculate_CAI(entries)
    df_GC = calculate_GC(entries)
    df_GC_skew = calculate_GCskew(entries)    
    df_melting_temp = calculcate_melting_temperature(entries)
    df_compressebility = calcultate_compressebility(entries)
    df_compressebility_frac = calcultate_compressebility_frac(entries)
    df_entropy = calculate_entropy(entries)
    
    to_concat = [df_CAI, df_melting_temp, df_GC, df_GC_skew, df_compressebility, df_compressebility_frac, df_entropy]
    if locations is not None:
        df_gene_length = calculcate_gene_lengths(entries, locations)
        to_concat.insert(2, df_gene_length)
    DF = pd.concat(to_concat, axis = 1)
    DF.index = names
    return DF

def Run(entries, scaling=False, verbose = False, splitted=False, skip_pep_freq=False):
    return RunAA(entries, scaling=scaling, verbose = verbose, splitted=splitted, skip_pep_freq=skip_pep_freq)
    
def RunAA(entries, scaling=False, verbose = False, splitted=False, skip_pep_freq=False):
    #global df_entropy, df_aah_scale_peps, df_prosite
    
    def msg(txt):
        if verbose:
            print(txt, file=sys.stderr)

    msg(f"Starting on {len(entries):_} sequences")
    reduced_entries = [(name,''.join([murphy_10_tab.get(x,x) for x in seq])) for name,seq in entries]

    names, sequences = list(zip(*entries))
    
    msg('Calculating prosite matches')
    df_prosite = prosite_matches(entries)

    msg('Calculating Shannon entropy')
    df_entropy = calculate_entropy(entries, is_protein=True)

    msg('Calculating AAH scale peptide occurences')
    df_aah_scale_peps = AAH_scale_patterns_score_all_peps(entries, splitted=splitted)

    msg('Calculating HTS occurences')
    df_HTS = calculate_HTS(entries, splitted=splitted)
    
    msg('Calculating %sbiopython ProteinAnalysis' % ('splitted ' if splitted else ''))
    df_biopython = biopython_proteinanalysis(entries, scaling=scaling, splitted=splitted)

    # msg('Calculating CTD')
    # df_CTD = calculate_CTD(entries)

    if not skip_pep_freq:
        msg("Calculating dipeptide_frequencies")
        df_dipeptide_frequencies = dipeptide_frequencies(entries, prefix="DIPEP")

        reduced_alphabet = ''.join(set(murphy_10_tab.values()))
        msg("Calculating reduced dipeptide_frequencies")
        df_reduced_dipeptide_frequencies = dipeptide_frequencies(reduced_entries, alphabet = reduced_alphabet, prefix="RED_DIPEP")

        msg("Calculating reduced tripeptide_frequencies")
        df_reduced_tripeptide_frequencies = tripeptide_frequencies(reduced_entries, alphabet = reduced_alphabet, prefix="RED_TRIPEP")
        msg("Concatenating")
        DF = pd.concat([
            df_prosite,
            df_entropy,
            df_aah_scale_peps,
            df_HTS,
            df_biopython,
#            df_CTD,
            df_dipeptide_frequencies,
            df_reduced_dipeptide_frequencies,
            df_reduced_tripeptide_frequencies,
        ], axis = 1)
    else:
        msg("Concatenating")
        DF = pd.concat([
            df_prosite,
            df_entropy,
            df_aah_scale_peps,
            df_HTS,
            df_biopython,
 #           df_CTD,
        ], axis = 1)

    #DF.insert(0, 'md5', [hashlib.md5(name.encode('utf-8')).hexdigest() for name, seq in entries])
    DF.insert(0, 'md5', [hashlib.md5(seq.encode('utf-8')).hexdigest() for name, seq in entries])
    DF.index = names
    
    return DF


def dna_aa_df_to_features(var, skip_pep_freq=False):
    if type(var) == pd.DataFrame:
        df = var
    else:
        df = pd.read_parquet(var)
    df.dna = df.dna.str.upper()
    coding_capacity, length, Nratio, GC = df['coding_capacity length Nratio GC'.split()].iloc[0].tolist()
    DF, DF_AA, DF_DNA = RunAAandDNA(df['gene dna'.split()].values, df['gene aa'.split()].values, locations=None, skip_pep_freq=skip_pep_freq)
    DF['genome_coding_capacity'] = coding_capacity
    DF['genome_length'] = length
    DF['genome_GC'] = GC
    DF['genome_Nratio'] = Nratio
    return DF

def call_genes(file, outdir='.', overwrite=False):
    import pyrodigal
    seqs_df_file = f"{outdir}/{os.path.basename(file).split('.fasta')[0]}.dna_aa.pa"
    if os.path.exists(seqs_df_file) and not overwrite: return pd.read_parquet(seqs_df_file)
    ic("working on", file)
    handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
    entries = list(SimpleFastaParser(handle))
    name = entries[0][0]
    meta = True
    p = pyrodigal.GeneFinder(meta=meta)
    length = sum([len(seq) for name, seq in entries])
    Nratio = sum([seq.count('N') for name, seq in entries])/length
    genes = []
    for name, seq in entries:
        name = name.split()[0]
        for i, gene in enumerate(p.find_genes(seq)):
            start = gene.begin
            end = gene.end
            strand = gene.strand
            dna_sequence = str(Seq(seq[(start-1):end]).reverse_complement()) if strand == -1 else str(seq)[(start-1):end]
            aa_sequence = create_validate_amino_acids_seq(str(gene.translate()).strip('*'))
            gene_name = f'{name}_{i+1}_{start}_{end}'
            genes.append([gene_name, start, end, strand, dna_sequence, aa_sequence])

    df = pd.DataFrame(genes, columns = 'gene start end strand dna aa'.split())
    coding_capacity = df.aa.str.len().sum()*3/length
    df['coding_capacity'] = coding_capacity
    df['Nratio'] = Nratio
    df['length'] = length
    df['GC'] = GC(''.join([seq for name, seq in entries]))
    ic("Saving to", seqs_df_file)
    df.to_parquet(seqs_df_file)
    return df

def create_triplets_from_fna(file, outfile=None, outdir='.', genome = None, min_size=5, dont_use = None, return_df=False):
    if not return_df and os.path.exists(outfile): return
    if not dont_use: dont_use = []

    df1 = call_genes(file, outdir=outdir)
    if df1.empty or df1.Nratio.mean() >= 1: return
    dfF = dna_aa_df_to_features(df1, skip_pep_freq=True)
    df = df1.merge(dfF, left_on='gene', right_index=True)

    if len(df) < min_size: return
    
    df = df.drop('AAseq md5'.split(), axis=1, errors='ignore')
    df['tmp_acc'] = df.gene
    cols = [x for x in df.columns if x.find('PEP') == -1 and x.find('PROSITE') == -1]
    df2 = df[cols]
    df1 = df2.shift(1).rename({x:f"{x}-LEFT" for x in df.columns.tolist()}, axis=1)
    df3 = df2.shift(-1).rename({x:f"{x}-RIGHT" for x in df.columns.tolist()}, axis=1)
    df2.columns = [f"{x}-MID" for x in df2.columns]
    dfT = pd.concat([df1, df2, df3], axis=1)
    dfT = dfT.iloc[1:-1]
    dfT = dfT[~(dfT['tmp_acc-LEFT'].isin(dont_use) | dfT['tmp_acc-MID'].isin(dont_use) | dfT['tmp_acc-RIGHT'].isin(dont_use))]
    cols = sorted(dfT.columns)
    dfT = dfT[cols]
    dfT.columns = [x.replace(':','_') for x in dfT.columns]
    #dfT.index = dfT['tmp_acc-MID']
    if return_df: return dfT
    
    dfT.to_parquet(outfile)
    
