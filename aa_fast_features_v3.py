#!/usr/bin/env python3
# Created: Mon Jul  9 20:28:32 2018
# Last changed: Time-stamp: <Last changed 2024-11-07 21:26:19 by Thomas Sicheritz-Pontén, thomas>

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
    rx = re.compile('([0-9]+)\.\.([0-9]+)')
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

dfCTD = pd.read_csv(io.StringIO("""Category,Group1,Group2,Group3
Hydrophobicity_ARGP820101,QSTNGDE,RAHCKMV,LYPFIW
Hydrophobicity_CASG920101,KDEQPSRNTG,AHYMLV,FIWC
Hydrophobicity_ENGD860101,RDKENQHYP,SGTAW,CVLIMF
Hydrophobicity_FASG890101,KERSQD,NTPG,AYHWVMFLIC
Hydrophobicity_PONP930101,KPDESNQT,GRHA,YMFWLCVI
Hydrophobicity_PRAM900101,RKEDQN,GASTPHY,CLVIMFW
Hydrophobicity_ZIMJ680101,QNGSWTDERA,HMCKV,LPFYI
Normalized van der Waals Volume,GASTPDC,NVEQIL,MHKFRYW
Polarity,LIFWCMVY,PATGS,HQRKNED
Polarizability,GASDT,CPNVEQIL,KMHFRYW
Charge,KR,ANCQGHILMFPSTWYV,DE
Secondary structure,EALMQKRH,VIYCWFT,GNPSD
Solvent accessibility,ALFCGIVW,RKQEND,MSPTHY"""))
CTDgroup1 = {dfCTD['Category'][i]: dfCTD['Group1'][i] for i in range(dfCTD.shape[0])}
CTDgroup2 = {dfCTD['Category'][i]: dfCTD['Group2'][i] for i in range(dfCTD.shape[0])}
CTDgroup3 = {dfCTD['Category'][i]: dfCTD['Group3'][i] for i in range(dfCTD.shape[0])}
CTDcategories = dfCTD.Category.tolist()
CTDdesc = [cat+'-G{}'.format(i) for cat in CTDcategories for i in range(1,4)]

def calculate_CTD(entries):
    # compute CTD composition
    L = []

    for name, seq in entries:
        tmp = [name]
        for cat in CTDcategories:
            g1 = sum([seq.count(aa) for aa in CTDgroup1[cat]])/len(seq)
            g2 = sum([seq.count(aa) for aa in CTDgroup2[cat]])/len(seq)
            g3 = sum([seq.count(aa) for aa in CTDgroup3[cat]])/len(seq)
            tmp.extend([g1, g2, g3])
        L.append(tmp)

    df = pd.DataFrame(L, columns = ['name'] + CTDdesc).set_index('name')
    return df


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

def seq2murphy10(seq):
    return ''.join([murphy_10_tab.get(x,x) for x in seq])

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

def run_function_splitted(seq, function, **args):
    Nseq, Mseq, Cseq = split_protein(seq)
    d = {}

    dW = function(seq, args)
    d.update(dW)
    
    dN = function(Nseq, args)
    dN = dict(zip(['%s_Nterm' % x for x in dN.keys()], list(dN.values())))
    d.update(dN)
    
    dM = function(Mseq, args)
    dM = dict(zip(['%s_Mregion' % x for x in dM.keys()], list(dM.values())))    
    d.update(dM)
    
    dC = function(Cseq, args)
    dC = dict(zip(['%s_Cterm' % x for x in dC.keys()], list(dC.values())))    
    d.update(dC)

    return d
    
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

def RunDNAandAA(DNA_entries, AA_entries, locations, verbose = False):
    return RunAAandDNA(DNA_entries, AA_entries, locations, verbose)

def RunReads(entries, verbose=False):
    names, sequences = list(zip(*entries))
    df_CAI = calculate_CAI(entries)
    df_GC = calculate_GC(entries)
    df_GC_skew = calculate_GCskew(entries)    
    df_melting_temp = calculcate_melting_temperature(entries)
    df_compressebility = calcultate_compressebility(entries)
    df_compressebility_frac = calcultate_compressebility_frac(entries)
    df_entropy = calculate_entropy(entries)

    DF = pd.concat([
        df_CAI,
        df_melting_temp,
        df_GC,
        df_GC_skew,
        df_compressebility,
        df_compressebility_frac,
        df_entropy,
        ], axis = 1)
    DF.index = names
    return DF

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

def seq2features(entries, min_length = 10, scaling=False, verbose = False):
    # entries is a list of (name, sequence) tuples
    # unpack as: names, sequences = zip(*entries)
    entries = [(x[0], create_validate_amino_acids_seq(x[1])) for x in entries]
    entries = [x for x in entries if x[0].find('pseudogene') == -1]
    entries = [x for x in entries if len(x[1]) >= min_length]
    lengths = [len(x[1]) for x in entries]

    if len(entries) == 0: return None
    if verbose:
        print("Converting %d sequences - min=%d max=%d" % (len(entries), min(lengths), max(lengths)))
    return RunAA(entries, scaling=scaling, verbose=verbose)

def seq_to_AA_DNA_features(AA_entries, DNA_entries, min_AA_length = 10, scaling=False, verbose = False):
    # entries is a list of (name, sequence) tuples
    # unpack as: names, sequences = zip(*entries)
    entries = [(x[0], create_validate_amino_acids_seq(x[1])) for x in entries]
    entries = [x for x in entries if x[0].find('pseudogene') == -1]
    entries = [x for x in entries if len(x[1]) >= min_length]
    lengths = [len(x[1]) for x in entries]
    if verbose:
        print("Converting %d sequences - min=%d max=%d" % (len(entries), min(lengths), max(lengths)))
    return RunAA(entries, scaling=scaling, verbose=verbose)

def seqfile2features(file, save = True, scaling=False, verbose = False, hdf_file = None, min_length = 10):
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    if not hdf_file:
        hdf_file = file + '.aa_features.hdf'
        
    if os.path.exists(hdf_file):
        df = pd.read_hdf(hdf_file)
    else:
        handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
        entries = list(SimpleFastaParser(handle))

    if len(entries) == 0:
        print("No sequences ... exiting", file=sys.stderr)
        return None
        
        df = seq2features(entries, scaling=scaling, verbose = verbose, min_length=min_length)
        if save:
            df.to_hdf(hdf_file, key = "df", complevel=1)
    return df

def DNA2features(file, min_length = 100, verbose=False):
    handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
    for name, seq in SimpleFastaParser(handle):
        if len(seq) < min_length :continue
        pass

def AA2fixed_rep(entries, return_cols = False):
    from scipy.stats import variation
    
    df = seq2features(entries, scaling=False)
    
    L = []
    L.append(df.count())
    L.append(df.mean())
    L.append(df.std())
    L.append(df.min())
    for q in [0.25, 0.5, 0.75]:
        L.append(df.quantile(q))
    L.append(df.max())
    
    ds = pd.concat(L, axis=1)
    ds.columns = 'count  mean   std   min   25%   50%   75%   max'.split()

    ds['kurtosis'] = df.kurtosis()
    ds['skew'] = df.skew()
    ds['mad'] = df.mad()
    ds['var'] = df.var()
    ds['variation'] = pd.Series(variation(df)).fillna(0).values

    s = pd.Series(np.concatenate([ds.loc[x].values for x in ds.index]))

    if return_cols:
        calcs = list(ds.columns)
        features = list(ds.index)
        cols = []
        [cols.extend([f +':' + c for c in calcs]) for f in features]
        return s, cols
        
    return s
 
def gbk2genes(gbfile, min_length=10):
    base = os.path.basename(gbfile).split('gbk')[0]
    handle = gzip.open(gbfile, 'rt') if gbfile.endswith('.gz') else open(gbfile)
    recordIter = SeqIO.parse(handle, 'genbank')
    genome_seq = None
    for n, record in enumerate(recordIter):
        accession = record.annotations['accessions'][-1] if 'accessions' in record.annotations else base
        organism = record.annotations['organism'] if 'organism' in record.annotations else base
        taxonomy = ';'.join(record.annotations.get('taxonomy',[]))
        year = int(record.annotations['date'][-4:]) if 'date' in record.annotations else None 
        names, dnas, aas, starts, ends, strands = [], [], [], [], [], []
        genome_seq = str(record.seq)
        for feature in record.features:
            if feature.type == 'CDS':
                fq = feature.qualifiers
                if not 'translation' in fq: continue
                translation = fq['translation'][0]
                if len(translation) < min_length: continue
                cds = str(feature.location.extract(record).seq)
                dnas.append(cds)
                aas.append(translation)
                starts.append(feature.location.start)
                ends.append(feature.location.end)
                strands.append(feature.location.strand)
                if 'protein_id' in fq:
                    name = fq['protein_id'][0]
                else:
                    name = f'{accession}_{i+1}_{start}_{end}'
                names.append(name)

    seq_length = len(genome_seq)
    df = pd.DataFrame(list(zip(*[names, starts, ends, strands, dnas, aas])), columns = 'gene start end strand dna aa'.split())

    df.aa = [create_validate_amino_acids_seq(aa_seq) for aa_seq in df.aa.tolist()]

    df['coding_capacity'] = df.aa.str.len().sum()*3/seq_length 
    df['Nratio'] = genome_seq.count('N')/seq_length
    df['length'] = seq_length
    df['GC'] = GC(genome_seq)
    df.attrs.update({'name':(base, organism), 'organism':organism, 'taxonomy':taxonomy, 'file':gbfile, 'accession':accession})
    if year: df['year'] = year
    return df

 
def fasta2genes(file, prefix='orf', min_length = 10, dna_entries=None, aa_entries=None, locations=None):
    ic("working on", file)
    base = os.path.basename(file).split('.fasta')[0].split('.fna')[0]
    handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
    entries = list(SimpleFastaParser(handle))
    if len(entries) == 0: return None
    name = entries[0][0]
    seq = ''.join([x[1] for x in entries]).upper()
    seq_length = sum([len(x[1]) for x in entries])
    Nratio = seq.count('N')/seq_length
    genome_seq = seq
    if not dna_entries:
        import pyrodigal

        meta = len(seq) < 100_000
        p = pyrodigal.GeneFinder(meta=meta)
        ic(f'running prdigal on {len(seq):_}nt')
        if not meta: p.train(seq)
        a = p.find_genes(seq)
        genes = []

        #dna_entries. aa_entries, locations = [], [], []
        for i, gene in enumerate(p.find_genes(seq)):
            start = gene.begin
            end = gene.end
            strand = gene.strand
            dna_sequence = str(Seq(seq[(start-1):end]).reverse_complement()) if strand == -1 else str(seq)[(start-1):end]
            aa_sequence = str(gene.translate()).strip('*')
            if len(aa_sequence) < min_length: continue
            ann = 'unknown'
            gene_name = f'{prefix}_{i+1}_{start}_{end}'
            genes.append([gene_name, start, end, strand, dna_sequence, aa_sequence])
        del p
        df = pd.DataFrame(genes, columns = 'gene start end strand dna aa'.split())
    else:
        names = [x[0] for x in dna_entries]
        starts = [x[0] for x in locations]
        ends = [x[1] for x in locations]
        strands = [x[2] for x in locations]
        dnas = [x[1] for x in dna_entries]
        aas = [x[1] for x in aa_entries]
        df = pd.DataFrame(list(zip(*[names, start, ends, strands, dnas, aas])), columns = 'gene start end strand dna aa'.split())

    df.aa = [create_validate_amino_acids_seq(aa_seq) for aa_seq in df.aa.tolist()]

    df['coding_capacity'] = df.aa.str.len().sum()*3/seq_length 
    df['Nratio'] = Nratio
    df['length'] = seq_length
    df['GC'] = GC(''.join([seq for name, seq in entries]))
    df.attrs.update({'name':base, 'file':file})
    return df

def create_intergenic_features(df):
    regions = (df.start - df.end.shift(1)).iloc[1:].values.tolist()
    L = len(regions)
    d = {}
    d['mean_length_intergenic_region'] = np.mean(regions)
    d['mean_length_pos_intergenic_region'] = np.mean([x for x in regions if x>=0])
    d['mean_length_neg_intergenic_region'] = np.mean([x for x in regions if x<0])
    d['n_pos_intergenic_regions_frac'] = len([x for x in regions if x>0])/L
    d['n_pos_intergenic_regions_frac'] = len([x for x in regions if x<0])/L
    d['n_zero_intergenic_regions_frac'] = len([x for x in regions if x ==0])/L
    swapped = (((df.strand == df.strand.shift(1)) == 0)*1).sum()
    d['strandness'] = swapped / len(df)   # TODO encode mean length of proteins on each swapped strand
    return d
    
def create_desc_features(df, features=None):
    # TODO use featuretools to combine genome and genes
    if not features: features = df.columns
    d = {}
    for feature in features:
        d[f"mean_{feature}"] = df[feature].mean()
        d[f"max_{feature}"] = df[feature].max()
        d[f"min_{feature}"] = df[feature].min()
        d[f"std_{feature}"] = df[feature].std()
        d[f"var_{feature}"] = df[feature].var()
        d[f"kurtosis_{feature}"] = df[feature].kurtosis()
        d[f"skew_{feature}"] = df[feature].skew()
        d[f"prod_{feature}"] = df[feature].prod()
    return d

def add_eng_base_features_inplace(df, features):
    for feature in features:
        df[f"cyclic_sin_{feature}"] = np.sin(2 * np.pi * df[feature])
        df[f"cyclic_cos_{feature}"] = np.cos(2 * np.pi * df[feature])


def create_genomic_features(file, dna_entries=None, aa_entries=None,
                            locations=None, file_format='auto', save=False, overwrite=False,
                            no_geneDF=False, skip_pep_freq=True, add_ranked=False, min_nr_genes =10):
    ic("working on", file)
    if file_format == 'auto': file_format = guess_file_format(file)
    if file_format == 'fasta':
        df = fasta2genes(file, dna_entries=dna_entries, aa_entries=aa_entries, locations=locations)
    elif file_format  == 'gbk': # uses the gbk translations instead of running prodigal
        df = gbk2genes(file)
    else:
        raise ValueError("Unknown file format")

    if len(df) < min_nr_genes:
        return pd.DataFrame() if no_geneDF else (pd.DataFrame(), pd.DataFrame())
    
    coding_capacity, length, Nratio, GC = df['coding_capacity length Nratio GC'.split()].iloc[0].tolist()
    DF, DF_AA, DF_DNA = RunAAandDNA(df['gene dna'.split()].values,
                                    df['gene aa'.split()].values,
                                    locations=df['start end strand'.split()].values, skip_pep_freq=skip_pep_freq)
    gene_feature_cols = [x for x in DF.select_dtypes(include='number') if x not in 'genome_size coding_capacity Nratio length GC'.split()]
    add_eng_base_features_inplace(DF, gene_feature_cols)
    genome_cols = ['coding_capacity', 'length', 'GC']
    d = df[genome_cols].iloc[0].to_dict()
    d['name'] = df.attrs.get('name', 'unknown_name')
    d.update(create_intergenic_features(df))
    gene_feature_cyclic_cols = [x for x in DF.columns if x.startswith('cyclic_')]
    features = gene_feature_cols + gene_feature_cyclic_cols
    d.update(create_desc_features(DF, features))
    genomic_features = df.from_dict(d, orient='index').T.set_index('name')
    genomic_features = genomic_features.astype(np.float32)
    #genomic_features = genomic_features.astype(np.float32).round(4)
    genomic_features.length = genomic_features.length.astype(np.int32)
    genomic_features.attrs.update(df.attrs)

    if no_geneDF: return genomic_features

    DF['complete_genome_coding_capacity'] = coding_capacity
    DF['complete_genome_length'] = length
    DF['complete_genome_GC'] = GC
    DF['complete_genome_Nratio'] = Nratio
    DF.attrs.update(df.attrs)
    if add_ranked:
        DF_ranked = rank_proteins(DF)
        DF = DF.merge(DF_ranked, left_index=True, right_index=True)
    return genomic_features, DF

def rank_proteins(df):
    features = df.select_dtypes(include='number').columns.tolist()
    features = [x for x in features if (not x.startswith('cyclic_')) and (not x.startswith('complete_genome'))]
    df_ranked = df[features].rank(method='dense', ascending=False).astype(int)
    df_ranked.columns = df_ranked.columns + '_ranked'
    df_ranked_pct = df[features].rank(method='dense', ascending=False, pct=True)
    df_ranked_pct.columns = df_ranked.columns + '_ranked_pct'
    return pd.concat([df_ranked, df_ranked_pct], axis=1)

def create_genomic_features_for_files(files, skip_pep_freq=True, add_ranked=False):
    from tqdm import tqdm
    df = pd.concat([create_genomic_features(file, skip_pep_freq=skip_pep_freq, no_geneDF=True, add_ranked=add_ranked) for file in tqdm(files, desc='Creating genomic features')])
    return df


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


def create_features(file, overwrite=False, return_df=False):
    outfile = file.split('.fasta')[0] + '.dna_aa_features.pa'
    if not overwrite and os.path.exists(outfile): return
    df1 = call_genes(file)
    if df1.empty or df1.Nratio.mean() >= 1: return
    dfF = dna_aa_df_to_features(df1, skip_pep_freq=True)
    df = df1.merge(dfF, left_on='gene', right_index=True)
    if return_df:
        return df

    df.to_parquet(outfile)
        

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


if __name__ == '__main__':
    from optparse import OptionParser
    usage = "%prog [options] file (or - for stdin)\n"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default = 0)
    parser.add_option("-S", "--add_sequences", action="store_true", dest="add_sequences", default = False)
    parser.add_option("-L", "--add_locations", action="store_true", dest="add_locations", default = False)
    parser.add_option("-D", "--DNA", action="store_true", dest="input_is_DNA", default = False)
    parser.add_option("-G", "--genomic", action="store_true", default = False)
    parser.add_option("--DNA-only", action="store_true", dest="DNA_only", default = False)
    parser.add_option("--split", action="store_true", dest="split", default = False)        
    parser.add_option("-m", "--min_length", action="store", type="int", dest="min_length", default = 10)
    parser.add_option("-f", "--format", action="store", type = "choice", dest="format", choices = ("hdf","csv", "csv.gz", 'pa'), default = 'pa')
    parser.add_option("-P", "--skip_pep_freq", action="store_true", default = 0)
    parser.add_option("-o", "--outfile", action="store", type=str, default = '')
    parser.add_option("-T", "--triplets", action="store", type=str, default = False)

    
    (options, args) = parser.parse_args()
    verbose = options.verbose
    min_length = options.min_length
    input_is_DNA = options.input_is_DNA
    use_genomic = options.genomic
    splitted = options.split
    skip_pep_freq = options.skip_pep_freq
    save_file = options.outfile
    triplets = options.triplets
    
    fasta_file = args[0]
    save_file = save_file if len(save_file) >1 else os.path.basename(fasta_file).replace('.gz','.splitted' if splitted else '') + '.' + options.format
    if os.path.exists(save_file): sys.exit(0)
        
    timer_start = time.time()

    if triplets:
        outfile = save_file.replace('.pa', '.trp.pa')
        create_triplets_from_fna(fasta_file, outfile)
    if use_genomic:
        df = create_genomic_features_for_files(args, skip_pep_freq=skip_pep_freq)
        df.to_parquet(save_file)
        sys.exit(0)

    fid = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)

    failed_seqs = 0
    msg("Reading", fasta_file)
    handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)
    entries = list(SimpleFastaParser(handle))

    entries = [x for x in entries if x[0].find('pseudogene') == -1]
    
    if input_is_DNA:
        msg("Running on DNA")
        DNA_entries = []
        for n, (name, seq) in enumerate(entries):
            seq = seq.strip('*').upper()
            if len(seq) < min_length*3: continue

            DNA_entries.append((name, seq))

        msg("Failed sequences", len(entries) - len(DNA_entries))

        locations = [extract_location(x) for x in DNA_entries]
        DNA_entries = [(x[0], x[1]) for x in DNA_entries]
        
        if options.DNA_only:
            DF = RunDNA(DNA_entries, verbose=verbose)

        else:

            AA_entries = [(x[0], str(Seq(x[1][:(len(x[1])//3)*3]).translate().strip('*'))) for x in DNA_entries]
            AA_entries = [(x[0], create_validate_amino_acids_seq(x[1])) for x in AA_entries]

            DF, DF_AA, DF_DNA = RunDNAandAA(DNA_entries, AA_entries, locations, verbose = verbose)
            if options.add_sequences: DF.insert(0, 'DNAseq', [x[1] for x in DNA_entries])
            if options.add_locations:
                DF.insert(0, 'gene_start', [x[0] for x in locations])
                DF.insert(0, 'gene_end', [x[1] for x in locations])
                DF.insert(0, 'gene_strand', [x[2] for x in locations])
            
    else:
        msg("Running on AA")
        AA_entries = []
        for n, (name, seq) in enumerate(entries):
            seq = seq.strip('*').upper()
            if len(seq) < min_length: continue
            seq = create_validate_amino_acids_seq(seq)
            #name = fix_name(name)
            AA_entries.append((name, seq))

        msg("Failed sequences", len(entries) - len(AA_entries))
        DF = Run(AA_entries, verbose=verbose, splitted=splitted, skip_pep_freq=skip_pep_freq)

    DF.index.name = 'name'
    if options.add_sequences: DF.insert(0, 'AAseq', [x[1] for x in AA_entries])
    runtime = time.time() - timer_start
    
    msg("Execution time for %d sequences: %s (%s per 1000 sequences)" % (len(entries), str(timedelta(seconds=runtime)), round(1000*runtime/len(entries),2)))
    #msg(DF.describe())

    msg("Saving to", save_file)
    if options.format == 'hdf':
        DF.to_hdf(save_file, key="df", complib='blosc:zstd', complevel=2)
    elif options.format.startswith('csv'):
        DF.to_csv(save_file, sep='\t', index=True)
    elif options.format.startswith('pa'):
        DF.to_parquet(save_file, index=True)

    print(biopython_proteinanalysis_seq.cache_info())
