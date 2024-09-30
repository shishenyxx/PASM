
import string
import os
import numpy as np

dna_comp = None
codon_dict = None

GENOME_DIR = os.environ['GENOME_DIR']
CODON_TABLE = GENOME_DIR + "/codon_table.txt"

def comp(seq_str):
    """complements the provided DNA sequence and returns it"""
    global dna_comp

    if dna_comp is None:
        dna_comp = string.maketrans("ATCGMRWSYKNatcgmrwsykn",
                                    "TAGCKYWSRMNtagckywsrmn")
    return seq_str.translate(dna_comp)


def revcomp(seq_str):
    """returns reverse complement of provided DNA sequence"""
    return comp(seq_str)[::-1]


def revcomp_nparray(vals):
    seqstr = from_nparray(vals)
    seqstr = revcomp(seqstr)
    return np.array([ord(x) for x in seqstr], dtype=np.uint8)
    

def from_nparray(vals):
    """converts a numpy array into a sequence string"""
    return "".join(chr(x) for x in vals)


def get_codon_dict():
    global codon_dict

    if codon_dict is not None:
        # codon dict already initialized
        return codon_dict

    # need to read and create new codon dict
    f = open(CODON_TABLE, "r")

    codon_dict = {}
    for line in f:
        words = line.rstrip().split(":")
        codon = words[0]
        aa = words[1]

        codon_dict[codon] = aa

    f.close()
    return codon_dict
