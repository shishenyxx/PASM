from subprocess import Popen, PIPE
import scipy.stats as stats
import sys
import re
import pandas as pd
import numpy as np
from scipy.stats import beta
from scipy import stats


def translate_bases(ref, depth, match):
    ref_lc = ref.lower()
    bases = ['A','a','T','t','C','c','G','g', "I", "i", "D", "d"]
    pos_n = {}
    count = {}
    for base in bases:
        count[base] = 0
    for base in bases:
        pos_n[base] = []
    I = 0; D = 0
    in_base= 0; del_base = 0; #mis_base = 0
    i = 0
    pos_q = 0
    end_reads = {}                                                                                                                           
    while(i < len(match)):                                                                                                                   
        if match[i] == '^':                                                                                                                  
            i += 2                                                                                                                           
        elif match[i] == "$":                                                                                                                
            i += 1                                                                                                                           
        elif match[i] == ".":                                                                                                                
            count[ref] += 1                                                                                                                  
            pos_n.setdefault(ref,[]).append(pos_q)                                                                                           
            pos_q += 1                                                                                                                       
            i += 1                                                                                                                           
        elif match[i] == ",":                                                                                                                
            count[ref_lc] += 1                                                                                                               
            pos_n.setdefault(ref_lc,[]).append(pos_q)                                                                                        
            pos_q += 1                                                                                                                       
            i += 1                                                                                                                           
        elif match[i] in bases:                                                                                                              
            count[match[i]] += 1
            pos_n.setdefault(match[i],[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] == '+':                                                                                                                
            I += 1                                                                                                                           
            if re.match(r"\d",match[i+2]):
                if match[i+3].isupper():
                    count["I"] += 1
                else:
                    count["i"] +=1
                c = 10*int(match[i+1]) + int(match[i+2])
                i += c+3
                pos_q += c+3
            else:
                if match[i+2].isupper():
                    count["I"] += 1
                else:
                    count["i"] +=1
                c = int(match[i+1])
                i += c+2
                pos_q += c+2
            in_base += c
        elif match[i] == '-':
            D += 1
            if re.match(r'\d',match[i+2]):
                if match[i+3].isupper():
                    count["D"] += 1
                else:
                    count["d"] +=1
                c = 10*int(match[i+1])+ int(match[i+2])
                i += c+3
            else:
                if match[i+2].isupper():
                    count["D"] += 1
                else:
                    count["d"] +=1
                c = int(match[i+1])
                i += c+2
            del_base += c
        elif match[i] == "*":                                                                                                                
            D += 1                                                                                                                           
            i += 1
        else:
            sys.stderr.write("ERROR: " + str(i) + " " + str(match[i]) + "\n")
            sys.exit(0)
    return count, pos_n



def wilson_binom_interval(success, total, alpha = 0.05):
    q_ = success / total
    crit = stats.norm.isf( alpha / 2.)
    crit2 = crit**2
    denom = 1 + crit2 / total
    center = (q_ + crit2 / (2 * total)) / denom
    dist = crit * np.sqrt(q_ * (1. - q_) / total + crit2 / (4. * total**2))
    dist /= denom
    ci_low = center - dist
    ci_upp = center + dist
    return ci_low, ci_upp



def main(argv):
    if len(argv) != 6:
        sys.stderr.write("usage: " + argv[0] + "<chr> <pos> <ref> <alt> <bam_file>\n")
        sys.exit(2)
    gt_chrom = argv[1]
    gt_pos = argv[2]
    gt_ref = argv[3]
    gt_alt = argv[4]
    bam_file = argv[5]

    query_position = str(gt_chrom) + ":" + str(gt_pos) + "-" + str(gt_pos)
    output, error = Popen(["samtools", "mpileup", "-r", query_position, \
           "-f", "/tscc/projects/ps-gleesonlab7/gleeson3/reference_fasta/human_g1k_v37_decoy.fasta",\
           "-Q 13 -q0 -AB -d5000000 ", \
           bam_file],\
           stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()
    output = output.decode()
    items = output.rstrip().split("\t")
    pos = int(items[1])
    ref = items[2].upper()
    depth = int(items[3])                                                                                                                    
    if depth > 0:                                                                                                                           
        match = items[4]
        quality = items[5]                                                                                                                  
        count, pos_n = translate_bases(ref, depth, match)
        num_ref = count[ref] + count[ref.lower()]
        if len(gt_ref) == len(gt_alt):
            num_alt = count[gt_alt] + count[gt_alt.lower()]
            oddsratio, pvalue = stats.fisher_exact([[count[gt_ref], count[gt_ref.lower()]], [count[gt_alt], count[gt_alt.lower()]]])
        elif len(gt_ref) > len(gt_alt): #deletion
            num_alt = count["D"] + count["d"]
            oddsratio, pvalue = stats.fisher_exact([[count[gt_ref], count[gt_ref.lower()]], [count["D"], count["d"]]])
        elif len(gt_ref) < len(gt_alt): #insertion
            num_alt = count["I"] + count["i"]
            oddsratio, pvalue = stats.fisher_exact([[count[gt_ref], count[gt_ref.lower()]], [count["I"], count["i"]]])
        ci_lower, ci_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
        a_count = count["A"] + count["a"]
        c_count = count["C"] + count["c"]
        g_count = count["G"] + count["g"]
        t_count = count["T"] + count["t"]
        i_count = count["I"] + count["i"]
        d_count = count["D"] + count["d"]
#        print("\t".join(["A", "C", "G", "T", "Insertion", "Deletion", "MAF", "CI_lower", "CI_upper", "Fisher_pvalue"]))
        print("\t".join(map(str, [a_count, c_count, g_count, t_count, i_count, d_count, num_alt/(num_ref+num_alt), ci_lower, ci_upper, pvalue])))
main(sys.argv)
