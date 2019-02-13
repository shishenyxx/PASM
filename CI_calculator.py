#!/usr/bin/python

import sys
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import re

def CI_calculator(f):
    file = ''
    mh=importr("yyxMosaicHunter")
    bases = ['A','a','T','t','C','c','G','g']
#    for chunk in f.chunks():
#        for line in chunk.split('\n'):
    for line in f:
        if True:
            if len(line.split()) > 0:
                items = line.strip().split("\t")
                chrom = items[0]
                pos = items[1]
                ref = items[2].upper()
                ref_lc = ref.lower()
                depth = items[3]
                ids = items[6]
                if depth > 0:
                    match = items[4]
                    quality = items[5]
                    pos_n = {}
                    m = 0
                    count = {}
                    for base in bases:
                        count[base] = 0
                    m_l = 0
                    I = 0; D = 0
                    o = 0; a = 0
                    star = 0
                    in_base = 0; del_base = 0; mis_base = 0

                    l = len(match)
                    i = 0
                    pos_q = 0
                    end_reads = {}
                    while(i < l):
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
                                c = 10*match[i+1]+match[i+2]
                                i += c+4
                                pos_q += c+4
                            else:
                                c = match[i+1]
                                i += c+3
                                pos_q += c+3
                            in_base += c
                        elif match[i] == '-':
                            D += 1
                            if re.match(r'\d',match[i+2]):
                                c = 10*match[i+1]+match[i+2]
                                i += c+4
                            else:
                                c = match[i+1]
                                i += c+3
                            del_base += c
                        else:
                            i += 1

                    allbases = {}
                    for base in ['A','T','C','G']:
                        count.setdefault(base, 0)
                        count.setdefault(base.lower(),0)
                        allbases[base] = count[base]+count[base.lower()]  

                    keys = sorted(allbases,key=lambda k:allbases[k],reverse=True)
                    if keys[0] == ref:
                        alt1 = keys[1]
                    else:
                        alt1 = keys[0]
                    ref_lc = ref.lower()
                    alt1_lc = alt1.lower()
                    quality_ref = ""; quality_alt = "";
                    quality_score_ref = []; quality_score_alt = []
                    for p in sorted(pos_n[ref]+pos_n[ref_lc]):
                        quality_ref += quality[p]
                        quality_score = ord(quality[p])-33
                        quality_score_ref.append(quality_score)
                
                    for p in pos_n[alt1]+pos_n[alt1_lc]:
                        quality_ref += quality[p]
                        quality_score = ord(quality[p])-33
                        quality_score_alt.append(quality_score)
            
                    if quality_alt == "":
                        quality_alt = "."

                    quality_ref = re.sub(r'!','',quality_ref)
                    quality_ref = re.sub(r'\\','',quality_ref)
                    quality_alt = re.sub(r'!','',quality_ref)
                    quality_alt = re.sub(r'\\','',quality_ref)
                    rcode = 'quality_ref <- "%s"' %(quality_ref, )
                    robjects.r(rcode)
                    rcode = 'quality_alt <- "%s"' %(quality_alt, )
                    robjects.r(rcode)
                    rcode = 'yyx_get_credible_interval(yyx_wrapped_mosaic_hunter_for_one_site(quality_ref, quality_alt)$likelihood_fun, c(0,1), 0.95)'
                    cred_int = robjects.r(rcode)
                    print cred_int[0],cred_int[1]
            else:
                continue

if __name__=="__main__":
    CI_calculator(open(sys.argv[1]))
