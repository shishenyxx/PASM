#!/bin/env python

import numpy as np

import genome.coord
from genome.info import Info

import sys

import re


GTYPE_UNKNOWN = -1

VCF_HEADER_NAMES = ['CHROM', "POS", "ID", "REF", "ALT", 
                    "QUAL", "FILTER", "INFO"]


class VCFError(ValueError):
    """defines an exception for parsing VCF files"""
    pass


class VCFReader(object):
    """A class for iterating over lines in a VCF file. Takes a an
    iterable object (e.g. File, GzipFile) as a constructor argument.
    Parses the VCF headers and then parses a row each time next() is
    called."""
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file

        # store meta information 
        self.meta_lines = []

       # read "meta" lines from VCF file that start with "##"
        for line in self.vcf_file:
            if line.startswith("##"):
                self.meta_lines.append(line.rstrip())
            else:
                self.meta_lines.append(line.rstrip())
                header = line
                break
        
        # read header
        if not header.startswith("#CHROM"):
            raise VCFError("expected header line to start with '#CHROM'")

        words = header[1:].rstrip().split("\t")

        if len(words) < len(VCF_HEADER_NAMES):
            raise VCFError("expected header to have at least %d columns" %
                           len(VCF_HEADER_NAMES))

        # check header names
        i = 0
        for pair in zip(words, VCF_HEADER_NAMES):
            i += 1
            if pair[0] != pair[1]:
                raise VCFError("expected header col %d to be %s "
                               "not '%s'" % (i, pair[1], pair[0]))
        
        if len(words) > len(VCF_HEADER_NAMES):
            if words[8] != "FORMAT":
                raise VCFError("expected header col 8 to be FORMAT "
                               "not '%s'" % words[8])

            self.header_words = words

            self.sample_names = []
            # parse sample names
            for word in words[9:]:
                self.sample_names.append(word)
            
        else:
            self.sample_names = []
            self.header_words = []

        self.n_sample = len(self.sample_names)
        self.vcf_row = VCFRow(self.n_sample, self.meta_lines)


    def __iter__(self):
        """required to make this an iterator object. just returns self"""
        return self


    def __next__(self):
        """parses next row in file, repopulates VCFRow and returns
        it. Note that for speed, the same VCF row is used each time."""
        next_line = self.vcf_file.__next__()
        self.vcf_row.parse_line(next_line)
        return self.vcf_row
    


class VCFRow(Info):

    # PAR regions for hg19 coords
    x_par1 = 2699520
    x_par2 = 154930290
    y_par1 = 2649520
    y_par2 = 59034050
    
    def __init__(self, n_sample, meta_lines):
        self.chrom_name = None
        self.start = None
        self.snp_id = None
        self.ref_base = None
        self.alt_base = None
        self.alt_alleles = None
        self.qual = None
        self.filter = None
        self.info = {}
        self.format = None
        self.gtypes = None

        self.n_sample = n_sample
        self.meta_lines = meta_lines


    def parse_line(self, line):
        words = line.rstrip().split("\t")        
        self.chrom_name = words[0]
        self.start = int(words[1])
        self.snp_id = words[2]
        self.ref_base = words[3]
        self.alt_base = words[4]
        self.alt_alleles = words[4].split(",")
        self.qual = words[5]
        self.filter = words[6]

        self.parse_info(words[7], self.meta_lines)

        
        #if not self.chrom_name.startswith("chr"):
        #    self.chrom_name = "chr" + self.chrom_name

        # the format string gives the order of genotype info in the
        # remaining columns example:
        # GT:GQ:DP:BQ:MQ:AD:FA:VAQ:FET:FT:SS
        #
        # GT - genotype separater of '/' or '/' indicates 
        #      phased or unphased (e.g. 0/1, 1|0, or 1/2),
        # DP - read depth at this position for this sample (Integer)
        # GQ - conditional genotype quality, encoded as a phred quality 
        #      -10log_10p (genotype call is wrong, conditioned on 
        #      site being variant) (Integer)
        

        if len(words) > len(VCF_HEADER_NAMES):
            self.format = words[8]
            self.gtypes = words[9:]

            if len(self.gtypes) != self.n_sample:
                raise VCFError("number of genotypes does not "
                               "match number of samples")

                               
    def is_in_par(self):
        '''True/False indicator of being in the pseudoautosomal region'''
        
        if self.chrom_name in ['chrY', 'Y']:
            if (self.start < self.y_par1) or (self.start > self.y_par2):
                return True

        if self.chrom_name in ['chrX', 'X']:
            if (self.start < self.x_par1) or (self.start > self.x_par2):
                return True

        return False
        

    def get_field_strs(self, format_code):
        """Returns a list of n_sample elements, corresponding to the
        provided format code. E.g. to get the genotype strings for all
        individuals: genotypes = vcf_row.get_field_strs("GT")
        """
        format_words = self.format.split(":")

        if format_code not in format_words:
            return ['.' for x in self.gtypes]

        i = format_words.index(format_code)

        return ['.' if ((x == '.') or (x == "./."))
                else x.split(":")[i] for x in self.gtypes]

        

    def get_int16_field(self, format_code):
        """Returns an numpy int array that of n_sample elements,
        corresponding to the provided format code. E.g. to get the
        read depth for all individuals: 
        depths = vcf_row.get_int_field("DP")
        """
        format_words = self.format.split(":")
        i = format_words.index(format_code)

        return np.array([0 if ((x == '.') or (x == "./.") or x.split(":")[i] == ".")
                         else int(x.split(":")[i])
                         for x in self.gtypes], dtype=np.int16)
        


    def get_num_alt_alleles(self):
        """Parses genotypes from a VCFRow. Returns integer
        array, n_samples long, that gives number of copies of
        ALT allele for each sample (0, 1, or 2). 
        Elements are set to GTYPE_UNKNOWN=-1 for samples without
        genotype information."""
        gtype_strs = self.get_field_strs("GT")
        n = len(gtype_strs)
        gtypes = np.zeros(n, dtype=np.int8)
        
        for i in range(n):
            if gtype_strs[i] == '.' or gtype_strs[i] == "" or \
              gtype_strs[i] == "./." or gtype_strs[i] == ".|.":
                gtypes[i] = GTYPE_UNKNOWN
                continue

            words = re.split("\||/", gtype_strs[i])
            
            if words[0] != "0" and words[0] != ".":
                gtypes[i] += 1
            if len(words) > 1:
                if words[1] != "0" and words[1] != ".":
                    gtypes[i] += 1
            
        return gtypes


    def get_num_ref_alleles(self):
        """Parses genotypes from a VCFRow. Returns integer
        array, n_samples long, that gives number of copies of
        REF allele for each sample (0, 1, or 2). 
        Elements are set to GTYPE_UNKNOWN=-1 for samples without
        genotype information."""
        gtype_strs = self.get_field_strs("GT")
        n = len(gtype_strs)
        gtypes = np.zeros(n, dtype=np.int8)
        
        for i in range(n):
            if gtype_strs[i] == '.' or gtype_strs[i] == "" or \
              gtype_strs[i] == "./." or gtype_strs[i] == ".|.":
                gtypes[i] = GTYPE_UNKNOWN
                continue
            
            words = re.split("\||/", gtype_strs[i])
            
            if words[0] == "0":
                gtypes[i] += 1
            if len(words) > 1:
                if words[1] == "0":
                    gtypes[i] += 1
        
        return gtypes


    def get_num_ref_alt_alleles(self):
        """Parses genotypes from a VCFRow. Returns integer
        array, n_samples long, that gives number of copies of
        REF allele for each sample (0, 1, or 2). 
        Elements are set to GTYPE_UNKNOWN=-1 for samples without
        genotype information."""
        gtype_strs = self.get_field_strs("GT")
        n = len(gtype_strs)
        ref_gtypes = np.zeros(n, dtype=np.int8)
        alt_gtypes = np.zeros(n, dtype=np.int8)

        for i in range(n):
            if gtype_strs[i] == '.' or gtype_strs[i] == "" or \
              gtype_strs[i] == "./." or gtype_strs[i] == ".|.":
                ref_gtypes[i] = GTYPE_UNKNOWN
                alt_gtypes[i] = GTYPE_UNKNOWN
                continue

            words = re.split("\||/", gtype_strs[i])

            if words[0] == "0":
                ref_gtypes[i] += 1
            elif words[0] != ".":
                alt_gtypes[i] += 1

            if len(words) > 1:
                if words[1] == "0":
                    ref_gtypes[i] += 1
                elif words[0] != ".":
                    alt_gtypes[i] += 1


        return ref_gtypes, alt_gtypes


    
    def get_alleles(self):
        """Parses genotypes from a VCFRow, but keeps coding as
        separate alleles. Returns integer array, n_samples*2 long,
        coded as 0=reference, 1=alternate, -1=unknown."""
        gtype_strs = self.get_field_strs("GT")
        n = len(gtype_strs)
        gtypes = np.empty(n*2, dtype=np.int8)
        gtypes[:] = GTYPE_UNKNOWN

        for i in range(n):
            if gtype_strs[i] == '.':
                gtypes[i] = GTYPE_UNKNOWN
            else:                
                # split on '|' or on '/'
                words = re.split("\||/", gtype_strs[i])
                if len(words) != 2:
                    gtypes[2*i] = GTYPE_UNKNOWN
                    gtypes[2*i+1] = GTYPE_UNKNOWN
                else:
                    if words[0] != '.':
                        val = int(words[0])
                        if val == 0:
                            # first allele is reference
                            gtypes[2*i] = 0
                        elif val > 0:
                            # first allele is non-reference
                            gtypes[2*i] = 1
                    if words[1] != '.':
                        val = int(words[1])
                        if val == 0:
                            # second allele is reference
                            gtypes[2*i+1] = 0
                        elif val > 0:
                            # second allele is non-reference
                            gtypes[2*i+1] = 1

        return gtypes
        

    def __str__(self):
        '''Returns a string representation of the first 8 "
        "fields for this variant'''
        
        # info_str = ";".join([str(k) + "=" + str(v) for k,v in self.info_dict.items()])
        # if info_str == "":
        #     info_str = "."

        
        return("\t".join([self.chrom_name, str(self.start), 
                          self.snp_id, self.ref_base, 
                          self.alt_base, str(self.qual), 
                          self.filter, ".", 
                          self.format, "\t".join(self.gtypes)]))

