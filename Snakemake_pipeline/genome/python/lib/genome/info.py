

import sys
import numpy as np



class Info(object):
    '''Parses the VCF info field'''


    # define flags for incomplete transcripts
    flags = set(['cds_start_NF', 'cds_end_NF'])

    # clinvar flags
    clinvar_flags = set(["pathogenic", "likely_pathogenic",
                         "pathogenic&likely_pathogenic",
                         "likely_pathogenic&pathogenic"])

    clinvar_benign_flags = set(["benign", "likely_benign",
                                "benign&likely_benign", "benign&benign",
                                "benign&benign&likely_benign",
                                "benign&likely_benign&benign"])
                                
    
    # define the set of loss-of-function consequences
    lof_consequences = set(["transcript_ablation", "splice_donor_variant", \
        "splice_acceptor_variant", "stop_gained", "frameshift_variant",  \
        "start_lost", "initiator_codon_variant",
        "conserved_exon_terminus_variant"])
    
    # define the set of missense (or non loss-of-function) consequences
    missense_consequences = set(["stop_lost", \
        "inframe_insertion", "inframe_deletion", "missense_variant", \
        "transcript_amplification", "protein_altering_variant",
        "splice_region_variant"])
        

    def parse_info(self, info_str, meta_lines):
        '''Parses the INFO column from VCF files'''

        vep_cats = []
        for line in meta_lines:
            if "ID=CSQ" in line:
                vep_cats = line.split("=")[-1].split("Format: ")[1].\
                  replace('"', '').replace('>', '').split("|")
          
        
        self.info = {}
        
        for item in info_str.split(";"):
            if "=" in item:
                elements = item.split("=")
                key = elements[0]
                value = "=".join(elements[1:]) 
                #key, value = item.split("=")
            else:
                key, value = item, True
            self.info[key] = value

        self.info['CSQ'] = self.parse_vep_info(self.info, vep_cats)
        self.most_severe_dict, self.most_severe_gene, self.most_severe_consequence, self.most_severe_flags = \
          self.get_most_severe(vep_cats)
        self.vep_alt_allele = self.get_vep_alt_allele()
        self.exac_af = self.get_exac_af()
        self.gnomad_exome_af = self.get_gnomad_exome_af()
        self.gnomad_genome_af = self.get_gnomad_genome_af()
        self.our_af, self.our_hom_alt, self.our_unaff_hom_alt, self.our_het_alt = self.get_our_af()
        self.gme_af = self.get_gme_af()
        self.max_af = self.get_max_af()
        self.exac_mpc = self.get_exac_mpc()

        

    def parse_vep_info(self, info, vep_cats):
        '''Parses the CSQ element of Info columns. This is
        variant annotation output from VEP'''

        vep_dict = {}
        
        if 'CSQ' in info:
            vep_annotations = info['CSQ'].split(",")

            for ann in vep_annotations:
                ann_fields = ann.split("|")

                for i in range(len(ann_fields)):
                    key = vep_cats[i]
                    value = ann_fields[i]

                    if value == "":
                        value = "."

                    if key not in vep_dict:
                        vep_dict[key] = []

                    vep_dict[key].append(value)

        else:
            return None

        return vep_dict


    
    def get_most_severe(self, vep_cats):
        '''Using the pick flag from VEP, this function selects the most
        severe annotation for this variant'''

        # make an empty dict with all the vep cats
        most_severe_dict = {}
        for cat in vep_cats:
            most_severe_dict[cat] = "."
            
        
        vep_dict = self.info['CSQ']
        
        if vep_dict is None:
            return most_severe_dict, ".", ".", "."
        
        allele_list = vep_dict['Allele']
        ann_len = len(allele_list)
        
        for i in range(ann_len):

            if vep_dict['PICK'][i] == "1":
                # VEP selected this as the most severe annotation

                flags = []

                # check to make sure there is good transcript support
                if vep_dict['FLAGS'][i] != "":
                    flags.append(vep_dict['FLAGS'][i])
                    #return most_severe_dict, ".", "."
                    

                if 'TSL' in vep_dict:
                    if vep_dict['TSL'][i] != "" and vep_dict['TSL'][i] != "1":
                        flags.append("TSL" + str(vep_dict['TSL'][i]))
                        #return most_severe_dict, ".", "."

                
                # take all the annotations and put them in the most severe dict
                for k in vep_dict.keys():
                    if vep_dict[k][i] != "":

                        if k == "MutationTaster_pred":
                            # mutation taster has multple entries joined by &, just
                            # take the first one
                            most_severe_dict[k] = vep_dict[k][i].split("&")[0]

                        else:
                            most_severe_dict[k] = vep_dict[k][i]

                
                most_severe_gene = vep_dict['SYMBOL'][i]
                most_severe_consequence = vep_dict['Consequence'][i]

                if flags == []:
                    flags = "."
                elif len(flags) == 1:
                    flags = flags[0]
                else:
                    flags = ",".join(flags)

                return most_severe_dict, most_severe_gene, most_severe_consequence, flags

        return most_severe_dict, ".", ".", "."
            


    def get_vep_alt_allele(self):
        '''The alleles in VEP's output do not always match those in the
        alt column. Converts the list of alt alleles to a list of
        alt alleles outputted by VEP'''

        vep_dict = self.info['CSQ']

        if vep_dict is None:
            return None

        if self.alt_base == "*":
            # vep doesnt annotate * alleles
            return None

        if self.alt_base in vep_dict['Allele']:
            return self.alt_base

        if len(self.ref_base) > 1 and "-" in vep_dict['Allele']:
            # this is a deletion, vep represents them as -
            return "-"

        if len(self.alt_base) > 1 and self.alt_base[1:] in vep_dict['Allele']:
            # this is an insertion, vep does not include the initial base
            return self.alt_base[1:]

        if len(self.alt_base.split(",")) > 1:
            # we do not expect to see multi-allelic sites
            return None

        #raise ValueError("This allele (%s) is not in VEP "
        #                "annotation at %s:%d\n" %
        #                (self.alt_base, self.chrom_name, self.start))
        return None


    def get_exac_af(self):
        '''Get the ExAC AF for this allele'''

        # ExAC AF entries can look like the following:
        # G:3.819e-04
        # C:0.007143&A:0

        vep_dict = self.info['CSQ']
        if vep_dict is None:
            return -1

        if 'ExAC_MAF' not in vep_dict:
            return -1

        exac_af = vep_dict['ExAC_MAF'][0] # use AF from 1st entry
        if exac_af == "" or exac_af == ".":
            return -1
        
        exac_afs = exac_af.split("&")
        for af in exac_afs:
            af_fields = af.split(":")
            allele = af_fields[0]
            allele_af = float(af_fields[1])
            
            if allele == self.vep_alt_allele:
                return allele_af

        # sys.stderr.write("Allele (%s) not found in ExAC annotation at %s:%d\n" %
        #                  (self.vep_alt_allele, self.chrom_name, self.start))

        return -1

        


    def get_gnomad_exome_af(self):
        '''Get the gnomAD exome AF for this allele'''

        vep_dict = self.info['CSQ']
        if vep_dict is None:
            return -1

        if 'gnomAD_AF' not in vep_dict:
            return -1

        gnomad_af = vep_dict['gnomAD_AF'][0] # use AF from 1st entry
        if (gnomad_af == "") or (gnomad_af == ".") or ("&" in gnomad_af):
            return -1
        
        return float(gnomad_af)




    def get_gnomad_genome_af(self):
        '''Get the gnomAD exome AF for this allele'''

        vep_dict = self.info['CSQ']
        if vep_dict is None:
            return -1

        if 'gnomADg_AF' not in vep_dict:
            return -1

        gnomad_af = vep_dict['gnomADg_AF'][0] # use AF from 1st entry
        if (gnomad_af == "") or (gnomad_af == ".") or ("&" in gnomad_af):
            return -1
        
        return float(gnomad_af)


    

    def get_our_af(self):
        '''Get the AF and number of homozygous alts from our cohort located in info'''

        # these AFs can be annotated like
        # Our_AF=C:0.00629921259843:C:0.00629921259843:C:0.013671875:*:0.013671875
        # take the 1st AF that matches our alt allele


        if "Our_AF" not in self.info:
            return -1, -1, -1, -1


        af_fields = self.info["Our_AF"].split(":")
        hom_alt_fields = self.info["Our_HomAlt"].split(":")

        if "Our_HetAlt" not in self.info:
            het_alt_fields = [-1, -1]
        else:
            het_alt_fields = self.info["Our_HetAlt"].split(":")            
        
        if "Our_Unaff_HomAlt" not in self.info:
            unaff_hom_alt_fields = [-1, -1]
        else:
            unaff_hom_alt_fields = self.info["Our_Unaff_HomAlt"].split(":")            
            
        ref_base = af_fields[0].split(">")[0]
        alt_base = af_fields[0].split(">")[1]
        allele_af = float(af_fields[1])
        hom_alt = int(hom_alt_fields[1])
        het_alt = int(het_alt_fields[1])
        unaff_hom_alt = int(unaff_hom_alt_fields[1])

        if (ref_base == self.ref_base) and (alt_base == self.alt_base):
            return allele_af, hom_alt, unaff_hom_alt, het_alt

        if (self.ref_base == self.alt_base) and (ref_base == self.ref_base[0]) and \
          (alt_base == self.alt_base[0]):
            # this is a weird situation: example ref=TGT, alt=TGT but annotation
            # converts Our_AF to T>T
            return allele_af, hom_alt, unaff_hom_alt, het_alt

        if (alt_base == self.alt_base):
            # the ref alleles don't match, happens for some sites after liftover
            return allele_af, hom_alt, unaff_hom_alt, het_alt

        

          
        raise ValueError("Allele (%s) not found in Our_AF annotation at %s:%d\n" %
                         (self.alt_base, self.chrom_name, self.start))


    def get_gme_af(self):
        '''Get the AF from the GME cohort'''

        if "GME_AF" not in self.info:
            return -1

        if ":" in self.info['GME_AF']:
            # There are multiple AFs in this line, need to
            # fix this in the GME input file by adding alt base
            # info, just throw these out for now
            return -1

        return float(self.info['GME_AF'])


    

    def get_max_af(self):
        '''Gets allele frequencies from our cohort, the GME cohort and
        ExAC and returns the max allele frequency'''

        return max([self.our_af, self.gme_af, self.exac_af])

        

    def is_lof(self):
        '''checks if a variant has a loss-of-function consequence'''

        # non-coding and UTR variants are not LOFs
        if "non_coding" in self.most_severe_consequence:
            return False

        if "UTR" in self.most_severe_consequence:
            return False

        if "NMD_transcript_variant" in self.most_severe_consequence:
            return False

        # check that at least one conseqeunce (joined by &) is a LOF
        if self.most_severe_consequence in self.lof_consequences:
            return True

        if sum([csq in self.lof_consequences \
                for csq in self.most_severe_consequence.split("&")]) > 0:
            return True
        
        return False
    

    def is_loftee_lof(self):
        '''checks if variant is a high quality LOF from loftee'''

        lof_call = self.most_severe_dict['LoF']
        lof_filter = self.most_severe_dict['LoF_filter']

    
        if lof_call == ".":
            # this is not a LOF from loftee
            return False

        if lof_call == "HC":
            return True
        elif lof_filter == "NON_CAN_SPLICE_SURR":
            # this is a LC call with a non_can_splice_surr filter
            # current versions of loftee change this to a flag
            return True
        else:
            # the rest of LC calls
            return False

    
    def is_missense(self):
        '''checks if a variant has a missense consequence'''

        if self.is_lof():
            # check if its a LOF first
            return False

        if ("splice_donor" in self.most_severe_consequence) and \
          ("UTR" in self.most_severe_consequence):
            # we are calling this a missense for now
            return True
        if ("splice_acceptor" in self.most_severe_consequence) and \
          ("UTR" in self.most_severe_consequence):
            return True

        if "non_coding" in self.most_severe_consequence:
            # do not include non-coding transcripts
            return False

        if "NMD_transcript_variant" in self.most_severe_consequence:
            return False

        # check that at least one conseqeunce (joined by &) is a missense
        if self.most_severe_consequence in self.missense_consequences:
            return True

        if sum([csq in self.missense_consequences for csq in self.most_severe_consequence.split("&")]) > 0:
            return True
        
        return False
        


    def is_zina_high_missense(self):
        '''checks if variant is a highly deleterious missense variant, do not include
        clinvar missense variants here (they get their own category)'''

        if self.is_missense() is False:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_flags:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_benign_flags:
            return False

        if "UTR" in self.most_severe_consequence:
            return False
        
        if "splice_region_variant" in self.most_severe_consequence:
            return True

        if 'REVEL_score' in self.most_severe_dict:
            revel_score = self.most_severe_dict['REVEL_score']
            if revel_score != "." and float(revel_score) >= 0.4:
                return True

        if 'MutationTaster_pred' in self.most_severe_dict:
            mt_score = self.most_severe_dict['MutationTaster_pred']
            if mt_score == "D":
                return True
        
        return False


    def is_high_missense(self):
        '''checks if variant is a highly deleterious missense variant, do not include
        clinvar missense variants here (they get their own category)'''

        if self.is_missense() is False:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_flags:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_benign_flags:
            return False

        if "splice_region_variant" in self.most_severe_consequence:
            return True

        if 'REVEL_score' in self.most_severe_dict:
            revel_score = self.most_severe_dict['REVEL_score']
            if revel_score != "." and float(revel_score) >= 0.4:
                return True
        
        return False



    
    def is_clinvar_missense(self):
        '''checks if the missense variant is pathogenic in clinvar'''

        if self.is_missense() is False:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_flags:
            return True

        return False
        


    def is_med_missense(self):
        '''checks if variant is a medium deleterious missense variant, do not include
        clinvar missense variants here (they get their own category)'''

        if not self.is_missense():
            return False

        if self.is_high_missense():
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_benign_flags:
            return False

        if self.most_severe_dict['CLIN_SIG'] in self.clinvar_flags:
            return False

        if 'REVEL_score' in self.most_severe_dict:
            revel_score = self.most_severe_dict['REVEL_score']
            if revel_score != "." and float(revel_score) >= 0.3 and float(revel_score) < 0.4:
                return True

        return False


    def is_synonymous(self):
        '''checks if a variant has a synonymous consequence'''
    
        if "synonymous" in self.most_severe_consequence:
            return True
        return False
    

    def get_exac_mpc(self):
        '''Reads the exacMPC info field to get a MPC score for this site'''

        if 'exacMPC' not in self.info:
            mpc_score = -1
            return mpc_score

        # make dict of exacMPC scores for this site
        exac_base_change = self.info['exacMPC'].split(":")[0]
        exac_val = self.info['exacMPC'].split(":")[1]

        site_key = self.ref_base + ">" + self.alt_base
        if site_key != exac_base_change:
            raise ValueError("Incorrect base change for MPC score")

        if exac_val == "NA":
            mpc_score = -1
        else:
            mpc_score = float(exac_val)
            
        return mpc_score

