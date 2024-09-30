

import sys
import os



class Individual(object):
    '''Creates an object for each individual'''
    
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    def __init__(self, family_id, indv_id, dad_id, mom_id, sex, status,
                 gleeson_id):
        
        self.family_id = family_id
        self.indv_id = indv_id
        self.dad_id = dad_id
        self.mom_id = mom_id
        self.sex = sex
        self.status = status
        self.gleeson_id = gleeson_id
        self.children = []
        self.mother = None
        self.father = None

        
    def __repr__(self):
        return 'Person(family_id="{}", person_id="{}", dad_id="{}", ' \
            'mom_id="{}", sex="{}", status="{}", gleeson_id="{}")'.format(
            self.family_id, self.indv_id, self.dad_id, self.mom_id,
            self.sex, self.status, self.gleeson_id)

    
    def is_male(self):
        '''returns True/False for whether the person is male'''
        return self.sex in self.male_codes

    
    def is_female(self):
        '''returns True/False for whether the person is male'''
        return self.sex in self.female_codes

    
    def is_affected(self):
        """returns true or false for affected, rather than the string value
        """
        # change how the affected status is encoded. Current ped files
        # encode "1" for unaffected, and "2" for affected. Change this to
        # True/False values, and catch any unknown affected statuses.
        if self.status not in set(["1", "2"]):
            raise ValueError("unknown status: " + self.status + ", "
                             "should be 1: unaffected, 2: affected")
        
        return self.status == "2"



class Family(object):
    '''Creates a family object'''
    def __init__(self, family_id, members=None):
        self.family_id = family_id
    
        self.members = members
        if self.members is None:
            self.members = []


    def __repr__(self):
        return 'Family(family_id="{}", members="{}")'.format(
            self.family_id, self.members)


    def add_member(self, indv_id, dad_id, mom_id, sex, status):
        '''Adds an indvidual to a family'''
        
        indv = Individual(self.family_id, indv_id, dad_id, mom_id, sex, status)
        self.members.append(indv)

    def get_member_ids(self):
        '''Makes a list of indv_ids for family members'''
        return [indv.indv_id for indv in self.members]

    def get_member_gleeson_ids(self):
        '''Makes a list of Gleeson IDs for family members'''
        return [indv.gleeson_id for indv in self.members]

    def get_member_affected_status(self):
        '''Makes a list of affected status for family members'''
        return [str(int(indv.is_affected())) for indv in self.members]


        

def open_ped(path):
    '''Opens a ped file and groups individuals into families
    
    Returns a dictionary keyed on family_id of dictionaries of
    Individual objects keyed on indv_id'''

    if not os.path.exists(path):
        sys.exit("Path to ped file does not exist: %s\n" % path)

    families = {}

    f = open(path, "r")
    for line in f:
        words = line.rstrip().split("\t")
        family_id = words[0]
        indv_id = words[1]
        dad_id = words[2]
        mom_id = words[3]
        sex = words[4]
        status = words[5]

        # ped file can optionally have a 7th column with gleeson ID
        if len(words) == 6:
            gleeson_id = None
        else:
            gleeson_id = words[6]

        indv = Individual(family_id, indv_id, dad_id, mom_id,
                          sex, status, gleeson_id)

        if family_id not in families:
            families[family_id] = {}

        families[family_id][indv_id] = indv
        
    f.close()
    return families



def load_families(path):
    '''Creates a dict of family data from a PED file'''

    family_dicts = open_ped(path)

    families = {}
    for family_id, family_dict in family_dicts.items():

        members = []
        for indv_id in family_dict.keys():
            indv = family_dict[indv_id]

            # link family members together
            if indv.mom_id != '0':
                
                indv.mother = family_dict[indv.mom_id]
                family_dict[indv.mom_id].children.append(indv)

            if indv.dad_id != '0':
                indv.father = family_dict[indv.dad_id]
                family_dict[indv.dad_id].children.append(indv)

            members.append(indv)

        family = Family(family_id, members)
        families[family_id] = family
        #families.append(family)

    return families
