import os, sys, re
import tables
import numpy as np

from genome.track import Track
import genome.trackstat
from genome.chrom import Chromosome

DEFAULT_ASSEMBLY = "hg18"

class GenomeDB():
    """This object provides a simple database interface overtop
    a set of directories containing HDF5 files. This is loosely based on
    Michael Hoffman's Genomedata package (see
    http://bioinformatics.oxfordjournals.org/content/26/11/1458.abstract
    or http://noble.gs.washington.edu/proj/genomedata/).

    The base directory for the database is obtained from the GENOME_DB
    environment variable or can specified as a constructor
    argument. Each HDF5 file under the base directory is considered
    a "Track", which can be accessed using a genome.track.Track object.
    These are typically opened for reading using the get_track method
    and can be created using the create_track method.

    For convenience, the GenomeDB class also provides several methods
    for obtaining the list of chromosomes that are in the database."""

    def __init__(self, path=None, assembly=None):
        if path is None:
            if 'GENOME_DB' in os.environ:
                path = os.environ['GENOME_DB']
            else:
                raise ValueError("No database path specified and "
                                 "GENOME_DB environment variable not set")

        if not path.endswith("/"):
            # add trailing '/' if none exists
            path = path + "/"

        if not os.path.exists(path):
            raise IOError("database path does not exist: %s" % path)

        if assembly is None:
            if 'GENOME_ASSEMBLY' in os.environ:
                assembly = os.environ['GENOME_ASSEMBLY']
            else:
                sys.stderr.write("no assembly specified and GENOME_ASSEMBLY "
                                 "environment variable not set--using %s "
                                 "by default\n" % DEFAULT_ASSEMBLY)
                assembly = DEFAULT_ASSEMBLY

        # add assembly to the end of the path
        assembly_path = path + assembly + "/" 

        if not os.path.exists(assembly_path):
            raise IOError("assembly '%s' does not exist under "
                          "database %s" % (assembly, path))

        if not os.path.isdir(assembly_path):
            raise IOError("assembly path is not directory: %s" %
                          assembly_path)

        self.assembly = assembly
        self.path = assembly_path
        
    
    def __enter__(self):
        sys.stderr.write("GenomeDB connection established\n")
        return self


    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.write("Cleaning up GenomeDB\n")
        return False


    def get_track_path(self, track_name):
        """Returns the filesystem path to the HDF5 file with the
        given track name"""
        if track_name.startswith("/"):
            track_name = track_name[1:]
        if track_name.endswith(".h5"):
            track_name = track_name[:-3]
        track_path = self.path + track_name + ".h5"

        return track_path

        
    def has_track(self, track_name):
        """Returns True if a track with the specified name exists"""
        track_path = self.get_track_path(track_name)

        if os.path.exists(track_path):
            return True

        return False

        
    def open_track(self, track_name, mode="r"):
        """Returns an open Track of the specified name. By default the
        track is opened in read mode, but other modes can be
        specified."""
        track_path = self.get_track_path(track_name)

        if not os.path.exists(track_path):
            raise ValueError("track %s does not exist" % track_name)
        
        return Track(track_name, track_path, mode)



    def init_track(self, track, data_type=np.float32, dflt=None):
        """initializes a track by creating arrays for every chromosome
        and setting them to a default value. Unless specified, the default
        value is nan for floats, and 0 for ints/uints"""
        dt = np.dtype(data_type)

        zlib_filter = tables.Filters(complevel=1, complib="zlib")
        
        if dflt is None:
            if np.issubdtype(dt, np.float):
                dflt = np.nan
            elif np.issubdtype(dt, np.int):
                dflt = 0
            elif np.issubdtype(dt, np.unsignedinteger):
                dflt = 0
            else:
                raise ValueError("don't know how initialize track with datatype "
                                 "%s" % dt.name)

        for chrom in self.get_all_chromosomes():
            shape = [chrom.length]
            atom = tables.Atom.from_dtype(dt, dflt=dflt)
            track.h5f.createCArray(track.h5f.root, chrom.name, atom, shape, 
                                   filters=zlib_filter)
            
            

    def get_track_stat(self, track):
        """Returns a TrackStat object containing statistics for an
        entire track (mean, max, sum, etc.)"""
        return genome.trackstat.get_stats(self, track)


    def set_track_stat(self, track):
        """computes and sets track statistics on the provided track object"""
        genome.trackstat.set_stats(self, track)

        
    def list_tracks(self, subdir=None):
        """returns a list of existing track names"""
        if subdir is None:
            path = self.path
        else:
            path = self.path + "/" + subdir + "/"
        
        filenames = os.listdir(path)

        track_names = []

        for filename in filenames:
            if filename.startswith("."):
                continue
            if filename.endswith(".h5"):
                track_name = filename[:-3]
                if subdir:
                    track_names.append(subdir + "/" + track_name)
                else:
                    track_names.append(track_name)
            elif os.path.isdir(path + filename):
                # recursively add tracks to this one                
                new_subdir = subdir + "/" + filename if subdir else filename
                track_names.extend(self.list_tracks(subdir=new_subdir))
        return track_names


    def create_track(self, track_name):
        """Creates a new HDF5 file in write mode, and returns
        a Track object wrapped around it."""
        track_path = self.get_track_path(track_name)

        if os.path.exists(track_path):
            sys.stderr.write("Track exists, overwriting (%s)\n" % (track_path))
            # raise IOError("Could not create track '%s' because it "
            #               "already exists.\nYou must remove "
            #               "the file '%s' before this track can be created."
            #               % (track_name, track_path))

        # create parent directories as needed
        dir_names = track_name.split("/")[:-1]
        base_dir = self.path
        for dir_name in dir_names:
            base_dir = base_dir + "/" + dir_name
            if not os.path.exists(base_dir):
                os.mkdir(base_dir)

        return Track(name=track_name, path=track_path, mode="w")



    def get_chromosome_dict(self):
        """Returns a dictionary of all chromosomes in the database,
        keyed on chromosome name"""
        chrom_list = self.get_all_chromosomes()        
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom.name] = chrom

        return chrom_dict
        


    def get_chromosomes(self, get_rand=False, get_auto=True, get_sex=True,
                        get_x=True, get_y=False, get_hap=False,
                        get_mito=False):
        """Returns a filtered list of Chromosomes from the
        database. Optional flags specify the subset of Chromosomes
        that are returned. By default the 22 autosomes and chrX are
        retrieved (but chrY, the mitochondrial chromosome, alternate
        haplotypes, and 'random' chromosomes are not)"""
        chrom_list = []
        chrom = None

        chrom_track = self.open_track("chromosome")

        flags = {'is_rand' : get_rand,
                 'is_auto' : get_auto,
                 'is_sex' : get_sex,
                 'is_x' : get_x,
                 'is_y' : get_y,
                 'is_hap' : get_hap,
                 'is_mito' :  get_mito}

        conditions = []
        for k, v in flags.items():
            if not v:
                # don't get this chromosome type
                conditions.append("(%s == False)" % k)

        row_iter = None
        if len(conditions) > 0:
            query = " & ".join(conditions)
            row_iter = chrom_track.h5f.root.chromosome.where(query)
        else:
            row_iter = chrom_track.h5f.root.chromosome

        for row in row_iter:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        chrom_track.close()

        return chrom_list
    

    def get_all_chromosomes(self):
        """Returns an unfiltered list of all of the chromosomes in the
        database"""
        chrom_list = []
        chrom = None

        chrom_track = self.open_track("chromosome")

        for row in chrom_track.h5f.root.chromosome:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        chrom_track.close()

        return chrom_list
        

        


    def get_chromosomes_from_args(self, args):
        """Convenience function, returns a list of chromosome objects
        that are obtained from parsing command line args that may be
        in any one of the following forms 'chr2' or '3', or '1-22'"""
        if isinstance(args, str):
            # if we were given a string, put it in a list
            args = [args]
        
        if len(args) < 1:
            raise ValueError("expected at least one value, got 0")
        
        chrom_name_dict = self.get_chromosome_dict()

        chrom_id_dict = dict([(str(c.idnum), c) for c \
                              in chrom_name_dict.values()])

        chrom_list = []
        for arg in args:
            words = arg.split("-")
            if len(words) == 2:
                # try parsing argument as a range of chromosome
                # ids like 1-22
                try:
                    start = int(words[0])
                    end = int(words[1])

                    if start <= end:
                        vals = [str(x) for x in range(start, end+1)]
                    else:
                        vals = [arg]
                except:
                    vals = [arg]
            else:
                vals = [arg]

            for val in vals:
                if val in chrom_name_dict:
                    chrom_list.append(chrom_name_dict[val])
                elif val in chrom_id_dict:
                    chrom_list.append(chrom_id_dict[val])
                else:
                    raise ValueError("unknown chromosome %s" % val)

        return chrom_list



    def get_chromosome(self, name):
        """Retrieves a single chromosome by name"""
        chrom_dict = self.get_chromosome_dict()
        return chrom_dict[name]

