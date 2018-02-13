__author__ = 'pfeijao'

import glob
import os
import re
from model import Genome, Chromosome
import subprocess
import json

LINEAR_END_CHR = "$"
CIRCULAR_END_CHR = "@"

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def check_gene_duplicates(gene_list):
    seen = set()
    duplicates = []
    for gene in [abs(gene_with_sign) for gene_with_sign in gene_list]:
      if gene in seen:
        duplicates.append(gene)
      seen.add(gene)
    return duplicates

def open_genome_file(filename, as_list=False):
    """
    Opens a genome file in GRIMM format.
    Example:

    >Genome1
    #chr1 - optional
    1 +2 -3 4 $
    #chr2
    -5 6 $
    >Genome2
    #chr1
    1 2 3 4 5 7 $

    """
    #TODO: add option for genes as strings. I would have to transform them to numbers, and keep a name dictionary
    # to go back when necessary.
    if as_list:
        genomes = []
    else:
        genomes = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split(" ")[0]
                genome = Genome(name)
                if as_list:
                    genomes.append(genome)
                else:
                    genomes[name] = genome
            elif line.startswith("#"):
                # chromosome; ignore name after, each new line is a new chromosome.
                continue
            else:
                if line.endswith(LINEAR_END_CHR):
                    circular = False
                elif line.endswith(CIRCULAR_END_CHR):
                    circular = True
                else:
                    raise RuntimeError("Invalid genome file %s. Unrecognized line:\n%s" % (filename, line))
                genes = map(int, re.split("\s+",line[:-1].strip()))
                duplicates = check_gene_duplicates(genes)
                if len(duplicates) > 0:
                  raise RuntimeError("Duplicated gene(s) %s on genome %s!\nIn the current version, RINGO cannot deal with duplicated genes. Please remove extra copies to be able to run RINGO on these genomes." % (duplicates, name))
                genome.add_chromosome(Chromosome(genes, circular))

    return genomes


def open_copy_number_file(genomes, filename):
    current_idx = 0
    chr_idx = 0
    genome = None
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip().split(" ")[0]
                genome = genomes[current_idx]
                current_idx += 1
                assert genome.name == name
                chr_idx = 0
            elif line.startswith("#"):
                # chromosome; ignore name after, each new line is a new chromosome.
                continue
            else:
                if line.endswith(LINEAR_END_CHR):
                    circular = False
                elif line.endswith(CIRCULAR_END_CHR):
                    circular = True
                else:
                    raise RuntimeError("Invalid file %s. Unrecognized line:\n%s" % (filename, line))
                genome.chromosomes[chr_idx].copy_number = map(int, line[:-1].strip().split(" "))


def open_adjacencies_file(filename):
    """
    Open a  file in adjacencies format. Genome identifiers are FASTA-like ">name" lines, and
    each following line is an adjacency in the format (x,y), where x and y are integers.
    """
    genome_list = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
              name = line[1:].strip().split(" ")[0]
              genome = []
              genome_list[name] = genome
              continue
            # read adjacency:
            genome.append(eval(line))
    # convert to genomes:
    return {name:Genome.from_adjacency_list(name, adj_list) for name, adj_list in genome_list.iteritems()}



def write_genomes_to_file(genomes, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format.
    """
    if isinstance(genomes, dict):
        iterator = genomes.itervalues()
    elif isinstance(genomes, list):
        iterator = iter(genomes)
    with open(filename, "w") as f:
        for genome in iterator:
            f.write(">%s\n" % genome.name)
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(gene) for gene in chromosome.gene_order]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


def write_genomes_copy_number_to_file(genomes, filename, write_chr_line=True):
    """
    Write genomes in a file with GRIMM format.
    """
    if isinstance(genomes, dict):
        iterator = genomes.itervalues()
    elif isinstance(genomes, list):
        iterator = iter(genomes)
    with open(filename, "w") as f:
        for genome in iterator:
            f.write(">%s\n" % genome.name)
            for idx, chromosome in enumerate(genome.chromosomes):
                if write_chr_line:
                    f.write("# chr%d\n" % (idx + 1))
                f.write("%s %s\n" % (" ".join([str(gene) for gene in chromosome.copy_number]),
                                     CIRCULAR_END_CHR if chromosome.circular else LINEAR_END_CHR))


# I/O JSON parameters:
def __read_parameters(filename):
  with open(filename, 'r') as f:
      return json.load(f)

def __write_parameters(param, filename):
  with open(filename,"w") as f:
      json.dump(param.__dict__, f, sort_keys = True, indent = 4)

#def read_simulation_parameters(folder):
#  return __read_parameters(os.path.join(folder, cfg.sim_paramfile()))

def write_simulation_parameters(param, output):
  __write_parameters(param, os.path.join(output,"params.cfg"))
