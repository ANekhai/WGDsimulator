import pyximport
pyximport.install()
from utils import model


def test_single_chromosome_genome_creation():
    genome = model.Genome.identity(10, 1)
    assert len(genome.chromosomes) == 1
    assert len(genome.chromosomes[0].gene_order) == 10


def test_multiple_chromosome_genome_creation():
    genome = model.Genome.identity(10, 2)
    assert len(genome.chromosomes) == 2
    assert len(genome.chromosomes[0].gene_order) == len(genome.chromosomes[1].gene_order)

#TODO: Add tests for genomes that do not divide evenly


#TODO: Add tests for BPGraph features