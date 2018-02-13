#!/usr/bin/env python2
import pyximport

pyximport.install()
import argparse
import os
import random
import math
import numpy as np
from utils import file_ops, model


import ConfigParser

__author__ = 'pfeijao' # MODIFIED BY ANEKHAI


# noinspection PyClassHasNoInit
class EventType:
    all = ["rearrangement", "deletion", "insertion", "duplication"]
    REARRANGEMENT, DELETION, INSERTION, DUPLICATION = all


# noinspection PyClassHasNoInit
class RearrangementType:
    all = ["reversal", "translocation"]
    REVERSAL, TRANSLOCATION = all


class SimParameters:
    def __init__(self, wgd_type="", num_genes=0, num_chr=0, del_p=[0, 0], ins_p=[0, 0], indel_length=[0, 0], duplication_p=[0, 0],
                 duplication_length=[0, 0], events=[0, 0]):
        self.WGD_type = wgd_type
        self.num_genes = num_genes
        self.num_chr = num_chr
        self.num_ev = events
        self.deletion_p = del_p
        self.insertion_p = ins_p
        self.indel_length = indel_length
        self.duplication_p = duplication_p
        self.duplication_length = duplication_length

        # Rearrangement, Insertion and Deletion prob:
        self.rearrangement_p = [1 - del_p[i] - ins_p[i] - duplication_p[i] for i in range(2)]

        assert self.rearrangement_p[0] + self.insertion_p[0] + self.deletion_p [0]+ self.duplication_p[0] == 1 and \
               self.rearrangement_p[1] + self.insertion_p[1] + self.deletion_p[1] + self.duplication_p[1] == 1
        # Rearrangement probabilities: # TODO: include as parameter; as of now, only reversals and transloc.
        if num_chr > 1:
            self.reversal_p = .9
            self.translocation_p = .1
        else:
            self.reversal_p = 1
            self.translocation_p = 0


class Simulation:
    def __init__(self, folder, sim_parameters=None):
        self.sim_parameters = sim_parameters
        self.folder = folder


    @staticmethod
    def apply_random_reversal(genome):
        chromosome = np.random.choice(genome.chromosomes)
        bp = sorted(np.random.choice(len(chromosome.gene_order) + 1, 2))
        chromosome.gene_order[bp[0]:bp[1]] = reversed([-x for x in chromosome.gene_order[bp[0]:bp[1]]])
        if chromosome.copy_number is not None:
            chromosome.copy_number[bp[0]:bp[1]] = reversed(chromosome.copy_number[bp[0]:bp[1]])


    @staticmethod
    def apply_random_translocation(genome):
        chromosomes = np.random.choice(genome.chromosomes, 2, replace=False)
        bp1 = np.random.choice(len(chromosomes[0].gene_order))
        bp2 = np.random.choice(len(chromosomes[1].gene_order))
        chromosomes[0].gene_order[bp1:], chromosomes[1].gene_order[bp2:] = \
            chromosomes[1].gene_order[bp2:], chromosomes[0].gene_order[bp1:]
        if chromosomes[0].copy_number is not None and chromosomes[1].copy_number is not None:
            chromosomes[0].copy_number[bp1:], chromosomes[1].copy_number[bp2:] = \
                chromosomes[1].copy_number[bp2:], chromosomes[0].copy_number[bp1:]


    @staticmethod
    def apply_random_segmental_duplication(genome, duplication_length_range, current_copy_number):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(duplication_length_range)
        # update to match COSER sims; lenght is fixed at max;
        #length = max(duplication_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        # print "L", length
        # position:
        position = np.random.choice(range(bp+1) + range(bp + length, chromosome.length()))
        block = chromosome.gene_order[bp:bp + length]
        # update gene copy number
        copy_number_block = []
        for gene in block:
            current_copy_number[abs(gene)] += 1
            copy_number_block.append(current_copy_number[abs(gene)])
        # apply dup:
        chromosome.copy_number[position:position] = copy_number_block
        chromosome.gene_order[position:position] = block


    @staticmethod
    def apply_random_deletion(genome, deletion_length_range):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(deletion_length_range)
        if bp + length > chromosome.length():
            length = chromosome.length() - bp
        chromosome.gene_order[bp:bp + length] = []
        chromosome.copy_number[bp:bp + length] = []
        # remove chromosome if empty:
        if len(chromosome.gene_order) == 0:
            genome.chromosomes.remove(chromosome)


    @staticmethod
    def apply_random_insertion(genome, gene, insertion_length_range, current_copy_number):
        chromosome = np.random.choice(genome.chromosomes)
        bp = np.random.choice(chromosome.length())
        length = np.random.choice(insertion_length_range)
        block = range(gene, gene + length)
        chromosome.gene_order[bp:bp] = block
        if chromosome.copy_number is not None:
            for gene in block:
                current_copy_number[abs(gene)] = 1
            chromosome.copy_number[bp:bp] = [1] * length
        return gene + length


    @staticmethod
    def apply_random_rearrangement(param, genome):
        rearrangement = np.random.choice([RearrangementType.REVERSAL, RearrangementType.TRANSLOCATION], 1,
                                         p=[param.reversal_p, param.translocation_p])
        if rearrangement == RearrangementType.REVERSAL:
            Simulation.apply_random_reversal(genome)
        elif rearrangement == RearrangementType.TRANSLOCATION:
            Simulation.apply_random_translocation(genome)
        else:
            raise RuntimeError("Unknown rearrangement type.")


    @staticmethod
    def apply_random_events(param, genome, current_insertion_gene, current_copy_number=None, postwgd=1):
        insertion_length_range = xrange(1, param.indel_length[postwgd] + 1)
        deletion_length_range = xrange(1, param.indel_length[postwgd] + 1)
        duplication_length_range = xrange(1, param.duplication_length[postwgd] + 1)

        # choose events and apply:
        event_count = {event: 0 for event in EventType.all}
        events = np.random.choice(
            [EventType.REARRANGEMENT, EventType.INSERTION, EventType.DELETION, EventType.DUPLICATION], param.num_ev[postwgd],
            p=[param.rearrangement_p[postwgd], param.insertion_p[postwgd], param.deletion_p[postwgd], param.duplication_p[postwgd]])

        for event in events:  # number of events, can be weighted by 'scaling' parameters
            if event == EventType.REARRANGEMENT:
                Simulation.apply_random_rearrangement(param, genome)

            elif event == EventType.DELETION:
                Simulation.apply_random_deletion(genome, deletion_length_range)
            elif event == EventType.INSERTION:
                current_insertion_gene = Simulation.apply_random_insertion(genome, current_insertion_gene,
                                                               insertion_length_range, current_copy_number)
            elif event == EventType.DUPLICATION:
                Simulation.apply_random_segmental_duplication(genome, duplication_length_range, current_copy_number)

            else:
                raise RuntimeError("Unknown evolutionary event.")
            event_count[event] += 1
        return event_count, current_insertion_gene


    #TODO: MAKE THIS WORK FOR BOTH CIRCULAR AND LINEAR (RIGHT NOW ONLY WORKS FOR LINEAR)
    @staticmethod
    def apply_WGD(genome):
        dup_chromosomes = [chromosome.clone() for chromosome in genome.chromosomes]
        for chromosome in dup_chromosomes:
            genome.chromosomes.append(chromosome)


    def run_simulation(self):

        param = self.sim_parameters
        # insertion and deletions parameters:

        # current insertion genes: (new genes)
        current_insertion_gene = param.num_genes + 1
        current_copy_number = None  # will init at root
        idx = 1

        current_genome = model.Genome.identity(param.num_genes, param.num_chr)
        initial_genome = current_genome.clone()
        initial_genome.name = "B"

        for chromosome in current_genome.chromosomes:
            chromosome.copy_number = [1] * len(chromosome.gene_order)
        current_copy_number = current_genome.gene_count()

        if param.WGD_type != "DD":
            for chromosome in current_genome.chromosomes:
                chromosome.copy_number = [1] * len(chromosome.gene_order)
            current_copy_number = current_genome.gene_count()

            _, current_insertion_gene = \
                            Simulation.apply_random_events(param, current_genome, current_insertion_gene, current_copy_number, postwgd=0)

        prewgd_genome = current_genome.clone()
        prewgd_genome.name = "R"

        Simulation.apply_WGD(current_genome)

        for chromosome in current_genome.chromosomes:
            chromosome.copy_number = [1] * len(chromosome.gene_order)
        current_copy_number = current_genome.gene_count()

        _, current_insertion_gene = \
            Simulation.apply_random_events(param, current_genome, current_insertion_gene, current_copy_number)

        wgd_genome = current_genome.clone()
        wgd_genome.name = "A"

        if param.WGD_type == "DD":
            genomes = [prewgd_genome, wgd_genome]
        elif param.WGD_type == "GH":
            genomes = [initial_genome, prewgd_genome, wgd_genome]

        return genomes

    # Custom simulation on the model of COSER paper, with D DCJs,
    def save_simulation(self, genomes):
        # Output sim result:
        output = self.folder
        if not os.path.exists(output):
            os.makedirs(output)

        # Genomes:
        file_ops.write_genomes_to_file(genomes, os.path.join(output, "blocks.txt"))

        for genome in genomes:
            file_ops.write_genomes_to_file([genome], os.path.join(output, genome.name + ".txt"))

        # Save parameters:
        param = self.sim_parameters
        file_ops.write_simulation_parameters(param, output)


# Main: Generate simulation
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulates rearrangement evolution w/ a WGD event")
    parser.add_argument("-wgd", "--WGD_type", type=str, default="DD", help="Type of WGD event, either 'DD' for double distance or 'GH' for genome halving. ")
    parser.add_argument("-n", "--num_genes", type=int, default=100, help="Number of genes in the root genome.")
    parser.add_argument("-c", "--num_chr", type=int, default=5, help="Number of chromosomes in the root genome.")
    parser.add_argument("-ev", "--num_ev", type=int, default=10, help="Number of rearrangement events before WGD.")
    parser.add_argument("-o", "--output", type=str, default="sim", help="Name of the output folder.")
    parser.add_argument("-dp", "--deletion_p", type=float, default=0.0, help="Percentage of deletions before WGD, from 0 to 1.0")
    parser.add_argument("-ip", "--insertion_p", type=float, default=0.0, help="Percentage of insertions before WGD, from 0 to 1.0")
    parser.add_argument("-dl", "--duplication_length", type=int, default=5, help="Maximum length of duplication events pre WGD.")
    parser.add_argument("-dup_p", "--duplication_p", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0 pre WGD")
    parser.add_argument("-dl2", "--duplication_length2", type=int, default=5,
                        help="Maximum length of duplication event post WGD.")
    parser.add_argument("-dup_p2", "--duplication_p2", type=float, default=0.0,
                        help="Percentage of duplications, from 0 to 1.0 post WGD")
    parser.add_argument("-il", "--indel_length", type=int, default=5, help="Maximum size of indel event in genes (pre-WGD).")
    parser.add_argument("-ev2", "--num_ev2", type=int, default=10, help="Number of rearrangement events after WGD.")
    parser.add_argument("-dp2", "--deletion_p2", type=float, default=0.0,
                        help="Percentage of deletions after WGD, from 0 to 1.0")
    parser.add_argument("-ip2", "--insertion_p2", type=float, default=0.0,
                        help="Percentage of insertions after WGD, from 0 to 1.0")
    parser.add_argument("-il2", "--indel_length2", type=int, default=5,
                        help="Maximum size of indel event in genes after WGD.")
    param = parser.parse_args()

    # Simulation parameters:
    sim_par = SimParameters(wgd_type=param.WGD_type, num_genes=param.num_genes, num_chr=param.num_chr,
                            del_p=[param.deletion_p, param.deletion_p2], ins_p=[param.insertion_p, param.insertion_p2],
                            indel_length=[param.indel_length, param.indel_length2 ], duplication_p=[param.duplication_p, param.duplication_p2],
                            duplication_length=[param.duplication_length, param.duplication_length2], events=[param.num_ev, param.num_ev2])
    # start sim object;
    sim = Simulation(param.output, sim_par)

    genomes = sim.run_simulation()

    sim.save_simulation(genomes)