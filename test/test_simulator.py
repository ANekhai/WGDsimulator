from src import simulator

import pyximport
pyximport.install()
from utils import model, dcj

test_sim = simulator.Simulation(folder=None)


# TODO: FIGURE OUT WHETHER OR NOT COPY NUMBER IS CORRECT, MAY HAVE TO MODIFY SOME TESTS FOR THIS IN THE NEAR FUTURE

# UTILITY FUNCTIONS

def gen_genome(genes, chromosomes):
    genome = model.Genome.identity(genes, chromosomes)
    for chromosome in genome.chromosomes:
        chromosome.copy_number = [1] * len(chromosome.gene_order)
    return genome


def get_chr(genome, number):
    return genome.chromosomes[number].gene_order


def get_copy_number(genome, number):
    return genome.chromosomes[number].copy_number

#################################
# BASIC FUNCTIONALITY TESTS #####
#################################


# TODO: FIX THIS TEST ONCE COPY NUMBER IS FIXED IN SIMULATOR CODE
def test_wgd_on_single_chromosome():
    genome = gen_genome(3, 1)
    test_sim.apply_WGD(genome)
    assert str(genome.chromosomes[0]) == str(genome.chromosomes[1])

    # assert 2 * sum(old_get_chr(genome, 0].copy_number) == sum(get_chr(genome, 0].copy_number)
    # assert 2 * sum(old_get_chr(genome, 0].copy_number) == sum(get_chr(genome, 1].copy_number)


def test_wgd_on_multiple_chromosomes():
    genome = gen_genome(10, 2)
    test_sim.apply_WGD(genome)
    assert str(genome.chromosomes[0]) == str(genome.chromosomes[2])
    assert str(genome.chromosomes[1]) == str(genome.chromosomes[3])


# REARRANGEMENT TESTS
def test_reversal_on_single_chromosome():
    genome = gen_genome(10, 1)
    test_sim.apply_random_reversal(genome)

    assert len(get_chr(genome, 0)) == sum(get_copy_number(genome, 0))
    assert len([i for i in get_chr(genome, 0) if i < 0]) > 0


def test_reversal_on_multiple_chromosomes():
    genome = gen_genome(10, 2)
    num_genes = 5
    test_sim.apply_random_reversal(genome)
    assert len(get_chr(genome, 0)) == num_genes
    assert len(get_chr(genome, 1)) == num_genes
    assert len([i for i in get_chr(genome, 0) if i < 0]) > 0 or \
           len([j for j in get_chr(genome, 1) if j < 0]) > 0


# TODO: translocations can just exchange two chromosomes (no event occurs)
def test_translocation_on_genome():
    genome = gen_genome(10, 2)
    test_sim.apply_random_translocation(genome)

    assert len(get_chr(genome, 0)) + len(get_chr(genome, 1)) == 10


def test_segmental_duplication_within_chromosome():
    genome = gen_genome(10, 1)
    dup_length = 3
    test_sim.apply_random_segmental_duplication(genome, [dup_length], genome.gene_count(), 0)

    assert get_copy_number(genome, 0).count(2) <= 2 * dup_length
    assert len(get_chr(genome, 0)) > 10


def test_segmental_duplication_between_chromosomes():
    genome = gen_genome(10, 2)
    test_sim.apply_random_segmental_duplication(genome, [3], genome.gene_count(), 1)

    assert get_copy_number(genome, 0).count(2) == get_copy_number(genome, 1).count(2)
    assert abs(len(get_chr(genome, 0)) - len(get_chr(genome, 1))) > 0


# TODO: tests for cases when seg dup larger than chromosome


def test_deletions_on_single_chromosome():
    genome = gen_genome(10, 1)
    num_genes = 10
    test_sim.apply_random_deletion(genome, [3])

    assert len(get_chr(genome, 0)) < num_genes
    assert len(get_chr(genome, 0)) >= num_genes - 3


def test_deletions_on_multiple_chromosomes():
    genome = gen_genome(10, 2)
    test_sim.apply_random_deletion(genome, [3])

    assert len(get_chr(genome, 0)) != len(get_chr(genome, 1))
    assert len(get_chr(genome, 0)) <= 5 and len(get_chr(genome, 1)) <= 5


def test_insertions_on_single_chromosome():
    genome = gen_genome(10, 1)
    num_genes = 10
    curr_copy_num = genome.gene_count()
    test_sim.apply_random_insertion(genome, num_genes + 1, [3], curr_copy_num)

    assert len(get_chr(genome, 0)) == num_genes + 3
    assert len(curr_copy_num.keys()) > num_genes


def test_insertions_on_multiple_chromosomes():
    genome = gen_genome(10, 2)
    curr_copy_num = genome.gene_count()
    test_sim.apply_random_insertion(genome, 11, [3], curr_copy_num)

    assert len(get_chr(genome, 0)) != len(get_chr(genome, 1))
    assert len(get_chr(genome, 0)) == 8 or len(get_chr(genome, 1)) == 8
    assert len(curr_copy_num.keys()) > 10


def test_running_events_on_genome():
    genome = gen_genome(100, 5)
    curr_copy_num = genome.gene_count()
    param = simulator.SimParameters(num_genes=100, num_chr=5, events=[100, 0])

    event_count, _ = test_sim.apply_random_events(param, genome, 101, curr_copy_num, postwgd=0)

    assert len(curr_copy_num.keys()) == 100
    assert event_count["rearrangement"] == 100


def test_simulation_parameters_with_no_rearrangements():
    param = simulator.SimParameters(ins_p=[.4, 0], del_p=[.4, 0], duplication_p=[.2, 0])


#####################################
# PIPELINE TESTS ####################
#####################################


def test_pipeline_no_wgd_only_rearrangements():
    param = simulator.SimParameters(wgd_type=None, num_genes=40, num_chr=2, events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) == sum(genomes[1].gene_count().values())


def test_pipeline_no_wgd_only_insertions():
    param = simulator.SimParameters(wgd_type=None, num_genes=40, num_chr=2, ins_p=[1, 0], indel_length=[3, 0], events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) < len(genomes[1].gene_count().keys())


def test_pipeline_no_wgd_only_deletions():
    param = simulator.SimParameters(wgd_type=None, num_genes=40, num_chr=2, del_p=[1, 0], indel_length=[3, 0], events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert len(genomes[0].chromosomes) >= len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) > len(genomes[1].gene_count().keys())


def test_pipeline_no_wgd_with_only_segmental_duplication():
    param = simulator.SimParameters(wgd_type=None, num_genes=40, num_chr=2, duplication_p=[1, 0], duplication_length=[3, 0], events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) < sum(genomes[1].gene_count().values())


def test_pipeline_no_wgd_with_only_inter_chromosomal_segmental_duplication():
    param = simulator.SimParameters(wgd_type=None, num_genes=40, num_chr=2, duplication_p=[1, 0], duplication_length=[3, 0], events=[10, 0], inter_p=1)
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) < sum(genomes[1].gene_count().values())


def test_pipeline_double_distance_only_rearrangements():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert 2 * sum(genomes[0].gene_count().values()) == sum(genomes[1].gene_count().values())


def test_pipeline_double_distance_only_insertions():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, ins_p=[0, 1], indel_length=[0, 3], events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) < len(genomes[1].gene_count().keys())


#TODO: FIX STRANGE BUG IN CODE WHERE NO GENES ARE DELETED?
def test_pipeline_double_distance_only_deletions():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, del_p=[0, 1], indel_length=[0, 3],
                                    events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) >= len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) >= len(genomes[1].gene_count().keys())
    assert 2 * sum(genomes[0].gene_count().values()) > sum(genomes[1].gene_count().values())

def test_pipeline_double_distance_only_segmental_duplications():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, duplication_p=[0, 1],
                                    duplication_length=[0, 3], events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) < sum(genomes[1].gene_count().values())


def test_pipeline_double_distance_with_only_inter_chromosomal_segmental_duplication():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, duplication_p=[0, 1],
                                    duplication_length=[0, 3], events=[0, 10], inter_p=1)
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) < sum(genomes[1].gene_count().values())


def test_double_distance_events_before_wgd_event():
    param = simulator.SimParameters(wgd_type="DD", num_genes=40, num_chr=2, del_p=[.3, 0], ins_p=[.3, 0], indel_length=[3, 0], duplication_p=[.2, 0],
                                    duplication_length=[3, 0], events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)
    genomes = sim.run_simulation()

    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert 2 * sum(genomes[0].gene_count().values()) == sum(genomes[1].gene_count().values())
    assert str(genomes[0].chromosomes[0]) == str(genomes[1].chromosomes[2])


def test_genome_halving_only_rearrangements_pre_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[2].gene_count().keys())
    assert str(genomes[1].chromosomes[0]) == str(genomes[2].chromosomes[2])


def test_genome_halving_only_insertions_pre_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, ins_p=[1, 0], indel_length=[3, 0],
                                    events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) < len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) == len(genomes[2].gene_count().keys())


def test_genome_halving_only_deletions_pre_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, del_p=[1, 0], indel_length=[3, 0],
                                    events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) >= len(genomes[1].chromosomes)
    assert len(genomes[0].gene_count().keys()) > len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) == len(genomes[2].gene_count().keys())


def test_genome_halving_only_segmental_duplications_pre_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, duplication_p=[1, 0],
                                    duplication_length=[3, 0], events=[10, 0])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) == len(genomes[2].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) < sum(genomes[1].gene_count().values())
    assert 2 * sum(genomes[1].gene_count().values()) == sum(genomes[2].gene_count().values())


def test_genome_halving_only_rearrangements_post_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[2].gene_count().keys())
    assert str(genomes[0].chromosomes[0]) == str(genomes[1].chromosomes[0]) and \
           str(genomes[0].chromosomes[1]) == str(genomes[1].chromosomes[1])


def test_genome_halving_only_insertions_post_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, ins_p=[0, 1], indel_length=[0, 3],
                                    events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) < len(genomes[2].gene_count().keys())


def test_genome_halving_only_deletions_post_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, del_p=[0, 1], indel_length=[0, 3],
                                    events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) >= len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) >= len(genomes[2].gene_count().keys())
    assert 2 * sum(genomes[1].gene_count().values()) > sum(genomes[2].gene_count().values())


def test_genome_halving_only_segmental_duplications_post_wgd():
    param = simulator.SimParameters(wgd_type="GH", num_genes=40, num_chr=2, duplication_p=[0, 1],
                                    duplication_length=[0, 3], events=[0, 10])
    sim = simulator.Simulation(folder=None, sim_parameters=param)

    genomes = sim.run_simulation()

    assert 2 * len(genomes[0].chromosomes) == len(genomes[2].chromosomes)
    assert len(genomes[0].gene_count().keys()) == len(genomes[1].gene_count().keys())
    assert len(genomes[1].gene_count().keys()) == len(genomes[2].gene_count().keys())
    assert sum(genomes[0].gene_count().values()) == sum(genomes[1].gene_count().values())
    assert 2 * sum(genomes[1].gene_count().values()) < sum(genomes[2].gene_count().values())

