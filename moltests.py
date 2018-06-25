from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.mllib.linalg.distributed import RowMatrix
from pyspark.mllib.linalg import Vectors
from pyspark.sql import SparkSession, Row
from pyspark.sql import SQLContext

from helpers import *
from helpers import *

EVALUATION = False
FLATMAPMETHOD = True
CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds82.smi'
#dataFile = '../mols/merged/Reninmerged.smi'

def main():
    sc = start_spark()

    compounds = load_data(sc, dataFile)
    #example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)

    fingerprints = select_fingerprints(compounds)

    fingerprints1 = fingerprints.partitionBy(4).cache()

    bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())

    #neighbours = calculate_neighbours(fingerprints1, bc_fingerprints).cache()

    #mol_neighbour_count, neighbours_with_counts = get_neighbours_counts(neighbours)

    neighbour_similarities = calculate_neighbours_similarities(fingerprints1, bc_fingerprints).cache()

    #neighbour_similarities.foreach(output)
    print("----------------------------")
    molecules_appearances = neighbour_similarities.map(lambda (mol_1,mol_2): (mol_1,1))
    #molecules_appearances.foreach(output)
    mol_neighbour_count = molecules_appearances.reduceByKey(lambda a, b: a + b)
    print ("-----------------------------------")
    #mol_neighbour_count.foreach(output)
    bc_mol_neighbour_count = sc.broadcast(mol_neighbour_count.collectAsMap())

    complete_similarities_fp = neighbour_similarities.map(lambda (a, b): (a, convert_single_neighbour_dict(b, bc_mol_neighbour_count.value),convert_owner_dict(a, bc_mol_neighbour_count.value)))
    complete_similarities_fp.foreach(output)
    complete_similarities_fp = complete_similarities_fp.filter(lambda (a,b,c): c < b[1] or (c == b[1] and a >= b[0])).map(lambda (a,b,c): (a,set([b])))
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    complete_similarities_fp.foreach(output)
    neighbours = complete_similarities_fp.reduceByKey(lambda a, b: a.union(b))
    neighbours = neighbours.map(lambda (mol, nbrs): (mol, sorted(nbrs, key=lambda x: (-x[1], x[0])))).cache()
    cluster_assignment = neighbours.map(lambda (a, b): (a,b,convert_owner_dict(a,bc_mol_neighbour_count.value),-1))
    cluster_assignment.foreach(output)
    #print(bc_mol_neighbour_count.value)

    #cluster_assignment = neighbours_with_counts.map(lambda (a,b,c): (a, convert_neighbours_dict(b, bc_mol_neighbour_count.value), c, -1))
    #cluster_assignment.foreach(output)
    print ("------------------------------------------")
    #assign cluster

    bc_invalid_clusters = sc.broadcast([])
    #old_cluster_count = 0
    #need_update = True
    bc_valid_condition = sc.broadcast(True)
    iteration = 0
    while bc_valid_condition.value:
    #for x in range(1,10):

        print ("----------------------------------ITERATION ",iteration ,"----------------------------------------------")
        cluster_assignment = cluster_assignment.map(lambda (mol_id,nb_count_tuple,mol_nbr_count,cluster_assigned) :
                                                    (mol_id,nb_count_tuple,mol_nbr_count,assign_cluster(mol_id, mol_nbr_count, nb_count_tuple, bc_invalid_clusters.value, cluster_assigned))).cache()

        #cluster_assignment.foreach(output)

        #valid = cluster_assignment.filter(lambda (a,b,c,d): a == d)
        valid = cluster_assignment.filter(lambda (a, b, c, d): a == d).map(lambda (a, b, c, d): a)

        bc_valids = sc.broadcast(valid.collect)
        #valid.foreach(output)
        print("Number of valids:", valid.count())

        invalid_clusters = cluster_assignment.filter(lambda (a,b,c,d): a == d)\
            .map(lambda (a,b,c,d): [mol[0] if mol[0] != a else None for mol in b])\
            .flatMap(lambda list: list)\
            .filter(lambda a : a is not None)

        invalid_clusters = invalid_clusters.map(lambda a:(a,1))

        #invalid_clusters.distinct().foreach(output)

        #dist_invalids = invalid_clusters.distinct().collect()
        dist_invalids = invalid_clusters.collectAsMap()

        print("Invalids:", dist_invalids.__len__())
        bc_invalid_clusters = sc.broadcast(dist_invalids)

        if dist_invalids.__len__() == 0:
            bc_valid_condition = sc.broadcast(False)
        else:
            bc_valid_condition = sc.broadcast(True)

        cluster_assignment = cluster_assignment.filter(lambda (mol_id,nbr_counts,mol_count, cluster_assigned): mol_id not in bc_invalid_clusters.value).cache()
        #cluster_assignment.foreach(output)
        cluster_assignment = cluster_assignment.map(lambda (mol_id, nbr_counts, mol_count, cluster_assigned): (mol_id, remove_invalid_nbrs_dict(nbr_counts, bc_invalid_clusters.value, mol_id == cluster_assigned), mol_count, cluster_assigned))
        iteration+=1

    if True:
        #remove tuple of neighbours, count
        valid_clusters = cluster_assignment.map(
            lambda (mol_id, nbr_counts, mol_count, cluster_assigned): (mol_id, 1))

        #valid_clusters.foreach(output)
        bc_valid_clusters = sc.broadcast(valid_clusters.collectAsMap())
        complete_clusters = neighbours.filter(lambda(mol_id, nbrs):mol_id in bc_valid_clusters.value)

        print("Number of clusters:", complete_clusters.count())

        cluster_combs = complete_clusters.cartesian(complete_clusters)\
            .filter(lambda (a,b): (len(a[1]) == len(b[1]) and a[0] >= b[0]) or (len(a[1]) < len(b[1])))
        #cluster_combs.foreach(output)
        cluster_combs = cluster_combs.map(lambda (a,b): (a[0], a[1].difference(b[1])) if a[0] != b[0] else (a[0], a[1]))

        cluster_combs = cluster_combs.reduceByKey(lambda a,b: a.intersection(b))
        #cluster_combs.foreach(output)
        print("Number of clusters:", valid_clusters.count())
        print ("--------------------------------------")
        #cluster_combs.foreach(output)
        #if EVALUATION:
        #    cluster_combs = cluster_combs.sortBy(lambda (a,b): a)
        #    cluster_combs1 = cluster_combs.map(lambda (a,b): (sort_list(list(b))))
        #RESULT
        #cluster_combs1.foreach(output)
        #output_results(sc, cluster_combs, compounds, EVALUATION)

def main_old():
    sc = start_spark()

    compounds = load_data(sc, dataFile)
    #example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)

    fingerprints = select_fingerprints(compounds)

    fingerprints1 = fingerprints.partitionBy(4).cache()

    bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())

    neighbours = calculate_neighbours(fingerprints1, bc_fingerprints).cache()

    mol_neighbour_count, neighbours_with_counts = get_neighbours_counts(neighbours)

    bc_mol_neighbour_count = sc.broadcast(mol_neighbour_count.collectAsMap())

    #print(bc_mol_neighbour_count.value)

    #neighbours.foreach(output)
    cluster_assignment = neighbours_with_counts.map(lambda (a,b,c): (a, convert_neighbours_dict(b, bc_mol_neighbour_count.value), c, -1))
    #cluster_assignment.foreach(output)
    print ("------------------------------------------")
    #assign cluster

    bc_invalid_clusters = sc.broadcast([])
    #old_cluster_count = 0
    #need_update = True
    bc_valid_condition = sc.broadcast(True)
    iteration = 0
    while bc_valid_condition.value:
    #for x in range(1,10):
        print ("----------------------------------ITERATION ",iteration ,"----------------------------------------------")
        cluster_assignment = cluster_assignment.map(lambda (mol_id,nb_count_tuple,mol_nbr_count,cluster_assigned) :
                                                    (mol_id,nb_count_tuple,mol_nbr_count,assign_cluster(mol_id, mol_nbr_count, nb_count_tuple, bc_invalid_clusters.value, cluster_assigned))).cache()

        #cluster_assignment.foreach(output)

        valid = cluster_assignment.filter(lambda (a,b,c,d): a == d)
        #valid.foreach(output)
        print("Number of valids:", valid.count())

        invalid_clusters = cluster_assignment.filter(lambda (a,b,c,d): a == d)\
            .map(lambda (a,b,c,d): [mol[0] if mol[0] != a else None for mol in b])\
            .flatMap(lambda list: list)\
            .filter(lambda a : a is not None)

        invalid_clusters = invalid_clusters.map(lambda a:(a,1))

        #invalid_clusters.distinct().foreach(output)

        #dist_invalids = invalid_clusters.distinct().collect()
        dist_invalids = invalid_clusters.collectAsMap()

        print("Invalids:", dist_invalids.__len__())
        bc_invalid_clusters = sc.broadcast(dist_invalids)

        if dist_invalids.__len__() == 0:
            bc_valid_condition = sc.broadcast(False)
        else:
            bc_valid_condition = sc.broadcast(True)

        cluster_assignment = cluster_assignment.filter(lambda (mol_id,nbr_counts,mol_count, cluster_assigned): mol_id not in bc_invalid_clusters.value).cache()
        #cluster_assignment.foreach(output)
        cluster_assignment = cluster_assignment.map(lambda (mol_id, nbr_counts, mol_count, cluster_assigned): (mol_id, remove_invalid_nbrs_dict(nbr_counts, bc_invalid_clusters.value, mol_id == cluster_assigned), mol_count, cluster_assigned))
        iteration+=1

    #remove tuple of neighbours, count
    #valid_clusters = cluster_assignment.map(lambda (mol_id,nbr_counts,mol_count, cluster_assigned): (mol_id, [x for x,y in nbr_counts]))
    valid_clusters = cluster_assignment.map(
        lambda (mol_id, nbr_counts, mol_count, cluster_assigned): (mol_id, 1))

    #valid_clusters.foreach(output)
    bc_valid_clusters = sc.broadcast(valid_clusters.collectAsMap())
    complete_clusters = neighbours.filter(lambda(mol_id, nbrs):mol_id in bc_valid_clusters.value)

    #complete_clusters.foreach(output)

    cluster_combs = complete_clusters.cartesian(complete_clusters)\
        .filter(lambda (a,b): (len(a[1]) == len(b[1]) and a[0] >= b[0]) or (len(a[1]) < len(b[1])))

    cluster_combs = cluster_combs.map(lambda (a,b): (a[0], set(a[1]).difference(set(b[1]))) if a[0] != b[0] else (a[0], set(a[1])))

    cluster_combs = cluster_combs.reduceByKey(lambda a,b: a.intersection(b))
    #cluster_combs.foreach(output)

    print("Number of clusters:", cluster_combs.count())
    print ("--------------------------------------")
    #cluster_combs.foreach(output)
    #if EVALUATION:
    #    cluster_combs = cluster_combs.sortBy(lambda (a,b): a)
    #    cluster_combs1 = cluster_combs.map(lambda (a,b): (sort_list(list(b))))
    #RESULT
    #cluster_combs1.foreach(output)
    #output_results(sc, cluster_combs, compounds, EVALUATION)


def start_spark():
    server_location = read_server()
    #conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
    conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
    sc = SparkContext(conf=conf)

    sc.addPyFile("helpers.py")
    return sc


def load_data(sc, dataFile):
    lines = sc.textFile(dataFile)

    smiles = lines.map(convertToBitVectorFP)

    print("Read All data, length: ", smiles.count())

    # remove file headers
    smiles = smiles.filter(lambda l: l[0] != "smiles" and l[0] != "SMILES")

    print("Filtered data, length: ", smiles.count())

    # assign index to each compound
    compounds = smiles.zipWithIndex()
    return compounds


def select_fingerprints(compounds):
    fingerprints = compounds.map(lambda (x, idx): (idx, x[2]))
    return fingerprints


def calculate_neighbours(fingerprints, bc_fingerprints):
    # create combinations for all the fingerprints, if condition is used to calculate only upper triangle of the similarity matrix
    cartesian_fp = fingerprints.flatMap(
        lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x >= y[0]])

    # keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
    similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1]))) \
        .filter(lambda (a, b, c): c >= SIMILARITY_THRESHOLD)

    if FLATMAPMETHOD:
        complete_similarities_fp = similarities_fp.flatMap(lambda (a, b, c): [(a, set([b])), (b, set([a]))])
    else:
        # a new rdd is created to get the lower trangle of the similarity matrix
        inverted_similarities = similarities_fp.map(lambda (a, b, c): (b, a, c))
        # the complete catrix is obtained through union
        complete_similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))

    # combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
    neighbours = complete_similarities_fp.reduceByKey(lambda a, b: a.union(b))

    return neighbours

def calculate_neighbours_similarities(fingerprints, bc_fingerprints):
    # create combinations for all the fingerprints, if condition is used to calculate only upper triangle of the similarity matrix
    cartesian_fp = fingerprints.flatMap(
        lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x >= y[0]])

    # keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
    similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1]))) \
        .filter(lambda (a, b, c): c >= SIMILARITY_THRESHOLD)

    if FLATMAPMETHOD:
        complete_similarities_fp = similarities_fp.flatMap(lambda (a, b, c): [(a,b)] if a == b else [(a, b), (b, a)])
    else:
        # a new rdd is created to get the lower trangle of the similarity matrix
        inverted_similarities = similarities_fp.map(lambda (a, b, c): (b, a, c))
        # the complete catrix is obtained through union
        complete_similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))

    return complete_similarities_fp


def get_neighbours_counts(neighbours):
    neighbours = neighbours.map(lambda (a,b):(a,b,len(b)))
    mol_neighbour_count = neighbours.map(lambda (a, b, c): (a, c))
    return mol_neighbour_count, neighbours


def output_results(sc, cluster_combs, compounds, evaluate):
    if evaluate:
        #cluster_combs = cluster_combs.sortBy(lambda (a,b): a)
        #cluster_combs1 = cluster_combs.map(lambda (a,b):(sort_list(list(b))))

        cluster_combs.collectAsMap()
        #invert id and molecules
        compounds = compounds.map(lambda (x,y): (y,x))
        bc_compounds = sc.broadcast(compounds.collectAsMap())
        print("Here")
        #from the clusters of indeces, get the actual molecule
        #ie. clusters of molecules
        cluster_mols = cluster_combs.map(lambda (cl_id,compound_ids): index_to_mol(compound_ids, compound_list=bc_compounds))
        print("Here")
        #assign ids to clusters
        cluster_groups = cluster_mols.zipWithIndex()
        print("Here")
        #flat the molecules in the format: Smiles string, Molecule Name, Cluster Ids
        cluster_groups = cluster_groups.flatMap(lambda (r, idx): [(Chem.MolToSmiles(cmp), cmp.GetProp("_Name"), idx) for cmp in r])
        print("Here")
        #Change molecule data to string
        cluster_groups = cluster_groups.map(lambda (r, name, idx): r + ' ' + name + ' ' + str(idx))

        #output to file
        print("Here")
        header = sc.parallelize(["smiles Name Cluster"])
        output_file = cluster_groups.union(header)
        #output_file.foreach(output)
        output_file.saveAsTextFile("../mols/resultsSpark/result")


if __name__ == '__main__':
    main_old()