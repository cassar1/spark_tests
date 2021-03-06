from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.mllib.linalg import Vectors
from pyspark.sql import SparkSession, Row
from pyspark.ml.feature import MinHashLSH
from pyspark.ml.linalg import Vectors
import pyspark.sql.functions as f
import sys
import time
import datetime

EVALUATION = False
FLATMAPMETHOD = True
CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
LSH_SIMILARITY = True
#dataFile = '../mols/compounds18.smi'
baseFile = '../mols/'
#baseFile = 'wasb://molecules-container@jurgenhdinsightstorage.blob.core.windows.net/'
#baseOutputFile = 'wasb://cluster-results@jurgenhdinsightstorage.blob.core.windows.net/'
executor_num = 1
#dataFile = '../mols/merged/Reninmerged.smi'

#conf = SparkConf().setAppName("DistributedClustering")
conf = SparkConf() \
    .setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf=conf)
spark = SparkSession.builder.config(conf = conf).getOrCreate()
sc.addPyFile("helpers.py")

from helpers import *

def main():
    print ("START ",datetime.datetime.now())
    compounds = load_data(sc, dataFile)
    #example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)
    compounds.partitionBy(executor_num)
    fingerprints = select_fingerprints(compounds)

    fingerprints1 = fingerprints
    #fingerprints1 = fingerprints.partitionBy(executor_num).cache()

    #fingerprints1.foreach(output)
    bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())
    if not LSH_SIMILARITY:
        neighbour_similarities = calculate_neighbours_similarities(fingerprints1, bc_fingerprints).cache()
    else:
        #vector_fps = fingerprints_to_vectors(fingerprints1)
        #bc_fingerprints = sc.broadcast(vector_fps.collectAsMap())
        #neighbour_similarities = calculate_LSH_neighbours(fingerprints1).cache()
        neighbour_similarities = calculate_LSH_custom(bc_fingerprints, fingerprints).cache()

    #neighbour_similarities.foreach(output)
    #print (neighbour_similarities.count())
    # for each molecule set (mol_id,1) to then reduce and finally get the final number
    molecules_appearances = neighbour_similarities.map(lambda (mol_1,mol_2): (mol_1,1))

    mol_neighbour_count = molecules_appearances.reduceByKey(lambda a, b: a + b)

    bc_mol_neighbour_count = sc.broadcast(mol_neighbour_count.collectAsMap())
    #counts = neighbour_similarities.countByKey()

    # change from (mol_id, nbr_id) to (mol_id, (nbr_id, nbr_count), nbr_count)
    complete_similarities_fp = neighbour_similarities.map(lambda (a, b): (a, convert_single_neighbour_dict(b, bc_mol_neighbour_count.value),convert_owner_dict(a, bc_mol_neighbour_count.value)))


#or (mol_nbr_count == nbr[1] and mol_id >= nbr[0])

    #this line is used to remove neighbours with nbr count smaller than the molecule
    #subset_similarities_fp = complete_similarities_fp.filter(lambda (mol_id, nbr, mol_nbr_count): mol_nbr_count <= nbr[1])\
    #    .map(lambda (a,b,c): (a,set([b])))

    subset_similarities_fp = complete_similarities_fp.map(lambda (a, b, c): (a, set([b])))


    neighbours = subset_similarities_fp.reduceByKey(lambda a, b: a.union(b))

    neighbours = neighbours.map(lambda (mol, nbrs): (mol, sorted(nbrs, key=lambda x: (-x[1], x[0])))).cache()
    #neighbours.foreach(output)
    # return (mol_id, neighbours, nbr_count of mol_id, -1)
    cluster_assignment = neighbours.map(lambda (a, b): (a,b,convert_owner_dict(a, bc_mol_neighbour_count.value),-1))

    # assign cluster
    bc_invalid_clusters = sc.broadcast([])

    bc_valid_condition = sc.broadcast(True)
    iteration = 0
    total_valids = sc.parallelize([])
    while bc_valid_condition.value:
        print ("----------------------------------ITERATION ",iteration ,"----------------------------------------------")
        cluster_assignment = cluster_assignment.map(lambda (mol_id,nb_count_tuple,mol_nbr_count,cluster_assigned) :
                                                    (mol_id,nb_count_tuple,mol_nbr_count,assign_cluster(mol_id, mol_nbr_count, nb_count_tuple, bc_invalid_clusters.value, cluster_assigned))).cache()

        #print("----------------------------ASSIGNMENT--------------------------------")
        valid = cluster_assignment.filter(lambda (a, b, c, d): a == d).map(lambda (a, b, c, d): (a, 1))
        total_valids = total_valids.union(valid)

        bc_valids = sc.broadcast(valid.collectAsMap())


        #invalid_clusters = complete_similarities_fp.filter(lambda (mol, nbr, count): mol in bc_valids.value) \
        #    .map(lambda (a, nbr, c): nbr[0] if nbr[0] != a else None) \
        #    .filter(lambda a : a is not None)

        #cluster_assignment.foreach(output)
        #complete_similarities_fp.foreach(output)

        invalid_clusters = cluster_assignment.filter(lambda (mol, nbr, count, assignment): mol in bc_valids.value) \
            .map(lambda (a, nbr, c, assign): [i[0] for i in nbr if i[0] != a]) \
            .flatMap(lambda ls: [l for l in ls])\
            .filter(lambda a: a is not None)

        #print("VALIDS", bc_valids.value)
        print("-----------------")
        #invalid_clusters.foreach(output)
        invalid_clusters = invalid_clusters.map(lambda a:(a,1))

        dist_invalids = invalid_clusters.collectAsMap()#.distinct()


        bc_invalid_clusters = sc.broadcast(dist_invalids)

        #print("ALL DICT", bc_invalid_clusters.value)
        print("NUMBER INVALIDS ", len(bc_invalid_clusters.value))

        if dist_invalids.__len__() == 0:
            bc_valid_condition = sc.broadcast(False)
        else:
            bc_valid_condition = sc.broadcast(True)

        # remove molecule clusters that are invalid
        cluster_assignment = cluster_assignment.filter(lambda (mol_id, nbr_counts, mol_count, cluster_assigned): mol_id not in bc_invalid_clusters.value).cache()
        # remove molecules in clusters that are invalid
        cluster_assignment = cluster_assignment.map(lambda (mol_id, nbr_counts, mol_count, cluster_assigned): (mol_id, [x for x in nbr_counts if x[0] not in bc_invalid_clusters.value], mol_count, cluster_assigned)).cache()
        #remove clusters that are valid
        cluster_assignment = cluster_assignment.filter(lambda (mol_id, nbr_counts, mol_count, cluster_assigned): mol_id not in bc_valids.value).cache()
        iteration += 1
        #cluster_assignment.foreach(output)

    bc_valid_clusters = sc.broadcast(total_valids.collectAsMap())

    cluster_centers_single_neighbours = complete_similarities_fp.filter(lambda(mol_id, nbr, nbr_counts):mol_id in bc_valid_clusters.value)\
        .map(lambda (a, b, c): (a, set([b[0]])))
    cluster_neighbours = cluster_centers_single_neighbours.reduceByKey(lambda a, b: a.union(b))
    #cluster_neighbours.partitionBy(executor_num)


    cluster_combs = cluster_neighbours.cartesian(cluster_neighbours)\
        .filter(lambda (a,b): (len(a[1]) == len(b[1]) and a[0] >= b[0]) or (len(a[1]) < len(b[1])))

    cluster_combs = cluster_combs.map(lambda (a,b): (a[0], set(a[1]).difference(set(b[1]))) if a[0] != b[0] else (a[0], set(a[1])))

    cluster_combs = cluster_combs.reduceByKey(lambda a,b: a.intersection(b))
    print("End ", datetime.datetime.now())
    cluster_combs = cluster_combs.sortByKey()
    print("number of cluster", cluster_combs.count())
    #cluster_combs.foreach(output)

    #if EVALUATION:
    #    cluster_combs = cluster_combs.sortBy(lambda (a,b): a)
    #    cluster_combs1 = cluster_combs.map(lambda (a,b): (sort_list(list(b))))


    #cluster_combs = cluster_combs.map(lambda (cl_id, compound_ids): (cl_id, [c[0] for c in compound_ids]))
    output_results(sc, cluster_combs, compounds, EVALUATION)

def start_spark():
    #server_location = read_server()
    conf = SparkConf()\
        .setAppName("DistributedClustering")
        #.setMaster("local")\

    #conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
    sc = SparkContext(conf=conf)

    sc.addPyFile("helpers.py")

    #from helpers import *
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
    BULK = True

    if BULK:
        cartesian_fp = fingerprints.map(
            lambda y: (y, bc_fingerprints.value.values()[0:y[0]+1]))
        #cartesian_fp.foreach(output)
        # keep only those that have similarity >= threshold a = fp1, b = fp2, c = similarity
        similarities_fp = cartesian_fp.map(lambda (a, b): (a, b, calculate_tanimoto_bulk(a[1], b)))\
            .flatMap(lambda (a,b,c): [(a[0], idx, c[idx]) for idx, x in enumerate(b)]) \
            .filter(lambda (a, b, c): c >= SIMILARITY_THRESHOLD)
    else:
        cartesian_fp = fingerprints.flatMap(lambda y: [(y,(x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x >= y[0]])
        similarities_fp = cartesian_fp.map(lambda (a,b): (a[0], b[0], calculate_tanimoto(a[1], b[1]))).filter(lambda (a,b,c) : c >= SIMILARITY_THRESHOLD)

    #similarities_fp.foreach(output)
    if FLATMAPMETHOD:
        complete_similarities_fp = similarities_fp.flatMap(lambda (a, b, c): [(a,b)] if a == b else [(a, b), (b, a)])
    else:
        # a new rdd is created to get the lower triangle of the similarity matrix
        inverted_similarities = similarities_fp.map(lambda (a, b, c): (b, a, c))
        # the complete matrix is obtained through union
        complete_similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))

    return complete_similarities_fp

def fingerprints_to_vectors(fingerprints):
    print("fingerprints_to_Vectors")
    vector_fingerprints = fingerprints.map(lambda (idx, x): (
        idx, Vectors.sparse(len(x), [(index, value) for (index, value) in enumerate(x) if value != 0]))).cache()
    return vector_fingerprints


def calculate_LSH_custom(bc_fingerprints, fingerprints):
    num_hashes = 100

    permutations = get_permutations(1024, num_hashes)

    buckets_allowed = get_buckets_allowed(num_hashes)

    bucket = buckets_allowed[0]

    fingerprints_hashes = fingerprints.map(
        lambda (idx, fp): (idx, fp, hash_fingerprints(fp, idx, bucket, permutations))).cache()

    hash_list = fingerprints_hashes.flatMap(lambda (idx, fp, hashes): hashes[1])

    print(hash_list.count())
    complete_hash_list = hash_list.reduceByKey(lambda a, b: a.union(b))
    bc_hash_list = spark.sparkContext.broadcast(complete_hash_list.collectAsMap())

    potential_neighbours = fingerprints_hashes.map(
        lambda (idx, fp, hash): (idx, get_neighbours(hash[0], bc_hash_list.value, idx)))

    neighbours = potential_neighbours.flatMap(
        lambda (idx, nbrs): (get_threshold_neighbours_flat2(bc_fingerprints.value, idx, nbrs, SIMILARITY_THRESHOLD)))

    complete_neighbours = neighbours.flatMap(lambda (a, b): [(a, b)] if a == b else [(a, b), (b, a)])
    return complete_neighbours

def get_neighbours_counts(neighbours):
    neighbours = neighbours.map(lambda (a,b):(a,b,len(b)))
    mol_neighbour_count = neighbours.map(lambda (a, b, c): (a, c))
    return mol_neighbour_count, neighbours


def output_results(sc, cluster_combs, compounds, evaluate):
    if evaluate:
        #cluster_combs = cluster_combs.sortBy(lambda (a,b): a)
        #cluster_combs1 = cluster_combs.map(lambda (a,b):(sort_list(list(b))))

        #cluster_combs.collectAsMap()
        #invert id and molecules
        compounds = compounds.map(lambda (x,y): (y,x))


        bc_compounds = sc.broadcast(compounds.collectAsMap())
        #from the clusters of indeces, get the actual molecule
        #ie. clusters of molecules

        cluster_mols = cluster_combs.map(lambda (cl_id,compound_ids): index_to_mol(compound_ids, compound_list=bc_compounds))

        #assign ids to clusters
        cluster_groups = cluster_mols.zipWithIndex()

        #flat the molecules in the format: Smiles string, Molecule Name, Cluster Ids

        cluster_groups = cluster_groups.flatMap(lambda (r, idx): [(Chem.MolToSmiles(cmp), cmp.GetProp("_Name"), idx) for cmp in r])

        #Change molecule data to string
        cluster_groups = cluster_groups.map(lambda (r, name, idx): r + ' ' + name + ' ' + str(idx))

        #output to file

        header = sc.parallelize(["smiles Name Cluster"])
        output_file = cluster_groups.union(header)

        #output_file.saveAsTextFile("../mols/resultsSpark/result")
        current_time_milli = int(round(time.time() * 1000))
        outputextension = str(current_time_milli)
        output_file.coalesce(1).saveAsTextFile(baseFile + "/output" + outputextension + "/result"+ outputextension)


if __name__ == '__main__':
    dataFile = baseFile + sys.argv[1] + '.smi'
    executor_num = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    print(executor_num , " , Executors")
    print("Data File ",dataFile)
    main()