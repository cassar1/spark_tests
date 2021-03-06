from __future__ import print_function
from pyspark import SparkConf, SparkContext
from helpers import *

EVALUATION = False
FLATMAPMETHOD = True

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds17.smi'

def main():
    sc = start_spark()

    compounds = load_data(sc, dataFile)
    #example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)

    fingerprints = select_fingerprints(compounds)

    fingerprints1 = fingerprints.partitionBy(6).cache()

    bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())

    neighbours = calculate_neighbours(fingerprints1, bc_fingerprints)

    neighbours.collect()
    neighbours.foreach(output)
    output_results(sc, neighbours, EVALUATION)

def start_spark():
    server_location = read_server()
    conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
    # conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
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


def output_results(sc, neighbours, output):
    if output:
        neighbours.collectAsMap()
        neighbours = neighbours.map(lambda (id, nbrs): (id, sorted(nbrs)))
        neighbours = neighbours.sortBy(lambda (id, nbrs): id)
        print("TOTAL NBRS", neighbours.count())

        #Change molecule data to string
        cluster_groups = neighbours.map(lambda (id, nbrs): str(id) + ' , ' + str(list(nbrs)))
        print("TOTAL NBRS", cluster_groups.count())
        #output to file

        header = sc.parallelize(["smiles Name Cluster"])
        output_file = cluster_groups.union(header)
        #output_file.foreach(output)
        output_file.saveAsTextFile("../mols/resultsSpark/resultNbrs")


if __name__ == '__main__':
    main()