from __future__ import print_function
from pyspark import SparkConf, SparkContext
from helpers import *

EVALUATION = False
FLATMAPMETHOD = True
server_location = read_server()
conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf = conf)

sc.addPyFile("helpers.py")

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds82.smi'
lines = sc.textFile(dataFile)

smiles = lines.map(convertToBitVectorFP)

print ("Read All data, length: ", smiles.count())

# remove file headers
smiles = smiles.filter(lambda l: l[0] != "smiles" and l[0] != "SMILES")

print ("Filtered data, length: ", smiles.count())

#assign index to each compound
compounds = smiles.zipWithIndex()
#example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)

fingerprints0 = compounds.map(lambda (x, idx): (idx, x[2]))

fingerprints1 = fingerprints0.partitionBy(6).cache()

bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())

# create combinations for all the fingerprints, if condition is used to calculate only upper triangle of the similarity matrix
cartesian_fp = fingerprints1.flatMap(lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x >= y[0]])
#cartesian_fp.foreach(output)

# keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))\
                              .filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD)



if FLATMAPMETHOD:
    complete_similarities_fp = similarities_fp.flatMap(lambda (a,b,c): [(a,set([b])),(b,set([a]))])
else:
# a new rdd is created to get the lower trangle of the similarity matrix
    inverted_similarities = similarities_fp.map(lambda (a,b,c): (b,a,c))

# the complete catrix is obtained through union
    complete_similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))
# the complete catrix is obtained through union
#similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))


#similarities_fp.foreach(output)

# sort by fp1

# combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
neighbours = complete_similarities_fp.reduceByKey(lambda a, b: a.union(b))#.map(lambda (a,b): (a,b,len(b)))

#.filter(lambda (a,b): len(b) > 0) \
    #\
    #.map(lambda (a,b): (a,b,len(b)))
neighbours.collectAsMap()
#neighbours.foreach(output)
if EVALUATION:
    neighbours.collectAsMap()
    neighbours = neighbours.map(lambda (id, nbrs): (id, sorted(nbrs)))
    neighbours = neighbours.sortBy(lambda (id, nbrs): id)
    print("TOTAL NBRS", neighbours.count())
    #invert id and molecules
    compounds = compounds.map(lambda (x,y): (y,x))
    bc_compounds = sc.broadcast(compounds.collectAsMap())


    #Change molecule data to string

    cluster_groups = neighbours.map(lambda (id, nbrs): str(id) + ' , ' + str(list(nbrs)))
    print("TOTAL NBRS", cluster_groups.count())
    #output to file

    header = sc.parallelize(["smiles Name Cluster"])
    output_file = cluster_groups.union(header)
    #output_file.foreach(output)
    output_file.saveAsTextFile("../mols/resultsSpark/resultNbrs")