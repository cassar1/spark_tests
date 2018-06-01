from __future__ import print_function
from pyspark import SparkConf, SparkContext

from helpers import *


server_location = read_server()
#conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf = conf)

sc.addPyFile("helpers.py")

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds190.smi'
lines = sc.textFile(dataFile)

smiles = lines.map(convertToBitVectorFP)

print ("Read All data, length: ", smiles.count())

# remove file headers
smiles = smiles.filter(lambda l: l[0] != "smiles" and l[0] != "SMILES")

print ("Filtered data, length: ", smiles.count())

#assign index to each compound
compounds = smiles.zipWithIndex()
#example entry in compounds: ((u'C[N+](C)(C)[C@@H]1[C@@H](O)O[C@H]2[C@@H](O)CO[C@H]21', <rdkit.Chem.rdchem.Mol object at 0x7f7dc0cd9050>, <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x7f7dc0cd90c0>, u'ZINC000000039940'), 0)

fingerprints1 = compounds.map(lambda (x, idx): (idx, x[2]))

bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())

# create combinations for all the fingerprints, filter is used to calculate only half of the similarity matrix
cartesian_fp = fingerprints1.flatMap(lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x > y[0]])
#cartesian_fp.foreach(output)

# keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))\
                              .filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD)

inverted_similarities = similarities_fp.map(lambda (a,b,c): (b,a,c))

similarities_fp = similarities_fp.union(inverted_similarities).map(lambda (a, b, c): (a, set([b])))

#similarities_fp.foreach(output)

# sort by fp1

# combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
neighbours = similarities_fp.reduceByKey(lambda a, b: a.union(b))\
    .filter(lambda (a,b): len(b) > 0)\
    .map(lambda (a,b): (a,b,len(b)))

neighbours.foreach(output)