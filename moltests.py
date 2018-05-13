from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.mllib.linalg.distributed import RowMatrix
from pyspark.mllib.linalg import Vectors
from pyspark.sql import SparkSession, Row
from pyspark.sql import SQLContext
from sklearn.utils.extmath import cartesian

from helpers import *


server_location = read_server()
#conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf = conf)

sc.addPyFile("helpers.py")

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3

lines = sc.textFile("mols/compounds5.smi")
smiles = lines.map(convertToBitVectorFP)

print ("Read All data, length: ", smiles.count())

# remove file headers
smiles = smiles.filter(lambda l: l[0] != "smiles")

print ("Filtered data, length: ", smiles.count())

# select only the fingerprint part, since convertToBitVectorFP returns an object
fingerprints = smiles.map(lambda x:x[2])

# Add index to fingerprints
fingerprints1 = fingerprints.zipWithIndex()

# convert index and value position
fingerprints1 = fingerprints1.map(lambda (x,y): (y,x))

bc_fingerprints = sc.broadcast(fingerprints1.collectAsMap())


# create combinations for all the fingerprints, filter is used to calculate only half of the similarity matrix
#cartesian_fp = fingerprints1.cartesian(fingerprints1).filter(lambda (a,b): a[0] >= b[0])
cartesian_fp = fingerprints1.flatMap(lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value if x <= y[0]] )
#print ("cartesian amount ", cartesian_fp.count())

# calculate similarities
# cartesian_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))

# keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))\
                              .filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD)\
                              .map(lambda (a, b, c): (a, set([b])))
# sort by fp1
#similarities_fp = similarities_fp.sortBy(lambda (a, b, c): a)

# change the second fingerprint to a list of fingerprints
#similarities_fp = similarities_fp.map(lambda (a, b, c): (a, set([b])))

# combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
neighbours = similarities_fp.reduceByKey(lambda a, b: a.union(b))\
    .filter(lambda (a,b): len(b) > 0)#\
    #.sortBy(lambda (a,b): (len(b), a), False)

#similarities_fp.collect()
#neighbours.collect()
#neighbours.foreach(output)
#overlapping_clusters = neighbours.map(lambda (a,b):b)\
#                                .zipWithIndex()\
#                                .map(lambda (x,y): (y,x))

#select first cluster as it would be lost in cartesian operation
#first_cluster = neighbours.filter(lambda (a,b): a == 0)

#overlapping_clusters.foreach(output)

#print ("--------------------------------------------")
# perform cartesian, with each cluster being grouped with clusters whose id is smaller (thus being higher ranked)
#.filter(lambda (a,b): a[0] > b[0])
#cluster_combs = neighbours.cartesian(neighbours).filter(lambda (a,b): a[0] > b[0])
cluster_combs = neighbours.cartesian(neighbours).filter(lambda (a,b): (len(a[1]) == len(b[1]) and a[0] >= b[0]) or (len(a[1]) < len(b[1])))
#sorted_cluster_combs = cluster_combs.sortBy(lambda (a,b): a[0], False)
#cluster_combs.foreach(output)

cluster_combs = cluster_combs.map(lambda (a,b):(a[0], a[1].difference(b[1]))  if a[0] != b[0] else (a[0], a[1]) )
#print ("----------------------------------------------")
#cluster_combs.foreach(output)

cluster_combs = cluster_combs.reduceByKey(lambda a, b: a.intersection(b))\
    .filter(lambda (a,b): len(b) > 0)

# combine all clusters, remove empty clusters and order by length descending
#cluster_combs = cluster_combs.union(first_cluster).filter(lambda (a,b): len(b) > 0).sortBy(lambda (a,b): (len(b),a),False)
                            #.filter(lambda (a,b): len(b) > 0)\

cluster_combs.foreach(output)