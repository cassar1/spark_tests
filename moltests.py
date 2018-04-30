from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.mllib.linalg.distributed import RowMatrix
from pyspark.mllib.linalg import Vectors
from pyspark.sql import SparkSession, Row
from pyspark.sql import SQLContext

from helpers import *


server_location = read_server()
conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf = conf)

sc.addPyFile("helpers.py")

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3

lines = sc.textFile("mols/compounds14.smi")
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


# create combinations for all the fingerprints, filter is used to calculate only half of the similarity matrix
cartesian_fp = fingerprints1.cartesian(fingerprints1).filter(lambda (a,b): a[0] >= b[0])
# calculate similarities
cartesian_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))
# keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
similarities_fp = cartesian_fp.filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD)
# sort by fp1
similarities_fp = similarities_fp.sortBy(lambda (a, b, c): a)

# change the second fingerprint to a list of fingerprints
similarities_fp = similarities_fp.map(lambda (a, b, c): (a, set([b])))

# combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
neighbours = similarities_fp.reduceByKey(lambda a, b: a.union(b))\
    .sortBy(lambda (a,b): len(b),False)

#similarities_fp.collect()
#neighbours.collect()
neighbours.foreach(output)
print ("--------------------------------------------")
overlapping_clusters = neighbours.map(lambda (a,b):b)\
                                .zipWithIndex()\
                                .map(lambda (x,y): (y,x))

print ("--------------------------------------------")
#select first cluster as it would be lost in cartesian operation
first_cluster = overlapping_clusters.filter(lambda (a,b): a == 0)

# perform cartesian, with each cluster being grouped with clusters whose id is smaller (thus beign higher ranked)
cluster_combs = overlapping_clusters.cartesian(overlapping_clusters).filter(lambda (a,b): a[0] > b[0])
cluster_combs = cluster_combs.map(lambda (a,b): (a[0],a[1].difference(b[1])))

cluster_combs = cluster_combs.reduceByKey(lambda a, b: a.intersection(b))

# combine all clusters, remove empty clusters and order by length descending
cluster_combs = cluster_combs.union(first_cluster)\
                            .filter(lambda (a,b): len(b) > 0)\
                            .sortBy(lambda (a,b): len(b),False)

cluster_combs.foreach(output)