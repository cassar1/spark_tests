from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.mllib.linalg.distributed import RowMatrix
from pyspark.mllib.linalg import Vectors
from pyspark.sql import SparkSession, Row
from pyspark.sql import SQLContext
from sklearn.utils.extmath import cartesian

from helpers import *

EVALUATION = True

server_location = read_server()
conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")
sc = SparkContext(conf = conf)

sc.addPyFile("helpers.py")

CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
#dataFile = '../mols/compounds5.smi'
dataFile = '../mols/merged/Reninmerged.smi'
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
cartesian_fp = fingerprints1.flatMap(lambda y: [(y, (x, bc_fingerprints.value[x])) for x in bc_fingerprints.value])
#if x != y[0]
# if x < y[0]
#cartesian_fp.foreach(output)

# keep only those that have similarity > threshold a = fp1, b = fp2, c = similarity
similarities_fp = cartesian_fp.map(lambda (a, b): (a[0], b[0], calculate_tanimoto(a[1], b[1])))\
                              .filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD)\
                              .map(lambda (a, b, c): (a, set([b])))
# sort by fp1

# combine the list of neighbours into key value, with key being the fingerprint, value being list of neighbours
neighbours = similarities_fp.reduceByKey(lambda a, b: a.union(b))\
    .filter(lambda (a,b): len(b) > 0)


#region Serial Part
#neighbours_sorted = neighbours.map(lambda (a,b): (a,sorted(b)))
#neighbours_sorted.sortBy(lambda (a,b): len(b), False)#.foreach(output)

#local_neighbours = neighbours_sorted.collect()

#res = []
#seen = [0] * len(local_neighbours)
#local_neighbours.sort(key=lambda x: len(x[1]), reverse=True)

'''for neighbours in local_neighbours:
    print (neighbours[0], neighbours[1])
    if seen[neighbours[0]]:
        continue
    tRes = []
    for nbr in neighbours[1]:
        #print (nbr)
        if not seen[nbr]:
            tRes.append(nbr)
            seen[nbr] = 1
    res.append(tRes)

cluster_combs = sc.parallelize(res)
cluster_combs = cluster_combs.zipWithIndex()
cluster_combs = cluster_combs.map(lambda (list,id):(id,list))'''
#endregion End Serial Part

cluster_combs = neighbours.cartesian(neighbours).filter(lambda (a,b): (len(a[1]) == len(b[1]) and a[0] <= b[0]) or (len(a[1]) < len(b[1])))

cluster_combs = cluster_combs.map(lambda (a, b): (a[0], a[1].difference(b[1]), b[0]) if a[0] != b[0] else (a[0], a[1], a[0]))
#cluster_combs.foreach(output)

potential_invalid = cluster_combs.filter(lambda (a, b, c): a not in b).map(lambda (a,b,c): (a,set([c])))#.map(lambda (a,b,c): a)
#potential_invalid.foreach(output)
#and (len(b) > 0)
print ("Potentially invalid clusters")
clusters_to_invalidate = potential_invalid.reduceByKey(lambda a, b: a.union(b))

#clusters_to_invalidate.foreach(output)


bc_clusters_to_invalidate = sc.broadcast(clusters_to_invalidate.collect())

actual_invalid_clusters = clusters_to_invalidate.map(lambda y: is_cluster_invalid(y[0], y[1], bc_clusters_to_invalidate.value))
actual_invalid_clusters = actual_invalid_clusters.filter(lambda a: a != -1)

print ("Actual Invalid Clusters")
#actual_invalid_clusters.foreach(output)
bc_clusters_to_invalidate = sc.broadcast(actual_invalid_clusters.collect())

cluster_combs_valid = cluster_combs.filter(lambda (a, b, c) : (a not in bc_clusters_to_invalidate.value) and (c not in bc_clusters_to_invalidate.value))\
    .map(lambda (a,b,c) : (a,b))

cluster_combs = cluster_combs_valid.reduceByKey(lambda a, b: a.intersection(b))\
    .filter(lambda (a,b): len(b) > 0)
print ("-----")
#cluster_combs.sortBy(lambda (a,b): len(b), True).foreach(output)


if EVALUATION:
    cluster_combs.collectAsMap()
    #invert id and molecules
    compounds = compounds.map(lambda (x,y): (y,x))
    bc_compounds = sc.broadcast(compounds.collectAsMap())

    #from the clusters of indeces, get the actual molecule
    #ie. clusters of molecules
    cluster_mols = cluster_combs.map(lambda (cl_id,compound_ids): index_to_mol(compound_ids, compound_list=bc_compounds))
    #assign ids to clusters
    cluster_groups = cluster_mols.zipWithIndex()
    #cluster_groups.foreach(output)
    #flat the molecules in the format: Smiles string, Molecule Name, Cluster Ids
    cluster_groups = cluster_groups.flatMap(lambda (r, idx): [(Chem.MolToSmiles(cmp), cmp.GetProp("_Name"), idx) for cmp in r])

    #Change molecule data to string
    cluster_groups = cluster_groups.map(lambda (r, name, idx): r + ' ' + name + ' ' + str(idx))

    #output to file

    header = sc.parallelize(["smiles Name Cluster"])
    output_file = cluster_groups.union(header)
    #output_file.foreach(output)
    output_file.saveAsTextFile("../mols/resultsSpark/result")
