from __future__ import print_function
from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession
import time
from numpy import array
from pyspark.mllib.clustering import BisectingKMeans, BisectingKMeansModel
from sklearn.externals import joblib

EVALUATION = True
#dataFile = '../mols/compounds5.smi'
dataFile = '../fyp_comparisons/dataset/SmilesMerged/ABL1merged.smi'
#dataFile = '../mols/percentages/3910ABL1.smi'
baseFile = '../mols/resultsSerial/bisecting'
modelFile = '../mols/bisectingModel'
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
    compounds = load_data(sc, dataFile)
    compounds.partitionBy(executor_num)
    fingerprints = select_fingerprints(compounds).cache()

    fp_only = fingerprints.map(lambda (id, smi, fp, name): fp)

    for x in [1500,2000]:
        start_time = time.time()
        model = BisectingKMeans.train(fp_only, k = x)
        #print(model.clusterCenters)
        #print("Clusters " ,len(model.clusterCenters))


        cost = model.computeCost(fp_only)
        #model.save(sc, baseFile + '/btreemodel')
        print ("Bisecting "+ str(cost))

        #model.clusterCenters.foreach(lambda ctr : print("Cluster Center"))

        all_fps = fingerprints.collect()
        cluster_assignment = []
        end_time1 = time.time()
        print("Clustering Time taken ", x, end_time1 - start_time)
        for fp in all_fps:
            cluster_assignment.append('{} {} {}'.format(fp[1], fp[3], model.predict(fp[2])))
            #print ( "FP ", fp[0], " SMI: ", fp[1], " ", model.predict(fp[2]))

        end_time = time.time()
        print ("Total Time taken " , x , end_time - start_time)
        if EVALUATION:
            header = sc.parallelize(["smiles Name Cluster"])
            clusters = sc.parallelize(cluster_assignment)

            output_file = header.union(clusters)
            #output_file.foreach(output)
            # output_file.saveAsTextFile("../mols/resultsSpark/result")
            current_time_milli = int(round(time.time() * 1000))
            outputextension = str(current_time_milli)
            output_file.coalesce(1).saveAsTextFile(baseFile + "/output" + str(x) + "/result" + outputextension)

    sc.stop()

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
    #compounds.foreach(output)
    fingerprints = compounds.map(lambda (x, idx): (idx, x[0], x[2], x[3]))

    #fingerprints.foreach(output)
    fingerprints_arr = fingerprints.map(lambda (idx, smi, fp, name): (idx, smi, convert_to_arr(fp), name))

    return fingerprints_arr


def convert_to_arr(fp):
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

if __name__ == '__main__':
    main()