from pyspark.ml.feature import MinHashLSH
from pyspark.ml.linalg import Vectors
import pyspark.sql.functions as f
from pyspark import SparkConf, SparkContext, SQLContext
from pyspark.sql import SparkSession, Row
from helpers import *


server_location = read_server()
conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")

spark = SparkSession.builder.config(conf = conf).getOrCreate()
#sqlc = SQLContext()
spark.sparkContext.addPyFile("helpers.py")


DEBUG = False
CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds82.smi'
#dataFile = "wasb://fypblobcontainer@jurgenfypstorage.blob.core.windows.net/compounds18.smi"
#spark.conf.set("fs.azure.account.key.jurgenfypstorage.blob.core.windows.net",
#               "MzAJwHDAFmqYQmEFU7E7MmTmzLIvQrOQhLfHoL5T+VANMhKEcZRLe05PZp0EuCdceWIwVhiltjE5j6vzzyq9FQ==")


lines = spark.sparkContext.textFile(dataFile)

smiles = lines.map(convertToBitVectorFP)

print ("Read All data, length: ", smiles.count())

# remove file headers
smiles = smiles.filter(lambda l: l[0] != "smiles" and l[0] != "SMILES")

print ("Filtered data, length: ", smiles.count())

#assign index to each compound
compounds = smiles.zipWithIndex()

fingerprints = compounds.map(lambda (x, idx): (idx,  Vectors.dense([int(n) for n in x[2]]))).cache()

schemaFps = spark.createDataFrame(fingerprints, ["id", "features"]).cache()
#print(schemaFps.show(n=10))
fingerprints.unpersist()

mh = MinHashLSH(inputCol="features", outputCol="hashes", numHashTables=5)
model = mh.fit(schemaFps)

# Feature Transformation
print("The hashed dataset where hashed values are stored in the column 'hashes':")
#model_transform = model.transform(schemaFps)


# model_transform.show()
# Compute the locality sensitive hashes for the input rows, then perform approximate
# similarity join.
# We could avoid computing hashes by passing in the already-transformed dataset, e.g.
# `model.approxSimilarityJoin(transformedA, transformedB, 0.6)`
print("Approximately joining dfA and dfB on distance smaller than 0.3:")
# print(model_transform)
mol_similarity = model.approxSimilarityJoin(schemaFps, schemaFps, 1 - SIMILARITY_THRESHOLD, distCol="JaccardDistance") \
    .select(f.col("datasetA.id").alias("source_id"),
            f.col("datasetB.id").alias("target_id"),
            f.col("JaccardDistance"))

#TO READD
if DEBUG:
    combinations_sim = mol_similarity.groupby(mol_similarity.source_id).agg(f.sort_array(f.collect_set("target_id")).alias("nbrs"))
else:
    combinations_sim = mol_similarity.groupby(mol_similarity.source_id).agg(f.collect_set("target_id").alias("nbrs"))

#combinations_sim = combinations_sim.withColumn("nbr_count", f.size("nbrs"))

if DEBUG:
    combinations_sim.orderBy("source_id").show(20, False)
combinations_sim.show(10)
