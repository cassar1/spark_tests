from __future__ import print_function
from pyspark import SparkConf, SparkContext, SQLContext
from pyspark.sql import SparkSession, Row
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, DoubleType, BinaryType
from pyspark.sql.functions import monotonically_increasing_id
from helpers import *
import pyspark.sql.functions as f
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem

EVALUATION = False

server_location = read_server()
conf = SparkConf().setMaster("local").setAppName("MoleculesTests")
#conf = SparkConf().setMaster(server_location).setAppName("MoleculesTests")

ss = SparkSession.builder.config(conf = conf).getOrCreate()
#sqlc = SQLContext()
ss.sparkContext.addPyFile("helpers.py")

def mapper(line):
    fields = line.split(' ')
    return Row(name=str(fields[0]))

def calculate_tanimoto(smiles1,smiles2):
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
        similarity = DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity)

        return similarity
    except Exception as e:
        print (str(e))
        print ("Error Smiles1", smiles1, " 2", smiles2)


CREATE_VECTORS = False
SIMILARITY_THRESHOLD = 0.3
dataFile = '../mols/compounds18.smi'
lines = ss.sparkContext.textFile(dataFile)

smiles = lines.map(mapper)
print ("Here1")
schemaSmiles = ss.createDataFrame(smiles).cache()
schemaSmiles.createOrReplaceTempView("smiles")
valid_smiles = ss.sql("SELECT * FROM smiles WHERE name != 'smiles'")
valid_smiles_id = valid_smiles.select("*").withColumn("id", monotonically_increasing_id())
print ("Here2")
#print(valid_smiles.show(n=10))

combinations = valid_smiles_id.alias("source").join(valid_smiles_id.alias("target") )\
    .where("source.Id <= target.Id")\
    .select(f.col("source.Id").alias("source_id"), f.col("source.Name").alias("source_smile"), f.col("target.Id").alias("target_id"),f.col("target.Name").alias("target_smile"))
print(combinations.show(n=10))
print ("Here3")

combinations_rdd = combinations.rdd.map(tuple)
#combinations_rdd.foreach(output)
similarities_fp = combinations_rdd.map(lambda (source_id, source_smiles,target_id,target_smiles): (source_id, target_id, calculate_tanimoto(source_smiles, target_smiles)))\
                              .filter(lambda (a,b,c): c >= SIMILARITY_THRESHOLD).cache()

#udf_tanimoto = f.udf(calculate_tanimoto, DoubleType())
#combinations_sim = combinations.withColumn("tanimoto",f.lit(calculate_tanimoto(combinations.source_smile,combinations.target_smile)))
#print(combinations_sim.show(n=10))
#combinations_sim = combinations_sim.filter(combinations_sim.tanimoto > SIMILARITY_THRESHOLD)
#print(combinations_sim.show(n=10))
print ("Here4")

schema = StructType([StructField("source",IntegerType(), False),StructField("target",IntegerType(), False),StructField("tanimoto",StringType(), False)])

combinations_sim = ss.createDataFrame(similarities_fp,schema=schema).cache()
print(combinations_sim.show(n=10))
#combinations_sim = combinations_sim.groupby(combinations_sim.source_id).agg(f.collect_set("target_id"))

#print(combinations_sim.show(n=10))

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

    header = ss.sparkContext.parallelize(["smiles Name Cluster"])
    output_file = cluster_groups.union(header)
    #output_file.foreach(output)
    output_file.saveAsTextFile("../mols/resultsSpark/resultNbrs")