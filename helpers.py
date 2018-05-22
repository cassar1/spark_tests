from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import SmilesWriter

def convertToFP(line):
    elements = line.split()
    mol = ""
    fingerprint = ""
    if(elements[0] != 'smiles'):
        mol = Chem.MolFromSmiles(elements[0])
        fingerprint = AllChem.GetMorganFingerprint(mol,2)
    return (elements[0], mol, fingerprint)


# Convert to bit fingerprints
def convertToBitVectorFP(line):
    elements = line.split()
    mol = ""
    fingerprint = ""
    if(elements[0] != 'smiles'):
        mol = Chem.MolFromSmiles(elements[0])
        #mol.SetProp("_Name",str(elements[1]))
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=1024)
        #base64Fp = fingerprint.ToBitString()
    return (elements[0], mol, fingerprint,elements[1])


def read_server():
    with open("serverPath.txt") as f:
        content = f.readline()
        return content


# Custom print method
def output(x):
    print (x)

def output_cluster(x):
    print x[0], ','.join(str(s) for s in x[1])


def calculate_tanimoto(fp1,fp2):
    return DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity)

def toCSVLine(data):
    return ','.join(str(d) for d in data)

def index_to_mol(indeces, compound_list):
    compounds = []
    for idx in indeces:
        #get mol item from list
        compound = compound_list.value[idx][1]
        compound.SetProp("_Name",str(compound_list.value[idx][3]))
        #print (compound_list.value[idx][3])
        compounds.append(compound)
    return compounds

'''def output_cluster_results(cluster_id, cluster):
    writer = SmilesWriter('../mols/resultsSerial/resultsclsmi.smi')
    writer.SetProps(['Cluster'])

    #cluster_id = 0
    for mol in cluster:
        mol.SetProp('Cluster', str(cluster_id))
        writer.write(mol)
    #cluster_id += 1
    writer.close()'''

def convert_to_smiles(cluster):
    print "x"

def persistDataToFile(data):
    file = open("../mols/resultsSpark/result.smi","w+")

    #file.write("%s" % data)
    for item in data:
        file.write("%s\n" % item)
    file.close()