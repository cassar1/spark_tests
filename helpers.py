from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem


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
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=1024)
        #base64Fp = fingerprint.ToBitString()
    return (elements[0], mol, fingerprint)


# Custom print method
def output(x):
    print (x)

def output_cluster(x):
    print x[0], ','.join(str(s) for s in x[1])


def calculate_tanimoto(fp1,fp2):
    return DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity)

def toCSVLine(data):
    return ','.join(str(d) for d in data)

def persistDataToFile(data):
    file = open("mols/result.csv","w+")

    #file.write("%s" % data)
    for item in data:
        file.write("%s\n" % item)
    file.close()