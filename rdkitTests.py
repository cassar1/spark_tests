from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def createFingerprint(mols, counts= False):
    fingerprints = []
    if not counts:
        #for mol in mols:
            #fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(smiles, 2, 1024))
        fingerprints = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols]
        return fingerprints
    else:
        info = {}
        fingerprints = [AllChem.GetMorganFingerprint(m,2,bitInfo=info) for m in mols]
        return fingerprints

def read_fingerprints(input_file):
    #molecules = Chem.SDMolSupplier(input_file)
    molecules = Chem.SmilesMolSupplier(input_file, titleLine=True, delimiter='\t')
    mols = []
    count = -1
    with open(input_file) as f:
        for line in f:
            count = count + 1
            if count == 0:
                continue
            sml = line.split(' ')[0]
            mols.append(Chem.MolFromSmiles(sml))

    moleculeFingerprints = createFingerprint(mols)
    return moleculeFingerprints

#adapted from RDKIT BUTINA Algorithm
def get_neighbours_list(fps, threshold):
    #start_time = time.time()
    n_points = len(fps)
    neighbour_list = [None] * n_points
    for i in range(n_points):
        neighbour_list[i] = []

    for i in range(n_points):
        #print "chemical", i
        for j in range(i):
            if i != j:
                distance = DataStructs.FingerprintSimilarity(fps[i], fps[j], metric= DataStructs.TanimotoSimilarity)

                if (distance >= threshold):
                    print "comparing chemicals ", i ," to chemical", j , "distance ", distance
                    neighbour_list[i].append(j)
                    neighbour_list[j].append(i)
    # sort list of neighbours by num neighbours
    tLists = [(len(y), x) for x, y in enumerate(neighbour_list)]
    tLists.sort(reverse=True)
    #print "time taken to calculate ", n_points ," : ", time.time() - start_time
    return tLists, neighbour_list


fingerprints = read_fingerprints("mols/compounds5.smi")
tuple_list, neighbours_list = get_neighbours_list(fingerprints, 0.1)

for tuple in tuple_list:
    print tuple