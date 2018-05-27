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
    if(elements[0] != 'smiles' and elements[0] != 'SMILES'):
        mol = Chem.MolFromSmiles(elements[0])
        #mol.SetProp("_Name",str(elements[1]))
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,2, nBits=1024)
        #base64Fp = fingerprint.ToBitString()
    else:
        print elements[0]
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
    try:

        return DataStructs.FingerprintSimilarity(fp1,fp2, metric=DataStructs.TanimotoSimilarity)
    except:
        print "Error fp", fp1, " 2", fp2

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


#region Recursion

def is_cluster_invalid(mol_id, potential_parents, complete_list):
    invalid = cluster_invalid(potential_parents, complete_list)
    if invalid:
        return mol_id
    else:
        return -1

# dataline = invalid molecule, neighbours, molecule invalidating the cluster
# we need to check whether the molecule invalidating the cluster is valid itself
# True cluster invalid
# False Cluster Valid

def cluster_invalid(potential_parents, complete_list):
    # has_valid_parent = has_valid_parent(potential_parents, complete_list)
    for parent in potential_parents:
        # print "Searching potential Parent ", parent
        mol_parents = get_mol_parents(parent, complete_list)
        # print "Parents for molecule ", parent
        if mol_parents is None:
            # print "returning count ", count
            return True

        parent_invalid = cluster_invalid(mol_parents[1], complete_list)

        # if any parent is valid, then this cluster is invalid
        if parent_invalid is False:
            return True
    return False


def get_mol_parents(mol_id, complete_list):
    for row in complete_list:
        if row[0] == mol_id:
            # print "Found Grandparents", row
            return row
    return None

#endregion Recursion

#region count based
def convert_neighbours(neighbours, complete_list):
    neighbour_counts = []
    for nbr in neighbours:
        for row in complete_list:
            if row[0] == nbr:
                #print ("Row ", row[0], " NBR ", nbr)
                neighbour_counts.append(row)
                break
    return neighbour_counts

def assign_cluster(this_id, this_count, neighbours, invalid_clusters):
    for nbr in neighbours:
        if nbr[0] not in invalid_clusters:
            if nbr[1] > this_count or (nbr[1] == this_count and this_id > nbr[0]):
                return nbr[0]
    return this_id

#endregion count based

#FOR TEST
def remove_el(nbrs, item_to_remove):
    nbrs.remove(item_to_remove)
    nbrs.sort()
    return nbrs

def sort_list(nbrs):
    nbrs.sort()
    return nbrs

def persistDataToFile(data):
    file = open("../mols/resultsSpark/result.smi","w+")

    #file.write("%s" % data)
    for item in data:
        file.write("%s\n" % item)
    file.close()