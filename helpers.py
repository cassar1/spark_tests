from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
from operator import itemgetter

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
        #print idx
        compound = compound_list.value[idx][1]
        compound.SetProp("_Name",str(compound_list.value[idx][3]))
        #print (compound_list.value[idx][3])
        compounds.append(compound)
    return compounds


#return list of tuples (molecule, count of neighbours)
#region count based
def convert_neighbours(neighbours, complete_list):
    neighbour_counts = []
    for nbr in neighbours:
        for row in complete_list:
            if row[0] == nbr:
                print ("Row ", row)
                neighbour_counts.append(row)
                break
    #neighbour_counts = sorted(neighbour_counts, key=itemgetter(1), reverse=True)
    neighbour_counts = sorted(neighbour_counts, key=lambda x: (-x[1], x[0]))
    #print ("Nbrs", neighbour_counts)
    return neighbour_counts

def convert_neighbours_dict(neighbours, complete_list):
    neighbour_counts = []
    for nbr in neighbours:
        #print ("Element, ",nbr," ", complete_list[nbr])
        neighbour_counts.append((nbr, complete_list[nbr]))
    neighbour_counts = sorted(neighbour_counts, key=lambda x: (-x[1], x[0]))
    return neighbour_counts


#a single neighbour is provided as a parameter, with complete list being a dictionary of mol_id and count
def convert_single_neighbour_dict(neighbour, complete_list):
    return (neighbour, complete_list[neighbour])
    #return set([(neighbour, complete_list[neighbour])])

def convert_owner_dict(owner, complete_list):
    return complete_list[owner]

def assign_cluster(this_id, this_count, neighbours, invalid_clusters, current_id):
    if this_id == current_id:
        return current_id
    #print neighbours
    for nbr in neighbours:
        if nbr[0] not in invalid_clusters:
            if nbr[1] > this_count or (nbr[1] == this_count and this_id > nbr[0]):
                return nbr[0]
            else:
                break
    return this_id

def remove_invalid_nbrs(neighbours, invalid_clusters, cluster_owner):
    valid_nbrs = []
    if not cluster_owner:
        for nbr in neighbours:
            if nbr[0] not in invalid_clusters:
                valid_nbrs.append(nbr)
    return valid_nbrs

def remove_invalid_nbrs_dict(neighbours, invalid_clusters, cluster_owner):
    valid_nbrs = []
    if not cluster_owner:
        for nbr in neighbours:
            if nbr[0] not in invalid_clusters:
                valid_nbrs.append(nbr)
    return valid_nbrs

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