from __future__ import division

import sys
import numpy as np

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

########################################################################
#                               LSH
########################################################################




def get_permutations(len_permutations=1024, num_permutations=100):
    return map(lambda _: np.random.permutation(len_permutations), range(num_permutations))


def get_buckets_allowed(min_hash):
    num_bits = (1024 - 1).bit_length()
    max_hashes_bucket = int(sys.maxint.bit_length() / num_bits)
    min_buckets = int(min_hash / max_hashes_bucket)

    bucket_sizes_allowed = []
    for i in range(min_buckets, 1000):
        if i is not 0 and not (min_hash % i):
            print "allowed buckets ", i
            hash_per_bucket = int(min_hash / i)

            if num_bits * hash_per_bucket > sys.maxint.bit_length():
                print 'numbers are too large to produce valid buckets'

            print "Hashes per bucket:", hash_per_bucket
            bucket_sizes_allowed.append(i)
    return bucket_sizes_allowed


def get_min_hash(fingerprint, permutations):
    # qfp_bits = [int(n) for n in AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)]
    qfp_bits = [int(n) for n in fingerprint]
    min_hash = []
    for perm in permutations:
        for idx, i in enumerate(perm):
            if qfp_bits[i]:
                min_hash.append(idx)
                break
    return min_hash


def hash_to_buckets(min_hash, num_buckets=10, nBits=1024):

    if len(min_hash) % num_buckets:
        raise Exception('number of buckets must be divisiable by the hash length')
    buckets = []
    hash_per_bucket = int(len(min_hash) / num_buckets)

    num_bits = (nBits - 1).bit_length()

    if num_bits * hash_per_bucket > sys.maxint.bit_length():
        raise Exception('numbers are too large to produce valid buckets')
    for b in range(num_buckets):

        buckets.append(reduce(lambda x, y: (x << num_bits) + y, min_hash[b:(b + hash_per_bucket)]))
    return buckets


def get_threshold_neighbours(bc_fingerprints, target_fp_idx, ids_to_compare, threshold):
    valid_neighbours = []

    for comparing_fp in ids_to_compare:
        if comparing_fp != target_fp_idx:
            target_fp = bc_fingerprints[target_fp_idx]
            nbr_fp = bc_fingerprints[comparing_fp]
            similarity = DataStructs.FingerprintSimilarity(target_fp, nbr_fp, metric=DataStructs.TanimotoSimilarity)
            if similarity > threshold:
                valid_neighbours.append(comparing_fp)
    return valid_neighbours


def hash_fingerprints(fingerprint, fp_index, num_buckets, permutations):
    hash_list = []
    # a dictionary with id of fingerprint and the list of hashes that it makes part of
    index_hashes = {}

    min_hash = get_min_hash(fingerprint, permutations)
    buckets = hash_to_buckets(min_hash, num_buckets)
    index_hashes[fp_index] = set()#[]
    for hash_bucket in buckets:
        #if hash_bucket not in hash_list:
            #hash_list[hash_bucket] = set()#[]

        hash_list.append((hash_bucket,set([fp_index])))
        index_hashes[fp_index].add(hash_bucket)

    return [index_hashes, hash_list]


def get_neighbours(index_hashes, hash_list, target_fp_idx):
    potential_neighbours = set()

    fp_hashes = index_hashes[target_fp_idx]

    for hash_id in fp_hashes:
        #for key in index_hashes:
        #    if hash_id in index_hashes[key]:
        #        potential_neighbours.add(key)
        potential_neighbours_in_hash = hash_list[hash_id]
        potential_neighbours.update(potential_neighbours_in_hash)

    #sorted_neighbours = sorted(potential_neighbours)
    return potential_neighbours

def get_threshold_neighbours_flat(bc_fingerprints, target_fp_idx, ids_to_compare, threshold):
    valid_neighbours = []

    for comparing_fp in ids_to_compare:
        if comparing_fp != target_fp_idx:
            target_fp = bc_fingerprints[target_fp_idx]
            nbr_fp = bc_fingerprints[comparing_fp]
            similarity = DataStructs.FingerprintSimilarity(target_fp, nbr_fp, metric=DataStructs.TanimotoSimilarity)
            if similarity > threshold:
                # print "fp:", comparing_fp, " sim: ", similarity
                valid_neighbours.append((target_fp_idx, comparing_fp))

    return valid_neighbours

def get_threshold_neighbours_flat2(bc_fingerprints, target_fp_idx, ids_to_compare, threshold):
    valid_neighbours = []
    nbr_fps = []
    for comparing_fp in ids_to_compare:
        nbr_fps.append(bc_fingerprints[comparing_fp])
    target_fp = bc_fingerprints[target_fp_idx]

    similarity_list = DataStructs.BulkTanimotoSimilarity(target_fp, nbr_fps)

    for idx, val in enumerate(ids_to_compare):
        if similarity_list[idx] > threshold:
            # print "fp:", comparing_fp, " sim: ", similarity
            valid_neighbours.append((target_fp_idx, val))

    return valid_neighbours


def get_neighbours_flat(target_fp_idx, ids_to_compare):
    valid_neighbours = []

    for comparing_fp in ids_to_compare:
        valid_neighbours.append((target_fp_idx, comparing_fp))

    return valid_neighbours