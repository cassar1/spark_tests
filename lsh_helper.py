from __future__ import division

import sys
import numpy as np
from rdkit import DataStructs


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

'''
def hash_fingerprints(fingerprint, fp_index, num_buckets, permutations):
    # a dictionary with id of fingerprint and the list of hashes that it makes part of
    index_hashes = {}

    min_hash = get_min_hash(fingerprint, permutations)
    buckets = hash_to_buckets(min_hash, num_buckets)
    index_hashes[fp_index] = set()#[]
    for hash_bucket in buckets:
        index_hashes[fp_index].add(hash_bucket)

    return index_hashes


def get_neighbours(index_hashes, target_fp_idx):
    potential_neighbours = set()

    fp_hashes = index_hashes[target_fp_idx]

    for hash_id in fp_hashes:
        for key in index_hashes:
            if key > target_fp_idx:
                if hash_id in index_hashes[key]:
                    potential_neighbours.add(key)

    #sorted_neighbours = sorted(potential_neighbours)
    return potential_neighbours
'''

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
    # print "Target fp:", target_fp_idx
    # print "Searching neighbours from ", len(ids_to_compare), " molecules"
    for comparing_fp in ids_to_compare:
        if comparing_fp != target_fp_idx:
            target_fp = bc_fingerprints[target_fp_idx]
            nbr_fp = bc_fingerprints[comparing_fp]
            similarity = DataStructs.FingerprintSimilarity(target_fp, nbr_fp, metric=DataStructs.TanimotoSimilarity)
            if similarity > threshold:
                # print "fp:", comparing_fp, " sim: ", similarity
                valid_neighbours.append((target_fp_idx, comparing_fp))

    return valid_neighbours