from .utils import Atom, Residue, ActiveSite
import numpy as np
import math
from math import factorial as fac

aa_hydro = {"ALA": 1.8, "ARG": -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5, 'GLU': -3.5, 'GLN': -3.5, 'GLY': -.4, 'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -.8, 'THR': -.7, 'TRP': -.9, 'TYR': -1.3, 'VAL': 4.2}
aa_mw = {"ALA": 89.094, "ARG": 174.203, 'ASN': 132.119, 'ASP': 133.104, 'CYS': 121.54, 'GLU': 147.131, 'GLN': 146.146, 'GLY': 75.067, 'HIS': 155.156, 'ILE': 131.175, 'LEU': 131.175, 'LYS': 146.189, 'MET': 149.208, 'PHE': 165.192, 'PRO': 115.132, 'SER': 105.093, 'THR': 119.119, 'TRP': 204.228, 'TYR': 181.191, 'VAL': 117.148}
aa_charge = {"ALA": 0, "ARG": 1, 'ASN': 0, 'ASP': -1, 'CYS': 0, 'GLU': -1, 'GLN': 0, 'GLY': 0, 'HIS':.9, 'ILE': 0, 'LEU': 0, 'LYS': 1, 'MET': -1, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0}

def location(site):
    """
    determines the "location" of a site
    the location is a point in 3-space corresponding to:
     total molecular weight
     total hydropathy index
     total charge
    """
    global aa_hydro
    global aa_mw
    global aa_charge

    site_hydro = 0
    site_mw = 0
    site_charge = 0

    for item in site.residues:
        site_hydro += aa_hydro[item.type]
        site_mw += aa_mw[item.type]
        site_charge += aa_charge[item.type]
    return (site_hydro, site_mw, site_charge)


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    #DEAL WITH DUPLICATES UGH

    # Euclidean distance
    
    return math.sqrt(sum([(x - y)**2 for x, y in zip(site_a, site_b)]))

def avg_to_cent(cluster, centroid):
    """
    compute average distance to centroid for given cluster
    """

    dist = []
    for site in cluster:
        dist.append(compute_similarity(location(site), centroid))
    return np.mean(dist)


def dbi(clusters):
    """
    compute Davies Bouldin index for a cluster
    """
    n = len(clusters)
    db_index = 0
    for a in clusters:
        a_l = kmeans_centroid(a)
        tmp = []  # store values to calculate min or max
        for b in clusters:
            b_l = kmeans_centroid(b)
            if b != a:
                numerator = avg_to_cent(a, a_l) + avg_to_cent(b, b_l)

                if b_l != a_l:
                    denom = compute_similarity(a_l, b_l)
                else:
                    return 0

                tmp.append(numerator/denom)
        db_index += max(tmp)
    return db_index / n

def kmeans_centroid(data):
    """
    input set of sites, return centroid as tuple
    """
    tuples = []
    for site in data:
        tuples.append(location(site))
    return tuple(np.mean(tuples, 0))

def kmeans_stop(max_iter, iterations, prev_centroid, centroid):
    """
    stop conditions for kmeans algorithm
    """
    if iterations > max_iter:
        return True
    if prev_centroid == centroid:
        print("centroids the same: " + str(iterations))
        return True


def kmeans(sites, k, max_iter):

    #randomly pick k points to serve as the first clusters
    clusters = [set() for x in range(k)]
    assignment = np.random.randint(len(sites), size=k)
    prev_centroid = []
    centroids = []
    
    for index, i in enumerate(assignment):
        clusters[index].add(sites[i])
        centroids.append(location(sites[i]))


    iterations = 0
    while not kmeans_stop(max_iter, iterations, prev_centroid, centroids):
        iterations += 1

        prev_centroid = centroids
        clusters = [set() for x in range(k)]
        for site in sites:
            # assign each site to a cluster
            dist = []
            for center in centroids:
                dist.append(compute_similarity(location(site), center))
            val, idx = min((val, idx) for (idx, val) in enumerate(dist))
            clusters[idx].add(site)

        # calculate centroid of each cluster
        centroids = []
        for cluster in clusters:
            if cluster != set():
                centroids.append(kmeans_centroid(cluster))
            else:
                random = np.random.randint(len(sites))
                centroids.append(location(sites[random]))
    print(dbi(clusters))

    print(clusters)

    return clusters



def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    k = 5
    if k > len(active_sites):
        raise ValueError('k = %s is greater than number of sites = %s' % (str(k), str(len(active_sites))))
    clusters = kmeans(active_sites, k, 10)

    return [[x for x in cluster] for cluster in clusters]


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    dbis = []  # track DBI of each clustering to pick optimal
    master_clust = [[] for x in range(len(active_sites))]  # track cluster history
    clusters = []  # holds the current clustering scheme
    for site in active_sites:
        tmp = set()
        tmp.add(site)
        clusters.append(tmp)

    master_clust[0] = clusters[:] 
    # this holds the "history" of clusters
    # first entry is a list of sets each containing one active site
    # final entry is a list of one set containing all active sites
    # dbi is used to pick optimal clustering from this master list


    centroids = []
    for cluster in clusters:
        centroids.append(kmeans_centroid(cluster))

    track = 0
    while len(clusters) >= 2:
        track += 1
        # build diagonal matrix to track distances
        # values that are zero or below the diagonal are set to inf
        #  so there is no special logic needed to find the minimum
        distance = np.zeros((len(clusters), len(clusters)), dtype=float)
        distance[:] = np.inf  # set 
        for i in range(len(clusters) - 1):
            for j in range(i + 1, len(clusters)):
                distance[i][j] = compute_similarity(centroids[i], centroids[j])

        # find the minimum distance; these two clusters will be merged
        min_ind = np.unravel_index(np.argmin(distance), (len(clusters), len(clusters)))

        cluster1 = clusters[min_ind[0]]
        cluster2 = clusters[min_ind[1]]
        add_clust = set()

        # create a set of the new sites
        for item in cluster2:
            add_clust.add(item)
        for item in cluster1:
            add_clust.add(item)

        # add the new site, remove the old
        clusters.append(add_clust)
        clusters.remove(cluster1)
        clusters.remove(cluster2)

        # calculate the DBI to determine the optimal clustering to use
        if len(clusters) > 1:
            dbis.append(dbi(clusters))

        master_clust[track] = clusters[:]  # create a copy of the current clusters

    min_dbi  = min(enumerate(dbis), key=lambda x: x[1] if x[1] > 0 else float('inf'))

    # print out the optimal cluster

    return ([[[x for x in cluster] for cluster in clusters] for clusters in master_clust], 
        [[y for y in cluster] for cluster in master_clust[min_dbi[0]]])




def rand_index(p, h):
    """
    calculates adjusted rand index for two clusterings
    ARI is between -1 and 1
    """
    rand_array = np.zeros((len(p), len(h)), dtype=float)
    for i, partition in enumerate(p):
        for j, hierarchical in enumerate(h):
            rand_array[i][j] = len(set(partition) & set(hierarchical))

    index = 0

    tmp_a = 0
    tmp_b = 0

    n = sum([sum(subrow) for subrow in rand_array])

    for k in range(len(rand_array)):
        cols = 0  # column sums
        for l in range(len(rand_array)):
            index += choose(rand_array[i][j], 2)
            cols += rand_array[i][j]
        tmp_b += choose(sum(rand_array[i]), 2)  #row-wise sums
        tmp_a += choose(cols, 2)  #column-wise sums

    expected_ind = (tmp_a*tmp_b) / choose(n, 2)

    max_ind = (tmp_a + tmp_b) / 2

    adj_rand_ind = (index - expected_ind) / (max_ind - expected_ind)
    print('rand')
    print(adj_rand_ind)

def choose(n, k):
    """returns binomial coefficient"""
    return fac(n) // fac(k) //fac(n - k)




