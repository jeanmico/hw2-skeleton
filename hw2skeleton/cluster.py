from .utils import Atom, Residue, ActiveSite
import numpy as np
import math

aa_hydro = {"ALA": 1.8, "ARG": -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5, 'GLU': -3.5, 'GLN': -3.5, 'GLY': -.4, 'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -.8, 'THR': -.7, 'TRP': -.9, 'TYR': -1.3, 'VAL': 4.2}
aa_mw = {"ALA": 89.094, "ARG": 174.203, 'ASN': 132.119, 'ASP': 133.104, 'CYS': 121.54, 'GLU': 147.131, 'GLN': 146.146, 'GLY': 75.067, 'HIS': 155.156, 'ILE': 131.175, 'LEU': 131.175, 'LYS': 146.189, 'MET': 149.208, 'PHE': 165.192, 'PRO': 115.132, 'SER': 105.093, 'THR': 119.119, 'TRP': 204.228, 'TYR': 181.191, 'VAL': 117.148}
aa_charge = {"ALA": 0, "ARG": 1, 'ASN': 0, 'ASP': -1, 'CYS': 0, 'GLU': -1, 'GLN': 0, 'GLY': 0, 'HIS':.9, 'ILE': 0, 'LEU': 0, 'LYS': 1, 'MET': -1, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0}

def location(site):
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

def kmeans_centroid(data):
    # convert clusters to tuples
    tuples = []
    for site in data:
        tuples.append(location(site))
    return tuple(np.mean(tuples, 0))

def kmeans_stop(max_iter, iterations, prev_centroid, centroid):
    if iterations > max_iter:
        return True
    if prev_centroid == centroid:
        print("centroids the same")
        return True


def kmeans(sites, k, max_iter):

    #randomly pick 5 points to serve as the first clusters
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

        for site in sites:
            # assign each site to a cluster
            dist = []
            for center in centroids:
                dist.append(compute_similarity(location(site), center))
            val, idx = min((val, idx) for (idx, val) in enumerate(dist))
            clusters[idx].add(site)

        print(len(clusters[0]))
        print(len(clusters[1]))
        print(len(clusters[2]))
        print(len(clusters[3]))
        print(len(clusters[4]))

        # calculate centroid of each cluster
        centroids = []
        for cluster in clusters:
            a = kmeans_centroid(cluster)
            centroids.append(kmeans_centroid(cluster))





def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!
    similarity = np.zeros((len(active_sites), len(active_sites)), dtype=float)
    for i in range(len(active_sites)-1):
        for j in range(i+1, len(active_sites)):
            similarity[i][j] = compute_similarity(location(active_sites[i]), location(active_sites[j]))

    kmeans(active_sites, 5, 10)

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!
    

    return []
