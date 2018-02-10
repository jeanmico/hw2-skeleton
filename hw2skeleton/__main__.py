import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings, remove_duplicates
from .cluster import cluster_by_partitioning, cluster_hierarchically, rand_index

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H | -C] <pdb directory> <output file>")
    sys.exit(0)

print(sys.argv[2])

remove_duplicates(sys.argv[2])

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites)
    write_clustering(sys.argv[3], clustering)

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clusterings, final = cluster_hierarchically(active_sites)
    write_mult_clusterings(sys.argv[3], clusterings)

if sys.argv[1][0:2] == "-C":
	print("comparing clusters")
	clusterings, h = cluster_hierarchically(active_sites)
	p = cluster_by_partitioning(active_sites)
	rand_index(p, h)

