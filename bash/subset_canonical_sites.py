#!/bin/python3

import sys

if len(sys.argv) != 4:
    print("Usage: python3 subset_singleton_sites.py <singleton_file> <sites_file> <out_file>")
    sys.exit(1)

singleton_file = sys.argv[1]
sites_file = sys.argv[2]
out_file = sys.argv[3]

#1 read singleton file

singleton = set()

with open(singleton_file,"r") as f:
    for line in f:
        chrom, pos = line.strip().split()
        singleton.add((chrom, pos))

#2 read sites file

sites = {}

with open(sites_file,"r") as f:
    for line in f:
        fields = line.strip().split()
        chrom, pos = (fields[0], fields[1])
        sites[(chrom,pos)] = fields

#3 extract singletons in sites

singleton_sites = []

for key in singleton:
    if key in sites:
        singleton_sites.append(sites[key])

#4 verification step

if len(singleton) != len(singleton_sites):
    print("Error singleton_sites file has fewer line than the number of singleton")
else:
    print(f"Singleton sites extracted successfully for {singleton_file}")

#4 write output file

with open(out_file,"w") as out:
    for row in singleton_sites:
        out.write("\t".join(row) + "\n")


