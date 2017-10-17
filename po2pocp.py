#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser("calculate pairwise pocp values from proteinortho-results")
parser.add_argument("-i", "--input", action = "store", required = True, dest = "input_pofile", help = "input proteinorthofile")
parser.add_argument("-o", "--output", action = "store", dest = "output", default = "pairwise_pocp_table.tab", help = "output filename (default = \"pairwise_pocp_table.tab\")")
args = parser.parse_args()

#read in po_file
def read_pofile(pofilename):
	sys.stderr.write("\nreading pofile \"{}\"\n".format(pofilename))
	ignoretokens = [ "# Species", "Genes", "Alg.-Conn."]
	infile = open(pofilename, "r")
	firstline = True
	podict = {}
	for line in infile:
		linetokens = line.rstrip().split("\t")
		if firstline:
			for tindex in range(len(linetokens)):
				if linetokens[tindex] not in ignoretokens:
					podict[tindex] = {"name" : linetokens[tindex], "oglist" : [], "total_ogcount" : 0}
			firstline = False
			continue
		for tindex in range(len(linetokens)):
			if tindex in podict:
				if linetokens[tindex] == "*":
					podict[tindex]["oglist"].append(0)
				else:
					og_gene_count = len(linetokens[tindex].split(","))
					podict[tindex]["oglist"].append(og_gene_count)
					podict[tindex]["total_ogcount"] += og_gene_count
	infile.close()
	return podict

#compare all against all
def pairwise_all_vs_all(podict):
	sys.stderr.write("starting pairwise comparisons\n")
	pocp_dict = {}
	from itertools import combinations
	indices = sorted(podict.keys())
	for pair in combinations(indices, 2):
		pocp = compare_pair(pair, podict)
		pocp_dict[pair] = { "pair" : pair, "pocp" : pocp }
		pocp_dict[pair[::-1]] = { "pair" : pair, "pocp" : pocp }
	return indices, pocp_dict

#compare pair
def compare_pair(pair, podict):
	xorg = podict[pair[0]]
	yorg = podict[pair[1]]
	sys.stderr.write("comparing {x} vs. {y}\n".format(x = xorg["name"], y= yorg["name"]))
	assert len(xorg["oglist"]) == len(yorg["oglist"]), "oglist-length does not match"
	totalgenes = xorg["total_ogcount"] + yorg["total_ogcount"]
	matchgenes = 0
	for og in range(len(xorg["oglist"])):
		if xorg["oglist"][og] != 0 and yorg["oglist"][og] != 0:
			matchgenes += xorg["oglist"][og] + yorg["oglist"][og]
	pocp = float(matchgenes) / totalgenes
	return pocp

#write output
def write_output(outputname, pocp_dict, podict, indices):
	sys.stderr.write("\nwriting to {}\n".format(outputname))
	header = "pairwise_pocp_values\t" + "\t".join([ podict[i]["name"] for i in indices ]) + "\n"
	outlines = [header]
	for indexA in indices:
		outline = podict[indexA]["name"]
		for indexB in indices:
			print " --> {} vs. {} ".format(indexA, indexB)
			if indexA == indexB:
				outline += "\t{}".format(1)
			else:
				outline += "\t{}".format(pocp_dict[(indexA, indexB)]["pocp"])
		outlines.append(outline + "\n")
	outfile = open(outputname, "w")
	for outline in outlines:
		outfile.write(outline)
	outfile.close()

def main():
	podict = read_pofile(args.input_pofile)
	indices, pocp_dict = pairwise_all_vs_all(podict)
	write_output(args.output, pocp_dict, podict, indices)
	sys.stderr.write("\nfinished!\n")

main()
