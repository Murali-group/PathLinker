# tests all paths of length 3 from pathlinker and cytoscape

from sets import Set

plRead = open("TGF_beta_1000-paths-pathlinker.txt", "r")
csRead = open("TGF_beta_1000-paths-cytoscape.txt", "r")

plSet = Set([])
csSet = Set([])

for line in plRead:
	plSet.add(line)
	
for line in csRead:
	csSet.add(line)
	
print plSet == csSet