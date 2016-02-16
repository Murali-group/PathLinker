# tests all paths of length <= 3

from sets import Set

plRead = open("Wnt-paths-pathlinker.txt", "r")
csRead = open("Wnt-paths-cytoscape.txt", "r")

plSet = Set([])
csSet = Set([])

for line in plRead:
	plSet.add(line)
	
for line in csRead:
	csSet.add(line)
	
print plSet == csSet