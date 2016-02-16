# tests 16501 weighted paths of log-transformed TNFalpha pathway

plRead = open("TNFalpha-weighted-paths-pathlinker.txt", "r")
csRead = open("TNFalpha-weighted-paths-cytoscape.txt", "r")

# dict that stores path: weight
pl = {}

epsilon = 0.0001

def sameFloat(a, b):
	return abs(a - b) <= epsilon

for line in plRead:
	if line[0] == '#':
		continue
	words = line.split( )
	weight = float(words[1])
	path = words[2]
	pl[path] = weight
	
errors = 0
	
for line in csRead:
	words = line.split( )
	weight = float(words[0])
	path = words[1]
	if path not in pl:
		print path + " " + str(weight)
		errors += 1
	else:
		if not sameFloat(weight, pl[path]):
			print path + " " + str(weight)
			errors += 1
			
print errors == 0
	