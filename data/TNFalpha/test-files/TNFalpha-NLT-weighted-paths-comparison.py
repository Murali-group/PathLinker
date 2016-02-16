# tests 10341 weighted paths of TNFalpha pathway

plRead = open("TNFalpha-weighted-paths-pathlinker.txt", "r")
csRead = open("TNFalpha-weighted-paths-cytoscape.txt", "r")

# dict that stores path: weight
pl = {}

epsilon = 0.0001

def sameFloat(a, b):
	return abs(a - b) <= epsilon

for line in plRead:
	words = line.split( )
	weight = float(words[0])
	path = words[1]
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
	