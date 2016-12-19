import sys
amin = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
nucl = ['A','T','G','C']

sel = sys.argv[1]
if(sel == "n"):
	arr = nucl
else:
	arr = amin

for i in range(-1, len(arr), 1):
	line = ""
	for j in range(0, len(arr), 1):
		if(i == -1):
			if(j == 0):
				line = "    " + arr[j] + "  "
				continue
			else:
				line = line + " " + arr[j] + "  "
				continue
		if(j == 0):
			line = arr[i] + "  "
		if(i == j):
			line = line + " 1  "
		else:
			line = line + "-1  "
	print(line)

