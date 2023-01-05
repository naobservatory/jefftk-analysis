import sys

state = "want locus"

for line in sys.stdin:
  if state == "want locus":
    if line.startswith("LOCUS"):
      locus = line.split()[1]
      state = "want organism"
  elif state == "want organism":
    if line.startswith("  ORGANISM"):
      state = "want tree"
      tree = []
  elif state == "want tree":
    if line.startswith(" "):
      for x in line.split(";"):
        tree.append(x.strip().removesuffix("."))
    else:
      print ("%s\t%s" % (locus, "; ".join(tree)))
      state = "want locus"
      
    
      
    
