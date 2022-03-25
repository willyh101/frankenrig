function orthoOris = pd2oo(prefOris)

o1 = mod(prefOris - 90, 360);
o2 = mod(prefOris + 90, 360);
orthoOris = [o1',o2'];