import sys
import re
import struct

# usage: python ml22gsp.py infile.sdf infile.ml2 outfile.gsp

class Mapper:
    def __init__(self):
        self.count = 0
        self.dict = {}
    def get_int(self, atom):
        if self.dict.has_key(atom):
            return self.dict[atom]
        else:
            self.dict[atom] = self.count
            self.count += 1
            return self.dict[atom]

vdict = Mapper()
edict = Mapper()

in_sdf  = sys.argv[1]
in_ml2  = sys.argv[2]
out_gsp = sys.argv[3]

g = open(in_sdf, 'r')
isFirstLine = True
line = g.readline()
molname = ''
activity = {}
while line:
    if isFirstLine:
        molname = line.strip()
        isFirstLine = False
    elif re.match(r'^>  <activity>', line):
        a = g.readline().strip()
        activity[molname] = {'Active':1, 'Inactive':-1}[a]
    elif re.match(r'^\$\$\$\$', line):
        isFirstLine = True
    line = g.readline()

outf = open(out_gsp, 'w')
f = open(in_ml2, 'r')

num_mol = 0

line = f.readline()
while line:
    mols = []
    molname = f.readline().strip() # header line
    hline = f.readline()
    token = hline.split()
    num_of_v = int(token[0])
    num_of_e = int(token[1])
    line = f.readline()
    while line:
        if re.match(r'^@', line):
            break
        line = f.readline()
    for i in range(num_of_v):
        vline = f.readline()
        token = vline.split()
        atom = token[5].strip()
        val = vdict.get_int(atom)
        mols.append('v %d %d' % (i, val))
    f.readline()    
    for i in range(num_of_e):
        eline = f.readline()
        token = eline.split()
        e_from  = int(token[1])
        e_to    = int(token[2])
        e_label = token[3].strip()
        val = edict.get_int(e_label) 
        mols.append('e %d %d %d' % (e_from-1, e_to-1, val))
    print >> outf, 't # 0 %d %s' % (activity[molname], molname)
    for m in mols:
        print >> outf, m
    print >> outf, ''
    num_mol += 1
    line = f.readline()

print '%d mols read' % num_mol
print 'num of v:', vdict.count
print vdict.dict
print 'num of e:', edict.count
print edict.dict

