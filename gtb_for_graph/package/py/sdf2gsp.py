import sys
import re
import struct

# usage: python sdf2gsp.py infile.sdf outfile.gsp

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
out_gsp = sys.argv[2]

outf = open(out_gsp, 'w')
f = open(in_sdf, 'r')

num_mol = 0

line = f.readline()
while line:
    mols = []
    activity = 0
    molname = line.strip() # header line
    f.readline()           # comment line
    f.readline()           # blank line
    hline = f.readline()
    token = struct.unpack('3s3s', hline[:6])
    num_of_v = int(token[0])
    num_of_e = int(token[1])
    for i in range(num_of_v):
        vline = f.readline()
        token = struct.unpack('10s10s10s1s3s', vline[:34])
        atom = token[4].strip()
        val = vdict.get_int(atom)
        mols.append('v %d %d' % (i, val))
    for i in range(num_of_e):
        eline = f.readline()
        token = struct.unpack('3s3s3s', eline[:9])
        e_from  = int(token[0])
        e_to    = int(token[1])
        e_label = int(token[2])
        val = edict.get_int(e_label) 
        mols.append('e %d %d %d' % (e_from-1, e_to-1, val))
    mols.append('')
    line = f.readline()
    while line:
        if re.match(r'^>', line):
            s = f.readline().strip()
            activity = {'Active': 1, 'Inactive':-1}[s]
            print >> outf, 't # 0 %d %s' % (activity, molname)
            for t in mols:
                print >> outf, t
            num_mol += 1
        elif re.match(r'^\$\$\$', line):
            break
        line = f.readline()
    line = f.readline()

print '%d mols read' % num_mol
print 'num of v:', vdict.count
print vdict.dict
print 'num of e:', edict.count
print edict.dict

