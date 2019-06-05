import sys

args = sys.argv
fname = args[1]

with open(fname, "r") as f:
    for line in f:
        token = line.strip().split()
        if len(token) > 1:
            print(' '.join(token[1:]))

