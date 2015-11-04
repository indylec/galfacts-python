import numpy as np
import argparse

parser=argparse.ArgumentParser()

parser.add_argument("field",help="Name of the field to generate weights")
parser.add_argument("bad_chans",help="Bad channels in a given field", nargs='+')

args=parser.parse_args()

bad_chans=args.bad_chans

field=args.field

#print bad_chans, field

bad_chans=np.asarray(bad_chans)

bad_chans=bad_chans.astype(np.int)

weights=np.ones(376)

for chan in bad_chans:
    weights[chan-1]=0.

fileout=field.upper()+'_weights.txt'

np.savetxt(fileout, weights)
    



