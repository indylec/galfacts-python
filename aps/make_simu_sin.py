import numpy as np
import glob
import struct

beams=[]
for b in range(0,7):
    beams.append("beam"+str(b))

print "beams used: ",beams

sim_dir="/local/scratch/GALFACTS/binary_tests/sin_test3/"#"/researchdata/fhgfs/arecibo-scratch/simulated/"
infile="average.dat"  #change these as required

#TOD file has structure RA DEC AST I Q U V, with 4byte values

mjds=glob.glob('5*')
mjds.sort()
#iterate over mjds, beams
for day in mjds:
    for beam in beams:
        print "Converting day {0}, {1}".format(day,beam)
        datfile=sim_dir+day+"/"+beam+"/"+infile
        try: #sometimes the beam/day combo has no file
             with file(datfile,'rb+') as f:
                fstring=f.read()
                nosamples=int(struct.unpack_from('i',fstring)[0])#1st 4 bytes is int with number of samples
                print nosamples, " samples in this file"
                print "replacing q values"
                xvals=np.arange(nosamples)*8.*np.pi/nosamples
                sinqu=2.*(np.sin(xvals))**2
                for samp in range (nosamples):
                    f.seek(20+samp*28)#q values start at the 20th byte and recur every 28 bytes (4*7)
                    qstring=f.read(4)
                    qval=struct.unpack_from('f',qstring)[0]
                    qval*=sinqu[samp]
                    f.seek(20+samp*28)#go back to the place you want to write over
                    f.write(struct.pack('f',qval))
                print "replacing u values"
                for samp in range (nosamples):
                    f.seek(24+samp*28)
                    ustring=f.read(4)
                    uval=struct.unpack_from('f',ustring)[0]
                    uval*=sinqu[samp]
                    f.seek(24+samp*28)#go back to the place you want to write over
                    f.write(struct.pack('f',uval))
        except IOError:
            print datfile, " missing. Moving to next file."
