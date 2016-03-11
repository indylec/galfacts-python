import numpy as np
import glob
import struct

#theoretical noise value, sigma: Tsys/sqrt(deltanu*tint)

thnoise=26.5/np.sqrt(4.2E6*0.2)

beams=[]
for b in range(0,7):
    beams.append("beam"+str(b))

print "beams used: ",beams

sim_dir="/local/scratch/GALFACTS/noise_simulation/noise_test3/"#"/researchdata/fhgfs/arecibo-scratch/simulated/"
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
                noiseq=np.random.normal(scale=thnoise, size=nosamples)
                for samp in range (nosamples):
                    f.seek(20+samp*28)#q values start at the 20th byte and recur every 28 bytes (4*7)
                    f.write(struct.pack('f',noiseq[samp]))
                noiseq=0.
                print "replacing u values"
                noiseu=np.random.normal(scale=thnoise, size=nosamples)
                for samp in range (nosamples):
                    f.seek(24+samp*28)
                    f.write(struct.pack('f',noiseu[samp]))
                noiseu=0.
        except IOError:
            print datfile, " missing. Moving to next file."
                
