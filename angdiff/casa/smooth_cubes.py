<<<<<<< HEAD
<<<<<<< .merge_file_NUVMMY
#This is just a set of CASA commands to smooth a cube, the example here is S1.
#Not meant to be run out of the box, needs tweaking for whatever you're going to do with it.

#cube to import


#define the smoothing beam

galbeam={"major":"36arcmin","minor":"36arcmin","pa":"0deg"}

importfits("s1cube","s1cube.image","whichhdu=0")

imsmooth(imagename='s1cube.image',beam=galbeam,targetres=True,outfile='s1cubesmooth.image')

exportfits('s1cubesmooth.image','s1cubesmooth.fits')

#rinse and repeat 
=======
=======
>>>>>>> bfab1cee915774ee0026ff6e68d072827b950814
import os

print "Running imsmooth script, producing smooth galfacts Q & U cubes for comparison with DRAO maps"

galbeam={"major":"36arcmin","minor":"36arcmin","pa":"0deg"}

####################################################
########               S1                  #########
####################################################

print "smoothing S1"
<<<<<<< HEAD
>>>>>>> .merge_file_slh2vX
=======
>>>>>>> bfab1cee915774ee0026ff6e68d072827b950814
