#program to plot Q, U, angle, reconstructed angle and fit vs. lambda sq.


import numpy as np
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser()

parser.add_argument("q_in")
parser.add_argument("u_in")
parser.add_argument("synrm_in")
parser.add_argument("cuberm_in")
parser.add_argument("zero_angle")
parser.add_argument("weights_in")
parser.add_argument("pixel_coords",help="x and y of desired pixel to plot",nargs=2)

args=parser.parse_args()

weights=np.loadtxt(args.weights_in)
goodch=np.nonzero(weights)[0]

qhead=fits.getheader(args.q_in)

nu_size=qhead['NAXIS3']
dnu=qhead['CDELT3']
nuref=qhead['CRVAL3']

print "Computing l2"
    
c2=299792458.**2
    
dl2 = c2/(dnu**2)

nu = np.arange(nu_size)*dnu+nuref

l2 = 0.5 * c2 * ((nu - 0.5 * dnu) ** -2 + (nu + 0.5 * dnu) ** -2)
l2 = np.flipud(l2)

cchan=nu_size/2

x=int(args.pixel_coords[0])
y=int(args.pixel_coords[1])

print x,y


print "Opening Q and U"

#open Q and U and extract LOS

qcube=fits.getdata(args.q_in)

qlos=qcube[:,y,x]

ucube=fits.getdata(args.u_in)

ulos=ucube[:,y,x]

#open rmsyn map

synrm=fits.getdata(args.synrm_in)

#compute ``raw" angle and generate rm best-guess angle

angle=0.5*np.arctan2(ulos,qlos)

cangle=angle[cchan]

target=synrm[y,x]*(l2-l2[cchan])+cangle

npi=np.around((target-angle)/np.pi)

rm_angle=angle+np.pi*npi

#open cube fit rm and angle map and extract fit coeffs

cuberm=fits.getdata(args.cuberm_in)

cubeang=fits.getdata(args.zero_angle)

a=cuberm[y,x]
b=cubeang[y,x]

l2_good=l2[goodch]
q_good=qlos[goodch]
u_good=ulos[goodch]
angle_good=angle[goodch]
rm_angle_good=rm_angle[goodch]
angle_fit=a*l2_good+b
diff=target-angle


#plot everything
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
fig=plt.figure()
ax1=fig.add_subplot(111)
#ax1.plot(l2_good,q_good,'b:',label='Q')
#ax1.plot(l2_good,u_good,'g:',label='U')
#ax1.plot(l2_good,np.ones(len(goodch))*np.pi/2.,'k--')
#ax1.plot(l2_good,np.ones(len(goodch))*-np.pi/2.,'k--')
#ax1.plot(l2_good,diff[goodch],'c-',label='diff')
ax1.plot(l2_good,target[goodch],'b-',label='target')
#ax1.plot(l2_good,npi[goodch],'b-',label='npi')
ax1.plot(l2_good,angle_good,'gx',label='raw angle')
ax1.plot(l2_good,rm_angle_good,'r+',label='pred. angle')
ax1.plot(l2_good,angle_fit,'k-', label='linear fit')
plt.xlabel('$\lambda^2$')
plt.ylabel('pol. angle [rad]')

plt.legend(loc="upper left")
leg = plt.gca().get_legend()
ltext=leg.get_texts()
plt.setp(ltext,fontsize='small')

outfile=str(x)+'_'+str(y)+'_l2_plot.pdf'

plt.savefig('/local2/scratch/GALFACTS/rm_synthesis/rmsyn/s1/'+outfile,bbox_inches='tight')





