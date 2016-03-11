import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('font',family='serif')
rc('text', usetex=True)
import sys


s1='/Users/leclercq/galfacts/aps/rm/thesis_plot_S1_RMSF.txt'
s2='/Users/leclercq/galfacts/aps/rm/thesis_plot_S2_RMSF.txt'
s3='/Users/leclercq/galfacts/aps/rm/thesis_plot_S3_RMSF.txt'
s4='/Users/leclercq/galfacts/aps/rm/thesis_plot_S4_RMSF.txt'
n2='/Users/leclercq/galfacts/aps/rm/thesis_plot_N2_RMSF.txt'
n3='/Users/leclercq/galfacts/aps/rm/thesis_plot_N3_RMSF.txt'
n4='/Users/leclercq/galfacts/aps/rm/thesis_plot_N4_RMSF.txt'


phi_range,s1_rmtf_re,s1_rmtf_im=np.loadtxt(s1,unpack=True)
phi_range,s2_rmtf_re,s2_rmtf_im=np.loadtxt(s2,unpack=True)
phi_range,s3_rmtf_re,s3_rmtf_im=np.loadtxt(s3,unpack=True)
phi_range,s4_rmtf_re,s4_rmtf_im=np.loadtxt(s4,unpack=True)
phi_range,n2_rmtf_re,n2_rmtf_im=np.loadtxt(n2,unpack=True)
phi_range,n3_rmtf_re,n3_rmtf_im=np.loadtxt(n3,unpack=True)
phi_range,n4_rmtf_re,n4_rmtf_im=np.loadtxt(n4,unpack=True)

s1_rmtf_abs=np.sqrt(s1_rmtf_re**2+s1_rmtf_im**2)
s2_rmtf_abs=np.sqrt(s2_rmtf_re**2+s2_rmtf_im**2)
s3_rmtf_abs=np.sqrt(s3_rmtf_re**2+s3_rmtf_im**2)
s4_rmtf_abs=np.sqrt(s4_rmtf_re**2+s4_rmtf_im**2)
n2_rmtf_abs=np.sqrt(n2_rmtf_re**2+n2_rmtf_im**2)
n3_rmtf_abs=np.sqrt(n3_rmtf_re**2+n3_rmtf_im**2)
n4_rmtf_abs=np.sqrt(n4_rmtf_re**2+n4_rmtf_im**2)

fig=plt.figure(figsize=(10,6))

ax=fig.add_subplot(111)

ax.plot(phi_range,s1_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,s2_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,s3_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,s4_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,n2_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,n3_rmtf_abs,'k-',linewidth=2, alpha=0.5)
ax.plot(phi_range,n4_rmtf_abs,'k-',linewidth=2, alpha=0.5)

ax.set_ylabel('$|\mathrm{R}(\phi)|$',fontsize=24)
ax.set_xlabel('$\phi[\mathrm{rad}/\mathrm{m}^2]$', fontsize=26)
#ax.yaxis.set_label_position('right')
#ax.yaxis.tick_right()
ax.tick_params(labelsize=22)
ax.xaxis.set_major_locator(MaxNLocator(prune='both'))

fig.savefig("/Users/leclercq/galfacts/aps/rm/all_fields_RMSF_plot.pdf",dpi=200,bbox_inches="tight")
