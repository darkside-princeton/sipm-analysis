import numpy as np
import sipm.sipm as sipm
import sipm.dataset as ds
import sipm.beta as beta
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sipm.constants as const
import ROOT
from numpy.random import normal
plt.style.use('darkside')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

bias = [63, 65, 67, 69, 71]
i=1
ch=0
v=bias[i]
dset = ds.Dataset(path='', mode='', channels=range(4))
dset.summed_integral_pe = []
dset.fprompt = []
file = ROOT.TFile("root/gamma_1120_{}V_bottom.root".format(v), "read")
tree = file.Get("tree")
nev = 0
nev_cut1 = 0
for iev, ev in enumerate(tree):
    nev += 1
    cut = True
    for ch in range(4):
        dset.ch[ch].baseline_avg.append(ev.bsl_avg[ch])
        # dset.ch[ch].baseline_med.append(ev.bsl_med[ch])
        dset.ch[ch].baseline_std.append(ev.bsl_std[ch])
        # dset.ch[ch].acquisition_max.append(ev.acq_max[ch])
        # dset.ch[ch].acquisition_min.append(ev.acq_min[ch])
        cut = cut and ev.bsl_std[ch]<2.5
    if cut:
        nev_cut1 += 1
        dset.fprompt.append(ev.f_prompt)
        dset.summed_integral_pe.append(ev.sum_pe)
dset.nev = nev
dset.nev_cut1 = nev_cut1
print('Bottom {}V {} events loaded. Cut1 fraction {:.2f}%'.format(v, nev, (1-nev_cut1/nev)*100))

import matplotlib.colors as colors

nbinsx = 500
range_minx = -50
range_maxx = 1200
nbinsy = 500
range_miny = -0.1
range_maxy = 1

fp_cut = 0.1
threshold = 20

plt.figure(0,figsize=(3,3))
plt.hist2d(dset.summed_integral_pe, dset.fprompt, bins=[np.linspace(range_minx, range_maxx, nbinsx+1), np.linspace(range_miny, range_maxy, nbinsy+1)], norm = colors.LogNorm())
plt.plot([range_minx, range_maxx], [fp_cut]*2, 'k--', linewidth=1)
# plt.plot([threshold, threshold], [range_miny, range_maxy], 'r--', linewidth=1)
plt.minorticks_on()
plt.grid()
plt.xlabel(r'$N_{\rm PE}$')
plt.ylabel(r'$F_{\rm prompt}$')
plt.savefig('fprompt_bottom_65v.png')
