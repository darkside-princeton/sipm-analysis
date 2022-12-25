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

def chisquare_two_hist_new(ly, alpha, data_hist, data_hist_bins, data_hist_err, tree, p_dict, range_pe):
    xmin = data_hist_bins[0]
    xmax = data_hist_bins[-1]
    nbins = len(data_hist_bins)-1
    bin_width = (xmax-xmin)/nbins
#     print(np.shape(data_hist), np.shape(data_hist_bins), np.shape(data_hist_err))
    
    simulated_pe = []
    for i,ev in enumerate(tree):
        lyr = normal(ly, ly*alpha)
        simulated_pe.append(normal(lyr*ev.Edep, np.sqrt((1+p_dict)*lyr*ev.Edep)))
    hSimPE, hSimPE_bins = np.histogram(simulated_pe, bins=np.linspace(xmin,xmax,nbins+1))
    hSimPE_err = np.sqrt(hSimPE)
#     print(np.shape(hSimPE), np.shape(hSimPE_bins), np.shape(hSimPE_err))
    simulated_pe = []
    range_bin = [int((range_pe[0]-xmin)/bin_width), int((range_pe[1]-xmin)/bin_width)]
    norm_sim = np.sum(data_hist[range_bin[0]:range_bin[1]])/np.sum(hSimPE[range_bin[0]:range_bin[1]])
    hSimPE = hSimPE*norm_sim
    hSimPE_err = hSimPE_err*norm_sim
    chi_square = 0
    for i in range(range_bin[0], range_bin[1]):
        chi_square += (data_hist[i] - hSimPE[i])**2/((data_hist_err[i])**2 + (hSimPE_err[i])**2)
    dof = range_bin[1]-range_bin[0] - 4
    print(ly, alpha, chi_square, dof)
    return chi_square

from scipy.stats import multivariate_normal
from scipy.stats import norm

def gauss2d(xy, ly, fano, Sx, Sy, theta, c):
    sxx = Sx**2*np.cos(theta)**2+Sy**2*np.sin(theta)**2
    sxy = (Sx**2-Sy**2)*np.sin(theta)*np.cos(theta)
    syy = Sx**2*np.sin(theta)**2+Sy**2*np.cos(theta)**2
    return -2*np.log(2*np.pi*Sx*Sy*multivariate_normal(mean=[ly,fano], cov=[[sxx, sxy],[sxy, syy]]).pdf(xy))+c


bias = [63, 65, 67, 69, 71]

ds1120bot = []
for i,v in enumerate(bias):
    dset = ds.Dataset(path='', mode='', channels=range(4))
    dset.summed_integral_pe = []
    dset.fprompt = []
    file = ROOT.TFile("root/gamma_1120_{}V_bottom.root".format(v), "read")
    tree = file.Get("tree")
    nev = 0
    nev_cut = 0
    for iev, ev in enumerate(tree):
        nev += 1
        cut = True
        for ch in range(4):
            cut = cut and ev.bsl_std[ch]<2.5
        if cut and ev.f_prompt>0.1:
            nev_cut += 1
            dset.summed_integral_pe.append(ev.sum_pe)
    dset.nev = nev
    dset.nev_cut = nev_cut
    print('Bottom {}V {} events loaded. Cut fraction {:.2f}%'.format(v, nev, (1-nev_cut/nev)*100))
    ds1120bot.append(dset)
    
ds1118bot_bkg = []
for i,v in enumerate(bias):
    dset = ds.Dataset(path='', mode='', channels=range(4))
    dset.summed_integral_pe = []
    dset.fprompt = []
    file = ROOT.TFile("root/backgrounds_1118_{}V_bottom.root".format(v), "read")
    tree = file.Get("tree")
    nev = 0
    nev_cut = 0
    for iev, ev in enumerate(tree):
        nev += 1
        cut = True
        for ch in range(4):
            cut = cut and ev.bsl_std[ch]<2.5
        if cut and ev.f_prompt>0.1:
            nev_cut += 1
            dset.summed_integral_pe.append(ev.sum_pe)
    dset.nev = nev
    dset.nev_cut = nev_cut
    print('Bottom {}V {} events loaded. Cut fraction {:.2f}%'.format(v, nev, (1-nev_cut/nev)*100))
    ds1118bot_bkg.append(dset)

xmin = 0
xmax = 1200
nbins = 300
bin_width = (xmax-xmin)/nbins

hist_gamma_bot = []
hist_gamma_bot_err = []
hist_gamma_bot_bins = []
hist_bkg_bot = []
hist_bkg_bot_err = []
hist_bkg_bot_bins = []
hist_dif_bot = []
hist_dif_bot_err = []
hist_dif_bot_bins = []
source = ['Cs-137', 'Backgrounds']
color = ['r', 'g']
bkg_boundary = [900, 900, 900, 900, 900]
# scale = [3.3, 2.2, 1.7, 1.3, 1]
# plt.figure(0,figsize=(4,4))
for i in range(len(ds1120bot)):
    hg,hgx = np.histogram(ds1120bot[i].summed_integral_pe, bins=np.linspace(xmin,xmax,nbins+1))
    hist_gamma_bot_err.append(np.sqrt(hg))
    # hg = hg/dset.ch[0].cumulative_time/bin_width
    norm_ = np.sum(hg[int((bkg_boundary[i]-xmin)/bin_width):])
    hg = hg/norm_
    hist_gamma_bot_err[i] = hist_gamma_bot_err[i]/norm_
    hist_gamma_bot.append(hg)
    hist_gamma_bot_bins.append(hgx)
    # plt.errorbar(0.5*(hgx[1:]+hgx[:-1]), hg, yerr=hist_gamma_bot_err[-1], label='{}V {}'.format(bias[i], source[0]), ls='none', capsize=0.5, elinewidth=0.5, fmt='C{}s'.format(i), markersize=0.5)

    hb,hbx = np.histogram(ds1118bot_bkg[i].summed_integral_pe, bins=np.linspace(xmin,xmax,nbins+1))
    hist_bkg_bot_err.append(np.sqrt(hb))
    norm_ = np.sum(hb[int((bkg_boundary[i]-xmin)/bin_width):])
    hb = hb/norm_
    hist_bkg_bot_err[i] = hist_bkg_bot_err[i]/norm_
    hist_bkg_bot.append(hb)
    hist_bkg_bot_bins.append(hbx)
    # plt.errorbar(0.5*(hbx[1:]+hbx[:-1]), hb, yerr=hist_bkg_bot_err[-1], label='{}V {}'.format(bias[i], source[1]), ls='none', capsize=0.5, elinewidth=0.5, fmt='C{}^'.format(i), markersize=0.5)

# plt.figure(1,figsize=(4,4))
hist_dif_bot_bins = hist_gamma_bot_bins
for i in range(len(hist_gamma_bot)):
    hist_dif_bot.append(hist_gamma_bot[i]-hist_bkg_bot[i])
    hist_dif_bot_err.append(np.sqrt(hist_gamma_bot_err[i]**2 + hist_bkg_bot_err[i]**2))
    # plt.errorbar(0.5*(hist_dif_bot_bins[i][1:]+hist_dif_bot_bins[i][:-1]), hist_dif_bot[i]*10**i, yerr=hist_dif_bot_err[i]*10**i, label='{}V'.format(bias[i]), ls='none', capsize=0.5, elinewidth=0.5, fmt='C{}o'.format(i), markersize=0.5)
    
import csv
p_dict_top = []
p_dict_bot = []
syserr2_top = [0]*5
syserr2_bot = [0]*5
for i,volt in enumerate(bias):
    with open('calibration_1122_{}V.csv'.format(volt)) as f:
        r = csv.reader(f)
        p_top = 0
        p_bot = 0
        line_count = 0
        for row in r:
            if line_count>0:
                if line_count<=4:
                    p_top += float(row[3])
                    syserr2_top[i] += (float(row[6])/float(row[5]))**2 + (float(row[4])/(1-float(row[3])))**2
                else:
                    p_bot += float(row[3])
                    syserr2_bot[i] += (float(row[6])/float(row[5]))**2 + (float(row[4])/(1-float(row[3])))**2
            line_count += 1
        p_top /= 4
        p_bot /= 4
        p_dict_top.append(p_top)
        p_dict_bot.append(p_bot)
print('DiCT probability:')
print(p_dict_top)
print(p_dict_bot)
print('Systematics from calibration:')
print('Top {}'.format(np.sqrt(syserr2_top)))
print('Bottom {}'.format(np.sqrt(syserr2_bot)))

ly_guess = [0.83, 0.895, 0.97, 1.015, 1.07]
ly_bot = []
alpha_bot = []
for iv,volt in enumerate(bias):
    print('Running for bottom {}V'.format(volt))
    data_hist = hist_dif_bot[iv]
    data_hist_bins = hist_dif_bot_bins[iv]
    data_hist_err = hist_dif_bot_err[iv]

    file = ROOT.TFile("root/pu_lar_cs137_Edep.root", "read")
    tr = file.Get("trEdep")

    lys = np.linspace(ly_guess[iv]*0.985, ly_guess[iv]*1.015, 10)
    alphas = np.linspace(0.02, 0.11, 10)
    X,Y = np.meshgrid(lys, alphas)

    # Spectrum fit range
    norm_min = 250
    norm_max = 900
    chi2map = np.array([[chisquare_two_hist_new(ly_, alpha_, data_hist, data_hist_bins, data_hist_err, tr, p_dict_bot[iv], [norm_min,norm_max]) for ly_ in lys] for alpha_ in alphas])
    min_chi2 = np.min(chi2map)
    chi2map = chi2map - min_chi2

    # Fit chi-square to find minimum and errors
    popt,pcov = curve_fit(gauss2d, np.stack((X, Y),axis=-1).reshape(-1,2), chi2map.flatten(), p0=[ly_guess[iv], 0.08, 0.002, 0.006, 0.2, 150], maxfev=10000)

    # Print best fit values and 68% CL
    ly_fit = popt[0]
    alpha_fit = popt[1]
    ly_68 = 1.515*np.sqrt(popt[2]**2*np.cos(popt[4])**2+popt[3]**2*np.sin(popt[4])**2)
    alpha_68 = 1.515*np.sqrt(popt[2]**2*np.sin(popt[4])**2+popt[3]**2*np.cos(popt[4])**2)
    print('Ly={:.3f}+/-{:.3f} PE/keV  alpha={:.3f}+/-{:.3f} (68% CL)'.format(ly_fit, ly_68, alpha_fit, alpha_68))
    ly_bot.append([ly_fit, ly_68])
    alpha_bot.append([alpha_fit, alpha_68])
    
    simPE = []
    nev = 0
    for i,ev in enumerate(tr):
        for j in range(1):
            ly_random = normal(loc=ly_fit,scale=ly_fit*alpha_fit) 
            simPE.append(normal(ly_random*tr.Edep, np.sqrt((1+p_dict_bot[iv])*ly_random*tr.Edep)))
            nev += 1
    print('{} events simulated'.format(nev))
    file.Close()

    ymin = 1e-2
    ymax = 1e2
    xmin = data_hist_bins[0]
    xmax = data_hist_bins[-1]
    nbins = len(data_hist_bins)-1
    bin_width = (xmax-xmin)/nbins
    hSim, hSim_bins = np.histogram(simPE, bins=np.linspace(xmin,xmax,nbins+1))
    hSimErr = np.sqrt(hSim)
    norm_ = np.sum(data_hist[int((norm_min-xmin)/bin_width):int((norm_max-xmin)/bin_width)])/np.sum(hSim[int((norm_min-xmin)/bin_width):int((norm_max-xmin)/bin_width)])
    hSim = hSim*norm_
    hSimErr = hSimErr*norm_
    chi2_inplot = np.sum(((hSim-data_hist)**2/(hSimErr**2 + data_hist_err**2))[int((norm_min-xmin)/bin_width):int((norm_max-xmin)/bin_width)])
    plt.figure(iv,figsize=(3,3))
    plt.errorbar(0.5*(data_hist_bins[1:]+data_hist_bins[:-1]), data_hist, yerr=data_hist_err, fmt='o', label='Data', ls='none', capsize=0.5, elinewidth=0.5, markersize=0.5)
    plt.errorbar(0.5*(hSim_bins[1:]+hSim_bins[:-1]), hSim, yerr=hSimErr, label=r'Simulation $L_y={:.3f}$PE/keV $\alpha$={:.3f}'.format(ly_fit, alpha_fit), fmt='o', ls='none', capsize=0.5, elinewidth=0.5, markersize=0.5)
    plt.plot([norm_min, norm_min], [ymin, ymax], 'k--', linewidth=1)
    plt.plot([norm_max, norm_max], [ymin, ymax], 'k--', linewidth=1)
    plt.text(xmax*0.8*0.5, ymax/10, r'$\chi^2/\rm DoF={:.1f}/{:.0f}$'.format(chi2_inplot,(norm_max-norm_min)/bin_width-4), fontsize=6,bbox=dict(facecolor='white',edgecolor='grey',alpha=1))
    plt.grid()
    plt.minorticks_on()
    plt.yscale('log')
    plt.ylabel('Rate [A.U.]')
    plt.xlabel(r'$N_{\rm PE}$')
    plt.legend(loc='upper right', fontsize=6)
    # plt.title('11/20 Top chamber {}V'.format(volt))
    plt.xlim(xmin, xmax*0.8)
    plt.ylim(ymin, ymax)
    plt.savefig('cs137_bot_{}v_fit.png'.format(bias[iv]))

    # (Ly,F) range
    XX,YY = np.meshgrid(np.linspace(lys[0],lys[-1],100),np.linspace(alphas[0],alphas[-1],100))
    plt.figure(iv+10,figsize=(4,3))
    plt.contourf(X,Y,chi2map,50)
    plt.colorbar(label=r'$\Delta\chi^2$')
    cs = plt.contour(XX,YY,gauss2d(np.dstack((XX,YY)),*popt),[popt[5]-2*np.log(2*(1-norm.cdf(1,0,1))), popt[5]-2*np.log(2*(1-norm.cdf(2,0,1)))], colors=['tab:red', 'tab:orange'])
    plt.clabel(cs, fmt={cs.levels[0]:'68\% CL', cs.levels[1]:'95\% CL'},fontsize=7)
    plt.plot([ly_fit],[alpha_fit],'wX', markersize=5)
    plt.xlabel(r'$L_y$ [PE/keV]')
    plt.ylabel(r'$\alpha$')
    plt.minorticks_on()
    plt.savefig('chisq_bot_{}v.png'.format(bias[iv]))

print(ly_bot)
print(alpha_bot)