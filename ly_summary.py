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
import csv

bias = [63,65,67,69,71]

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

ly_theo = 0.89
ly_theo_low = 0.684
ly_theo_high = 1.062

ly_top = np.array([[0.1475808886427955, 0.0007715836820109614], [0.1585964260959682, 0.0007931781937859157], [0.1697453980806501, 0.0007166121754199277], [0.1759975775049456, 0.0007159596763145591], [0.18678717660479535, 0.0008976510532069087]])
ly_bot = np.array([[0.8312414554260226, 0.0026280867682118027], [0.8942028907171082, 0.0031565289942957534], [0.9686989293017733, 0.004094245145486982], [1.0137736017623442, 0.004774948643927789], [1.0698533447130885, 0.004912927411632563]])
alpha_top = np.array([[0.06424253995801037, 0.014418556315187906], [0.06624683135008762, 0.015518252189349272], [0.06445076988315324, 0.011803558194587609], [0.0675716607348285, 0.011781073498932036], [0.06713486045477822, 0.014007084409987634]])
alpha_bot = np.array([[0.07919821499844223, 0.005449666196894078], [0.08373716801593523, 0.006071637316807423], [0.07810559596660509, 0.007244255239888636], [0.07904752972983427, 0.007564209950160533], [0.07219723963673799, 0.007818329664193305]])
ly_top_toterr = np.sqrt(np.array(ly_top)[:,1]**2 + (np.array(ly_top)[:,0]**2*syserr2_top))
ly_bot_toterr = np.sqrt(np.array(ly_bot)[:,1]**2 + (np.array(ly_bot)[:,0]**2*syserr2_bot))

plt.figure(1,figsize=(3,3))
plt.errorbar(bias,np.array(ly_top)[:,0],yerr=ly_top_toterr, fmt='o', label='Tyvek', ls='none', elinewidth=2, capsize=2, markersize=2)
plt.errorbar(bias,np.array(ly_bot)[:,0],yerr=ly_bot_toterr, fmt='o', label='ESR', ls='none', elinewidth=2, capsize=2, markersize=2)
plt.xlabel(r'Bias Voltage [$\rm V$]')
plt.ylabel(r'Light Yield [$\rm PE/keV$]')
plt.legend()
plt.grid()
plt.minorticks_on()
plt.xlim(62,72)
plt.ylim(0,1.2)
plt.savefig('ly_top_bot.png')

plt.figure(2,figsize=(3,3))
ly_ratio = np.array(ly_top)[:,0]/np.array(ly_bot)[:,0]
ly_ratio_err = ly_ratio * np.sqrt((ly_bot_toterr/np.array(ly_bot)[:,0])**2+(ly_top_toterr/np.array(ly_top)[:,0])**2)
plt.errorbar(bias,ly_ratio,yerr=ly_ratio_err, fmt='o', ls='none', elinewidth=2, capsize=2, markersize=2, label='Experiment')
plt.fill_between(x=[62,72],y1=[ly_theo_low]*2, y2=[ly_theo_high]*2, facecolor='k',alpha=0.2)
plt.plot([62,72], [ly_theo]*2, 'k--', label='Theory')

lyr_mu = np.sum(ly_ratio/ly_ratio_err**2)/np.sum(1/ly_ratio_err**2)
lyr_sigma = 1/np.sqrt(np.sum(1/ly_ratio_err**2))
print('Average Tyvek-to-ESR Ly ratio = {:.4f}+/-{:.4f}'.format(lyr_mu, lyr_sigma))
variation = (np.max(ly_ratio)-np.min(ly_ratio))/2
print('max-min/2 = {:.4f}'.format(variation)) 
print('max-min/2/avg = {:.2f}%'.format(variation/lyr_mu*100)) 
# total_sys = np.sqrt(variation**2 + syserr_calib**2)
# print('Total systematics = {:.2f}%'.format(total_sys*100))
# print('Total systematics = {:.4f}PE/keV'.format(total_sys*lyr_mu))
plt.fill_between
plt.xlabel(r'Bias Voltage [$\rm V$]')
plt.ylabel(r'$L_{y}^{\rm Tyvek}/L_{y}^{\rm ESR}$')
plt.grid()
plt.minorticks_on()
plt.legend()
plt.xlim(62,72)
plt.ylim(0,1.4)
plt.savefig('ly_ratio.png')

plt.figure(3,figsize=(3,3))
plt.errorbar(bias,np.array(alpha_top)[:,0],yerr=np.array(alpha_top)[:,1], fmt='o', label='Tyvek', ls='none', elinewidth=2, capsize=2, markersize=2)
plt.errorbar(bias,np.array(alpha_bot)[:,0],yerr=np.array(alpha_bot)[:,1], fmt='o', label='ESR', ls='none', elinewidth=2, capsize=2, markersize=2)
plt.xlabel(r'Bias Voltage [$\rm V$]')
plt.ylabel(r'Non-uniformity $\alpha$')
plt.legend()
plt.grid()
plt.minorticks_on()
plt.xlim(62,72)
plt.ylim(0,0.14)
plt.savefig('alpha_top_bot.png')