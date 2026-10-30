import numpy as np
import matplotlib.pyplot as plt

UT = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/tables/mortality_UNIVERSAL.csv',skip_header=1,usecols=range(1,24),dtype='float',delimiter=',')

sals = []           # this will be the array of SAL values representing the percentiles in the original input UNIVERSAL table 
wlvs = []           # this will be the array of WLV values representing the percentiles in the original input UNIVERSAL table
ptiles = [0,0.01]
for p in range(5,100,5):
    ptiles.append(p/100)
ptiles.append(0.99)
ptiles.append(1)


            
sal_cols = range(8,17)
wlv_cols = range(17,26)

covs = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',skip_header=1,dtype='str',delimiter=',',usecols=[0])
grps = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',skip_header=1,dtype='int',delimiter=',',usecols=[5])
niche_pt = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',dtype='str',delimiter=',',usecols=sal_cols)[0]
niche_sal = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',skip_header=1,dtype='float',delimiter=',',usecols=sal_cols)
niche_wlv = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',skip_header=1,dtype='float',delimiter=',',usecols=wlv_cols)


niche_p = []
for t in niche_pt:
    niche_p.append(float(t.split('_')[1])/100.0)


fig0,ax0= plt.subplots()
ax0.set_xlim(0,0.5)
ax0.set_ylim(0,35)
ax0.set_xlabel('Water level variability (m)')
ax0.set_ylabel('Salinity (ppt)')
fig0.suptitle('Probability of Mortality: All species')


for nc in range(0,len(covs)):
    
    sals = []           # this will be the array of SAL values representing the percentiles in the original input UNIVERSAL table 
    wlvs = []           # this will be the array of WLV values representing the percentiles in the original input UNIVERSAL table

    cov = covs[nc]
    grp = grps[nc]


    if grp in [4,5,6,7,9,10,11,12,13]:
    
        for p in ptiles:
            sals.append(np.interp(p,niche_p,niche_sal[nc]))
            wlvs.append(np.interp(p,niche_p,niche_wlv[nc]))

        fig,ax= plt.subplots()
        fig.colorbar(ax.pcolor(wlvs,sals,UT,cmap='coolwarm'),orientation='horizontal')#Greys
        ax0.pcolor(wlvs,sals,UT,cmap='coolwarm')#Greys
        ax.set_xlim(0,0.5)
        ax.set_ylim(0,35)
        ax.set_xlabel('Water level variability (m)')
        ax.set_ylabel('Salinity (ppt)')
        fig.suptitle('Probability of Mortality: %s' % cov)
        plt.savefig('./tables/new/mortality_%s.png' % cov)

    with open('./tables/new/mortality_%s.csv' % cov,mode='w') as oc:
        line = 'mortality_%s' % cov
        for wlv in wlvs:
            line = '%s,%0.3f' % (line,wlv)
        oc.write('%s\n' % line)
        for r in range(0,len(sals)):
            line = '%0.4f' % sals[r]
            for c in range(0,len(wlvs)):
                line = '%s,%0.4f' % (line,UT[r][c])
            oc.write('%s\n' % line)
            
