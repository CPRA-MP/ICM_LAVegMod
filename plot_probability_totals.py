import numpy as np
import matplotlib.pyplot as plt

coverages = np.genfromtxt('C:/Eric/git/ICM_LAVegMod/LAVegMod_coverage_attributes.csv',delimiter=',',skip_header=1,dtype='str')
coverage_type = 'Emergent Marsh'

flotant_classes = [4,5,6,7]
emergent_classes = [10,11,12,13]
forest_classes = [8,9]
island_classes = [14]

total_est = np.zeros([39,24])
total_mort = np.zeros([39,24])


for row in coverages:
    sym = row[0]
    cov = int(row[5])
    dis = int(row[6])

    et = 'C:/Eric/git/ICM_LAVegMod/tables/establishment_%s.csv' % sym
    mt = 'C:/Eric/git/ICM_LAVegMod/tables/mortality_%s.csv' % sym

    if cov in emergent_classes:
        e = np.genfromtxt(et,delimiter=',',skip_header=1,usecols=range(1,25))
        m = np.genfromtxt(mt,delimiter=',',skip_header=1,usecols=range(1,25))
        
        try:
            dummy = wlv[0]
        except:
            wlv = np.genfromtxt(et,delimiter=',',usecols=range(1,25))[0]
            sal = np.genfromtxt(et,delimiter=',',skip_header=1,usecols=[0])

        for r in range(0,e.shape[0]):
            for c in range(0,e.shape[1]):
                if e[r][c] > 0:
                    total_est[r][c] = 1

        for r in range(0,m.shape[0]):
            for c in range(0,m.shape[1]):
                if m[r][c] > 0:
                    total_mort[r][c] = 1

fig,ax= plt.subplots()
ax.pcolor(wlv,sal,total_est)
ax.set_xlim(0,0.5)
ax.set_ylim(0,35)
ax.set_xlabel('Water level variability (m)')
ax.set_ylabel('Salinity (ppt)')
fig.suptitle('Non-zero Probability of Establishment: %s' % coverage_type)
