import numpy as np
import os

## Location of all input files:
path_asc = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_LAVegMod\Snedden_mapping_scheme' 
os.chdir(path_asc)

#for Y in years: # come back here when doing all 50 years
#    print('On year '+ str(Y))

 #   if Y <10:
 #       #read in the LAVegMod output coverage values
 #       LVMout = np.genfromtxt('MP2023_'+ S +'_'+ G+'_' + C+ '_' + U + '_' + V + '_SLA_O_0'+str(Y)+'_0'+str(Y)+'_V_vegty.asc+',skip_header=372, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623
 #       #read in the coverage names
 #       sp_names = np.genfromtxt('MP2023_'+ S +'_'+ G+'_' + C+ '_' + U + '_' + V + '_SLA_O_0'+str(Y)+'_0'+str(Y)+'_V_vegty.asc+',skip_header=371, skip_footer = 173898, delimiter=',', dtype='str') 
 #       LVMout = LVMout[LVMout[:, 0].argsort()]
  #  else: 
  #      LVMout = np.genfromtxt('MP2023_'+ S +'_'+ G+'_' + C+ '_' + U + '_' + V + '_SLA_O_' +str(Y)+'_'+str(Y)+'_V_vegty.asc+',skip_header=372, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623    
  #      sp_names = np.genfromtxt('MP2023_'+ S +'_'+ G+'_' + C+ '_' + U + '_' + V + '_SLA_O_' +str(Y)+'_'+str(Y)+'_V_vegty.asc+',skip_header=371, skip_footer = 173898, delimiter=',', dtype='str')
  #      LVMout = LVMout[LVMout[:, 0].argsort()]

LVMout = np.genfromtxt('MP2023_S07_G500_C000_U00_V00_SLA_O_03_03_V_vegty.asc+',skip_header=372, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623
LVMout = LVMout[LVMout[:, 0].argsort()]

sp_names = np.genfromtxt('MP2023_S00_G000_C000_U00_V00_SLA_I_00_00_V_vegty.asc+',skip_header=371, skip_footer = 173898, delimiter=',', dtype='str')
sp_names = (np.char.strip(sp_names,'"'))

LVMgrid = np.genfromtxt('MP2023_S00_G000_C000_U00_V00_SLA_I_00_00_V_vegty.asc+',skip_header=6, skip_footer = 173899, delimiter=' ', dtype=int) # skip the top portion of the asc+ files and only read in the column data starting on line 623

## Location of output: 
path = r'C:\Users\madel\Documents\UNO_CPRA_LAVegMod\One_species_test_Jan2022'


for species_col in range(1,len(sp_names)-15): #cut off cell ID and the last 15 columns, which are barrier island species and summary values
## Or set species_col to the species of interest and change the indentation below
    mapping_dict = {} # empty dictionary that you will fill with grid_id as keys 
    for row in LVMout:
            gid = row[0]
            val = row[species_col]
            mapping_dict[gid] = val            # you’ll want to ensure that the format for gid is the same type as gid_val in the script below – if you try to pass a string key and it’s looking for an integer key, you’ll get errors
    print(sp_names[species_col])
    
    asc_header = 'ncols 1052\nnrows 365\nxllcorner 404710\nyllcorner 3199480\ncellsize 480\nNODATA_value -9999\n'
    asc_out_file = r'%s\%s.asc' % (path, 'MP2023_S00_G000_C000_U00_V00_SLA_I_03_03_V_'+sp_names[species_col])
    
    with open(asc_out_file, mode='w') as outf:
        outf.write(asc_header)
        for row in LVMgrid:
            nc = 0           
            for col in row:
                gid_map = row[nc]
                if gid_map > 0:
                    gid_val = mapping_dict[gid_map] # here mapping dict is a dictionary with grid cell ID as the key and some value for each key - if the ASC grid has a no data cell (-9999) then it will not have a key in the dictionary and is written out as no data (-9999)
                else:
                    gid_val = gid_map
                if nc == 0:
                    rowout = '%s'  % gid_val
                else:
                    rowout = '%s %s' % (rowout,gid_val)
                nc += 1
            outf.write('%s\n' % rowout)