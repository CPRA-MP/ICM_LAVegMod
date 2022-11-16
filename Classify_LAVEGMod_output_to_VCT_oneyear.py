import numpy as np
import os
import sys

# Snedden Community Vegetation Types (integer values for SOM_species.csv)
#  1: Maidencane          
#  2: Three-Square        
#  3: Rosseau Cane       
#  4: Paspalum            
#  5: Wiregrass         
#  6: Bulltongue         
#  7: Needlerush        
#  8: Bulrush        
#  9: Brackish Mix      
# 10: Oystergrass    
# 11: Saltgrass  
# 12: majority barrier island 
# 13: majority non-SOM bottomland hardwood forest + swamp forest
# 14: No SOM species found (i.e. all water) 


S = 'S%02d' % int(sys.argv[1]) # 'S07'
G = 'G%03d' % int(sys.argv[2]) # 'G500'
years = range(1,52) #[2,12,22,32,42,52]

veg_dir = r'%s/%s/veg' % (S,G)

species_lookup = '/ocean/projects/bcs200002p/ewhite12/code/ICM_LAVegMod/cellxspecies.csv'
SOM_lookup     = '/ocean/projects/bcs200002p/ewhite12/code/ICM_LAVegMod/SOM_species.csv'
BMU_lookup     = '/ocean/projects/bcs200002p/ewhite12/code/ICM_LAVegMod/BMU_to_community.csv'

LVMout_file_excon = r'%s/MP2023_S00_G000_C000_U00_V00_SLA_I_00_00_V_vegty.asc+' % veg_dir
for y in years:
    print(' mapping Snedden Vegetation Community Types: Year %02d' % y)
    if y == 0:
        LVMout_file = LVMout_file_excon
    else:
        LVMout_file = r'%s/MP2023_%s_%s_C000_U00_V00_SLA_O_%02d_%02d_V_vegty.asc+' % (veg_dir,S,G,y,y)
    
    asc_out_file = r'%s/MP2023_%s_%s_C000_U00_V00_SLA_O_%02d_%02d_V_VCT.asc' % (veg_dir,S,G,y,y)
     
  
    LVMout = np.genfromtxt(LVMout_file,skip_header=372, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623
    LVMout = LVMout[LVMout[:, 0].argsort()]
    
    sp_names = np.genfromtxt(LVMout_file_excon,skip_header=371, skip_footer = 173898, delimiter=',', dtype='str')
    sp_names = (np.char.strip(sp_names,'"'))
    
    LVMgrid = np.genfromtxt('%s/veg_grid.asc' % veg_dir, skip_header=6, delimiter=' ', dtype=int)
    
    cellxspecies = np.genfromtxt(species_lookup,delimiter=',',dtype='float') 
    SOM_species = np.genfromtxt(SOM_lookup,delimiter=',',skip_header=1,dtype='str') 
    BMU_to_community = np.genfromtxt(BMU_lookup,delimiter=',',skip_header=1,dtype='int') 
    
    SOM_indices_of_missing_species = []
    SOM_match_indices = []
    veg_match_indices = []
    match = 0
    for j in range(0,len(SOM_species)):
        for i in range(0,len(sp_names)):
            if SOM_species[j] == sp_names[i]:
                match = 1
                veg_match_indices.append(i)
                SOM_match_indices.append(j)
        if match < 1:
            SOM_indices_of_missing_species.append(j)    
        match = 0
    
    
    
    cat = {}
    SOM_dist = {}
    BMU = {}
    
    for row in LVMout:
        cell_id = row[0]
        BI = sum(row[40:48]) #columns 40-47 (when cell ID = column 0) are all barrier island species (have to do one col more because python)
        if BI>= 0.5:
            cat[cell_id] = 13 #majority barrier island
            SOM_dist[cell_id] = 9999
            BMU[cell_id] = 9999
        elif sum(row[49:51])>=0.5: #column 49 and 50 (when cell ID = column 0) are % land bottomland hardhood forest and swamp forest
            cat[cell_id] = 12 #majority forest
            SOM_dist[cell_id] = 9999
            BMU[cell_id] = 9999
        else:        
            form_row = np.zeros(len(SOM_species)) #form_row = formatted row for SOM species
            form_row[SOM_match_indices] = row[veg_match_indices]
            form_row[10] = row[20] + row[15] #pahe2 + pahe2_float #hard coded watch out
            if sum(form_row)== 0:
                cat[cell_id] = 14 #majority something not included in VCT
                SOM_dist[cell_id] = 9999
                BMU[cell_id] = 9999
                continue 
            else:
                form_row = (form_row/sum(form_row))*100
                form_row[SOM_indices_of_missing_species] = 'Nan'
                dist = []
                for n in range(0,len(cellxspecies)):      
                    D = (form_row-cellxspecies[n])
                    DS = D*D
                    dist.append(np.sqrt(np.nansum(DS)))
                SOM_dist[cell_id] = min(dist)
                BMU[cell_id] = int(dist.index(SOM_dist[cell_id])) + 1     
                cat[cell_id] = BMU_to_community[int(BMU[cell_id])-1][1] 
    
    ## Location of output file:
    asc_header = 'ncols 1052\nnrows 365\nxllcorner 404710\nyllcorner 3199480\ncellsize 480\nNODATA_value 9999\n'
    
    ## Output file name: 
    with open(asc_out_file, mode='w') as outf:
        outf.write(asc_header)
        for row in LVMgrid:
            nc = 0           
            for col in row:
                gid_map = row[nc]
                if gid_map > 0 and gid_map<173898:
                    gid_val = cat[gid_map] # here mapping dict is a dictionary with grid cell ID as the key and some value for each key - if the ASC grid has a no data cell (-9999) then it will not have a key in the dictionary and is written out as no data (-9999)
                else:
                    gid_val = 9999
                if nc == 0:
                    rowout = '%s'  % gid_val
                else:
                    rowout = '%s %s' % (rowout,gid_val) #building the whole row
                nc += 1
            outf.write('%s\n' % rowout)
