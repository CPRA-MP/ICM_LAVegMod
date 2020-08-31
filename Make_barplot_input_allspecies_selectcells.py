#Generates a csv properly formatted for use in 'Barplots_allspecies_selectcells.py' 
import numpy as np
import os
import pandas

G = 'G028' ##INPUT
S = 'S03' ##INPUT
years = range(1,51) 

#specify all directories needed: 
path_cellID = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_LAVegMod\Rus_1_2_July2020'
path_asc = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_LAVegMod\G028\veg'
outdir = path_cellID

#read the cell IDs of where we want output 
os.chdir(path_cellID)
cell_ID =  np.genfromtxt('sample.csv',skip_header=1, delimiter=',', dtype='int') ###INPUT csv containing the cell IDs for the cells of interest 
cell_ID = cell_ID[:,0] #first column contains the cell IDs needed

#move to the folder with the LAVegMod asc+ output file
os.chdir(path_asc)

#build a dataframe with the following columns:
df = pandas.DataFrame(columns=['pro_no','S','year','coverage_code','cell_ID','value'])

#populate the dataframe:
m = 0
for Y in years: 
    print('On year '+ str(Y))
    if Y <10:
        #read in the LAVegMod output coverage values
        LVMout = np.genfromtxt('MP2023_'+ str(S) +'_'+ str(G)+'_C000_U00_V00_SLA_O_0'+str(Y)+'_0'+str(Y)+'_V_vegty.asc+',skip_header=622, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623
        #read in the coverage names
        sp_names = np.genfromtxt('MP2023_'+ str(S) +'_'+ str(G)+'_C000_U00_V00_SLA_O_0'+str(Y)+'_0'+str(Y)+'_V_vegty.asc+',skip_header=621, skip_footer = 187553, delimiter=',', dtype='str') 
    else: 
        LVMout = np.genfromtxt('MP2023_'+ str(S) +'_'+ str(G)+'_C000_U00_V00_SLA_O_'+str(Y)+'_'+str(Y)+'_V_vegty.asc+',skip_header=622, delimiter=',', dtype='float') # skip the top portion of the asc+ files and only read in the column data starting on line 623    
        sp_names = np.genfromtxt('MP2023_'+ str(S) +'_'+ str(G)+'_C000_U00_V00_SLA_O_'+str(Y)+'_'+str(Y)+'_V_vegty.asc+',skip_header=621, skip_footer = 187553, delimiter=',', dtype='str')
    for cell in cell_ID:
        coverages = LVMout[cell-1]
        for i in range(1,len(sp_names)):
            df.loc[m] = pandas.Series({'pro_no':G,'S':S,'year':Y,'coverage_code':sp_names[i],'cell_ID':cell,'value':coverages[i]})
            m += 1
 
#save the dataframe as a csv
out_name = ('\MP2023_'+G+'_'+S+'_'+'barplot_input')
df.to_csv(outdir+out_name+'.csv')