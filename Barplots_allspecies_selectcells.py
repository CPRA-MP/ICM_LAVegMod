####builds bar plots for one cell of LAVegMod for all 50 years

import matplotlib.pyplot as plt
import numpy as np


outdir = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_LAVegMod\flotant_cells_090920'
input_file = r'C:\Users\madel\Coastal Hydro Dropbox\Madeline Foster-Martinez\MRF_LAVegMod\MP2023_G028_S03_flotant090920.csv'

data_in = np.genfromtxt(input_file,skip_header=1,delimiter=',',dtype='str')

print(' - reading in data from %s' % input_file)

# columns in the input file:
#   Project number
#   Scenario
#   Year
#   Coverage code
#   Cell ID
#   Value (fraction of cell covered [0-1])

# Read input data tables into nested dictionaries of timeseries arrays

years = []
PSEHts = {}


for row in data_in:
    prj = row[1] #project ID, e.g. G028
    s = row[2] #scenario, e.g. 03
    yr = int(row[3]) #year
    cc = row[4] #coverage code, e.g. 'COES' or 'WATER'
    cc = cc.lstrip(' ') 
    cid = int(row[5]) #cell ID
    area = float(row[6]) #coverage value of that type [0-1]
     
    # the year and timeseries arrays will be ordered by the order of the input file...if not chronological in the input file, these arrays will be out of order
    # can fix this by adding one more layer to the dictionary that uses year as one more key...but this will complicate plotting functions below  
    # ignore for now since all files are chronologically structured
    
    # append year to timeseries dt array
    if yr not in years:
        years.append(yr)
    
    if prj not in PSEHts.keys():
        PSEHts[prj] = {}
        
    if s not in PSEHts[prj].keys():
        PSEHts[prj][s] = {}
        
    if cid not in PSEHts[prj][s].keys():
        PSEHts[prj][s][cid] = {}
 
    if cc not in PSEHts[prj][s][cid].keys() and cc not in ['FIBS', 'pL_BF', 'pL_SF', 'pL_FM', 'pL_IM', 'pL_BM', 'pL_SM','DEAD_Flt','CELLID']:                  
        PSEHts[prj][s][cid][cc] = []   
    if cc not in ['FIBS', 'pL_BF', 'pL_SF', 'pL_FM', 'pL_IM', 'pL_BM', 'pL_SM','DEAD_Flt','CELLID']:
        PSEHts[prj][s][cid][cc].append(area)
      
sp_names = list(PSEHts[prj][s][cid].keys())

for S in list(PSEHts[prj].keys()): #scenarios
    for P in list(PSEHts.keys()): #projects
        for Cl in list(PSEHts[P][S].keys()): # cell ids 
            print(' - plotting %s %s %s' % (S,P,Cl) )
            bars = []
            legtxt = []
            col = []
            bot = []
            for sp in range(0,len(sp_names)): 
                bars.append(PSEHts[P][S][Cl][sp_names[sp]])
                legtxt.append(sp_names[sp])
                if sp_names[sp] == 'WATER':
                    col.append([0,0,1,1]) # water = blue
                elif sp_names[sp] == 'NOTMOD':
                    col.append([192/255,192/255,192/255,1]) # not mod = gray
                elif sp_names[sp] == 'BAREGRND_OLD':
                    col.append([153/255,102/255,51/255,1]) # bareground old = brown
                elif sp_names[sp] == 'BAREGRND_NEW':
                    col.append([204/255,153/255,0,1]) # bareground new = brown
                elif sp_names[sp] == 'SAV':
                    col.append([153/255,204/255,255/255,1]) # SAV = light blue
                elif sp_names[sp] == 'PAHE2_Flt':
                    col.append([204/255,0/255,204/255,1]) # 
                elif sp_names[sp] == 'ELBA2_Flt':
                    col.append([1,51/255,1,1]) # 
                elif sp_names[sp] == 'BAREGRND_Flt':
                    col.append([255/255,204/255,255/255,1]) # 
                elif sp_names[sp] == 'NYAQ2':
                    col.append([76/255,0/255,153/255,1]) # 
                elif sp_names[sp] == 'SANI':
                    col.append([127/255,0/255,255/255,1]) # 
                elif sp_names[sp] == 'TADI2':
                    col.append([178/255,102/255,255/255,1]) # 
                elif sp_names[sp] == 'QULA3':
                    col.append([0/255,102/255,51/255,1]) # 
                elif sp_names[sp] == 'QULE':
                    col.append([0/255,153/255,76/255,1]) # 
                elif sp_names[sp] == 'QUNI':
                    col.append([0/255,204/255,102/255,1]) # 
                elif sp_names[sp] == 'QUTE':
                    col.append([0/255,255/255,128/255,1]) # 
                elif sp_names[sp] == 'QUVI':
                    col.append([51/255,255/255,153/255,1]) # 
                elif sp_names[sp] == 'ULAM':
                    col.append([153/255,205/255,154/255,1]) 
                elif sp_names[sp] == 'COES':
                    col.append([75/255,153/255,0/255,1]) # 
                elif sp_names[sp] == 'MOCE2':
                    col.append([102/255,204/255,0/255,1]) # 
                elif sp_names[sp] == 'PAHE2':
                    col.append([128/255,255/255,0/255,1]) # 
                elif sp_names[sp] == 'SALA2':
                    col.append([178/255,255/255,102/255,1]) # 
                elif sp_names[sp] == 'ZIMI':
                    col.append([229/255,255/255,204/255,1]) #                      
                elif sp_names[sp] == 'CLMA10':
                    col.append([204/255,204/255,0/255,1]) # 
                elif sp_names[sp] == 'ELCE':
                    col.append([255/255,255/255,0/255,1]) # 
                elif sp_names[sp] == 'POPU5':
                    col.append([255/255,255/255,153/255,1]) # 
                elif sp_names[sp] == 'SALA':
                    col.append([255/255,255/255,102/255,1]) # 
                elif sp_names[sp] == 'IVFR':
                    col.append([241/255,247/255,40/255,1]) # 
                elif sp_names[sp] == 'PAVA':
                    col.append([255/255,255/255,204/255,1]) # 
                elif sp_names[sp] == 'PHAU7':
                    col.append([243/255,249/255,83/255,1]) # 
                elif sp_names[sp] == 'SCCA11':
                    col.append([226/255,232/255,54/255,1]) # 
                elif sp_names[sp] == 'TYDO':
                    col.append([220/255,228/255,8/255,1]) #                     
                elif sp_names[sp] == 'SCAM6':
                    col.append([255/255,128/255,0/255,1])                    
                elif sp_names[sp] == 'SCRO5':
                    col.append([255/255,153/255,51/255,1])                                  
                elif sp_names[sp] == 'SPPA':
                    col.append([255/255,178/255,102/255,1])                    
                elif sp_names[sp] == 'SPCY':
                    col.append([255/255,229/255,204/255,1])
                elif sp_names[sp] == 'DISP':
                    col.append([153/255,0,0,1])
                elif sp_names[sp] == 'JURO':
                    col.append([204/255,0,0,1])
                elif sp_names[sp] == 'AVGE':
                    col.append([255/255,0,0,1])
                elif sp_names[sp] == 'SPAL':
                    col.append([255/255,102/255,102/255,1]) 
                elif sp_names[sp] == 'BAHABI':
                    col.append([102/255,0/255,51/255,1]) 
                elif sp_names[sp] == 'DISPBI':
                    col.append([153/255,0/255,76/255,1]) 
                elif sp_names[sp] == 'PAAM2':
                    col.append([204/255,0/255,102/255,1]) 
                elif sp_names[sp] == 'SOSE':
                    col.append([255/255,0/255,127/255,1]) 
                elif sp_names[sp] == 'SPPABI':
                    col.append([255/255,51/255,153/255,1]) 
                elif sp_names[sp] == 'SPVI3':
                    col.append([255/255,102/255,178/255,1]) 
                elif sp_names[sp] == 'STHE9':
                    col.append([255/255,153/255,204/255,1]) 
                elif sp_names[sp] == 'UNPA':
                    col.append([255/255,204/255,229/255,1]) 
                B = np.array(bars) #shape of B is (47,50), row for each coverage and column for each year
                bot.append(list(np.sum(B,0)))

            BB = np.where(B > 0, 1, B) #if the coverage value is greater than 0, set it to 1 
            BBs = sum(BB) #this sum gives the number of unique coverages every year
           
            wid = 0.75  #set column width
            m=0
            plt.figure(num=None, figsize=[7.5,4],dpi=300)
            if sum(bars[31])>0:
                i = 31
                pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                m+=1
                bot1 = np.array(bars[i])
            if sum(bars[28])>0:
                i = 28
                if m == 0:
                    pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                    m+=1
                    bot1 = np.array(bars[i])
                elif m==1:
                    pbar1 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot1)
                    m+=1
                    bot2 = np.array(bot1)+np.array(bars[i]) 
            if sum(bars[29])>0:
                i = 29
                if m == 0:
                    pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                    m+=1
                    bot1 = np.array(bars[i])
                elif m==1:
                    pbar1 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot1)
                    m+=1
                    bot2 = np.array(bot1)+np.array(bars[i])                   
                elif m==2:
                    pbar2 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot2)
                    m+=1
                    bot3 = np.array(bot2)+np.array(bars[i])    
            if sum(bars[30])>0:
                i = 30
                if m == 0:
                    pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                    m+=1
                    bot1 = np.array(bars[i])
                elif m==1:
                    pbar1 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot1)
                    m+=1
                    bot2 = np.array(bot1)+np.array(bars[i])                   
                elif m==2:
                    pbar2 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot2)
                    m+=1
                    bot3 = np.array(bot2)+np.array(bars[i])                   
                elif m==3:
                    pbar3 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot3)
                    m+=1
                    bot4 = np.array(bot3)+np.array(bars[i])                 
            if sum(bars[32])>0:
                i = 32
                if m == 0:
                    pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                    m+=1
                    bot1 = np.array(bars[i])
                elif m==1:
                    pbar1 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot1)
                    m+=1
                    bot2 = np.array(bot1)+np.array(bars[i])                   
                elif m==2:
                    pbar2 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot2)
                    m+=1
                    bot3 = np.array(bot2)+np.array(bars[i])                   
                elif m==3:
                    pbar3 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot3)
                    m+=1
                    bot4 = np.array(bot3)+np.array(bars[i])                     
                elif m==4:
                    pbar4 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot4)
                    m+=1
                    bot5 = np.array(bot4)+np.array(bars[i])                
            for i in range(0,len(sp_names)):    
                if i < 28 or i > 32:              
                    if sum(bars[i])>0:
                        if m == 0:
                            pbar0 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i])
                            m+=1
                            bot1 = np.array(bars[i])
                        elif m==1:
                            pbar1 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot1)
                            m+=1
                            bot2 = np.array(bot1)+np.array(bars[i])                   
                        elif m==2:
                            pbar2 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot2)
                            m+=1
                            bot3 = np.array(bot2)+np.array(bars[i])                   
                        elif m==3:
                            pbar3 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot3)
                            m+=1
                            bot4 = np.array(bot3)+np.array(bars[i])                     
                        elif m==4:
                            pbar4 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot4)
                            m+=1
                            bot5 = np.array(bot4)+np.array(bars[i])
                        elif m==5:
                            pbar5 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot5)
                            m+=1
                            bot6 = np.array(bot5)+np.array(bars[i])
                        elif m==6:
                            pbar6 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot6)
                            m+=1
                            bot7 = np.array(bot6)+np.array(bars[i])
                        elif m==7:
                            pbar7 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot7)
                            m+=1
                            bot8 = np.array(bot7)+np.array(bars[i])
                        elif m==8:
                            pbar8 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot8)
                            m+=1
                            bot9 = np.array(bot8)+np.array(bars[i])
                        elif m==9:
                            pbar9 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot9)
                            m+=1
                            bot10 = np.array(bot9)+np.array(bars[i])
                        elif m==10:
                            pbar10 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot10)
                            m+=1
                            bot11 = np.array(bot10)+np.array(bars[i])
                        elif m==11:
                            pbar11 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot11)
                            m+=1
                            bot12 = np.array(bot11)+np.array(bars[i])
                        elif m==12:
                            pbar12 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot12)
                            m+=1
                            bot13 = np.array(bot12)+np.array(bars[i])
                        elif m==13:
                            pbar13 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot13)
                            m+=1
                            bot14 = np.array(bot13)+np.array(bars[i])
                        elif m==14:
                            pbar14 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot14)
                            m+=1
                            bot15 = np.array(bot14)+np.array(bars[i])
                        elif m==15:
                            pbar15 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot15)
                            m+=1
                            bot16 = np.array(bot15)+np.array(bars[i])
                        elif m==16:
                            pbar16 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot16)
                            m+=1
                            bot17 = np.array(bot16)+np.array(bars[i])
                        elif m==17:
                            pbar17 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot17)
                            m+=1
                            bot18 = np.array(bot17)+np.array(bars[i])
                        elif m==18:
                            pbar18 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot18)
                            m+=1
                            bot19 = np.array(bot18)+np.array(bars[i])
                        elif m==19:
                            pbar19 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot19)
                            m+=1
                            bot20 = np.array(bot19)+np.array(bars[i])
                        elif m==20:
                            pbar20 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot20)
                            m+=1
                            bot21 = np.array(bot20)+np.array(bars[i])
                        elif m==21:
                            pbar21 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot21)
                            m+=1
                            bot22 = np.array(bot21)+np.array(bars[i])
                        elif m==22:
                            pbar22 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot22)
                            m+=1
                            bot23 = np.array(bot22)+np.array(bars[i])
                           
                        elif m==23:
                            pbar23 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot23)
                            m+=1
                            bot24 = np.array(bot23)+np.array(bars[i])
                        elif m==24:
                            print(i)
                            print(len(legtxt))
                            print(len(col))
                            print(len(bars))
                            pbar24 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot24)
                            m+=1
                            bot25 = np.array(bot24)+np.array(bars[i])
                        elif m==25:
                            pbar25 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot25)
                            m+=1
                            bot26 = np.array(bot25)+np.array(bars[i])
                        elif m==26:
                            pbar26 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot26)
                            m+=1
                            bot27 = np.array(bot26)+np.array(bars[i])
                        elif m==27:
                            pbar27 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot27)
                            m+=1
                            bot28 = np.array(bot27)+np.array(bars[i])
                        elif m==28:
                            pbar28 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot28)
                            m+=1
                            bot29 = np.array(bot28)+np.array(bars[i])
                        elif m==29:
                            pbar29 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot29)
                            m+=1
                            bot30 = np.array(bot29)+np.array(bars[i])
                        elif m==30:
                            pbar30 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot30)
                            m+=1
                            bot31 = np.array(bot30)+np.array(bars[i])
                        elif m==31:
                            pbar31 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot31)
                            m+=1
                            bot32 = np.array(bot31)+np.array(bars[i])
                        elif m==32:
                            pbar32 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot32)
                            m+=1
                            bot33 = np.array(bot32)+np.array(bars[i])
                        elif m==33:
                            pbar33 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot33)
                            m+=1
                            bot34 = np.array(bot33)+np.array(bars[i])
                        elif m==34:
                            pbar34 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot34)
                            m+=1
                            bot35 = np.array(bot34)+np.array(bars[i])
                        elif m==35:
                            pbar35 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot35)
                            m+=1
                            bot36 = np.array(bot35)+np.array(bars[i])
                        elif m==36:
                            pbar36 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot36)
                            m+=1
                            bot37 = np.array(bot36)+np.array(bars[i])
                        elif m==37:
                            pbar37 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot37)
                            m+=1
                            bot38 = np.array(bot37)+np.array(bars[i])
                        elif m==38:
                            pbar38 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot38)
                            m+=1
                            bot39 = np.array(bot38)+np.array(bars[i])
                        elif m==39:
                            pbar39 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot39)
                            m+=1
                            bot40 = np.array(bot39)+np.array(bars[i])
                        elif m==40:
                            pbar40 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot40)
                            m+=1
                            bot41 = np.array(bot40)+np.array(bars[i])  
                        elif m==41:
                            pbar41 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot41)
                            m+=1
                            bot42 = np.array(bot41)+np.array(bars[i])     
                        elif m==42:
                            pbar42 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot42)
                            m+=1
                            bot43 = np.array(bot42)+np.array(bars[i])
                        elif m==43:
                            pbar43 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot43)
                            m+=1
                            bot44 = np.array(bot43)+np.array(bars[i])
                        elif m==44:
                            pbar44 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot44)
                            m+=1
                            bot45 = np.array(bot44)+np.array(bars[i])
                        elif m==45:
                            pbar45 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot45)
                            m+=1
                            bot46 = np.array(bot45)+np.array(bars[i])
                        elif m==46:
                            pbar46 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot46)
                            m+=1
                            bot47 = np.array(bot46)+np.array(bars[i])                     
                        elif m==47:
                            pbar47 = plt.bar(years,bars[i],width=wid,label=legtxt[i],color=col[i],bottom=bot47)
                            m+=1
                            bot48 = np.array(bot47)+np.array(bars[i])     
             
            #set axes and chart titles
            x_txt = 'Year'
            y_txt = 'Coverage Type'
            title_txt = 'Coverage:  %s - %s - %s' % (S,P,Cl)
            #print the number of unique coverages at the top of each bar
            for r in range(0,50,2):
                plt.text(r+0.55,0.95,str(int(BBs[r])),fontsize=6)
            for r in range(1,51,2):
                plt.text(r+0.55,0.97,str(int(BBs[r])),fontsize=6) #vertically stagger every other one
            
            plt.legend(ncol=3,fontsize='small',bbox_to_anchor=(0.5,-0.2), loc='upper center', borderaxespad=0.00)
            plt.ylabel(y_txt)
            plt.xlabel(x_txt)
            plt.title(title_txt)
           
            plt.savefig(r'%s\%s_%s_coverage_timeseries_%s.png' % (outdir,S,P,Cl), bbox_inches='tight') 
            plt.show()
            plt.close()
