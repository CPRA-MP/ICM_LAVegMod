import numpy as np
import os
import sys
import subprocess
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib import cm
import rasterio as rio
from rasterio.plot import plotting_extent

# read in simulation and plot settings from command arguments passed into this .py script
s = int(sys.argv[1])                 # s = '7'
g = int(sys.argv[2])                 # g = '503'
year = int(sys.argv[3])         # year = 2015
startyear = int(sys.argv[4])    # startyear = 2015
#ftype = sys.argv[5]

overwrite = 0       # overwrite=0 will not overwrite files if they exist, overwrite=1 will overwrite pre-existing files
spinup_years = 2
elapsedyear = year - startyear + 1
footnote = ''
asc_grid_rows = 371
ngrid = 173898



#############################
##      Setup folders      ##
#############################
print('Begin VCT TIF mapping for S%02d G%03d - yr %s' % (s,g,year))
if overwrite == 1:
    print(' - File overwrite flag setting turned ON (1) - output files will be overwritten if they already exist.')
else:
    print(' - File overwrite flag setting turned OFF (0) - output files will NOT be overwritten if they already exist.')      

print('\nsetting up folders')
veg_fol          = 'S%02d/G%03d/veg' % (s,g)
in_fol          = 'S%02d/G%03d/geomorph/input' % (s,g)
out_fol         = 'S%02d/G%03d/geomorph/output' % (s,g)
xyz_fol         = '%s/xyz' % out_fol 
tif_fol         = '%s/tif' % out_fol
png_fol         = '%s/png' % out_fol

for fol in [xyz_fol,tif_fol,png_fol]:
    try:
        if os.path.isdir(fol) == False:
            os.mkdir(fol)
    except:
        print('could not build %s' % fol)


if 1<2: #try:
    GRIDasc_pth     = '%s/veg_grid.asc' % (veg_fol)
    VCTasc_pth      = '%s/MP2023_S%02d_G%03d_C000_U00_V00_SLA_O_%02d_%02d_V_VCT.asc' % (veg_fol,s,g,elapsedyear,elapsedyear)
    VCTtif_pth      = '%s/MP2023_S%02d_G%03d_C000_U00_V00_SLA_N_%02d_%02d_V_VCT.tif' % (tif_fol,s,g,elapsedyear,elapsedyear)
    LTtif_pth       = '%s/MP2023_S%02d_G%03d_C000_U00_V00_SLA_O_%02d_%02d_W_lndtyp.tif' % (tif_fol,s,g,elapsedyear,elapsedyear)
    LVgridtif_pth   = '%s/MP2023_S00_G000_C000_U00_V00_SLA_I_00_00_W_grid30.tif' % (in_fol)
    png_pth         = '%s/MP2023_S%02d_G%03d_C000_U00_V00_SLA_O_%02d_%02d_W_VCT.png' % (png_fol,s,g,elapsedyear,elapsedyear)

    png_title = 'Vegetation Community Type - S%02d - G%03d - Year %02d' % (s,g,elapsedyear-spinup_years)
    
    if os.path.isfile(VCTasc_pth) == False:
        cmdstr = ['unzip', '%s.zip' % VCTasc_pth, '-d', '%s' % veg_fol]
        cmdout = subprocess.check_output(cmdstr).decode()
        
    if os.path.isfile(GRIDasc_pth) == False:
        cmdstr = ['unzip', '%s.zip' % GRIDasc_pth, '-d', '%s' % veg_fol]
        cmdout = subprocess.check_output(cmdstr).decode()
        
    if os.path.isfile(LTtif_pth) == False:
        cmdstr = ['unzip', '%s.zip' % LTtif_pth, '-d', '%s' % tif_fol]
        cmdout = subprocess.check_output(cmdstr).decode()
        
    if os.path.isfile(LVgridtif_pth) == False:
        cmdstr = ['unzip', '%s.zip' % LVgridtif_pth, '-d', '%s' % in_fol]
        cmdout = subprocess.check_output(cmdstr).decode()
      



    ################################################
    ##      Check for old files for overwrite     ##
    ################################################    
    if os.path.isfile(VCTtif_pth) == True:
        build_TIF = overwrite
        print('\nTIF raster file already exists - will use overwrite flag setting (%d) - FFIBS ' % (overwrite))
    else:
        build_TIF = 1
    
    if os.path.isfile(png_pth) == True:
        mapPNG = overwrite
        print('\nPNG image file already exists - will use overwrite flag setting(%d) - FFIBS ' % (overwrite))
    else:
        mapPNG = 1
    
    
    ############################################
    ##           Read in TIF rasters          ##
    ############################################
    if build_TIF == True:
        print('\ncombining FFIBS and landtype files')
    
        # open and read landtype TIF raster with rasterio  - then filter for NoData
        with rio.open(LTtif_pth) as open_tif:
            LTtif = open_tif.read(1)
            LTtrans = open_tif.transform        # save transformation settings to properly set raster resolution and coordinates when writing output TIF
        LTcrs = 'EPSG:26915'                    # EPSG code for UTM Zone 15N projection
    
        # read in ICM-LAVegMod grid raster        
        with rio.open(LVgridtif_pth) as open_tif:
            LVgridtif = open_tif.read(1)
        LVgridtif = LVgridtif[:-20]       # input grid TIF for grid has 20 additional rows of NoData at the bottom of the raster compared to the LT tif from Morph
        
        if LVgridtif.shape != LTtif.shape:
            print('Grid and landtype rasters are not of the same shape')
            sys.exit()
    
        # read in VCT score for each LAVegMod grid cell
        VCT_asci = np.genfromtxt(VCTasc_pth,skip_header=6,dtype='int',delimiter=' ')
        GRID_asci = np.genfromtxt(GRIDasc_pth,skip_header=6,dtype='int',delimiter=' ')
    
        if VCT_asci.shape != GRID_asci.shape:
            print('Grid and VCT ASCI rasters are not of the same shape')
            sys.exit()
    
        
        VCT_d = {}
        
        
        rows,cols = GRID_asci.shape
        for row in range(0,rows):
            for col in range(0,cols):                
                cID = int(GRID_asci[row][col])               # grid cell ID must be integer since grid30.tif is an integer raster
                VCT_d[cID] = VCT_asci[row][col]
    
    
        ############################
        ##  VCT_30 data values  ##
        ############################
        ## -9999        :   either NoData in GridID or Landtype rasters
        ## 1            :   Maidencane          
        ## 2            :   Three-Square        
        ## 3            :   Rosseau Cane       
        ## 4            :   Paspalum            
        ## 5            :   Wiregrass         
        ## 6            :   Bulltongue         
        ## 7            :   Needlerush        
        ## 8            :   Bulrush        
        ## 9            :   Brackish Mix      
        ## 10           :   Oystergrass    
        ## 11           :   Saltgrass  
        ## 12           :   majority non-SOM bottomland hardwood forest + swamp forest
        ## 13           :   majority barrier island 
        ## 14           :   No SOM species found (i.e. all water) 
        ## 200          :   pixel is water
        ## 300          :   pixel is unvegetated mudflat
        ## 400          :   pixel is upland/developed land
        ## 500          :   pixel is flotant
        ## +9999        :   land pixel is of vegetated land landtype, but there is NoData in the FFIBS output data from LAVegMod
        ############################
        # map VCT score to 30-m grid when the Landtype = 1, otherwise remap Landtype to higher values
        
        VCT_30 = LVgridtif*0.0     # build zero array in equal dimensions of grid30 tif
        rows,cols = LVgridtif.shape
        for r in range(0,rows):
            if r in range(0,rows,int(rows/10)):
                print('processed %s of %s rows...' % (r,rows))
            for c in range(0,cols):
                gridID = LVgridtif[r][c]
                LT   = LTtif[r][c]
                if gridID > 0:
                    if LT > 0:
                        if LT == 1:
                            VCT = VCT_d[gridID]
                            if VCT > 0:            # FFIBS is NoData, but it is on vegetated landtype - set to new NoData value of +9999
                                vct2map = VCT
                            else:
                                vct2map = 9999
                        else:
                            vct2map = LT*100
                    else:
                        vct2map = -9999              # landtype is NoData, keep NoData value as -9999
                else:
                    vct2map = -9999                  # gridID is NoData, keep NoData value as -9999
                VCT_30[r][c] = vct2map
    
        # update grids to account for NoData
        VCT_30_na = VCT_30 #np.ma.masked_where(VCT_30 < -9990 ,VCT_30,copy=True)
        
        ##################################################################
        ##       Export FFIBS landtype combination raster as a TIF      ##
        ##################################################################
        print('\nSaving %s' % VCTtif_pth)
        with rio.open(VCTtif_pth,'w',dtype=rio.int16,count=1,driver='GTiff',height=rows,width=cols,crs=LTcrs,transform=LTtrans) as dest:
            dest.write(VCT_30.astype(rio.int16),1)
#except:
#    print('failed to convert and/or map %s' % ftype)
