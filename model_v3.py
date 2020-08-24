#!/usr/bin/env python

##\file
##\brief Top level file for the ecological modeling.
##\detail This file describes the model proper and the major
# model-specific coding elements used to build this model.
##\section Includes Included Libraries
# The file makes use of the following libraries
# - Standard Python Libraries
#   - copy
#   - exceptions
#   - itertools
#   - re
#   - StringIO
#   - sys
#   - time
# - Third Party Libraries
#   - pandas version 0.13.1
#   - xlrd
# - Model defined libraries
#   - config
#   - event
#   - function
#   - landscape

# Diagnostic modules
#from memory_profiler import profile
#from pympler.asizeof import asizeof

# STD Python modules
import collections
import copy
from builtins import Exception as exceptions
#import exceptions
import glob
import itertools
import math
import os
import re
import sys
import time

# Third party modules
import numpy
import pandas
import xlrd

# This model's modules
import config
import event
import function
import landscape
import plantingmodel


##\class SpeciesModel
##\brief Base class for all species models.
##\detail This is a generic base class defining the functionality
# of a single species model for a single location in the landscape.
# All species models should inherit
# from this class to make sure they get the basic elements.
##\role{Ecology/Machinery}
class SpeciesModel(object):
    ##\brief A reference to the model Params object
    params   = None
    ##\brief A reference to the model DispersalModel object
    dspModel = None

    ## Class constructor.
    #\param [in] self      The object's reference to itself.
    #\param [in] index     The species index.
    #\param [in] abr       The species USDA alpha-numeric vegetation code.
    #\param [in] name      The full species name.
    #\param [in] modelType The type of a model.
    #\param [in] habitat   The habitat provided by the species (salinity regime)
    #\param [in] dspClass  The dispersal class of the species (0,1,2,3=high)
    #\param [in] ffibsScore The FFIBS (forested, fresh, intermediate, brackish, saline) score for that species
    #\param [in] cover     The fraction of model cell that is covered by a species.
    #\param [in] loc       The x,y location of a cell. Given as the row,column index with in the model grid.
    def __init__(self, index=0, name='', abr='', modelType='', habitat='', dspClass=0, ffibsScore=0, cover=0, loc=0):
        ##\brief Species index
        ##\details An integer that is a unique numerical species index.
        ## - Type: integer
        ## - Range: [0, inf)
        #self.index     = index

        ##\brief The long name for a species
        ##\details A string used to make output pretty
        ## - Type: string
        #self.name      = name

        ##\brief Species code
        ##\details A string abbreviation for the species. The value is the species' USDA alpha-numeric vegetation code.
        # This is an important member value in the model. It is used as the key in the dictionary that stores the list
        # of species in a patch.
        ## - Type: string
        self.abr       = abr

        ##\brief The model type
        ##\details A string containing the model type.
        ##\n
        ## This member should one of the following values (as a string):\n
        # - BarrierIslandModel
        # - EmergentWetlandModel
        # - BottomlandHardwoodForestModel
        # - SAVModel
        # - NullModel_Coverage
        # - NullModel
        ##\n
        self.modelType = modelType

        ##\brief The habitat type
        ##\details A string containing the habitat type.
        ##\n
        ## This member should one of the following values (as a string):\n
        # - Fresh
        # - Intermediate
        # - Brackish
        # - Saline
        # - NA
        ##\n
        self.habitat = habitat

        ##\brief The dispersal class of a species
        ##\details Dispersal classes describe far a species can establish, 1 is from one cell away, 2 is from two cells away,  3 is from anywhere, and 0 is for null
        ## - Type: float
        ## - Range: [0, 3]
        self.dspClass     = dspClass

        ##\brief The FIBS score of a species corresponds to its salinity regime: Fresh = 0.25, Intermediate = 1.5, Brackish = 7.15 or 11.5, Saline = 17.5 or 24, null coverages are nan
        ##\details FIBS score
        ## - Type: float
        ## - Range: [0.25, 24]
        self.ffibsScore     = ffibsScore

        ##\brief The area occupied by a species
        ##\details This is the fraction of a patch covered by a species.
        ## - Type: float
        ## - Range: [0, inf)
        self.cover     = cover

        ##\brief Species Location
        ##\details This is the location of a species. The location is stored as a python tuple
        # that contains the row and column index of the landscape where the species is present
        ## - Type: tuple
        self.loc       = loc

    ##\brief Object configuration.
    #\param self The object's reference to itself.
    #\param index The species index.
    #\param abr The species USDA alpha-numeric vegetation code.
    #\param name The full species name.
    #\param modelType The type of a model.
    #\param habitat The habitat provided by the species (salinity regime).
    #\param dspClass The dispersal class the species belongs (low, med, high).
    #\param ffibsScore The FFIBS score of the species.
    #\param cover The fraction of model cell that is covered by a species.
    #\param loc The x,y location of a cell. Given as the row,column index with in the model grid.
    ##\detail This function is used to configure an object after it has been created.
    # this is used SpeciesModelList when it is creating the list of species from
    # the species growth/senescence tables and by PatchModel when creating the
    # individual species model for each location.
    def config(self, index=0, name='', abr='', modelType='', habitat='', dspClass=0, ffibsScore=0, cover=0, loc=0):
        #self.index     = index
        #self.name      = name
        self.abr       = abr
        self.modelType = modelType
        self.habitat   = habitat
        self.dspClass = dspClass
        self.ffibsScore = ffibsScore
        self.cover     = cover
        #self.loc       = loc

    ##\brief Computes the probability of senescence.
    ##\detail This is a generic function that returns the probability of plants senesencing.
    # This function should be redefined by each class that inherits form SpeciesModel.
    def senescence(self, loc):
        return 0.0

    ##\brief Computes the probability of establishment
    ##\details This is a generic function that returns the probability of plants growing.
    # This function should be redefined by each class that inherits from SpeciesModel.
    def growth(self, loc):
        return 0.0

    ##\brief Computes the probability of establishment for high dispersal species
    ##\details This is a generic function that returns the probability of plants growing.
    # This function should be redefined by each class that inherits from SpeciesModel.
    def spread(self, loc):
        return 0.0


    def __float__(self):
        return -3.0
        #return self.cover

    #def __str__(self):
    #    ret = ''
    #    #ret =  str(self.index) + ', '
    #    #ret += self.name            + ', '
    #    ret += self.abr             + ', '
    #    ret += self.modelType       + ', '
    #    ret += str(self.cover)      + ', '
    #    ret += str(self.loc)
    #    return ret

##\class NullModel
##\brief A place holder.
##\role{Machinery}
##\detail A place holder for items in the output that are not coverage values (e.g. FFIBS).
## Having this class makes the rest of the model code simpler.
class NullModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat = '', dspClass=0, ffibsScore=0, cover=0, loc=0):
        SpeciesModel.__init__(self, index, name, abr, 'NullModel', habitat, dspClass, ffibsScore, cover, loc)

    ##\brief
    ##\detail
    def senescence(self, loc):
        return 0.0

    ##\brief
    ##\detail
    def growth(self, loc):
        return 0.0

    def spread(self,loc):
        return 0.0
    
##\class NullModel_Coverage
##\brief A place holder.
##\role{Machinery}
##\detail A place holder for coverage types that
# do not really do anything but do need to be counted when adding up coverage (e.g. water).
# Having this class makes the rest of the model code simpler.
class NullModel_Coverage(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat = '', dspClass=0, ffibsScore=0, cover=0, loc=0):
        SpeciesModel.__init__(self, index, name, abr, 'NullModel_Coverage', habitat, dspClass, ffibsScore, cover, loc)

    ##\brief
    ##\detail
    def senescence(self, loc):
        return 0.0

    ##\brief
    ##\detail
    def growth(self, loc):
        return 0.0

    def spread(self,loc):
        return 0.0

##\class BottomlandHardwoodForestModel

##\brief This class handles the ecology for the bottomland hardwood forest species.
##\role{Ecology}
##\detail bottomland hardwood forest types change state based on the
# height of the habitat above mean water surface.
class BottomlandHardwoodForestModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat='', dspClass=0, ffibsScore=0, cover=0, loc=0, Pdata=None, Ddata=None):
        SpeciesModel.__init__(self, index, name, abr, 'BottomlandHardwoodForestModel', habitat, dspClass, ffibsScore, cover, loc)
        self.P = function.Function(Pdata['elvValue'], Pdata['rate'])
        self.D = function.Function(Ddata['elvValue'], Ddata['rate'])

    ## Senescence function
    # This function determines the probability that a bottomland hardwood forest species
    # will experience senescence. The probability of senescence is a function
    # of the elevation of the local habitat above the mean water surface. The
    # exact functional relationship between height above water and the probability
    # of senescence is determined by the senescence table for each species.
    # The senescence tables are defined in a model input file.
    def senescence(self, loc):
        try:
            elv = SpeciesModel.params.heightAboveWater[loc]
            sal = SpeciesModel.params.meanSal[loc]
            bi  = SpeciesModel.params.biEstCond[loc]

            if bi or sal > 1.0:
                return 1.0

            return self.D[elv]

        except RuntimeError as error:
            msg = 'BottomlandHardwoodForestModel.senescence(): Error: location out of range. Additional info follows\n'
            msg += str(error)
            raise RuntimeError(msg)

    ## Growth function
    # This function determines the probability that a bottomland hardwood forest species will become
    # established at the current location. The probability of a bottomland hardwood forest species becoming
    # established is depended on three major factors. First, the salinity must be below 1.0 ppt.
    # Second, there must be a period of 14 days with no flooding followed by a period of 14 with
    # water depths below 14 cm. Finally, the height of the habitat above mean water level determines
    # the final probability.


    # I'm not entirely happy with the current state of this function because the computation
    # of the basic establishment conditions (14 day no flood, 14 days water depth < 14 cm) is
    # handled external to the model. The problem is that the current approach divides
    # the responsibilities for representing the ecology of upland forest species (now bottomland hardwood forest).
    # An external program determines the establishment conditions while the model proper handles
    # the final computation of the probability of establishment. I would like to fix this.
    def growth(self, loc):
        try:
            elv = SpeciesModel.params.heightAboveWater[loc]
            sal = SpeciesModel.params.meanSal[loc]
            est = SpeciesModel.params.treeEstCond[loc]
            bi  = SpeciesModel.params.biEstCond[loc]
            dsp = SpeciesModel.dspModel[loc][self.abr]

            if bi or sal > 1.0:
                return 0

            return self.P[elv] * dsp * est
        except RuntimeError as error:
            msg = 'BottomlandHardwoodForestModel.growth(): Error: location out of range. Additional info follows\n'
            msg += str(error)
            raise RuntimeError(msg)

    def spread(self,loc):
        if self.dspClass == 3:
            try:
                elv = SpeciesModel.params.heightAboveWater[loc]
                sal = SpeciesModel.params.meanSal[loc]
                est = SpeciesModel.params.treeEstCond[loc]
                bi  = SpeciesModel.params.biEstCond[loc]

                if bi or sal > 1.0:
                    return 0

                return self.P[elv]*est


            except RuntimeError as error:
                msg = 'BottomlandHardwoodForestModel.growth(): Error: location out of range. Additional info follows\n'
                msg += str(error)
                raise RuntimeError(msg)
        else:
            return 0.0

##\class EmergentWetlandModel
##\brief This class handles the ecology for the emergent wetland species.
##\role{Ecology}
##\detail

class EmergentWetlandModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat = '', dspClass=0, ffibsScore=0, cover=0, loc=0, Pdata=None, Ddata=None):
        SpeciesModel.__init__(self, index, name, abr, 'EmergentWetlandModel', habitat, dspClass, ffibsScore, cover, loc)
        self.P     = function.Function2DFast(Pdata['waValue'], Pdata['salValue'], Pdata['rate'])
        self.D     = function.Function2DFast(Ddata['waValue'], Ddata['salValue'], Ddata['rate'])

    ## Senescence function
    # This function determines the probability that an emergent wetland species
    # will experience senescence. The probability of senescence is a function
    # of variation in water stage height (computed as the standard deviation of stage)
    # and salinity. The exact relationship between the probability of senescence and
    # the environmental factors is describe the by scenescence table for each species.
    # The senescence tables are defined in a model input file.
    def senescence(self, loc):
        waveAmp = SpeciesModel.params.waveAmp[loc]
        meanSal = SpeciesModel.params.meanSal[loc]
        elv     = SpeciesModel.params.heightAboveWater[loc]
        bi      = SpeciesModel.params.biEstCond[loc]

        #if bi or elv > SpeciesModel.params.elevationThreshold:
            #return 1.0

        if bi: return 1.0

        return self.D[waveAmp,meanSal]

   ## Growth function
    # This function determines the probability that an emergent wetland species will become
    # established at the current location. The probability of an wetland species becoming
    # established is depended on
    # three major factors. First, the salinity must be below 1.0 ppt.
    # Second, there must be a period of 14 days with no flooding followed by a period of 14 with
    # water depths below 14 cm. Finally, the height of the habitat above mean water level determines
    # the final probability.
    def growth(self, loc):
        waveAmp = SpeciesModel.params.waveAmp[loc]
        meanSal = SpeciesModel.params.meanSal[loc]
        elv     = SpeciesModel.params.heightAboveWater[loc]
        bi      = SpeciesModel.params.biEstCond[loc]
        dsp     = SpeciesModel.dspModel[loc][self.abr]


#No_elev_threshold#        if bi or elv > SpeciesModel.params.elevationThreshold:
        #if bi > SpeciesModel.params.elevationThreshold:

        if bi: return 0.0

        return self.P[waveAmp, meanSal] * dsp

    def spread(self,loc):
        if self.dspClass == 3:
            waveAmp = SpeciesModel.params.waveAmp[loc]
            meanSal = SpeciesModel.params.meanSal[loc]
            elv     = SpeciesModel.params.heightAboveWater[loc]
            bi      = SpeciesModel.params.biEstCond[loc]

            if bi: return 0.0

            return self.P[waveAmp, meanSal]
        else:
            return 0.0


class SwampForestModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat = '', dspClass=0, ffibsScore=0, cover=0, loc=0, Pdata=None, Ddata=None):
        SpeciesModel.__init__(self, index, name, abr, 'SwampForestModel', habitat, dspClass, ffibsScore, cover, loc)
        self.P     = function.Function2DFast(Pdata['waValue'], Pdata['salValue'], Pdata['rate'])
        self.D     = function.Function2DFast(Ddata['waValue'], Ddata['salValue'], Ddata['rate'])

    ## Senescence function
    # This function determines the probability that an emergent wetland species
    # will experience senescence. The probability of senescence is a function
    # of variation in water stage height (computed as the standard deviation of stage)
    # and salinity. The exact relationship between the probability of senescence and
    # the environmental factors is describe the by scenescence table for each species.
    # The senescence tables are defined in a model input file.
    def senescence(self, loc):
        waveAmp = SpeciesModel.params.waveAmp[loc]
        meanSal = SpeciesModel.params.meanSal[loc]
        bi      = SpeciesModel.params.biEstCond[loc]

        if bi: return 1.0

        return self.D[waveAmp,meanSal]

   ## Growth function
    # This function determines the probability that an emergent wetland species will become
    # established at the current location. The probability of an wetland species becoming
    # established is depended on
    # three major factors. First, the salinity must be below 1.0 ppt.
    # Second, there must be a period of 14 days with no flooding followed by a period of 14 with
    # water depths below 14 cm. Finally, the height of the habitat above mean water level determines
    # the final probability.
    def growth(self, loc):
        waveAmp = SpeciesModel.params.waveAmp[loc]
        meanSal = SpeciesModel.params.meanSal[loc]
        est     = SpeciesModel.params.treeEstCond[loc]
        bi      = SpeciesModel.params.biEstCond[loc]
        dsp     = SpeciesModel.dspModel[loc][self.abr]

        if bi: return 0.0

        return self.P[waveAmp, meanSal] * dsp * est

    def spread(self,loc):
        if self.dspClass == 3:
            waveAmp = SpeciesModel.params.waveAmp[loc]
            meanSal = SpeciesModel.params.meanSal[loc]
            est     = SpeciesModel.params.treeEstCond[loc]
            bi      = SpeciesModel.params.biEstCond[loc]

            if bi: return 0.0

            return self.P[waveAmp, meanSal] * est
        else:
            return 0.0

##\class SAVModel
##\brief This class handles the ecology for the SAV species.
##\role{Ecology}
##\detail
#
class SAVModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat = '', dspClass=0, ffibsScore=0, cover=0, loc=0, SAVData=None):
        SpeciesModel.__init__(self, index, name, abr, 'SAVModel', habitat, dspClass,ffibsScore, cover, loc)

        if SAVData == None:
            errorMessage = 'SAVModel: Error: No SAVData object passed to constructor (i.e. __init__() )\n'
            raise RuntimeError(errorMessage)

        self.betaIntercept = SAVData['Intercept']
        self.betaTemp      = SAVData['Temp']
        self.betaSal       = SAVData['Sal']
        self.betaDepth     = SAVData['Depth']

    def senescence(self, loc):
        return 0;

    def growth(self, loc):
        d = SpeciesModel.params.smDepth[loc]
        s = SpeciesModel.params.smSal[loc]
        t = SpeciesModel.params.smTemp[loc]
        bi = SpeciesModel.params.biEstCond[loc]

        if bi:
            return 0.0

        return min( 1.0, max(0.0, self.betaIntercept + self.betaTemp*t + self.betaSal*s + self.betaDepth*d   )  )

##\class BarrierIslandModel
##\brief This class handles the ecology for the barrier island species.
##\role{Ecology}
##\detail
#
class BarrierIslandModel(SpeciesModel):
    def __init__(self, index=0, name='', abr='', habitat='', dspClass=0, ffibsScore=0, cover=0, loc=0, Pdata=None, Ddata=None):
        SpeciesModel.__init__(self, index, name, abr, 'BarrierIslandModel', habitat, dspClass, ffibsScore, cover, loc)
        self.P = function.Function(Pdata['elvValue'], Pdata['rate'])
        self.D = function.Function(Ddata['elvValue'], Ddata['rate'])

    def senescence(self, loc):
        elv = SpeciesModel.params.biHeightAboveWater[loc]
        bi  = SpeciesModel.params.biEstCond[loc]

        if not bi:
            return 1.0

        return self.D[elv]

    def growth(self, loc):
        elv = SpeciesModel.params.biHeightAboveWater[loc]
        bi  = SpeciesModel.params.biEstCond[loc]
        #dsp = SpeciesModel.dspModel[self.loc][self.abr]

        if not bi:
            return 0.0

        return self.P[elv]# * dsp

##\class
##\brief
##\role{Ecology}
##\detail
#
class FloatingMarshModel(EmergentWetlandModel):
    def __init__(self, index=0, name='', abr='', habitat ='', dspClass=0, ffibsScore=0, cover=0, loc=0, Pdata=None, Ddata=None):
        EmergentWetlandModel.__init__(self, index, name, abr, habitat, dspClass, ffibsScore, cover, loc, Pdata, Ddata)
        self.modelType = 'FloatingMarshModel'
        self.transList = None

##\brief Growth/Senescence table storage
##\detail This class stores the growth and senescence tables
# for each species. The class is structured as a dictionary.
# Each species is identified by a unique abbreviation. For
# example, _Phragmities australis_ (Roseau cane) is represented
# by the abbreviation PAHU7. The abbreviations are based on the
# USDA species abbreviations. The abbreviations for each of the species
# used in this model are defined in the establishment table file.
#
# Each unique abbreviation is associated with the growth and senescence
# tables. These tables are defined in the establishment and senescence file.
#
# This class' job is to make the data from the establishment and senescence
# tables available to the other model components.
#
class SpeciesModelList(dict):
    ##\brief Class constructor
    def __init__(self):
        dict.__init__(self)

    ##\brief Load data from establishment and senescence files.
    ##@param [in] estFilename The name of the file containing the establishment tables.
    ##@param [in] mortFilename The name of the file containing the senescence tables.
    ##\details This member function loads the establishment and senescence tables
    # from the filenames passed in as arguments. The files are expected to be
    # MSExcel spreadsheets stored in the XLSX format. The data in each file
    # should be divided into a set of separate tabbed sheets. One sheet should
    # be named "VegTypeNames". This sheet should list all the species to include
    # in the model. There should also be one table for each of the species included in the
    # model. There are three exceptions to this: SAV, the bottomland hardwood forest species and the
    # barrier island species. SAV parameters are read from the model configuration file.
    # The barrier island species and the bottomland hardwood forest species are each represented by a single
    # sheet.
    #
    # The function performs the following steps:\n
    def config(self, estFilename, mortFilename):
        errorMessage = ''

        ###################################################
        ##
        # - __Step 1__: Read the VegTypeNames sheet from the
        #   establishment and senescence files.
        #
        #   The VegTypeNames file is where we get the
        #   names of the species in the model. If this
        #   sheet does not exist in the Excel files,
        #   or if the files do not exist, raise and
        #   exception.
        ###################################################

        try:
            estData = pandas.read_excel(estFilename, 'VegTypeNames', index_col=None)
        except IOError as error:
            errorMessage += 'SpeciesModelList: Error: Could not open file for reading: ' + str(estFilename) + '\n'
            errorMessage += 'SpeciesModelList: Error: Additional error info: ' + str(error) + '\n';
        except xlrd.biffh.XLRDError as error:
            errorMessage += 'SpeciesModelList: Error: While opening : ' + str(estFilename) + '\n'
            errorMessage += 'SpeciesModelList: Error: Additional error info: ' + str(error) + '\n';

        try:
            mortData = pandas.read_excel(mortFilename, 'VegTypeNames', index_col=None)
        except IOError as error:
            errorMessage += 'SpeciesModelList: Error: Could not open file for reading: ' + str(estFilename) + '\n' #shouldn't this be mortFilename?
            errorMessage += 'SpeciesModelList: Error: Additional error info: ' + str(error) + '\n';
        except xlrd.biffh.XLRDError as error:
            errorMessage += 'SpeciesModelList: Error: While opening : ' + str(estFilename) + '\n'
            errorMessage += 'SpeciesModelList: Error: Additional error info: ' + str(error) + '\n';

        if len(errorMessage):
            raise RuntimeError(errorMessage);

        ###################################################
        ##
        # - __Step 2__: Get the list of species from the 'Symbol'
        #   column. Raise exceptions as needed.
        #
        ###################################################
        try:
            speciesList = estData['Symbol']
        except KeyError as error:
            errorMessage += 'SpeciesModelList: Error: Symbol column not defined in establishment tables \n'
            raise RuntimeError(errorMessage);

        ###################################################
        ##
        # - __Step 3__: Read the establishment and senescence tables.
        #
        #   For each species, configure the P & D
        #   Functions for each species and push everything
        #   into a dictionary (dict).
        #
        ###################################################

        for spSymbol in speciesList:

            spInfo      = estData[ estData['Symbol'] == spSymbol ].to_dict('record')[0]
            spID        = spInfo['ID']
            spName      = spInfo['Common Name']
            spModelType = spInfo['ModelType']
            spHabitat   = spInfo['Habitat']
            spDspClass  = spInfo['Dispersal Class']
            spFFIBSScore = spInfo['FFIBS score']

            print(('SpeciesModelList: Msg: Configuring model for species ' + str(spSymbol) + ', model type is ' + spModelType + '.'))

            ###########################################
            # Bottomland Hardwood Forest Model
            #
            ###########################################
            if ( spModelType == 'BottomlandHardwoodForestModel'):
                reader = function.ReadVegTable1D()
                Pdata = None
                Ddata = None

                try:
                    Pdata = reader.read(estFilename, 'BottomlandHardwoodForest', spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'


                try:
                    Ddata = reader.read(mortFilename, 'BottomlandHardwoodForest', spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                spModel = BottomlandHardwoodForestModel(index=spID, name=spName, abr=spSymbol, habitat = spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1, Pdata=Pdata, Ddata=Ddata)
                self.__setitem__(spSymbol, spModel)

            ###########################################
            # Emergent Wetland Model
            #
            ###########################################
            elif (spModelType == 'EmergentWetlandModel'):
                reader = function.ReadVegTable2D()
                Pdata = None
                Ddata = None
                try:
                    Pdata = reader.read(estFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n';
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                try:
                    Ddata = reader.read(mortFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                spModel = EmergentWetlandModel(index=spID, name=spName, abr=spSymbol, habitat = spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1, Pdata=Pdata, Ddata=Ddata)
                self.__setitem__(spSymbol, spModel)


            ###########################################
            # Swamp Forest Model
            #
            ###########################################
            elif (spModelType == 'SwampForestModel'):
                reader = function.ReadVegTable2D()
                Pdata  = None
                Ddata  = None
                try:
                    Pdata = reader.read(estFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n';
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                try:
                    Ddata = reader.read(mortFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                spModel = SwampForestModel(index=spID, name=spName, abr=spSymbol, habitat = spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1, Pdata=Pdata, Ddata=Ddata)
                self.__setitem__(spSymbol, spModel)


            ###########################################
            # Floating Marsh Model
            #
            ###########################################
            elif (spModelType == 'FloatingMarshModel'):
                reader = function.ReadVegTable2D()
                Pdata = None
                Ddata = None
                try:
                    Pdata = reader.read(estFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n';
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                try:
                    Ddata = reader.read(mortFilename, spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                spModel = FloatingMarshModel(index=spID, name=spName, abr=spSymbol, habitat = spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1, Pdata=Pdata, Ddata=Ddata)
                self.__setitem__(spSymbol, spModel)



            ###########################################
            # Submerged Aquatic Vegetation Model
            #
            ###########################################
            elif ( spModelType == 'SAVModel'):
                savData = None
                try:
                    savData = pandas.read_excel(estFilename, spSymbol, index_col=0 ).to_dict('dict')['Value']
                except IOError as error:
                    errorMessage += 'SpeciesModelList: Error : Could not open file for reading : ' + estFilename + '\n'
                    errorMessage += 'SpeciesModelList: Extra error info: ' + str(error) + '\n'
                except xlrd.biffh.XLRDError as error:
                    errorMessage += 'SpeciesModelList: Error : Could not find sheet name ' + spSymbol + ' in file ' + estFilename + '\n'
                    errorMessage += 'SpeciesModelList: Extra error info: ' + str(error) + '\n'
                except KeyError as error:
                    errorMessage += 'SpeciesModelList: Error : Dictionary key error. SAV table probably has wrong column headings\n'
                    errorMessage += 'SpeciesModelList: Error : '

                spModel = SAVModel(index=spID, name=spName, abr=spSymbol, habitat=spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1, SAVData=savData)
                self.__setitem__(spSymbol, spModel)

            ###########################################
            # Barrier Island Model
            #
            ###########################################
            elif ( spModelType == 'BarrierIslandModel'):
                reader = function.ReadVegTable1D()
                Pdata = None
                Ddata = None

                try:
                    Pdata = reader.read(estFilename, 'BarrierIsland', spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'


                try:
                    Ddata = reader.read(mortFilename, 'BarrierIsland', spSymbol)
                except RuntimeError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'
                except ValueError as error:
                    errorMessage += str(error) + '\n'
                    errorMessage += 'SpeciesModelList: Error: error while working on species ' + str(spSymbol) + '\n'

                spModel = BarrierIslandModel(index=spID, name=spName, abr=spSymbol, habitat=spHabitat, dspClass=spDspClass,ffibsScore=spFFIBSScore, cover=0, loc=-1, Pdata=Pdata, Ddata=Ddata)
                self.__setitem__(spSymbol, spModel)

            ###########################################
            # Null Model
            #
            ###########################################
            elif ( spModelType == 'NullModel'):
                spModel = NullModel(index=spID, name=spName, abr=spSymbol, habitat=spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1)
                self.__setitem__(spSymbol, spModel)
                
            ###########################################
            # Null Model Coverage
            #
            ###########################################
            elif ( spModelType == 'NullModel_Coverage'):
                spModel = NullModel_Coverage(index=spID, name=spName, abr=spSymbol, habitat=spHabitat, dspClass=spDspClass, ffibsScore=spFFIBSScore, cover=0, loc=-1)
                self.__setitem__(spSymbol, spModel)

            ###########################################
            # Error State
            #
            ###########################################
            else:
                errorMessage += 'SpeciesModelList: Error: Unknown model type specified : ' + spModelType + '\n'

        if len(errorMessage):
            raise RuntimeError(errorMessage)

##\class Params
##\brief This class is the interface between the outside world and the model.
##\role{Machinery}
##\detail This class is an interface between the internal components of the model
# and the outside world. That is, everything that is read in from outside the
# model is stored in an object of class Params. Params then makes all of this
# information available to the other model components.
#
# It is important to note that
# there is only meant to be ONE object derived from this class, and that object
# is owned by the object of class Model.
#
# Params handles reading values from the model configuration file. Params
# also manages all of the model IO streams, as well as the landscape.Landscape
# objects that hold the spatial data files used by the model, including the
# hydrology, salinity, elevation, initial conditions and other spatial data sets. Params also
# holds a special object of type model.SpeciesModelList. This class is
# used to store the the growth and senescence tables for each species along with
# any other species specific information.
class Params(object):
    ##\brief Class constructor
    def __init__(self):
        ##\brief Holds the configuration information
        self.configDict     = config.Config()
        ##\brief Model start year
        self.startYear      = 0
        ##\brief Model end year
        self.endYear        = 0
        ##\brief List of years to read data from the wetland morph output file.
        self.wetlandMorphYears = []
        ##\brief List of model input filenames.
        #
        self.inputStrm      = {'WaveAmplitudeFile'               :None, \
                              'MeanSalinityFile'                 :None, \
                              'SummerMeanWaterDepthFile'         :None, \
                              'SummerMeanSalinityFile'           :None, \
                              'SummerMeanTempFile'               :None, \
                              'HeightAboveWaterFile'             :None, \
                              'BarrierIslandEstCondFile'         :None, \
                              'BarrierIslandHeightAboveWaterFile':None, \
                              'WetlandMorphLandWaterFile'        :None, \
                              'AcuteSalinityStressFile'         :None, \
                              'TreeEstCondFile'                  :None }

        ##\brief List of model output filenames.
        #
        self.outputStrm      = { 'OutputFile'            :None}

        ##\brief The annual wave amplitude map
        self.waveAmp           = landscape.Landscape()
        ##\brief The annual mean salinity map
        self.meanSal           = landscape.Landscape()
        ##\brief The summer mean water depth map
        self.smDepth           = landscape.Landscape()
        ##\brief The summer mean water salinity map
        self.smSal             = landscape.Landscape()
        ##\brief The summer mean salinity map
        self.smTemp            = landscape.Landscape()
        ##\brief The barrier island establishment condition map
        self.biEstCond         = landscape.Landscape()
        ##\brief The barrier height above water file
        self.biHeightAboveWater = landscape.Landscape()
        ##\brief The map of heights above water map
        self.heightAboveWater   = landscape.Landscape()
        ##\brief The land/water map
        self.landWater         = landscape.Landscape()
        ##\brief The map of accutue salinity stress
        self.acuteSal         = landscape.Landscape()
        ##\brief The tree establishment conditions map
        self.treeEstCond       = landscape.Landscape() # TreeEstCondFile

        ##\brief The model initial conditions.
        #
        self.initCond        = landscape.LandscapePlus() # InitialConditionFile

        ##\brief The species model list (SpeciesModelList)
        #
        self.spModelList     = SpeciesModelList()

        ##\brief The elevation threshold for tree establishment (soon to be deprecated).
        #
        self.elevationThreshold = 0.1525

        self.yearFormatString   = '{:02d}'

    ##\brief Open input/output files
    ##\param key The keyword name for the file as it appears in the configuration file.
    ##\param mode Should the file be opened for reading or writing.
    ##\details This is a convenience function. It checks to make sure the
    # file actually exists and opens if it does. If any problems are encountered,
    # such as file not existing or cannot be opened, an exception is raised.
    def open_file(self, key, mode):
        try:
            filename = self.configDict[key];
        except KeyError as error:
            errorMessage  = 'Params: Error: ' + str(key) + ' is not defined in the configuration file.\n'
            errorMessage += 'Params: Error: Additional error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage);

        try:
            strm = open(filename, mode)
        except IOError as error:
            errorMessage  = 'Params: Error: Could not open the file ' + filename + (' for reading' if mode == 'r' else ' for writing') +'\n'
            errorMessage += 'Params: Error: Additional error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage);

        return strm

    ##\brief Read the model configuration information.
    ##\param argv A object containing the name of the configuration file.
    ##\details This member function reads in the information
    # contained in the model configuration file and sets up
    # a number of program elements based on this information
    #
    # The function performs the follwing steps:
    #
    def config(self, argv):

        errorMessage = ''


        ###################################################
        ##
        # - __Step 1__: Read in the configuration data from the
        #   configuration file.
        #
        ###################################################
        try:
            self.configDict.config(argv);
        except RuntimeError as error:
            errorMessage = 'Params: Error: An error occurred reading the configuration file.\n' + str(error)
            raise RuntimeError(errorMessage)

        ###################################################
        ##
        # - __Step 1.5__: Get/set some detail assets for the model
        #
        ###################################################
        if 'YearFormatString' in self.configDict:
            self.yearFormatString = self.configDict['YearFormatString']
            print(('Params: Msg: Using user defined year format string : ' + self.yearFormatString))

        ###################################################
        ##
        # - __Step 2__: Get the start end end year of the simulation
        #
        ###################################################
        try:
            self.startYear = int( self.configDict['StartYear'] )
            print(('Params: Msg: StartYear = ' + str(self.startYear)))
        except KeyError as error:
            errorMessage += 'Params: Error: StartYear is not defined in the configuration file\n'

        try:
            self.endYear = int( self.configDict['EndYear'])
            print(('Params: Msg: EndYear = ' + str(self.endYear)))
        except KeyError as error:
            errorMessage += 'Params: Error: EndYear is not defined in the configuration file\n'

        ###################################################
        ##
        # - __Step 3__: Read the initial conditions file
        #
        ###################################################
        try:
            print(('Params: Msg: Reading initial conditions from ' + self.configDict['InitialConditionFile']))
            initCondStrm = self.open_file('InitialConditionFile','r')
            reader = landscape.ReadASCIIGridPlus()
            reader.read(initCondStrm, self.initCond)
            initCondStrm.close();

            # Check to see if the input file contains a DEAD_Flt class.
            # If there isn't one, the add it to the initial conditions.
            if not('DEAD_Flt' in self.initCond.table):
                print ('Params: Msg: Adding DEAD_Flt class to the initial conditions because the type was not defined.')
                self.initCond.table['DEAD_Flt'] = 0.0
            if not('BAREGRND_Flt' in self.initCond.table):
                print ('Params: Msg: Adding BAREGRND_Flt class to the initial conditions because the type was not defined.')
                self.initCond.table['BAREGRND_Flt'] = 0.0
            if not('BAREGRND_NEW' in self.initCond.table):
                print ('Params: Msg: Adding BAREGRND_NEW class to the initial conditions because the type was not defined.')
                self.initCond.table['BAREGRND_NEW'] = 0.0
            if not('BAREGRND_OLD' in self.initCond.table):
                print ('Params: Msg: Adding BAREGRND_OLD class to the initial conditions because the type was not defined.')
                self.initCond.table['BAREGRND_OLD'] = 0.0
            if not('FFIBS' in self.initCond.table):
                print ('Params: Msg: Adding FFIBS class to the initial conditions because the type was not defined.')
                self.initCond.table['FFIBS'] = -9999 #numpy.nan #0.0 #numpy.nan
            if not('pL_BF' in self.initCond.table):
                print ('Params: Msg: Adding pl_BF class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_BF'] = 0.0
            if not('pL_SF' in self.initCond.table):
                print ('Params: Msg: Adding pl_SF class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_SF'] = 0.0
            if not('pL_FM' in self.initCond.table):
                print ('Params: Msg: Adding pl_FM class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_FM'] = 0.0
            if not('pL_IM' in self.initCond.table):
                print ('Params: Msg: Adding pl_IM class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_IM'] = 0.0
            if not('pL_BM' in self.initCond.table):
                print ('Params: Msg: Adding pl_BM class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_BM'] = 0.0
            if not('pL_SM' in self.initCond.table):
                print ('Params: Msg: Adding pl_SM class to the initial conditions because the type was not defined.')
                self.initCond.table['pL_SM'] = 0.0


            #at this point, the initiCond.table is correct (values and nans only for FFIBS)

        except RuntimeError as error:
            errorMessage += str(error)

        #self.make_lon_lat_file(self.initCond)

        ###################################################
        ##
        # - __Step 4__: Configure the SpeciesModelList object
        #
        ###################################################
        try:
            estFilename = self.configDict['EstFilename']
        except KeyError as error:
            errorMessage += 'Params: Error: KeyError in EstFilename : ' + str(error) + '\n'
        except RuntimeError as error:
            errorMessage += 'Params: Error: RuntimeError in EstFilename : ' + str(error)

        try:
            mortFilename = self.configDict['MortFilename']
        except KeyError as error:
            errorMessage += 'Params: Error: KeyError in MortFilename : ' + str(error) + '\n'
        except RuntimeError as error:
            errorMessage += 'Params: Error: RuntimeError in MortFilename : ' + str(error)

        try:
            self.spModelList.config(estFilename, mortFilename)
        except KeyError as error:
            errorMessage += 'Params: Error: KeyError in spModelList.config(): ' + str(error) + '\n'
        except RuntimeError as error:
            errorMessage += 'Params: Error: RuntimeError in spModelList.config(): ' + str(error)

        ###################################################
        ##
        # - __Step 5__: Open all the input files for reading.
        #
        ###################################################
        for key in iter(self.inputStrm.keys()):
            try:
                self.inputStrm[key] = self.open_file(key, 'r')
                print(('Params: Msg: Opened file for reading: ' + str(key) + ' = ' + self.configDict[key]))
            except RuntimeError as error:
                errorMessage += str(error)

        ###################################################
        ##
        # - __Step 5.5__: Get the years to read the wetland morph data
        #
        ###################################################
        try:
            yearsString = self.configDict['WetlandMorphYears']
            yearList    = [ x for x in yearsString.split(',')]
            for year in yearList:
                rangeMatch = re.match('[0-9]+:[0-9]+', year)
                if rangeMatch is not None:
                    yearList.remove(year)
                    start,end = year.split(':')
                    for y in range( int(start),int(end)+1 ):
                        yearList.append(str(y))
            self.wetlandMorphYears = [ int(x) for x in yearList ]
            print(('Params: Msg: Wetland Morph land/water read years = ' + str(self.wetlandMorphYears)))
        except KeyError as error:
            errorMessage += 'Params: Error: WetlandMorphYears not defined in the configuration file\n'


        ###################################################
        ##
        # - __Step 6__: Open all the output files for writing.
        #
        ###################################################

        # Steps 6.1 to 6.3: Set up the single year files.
        # Step 6.1: Get the filename template for the single year files.
        try:
            outputTemplate = None
            outputTemplate = self.configDict['OutputTemplate']
            print(('Params: Msg: OutputTemplate = ' + self.configDict['OutputTemplate']))
        except KeyError as error:
            errorMessage += 'Params: Error: OutputTemplate is not defined in the configuration file\n'
            errorMessage += 'Params: Error: Additional error info: ' + str(error)

        # Step 6.2: Get the years we want output stored in separate files
        try:
            outputYear = None
            outputYear = self.configDict['OutputYears']
            outputYear = [ int(x) for x in outputYear.split(',')]
            print(('Params: Msg: OutputYears = ' + str(outputYear)))
        except KeyError as error:
            errorMessage += 'Params: Error: OutputYears is not defined in the configuration file\n'
            errorMessage += 'Params: Error: Additional error info: ' + str(error)

        # Step 6.3: Build the file names for the separate files.
        if outputTemplate != None and outputYear != None:
            for year in outputYear:
                filename = re.sub(r'<YEAR>', self.yearFormatString.format(year), outputTemplate);
                self.configDict['__Single_'+str(year)] = filename
                self.outputStrm['__Single_'+str(year)] = None


        # Step 6.4: Open all the output files.
        self.outputStrm.pop('OutputFile')
        for key in iter(self.outputStrm.keys()):
            try:
                self.outputStrm[key] = self.open_file(key, 'w')
                print(('Params: Msg: Output file opened ' + str(self.configDict[key])))
            except RuntimeError as error:
                errorMessage += str(error)
                errorMessage += 'Params: Error: This is very odd\n'

        ###################################################
        ##
        # - __Step 7__: Get the files for the planting models read
        #
        ###################################################

        # try:
        #     self.plantingDict = dict()
        #     self.plantingPath = self.configDict['PlantingPath']
        #     shapefileList = glob.glob( os.path.join(self.plantingPath,'*.shp') )
        #     if ( len(shapefileList) == 0 ):
        #         raise exceptions.Warning('Params: Msg: No shapefiles found in PlantingPath = ' + self.plantingPath)
        #
        #     driver = ogr.GetDriverByName('ESRI Shapefile')
        #     for file in shapefileList:
        #         dataSrc = driver.Open(file, 0)
        #         if dataSrc == None:
        #             raise exceptions.Warning('Params: Msg: Could not open shapefile named ' + file)
        #
        #         print 'Params: Msg: Reading planting shapefile data from ' + file
        #         self.plantingDict[file] = dataSrc
        #
        # except exceptions.KeyError:
        #     print 'Params: Msg: No Plantings files found. Moving on.'
        # except exceptions.Warning as warning:
        #     print warning

        try:
            self.plantingsStrm = self.open_file('PlantingTextFile','r')
            print(('Params: Msg: Reading planting infomration from ' + self.configDict['PlantingTextFile']))
        except RuntimeError as error:
            print ('Params: Warning: Plantings will not be used in this run. More information follows.')
            print (error)
            self.plantingsStrm = None


        ###################################################
        ##
        # - __Step 8__: If any errors have occurred, raise an
        # exception
        #
        ###################################################
        if len(errorMessage):
            raise RuntimeError(errorMessage)

        ###################################################
        ##
        # - __Step 9__: Give all of the landscape objects their
        #  correct size
        ###################################################
        print ('Params: Msg: Giving all data layers their proper size')
        reader = landscape.ReadASCIIGrid()
        reader.read(self.inputStrm['WaveAmplitudeFile'],        self.waveAmp)
        reader.read(self.inputStrm['MeanSalinityFile'],         self.meanSal)
        reader.read(self.inputStrm['SummerMeanWaterDepthFile'], self.smDepth)
        reader.read(self.inputStrm['SummerMeanSalinityFile'],   self.smDepth)
        reader.read(self.inputStrm['SummerMeanTempFile'],       self.smTemp)
        reader.read(self.inputStrm['HeightAboveWaterFile'],     self.heightAboveWater)
        reader.read(self.inputStrm['BarrierIslandEstCondFile'], self.biEstCond)
        reader.read(self.inputStrm['BarrierIslandHeightAboveWaterFile'], self.biHeightAboveWater)
        reader.read(self.inputStrm['WetlandMorphLandWaterFile'],self.landWater)
        reader.read(self.inputStrm['AcuteSalinityStressFile'], self.acuteSal)
        reader.read(self.inputStrm['TreeEstCondFile'],          self.treeEstCond)

        print ('Params: Msg: Rewinding all the input streams')
        for strm in iter(self.inputStrm.values()):
            strm.seek(0,0)

        if 'XFile' in self.configDict:
            print ('Params: Msg: Access to XFile denied')

    ##\brief Actions take at the end of the simulation.
    ##\details The primary responsibility of this function
    # is to close all of the open file streams. It is
    # called at the end of the simulation as the program
    # is getting ready to exit.
    def done(self):
        for iter in iter(self.inputStrm.values()):
            iter.close();

        for iter in iter(self.outputStrm.values()):
            iter.close();

    def make_lon_lat_file(self, initCond):
        print ('Params: Msg: Building locations data')
        nrow = int( initCond.nrow )
        ncol = int( initCond.ncol )


        ret = pandas.DataFrame(index=list(range(0,nrow*ncol)), columns=['CellID','row','col','x','y'])
        count = 0
        for r,c in itertools.product( list(range(0,nrow)), list(range(0,ncol))):
            if initCond.has_data_at(r,c):
                cellID = initCond.data[(r,c)]
                lon    = initCond.xllcorner + (c + 0.5) * initCond.cellsize
                lat    = initCond.yllcorner + ((nrow-1) - r + 0.5) * initCond.cellsize
                ret['CellID'][count] = cellID
                ret['row']   [count] = r
                ret['col']   [count] = c
                ret['x']     [count] = lon
                ret['y']     [count] = lat
                count += 1

        print ('Params: Msg: Writting locations data')
        ret.to_csv('locations.csv')


    def __del__(self):
        pass

##\class PatchModel
##\brief This class handles the plant community dynamics at a single location
##\role {Ecology/Machinery}
##\detail This class handles the ecology of a single location in the landscape.
# The information for each location is stored as a Python dictionary. They
# keys for the dictionary are the abbreviated names, stored as a string, for
# each of the species. The value for each entry in the dictionary is an object
# instantiated from a class that inherits from SpeciesModel. That is, each value
# is an object of one of the following classes: BottomlandHardwoodForestModel, EmergentWetlandModel,
# SAVModel, BarierIslandModel, FloatingMarshModel, NullModel_Coverage or NullModel.
class PatchModel(dict):
    def __init__(self):
        dict.__init__(self)
        self.params = None

    def config(self, params):
        self.params = params

    def update(self, loc, spCoverList, spModelList, firstyear):
        occupied         = 0.0
        unoccupied       = 0.0
        lost             = 0.0
        growthLikelihood = 0.0
        spreadLikelihood = 0.0

        #######################################################
        # WARNING: The order of steps here is important. Additional
        #          care is required here because you can change
        #          the order without actually causing the model to break.
        #          Changing the order will change the way the model simulates
        #          the ecology the plants. Unless you really know what
        #          you are doing, do not change the order. Even if
        #          you think you know what you are doing, think about
        #          what you are going to do carefully. Perhaps you should
        #          go take and nap and think things over before you make
        #          changes here. If you break things, it will be on your head.
        #
        #######################################################

        #######################################################
        # Step 1: Work out the increase and decrease in cover
        #         for all types except, SAV, WATER and
        #         FloatingMarshModel. Those types have to be
        #         handled in a special way.
        #
        #######################################################


        # Step 1.05: Allow "weedy" species to establish on new bareground. If none can
        # establish, then convert it to old bareground.
        
        if spCoverList['BAREGRND_NEW']>0.0:
 #           for spName, spModel in itertools.filterfalse(lambda kv: kv[0] == 'BAREGRND_NEW' or kv[0] == 'BAREGRND_OLD' or kv[0] == 'SAV' or kv[0] == 'WATER' or kv[0] == 'FFIBS' or kv[1].modelType == 'FloatingMarshModel', list(spModelList.items())):
    #technicaly not needed because spread will just be 0 for all the others... 
            for spName, spModel in itertools.filterfalse(lambda kv: kv[1].modelType == 'NullModel' or kv[1] == 'NullModel_Coverage' or kv[1].modelType == 'FloatingMarshModel', list(spModelList.items())):             
                spreadLikelihood    += spModel.spread(loc)
            if spreadLikelihood:
                for spName, spModel in itertools.filterfalse(lambda kv: kv[0] == 'BAREGRND_NEW' or kv[0]== 'BAREGRND_OLD' or kv[0]=='SAV' or kv[0]=='WATER' or kv[0] == 'FFIBS' or kv[1].modelType=='FloatingMarshModel', list(spModelList.items())):
                    spread               = spModel.spread(loc)/spreadLikelihood
                    spCoverList[spName] += spread * spCoverList['BAREGRND_NEW']
            else:
                spCoverList['BAREGRND_OLD'] += spCoverList['BAREGRND_NEW']

        spCoverList['BAREGRND_NEW']=0.0 #it either already was 0.0, or it was just reassigned to BG OLD or vegetation (meaning it now needs to be 0.0)
        

        # Step 1.1: Work out the loss of cover for all species.
        # Skip the Floating Marsh types, because they are handled separately.
        # Skip the SAV and WATER types because they do not "senesce" in the same way
        # as the other types, and because their growth is not computed the same way
        # as the other types.

        spreadLikelihood = 0.0

        for spName, spModel in itertools.filterfalse(lambda kv: kv[0] == 'BAREGRND_NEW' or kv[0] == 'BAREGRND_OLD' or kv[0] == 'SAV' or kv[0] == 'WATER' or kv[0] == 'FFIBS' or kv[1].modelType == 'FloatingMarshModel', list(spModelList.items())):
            cover                = spCoverList[spName]
            death                = spModel.senescence(loc)
            #occupied            += cover
            lost                += death * cover
            spCoverList[spName] -= death * cover
            growthLikelihood    += spModel.growth(loc)
            spreadLikelihood    += spModel.spread(loc)

            growth = 999
            spread = 999


        unoccupied += lost

        spCoverList['BAREGRND_NEW'] = unoccupied #Newly dead areas become new

        # step 1.15: Work out the total area that is currently unoccupied, i.e. area in BAREGRND_OLD and newly dead area (same as BAREGROUND_NEW)
        # This is the area available for a terrestrial species to claim via establishment
        unoccupied       += spCoverList['BAREGRND_OLD']
        unoccupied_0 = unoccupied


        # Step 1.2: Work out the gain in cover for all species.
        # Again, skip SAV, WATER and the floating marsh types because they do not
        # operate the same way as the other types.
        if growthLikelihood:
            for spName, spModel in itertools.filterfalse(lambda kv: kv[0] == 'BAREGRND_NEW' or kv[0]== 'BAREGRND_OLD' or kv[0]=='SAV' or kv[0]=='WATER' or kv[1].modelType == 'NullModel' or kv[1].modelType=='FloatingMarshModel', list(spModelList.items())):
                growth               = spModel.growth(loc)/growthLikelihood
                spCoverList[spName] += growth * unoccupied_0
                unoccupied          -= growth * unoccupied_0 #Should always end up at 0

        else:
            if spreadLikelihood and firstyear == 0:
                for spName, spModel in itertools.filterfalse(lambda kv: kv[0] == 'BAREGRND_NEW' or kv[0]== 'BAREGRND_OLD' or kv[0]=='SAV' or kv[0]=='WATER' or kv[1].modelType == 'NullModel' or kv[1].modelType=='FloatingMarshModel', list(spModelList.items())):
                    spread               = spModel.spread(loc)/spreadLikelihood
                    spCoverList[spName] += spread * unoccupied_0
                    unoccupied          -= spread * unoccupied_0 #Should always end up at 0


        if growthLikelihood:
            spCoverList['BAREGRND_NEW']=unoccupied #should be the same as setting to 0
            spCoverList['BAREGRND_OLD']=unoccupied #should be the same as setting to 0
        elif spreadLikelihood and firstyear==0:
            spCoverList['BAREGRND_NEW']=unoccupied #should be the same as setting to 0
            spCoverList['BAREGRND_OLD']=unoccupied #should be the same as setting to 0


        #######################################################
        # Step 2: Work out the change in cover for the floating
        #         types.
        #
        #######################################################

        # Step 2.1: Compute the area lost per floating species and the total area lost by all floating species
        growthLikelihood    = 0.0
        deadThinFloating    = 0.0
        deadThickFloating   = 0.0

        # For thin mat, the loss is:
        spModel = spModelList['ELBA2_Flt']
        cover = spCoverList['ELBA2_Flt']
        death = spModel.senescence(loc)
        spCoverList['ELBA2_Flt']    -= death * cover
        deadThinFloating += death*cover
        growthLikelihood       += spModel.growth(loc)

        #For thick mat, the loss is:
        spModel = spModelList['PAHE2_Flt']
        cover = spCoverList['PAHE2_Flt']
        death = spModel.senescence(loc)
        spCoverList['PAHE2_Flt']    -= death * cover
        deadThickFloating += death*cover
        growthLikelihood       += spModel.growth(loc)

        # Step 2.2: Compute the area gained by each floating species.
        #If there can be growth, it can be on dead thin mat, dead thick mat, or bareground float

        if growthLikelihood:
            deadFloating_0 = deadThickFloating + deadThinFloating + spCoverList['BAREGRND_Flt']
            spCoverList['BAREGRND_Flt'] = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'FloatingMarshModel', iter(spModelList.items())):
                growth = spModel.growth(loc)/growthLikelihood
                spCoverList[spName]    += growth * deadFloating_0
        else:
            spCoverList['WATER']    += deadThinFloating + spCoverList['BAREGRND_Flt']
            spCoverList['DEAD_Flt'] += deadThinFloating + spCoverList['BAREGRND_Flt'] #passed to Morph, where it is set as water 1 m deep
            spCoverList['BAREGRND_Flt']    = deadThickFloating


        #######################################################
        # Step 3: Work out the change in cover for SAV and WATER
        #
        #######################################################

        SAVModel             = spModelList['SAV']
        totalWater           = spCoverList['SAV'] + spCoverList['WATER']
        spCoverList['SAV']   = SAVModel.growth(loc) * totalWater
        spCoverList['WATER'] = totalWater - spCoverList['SAV']

        #######################################################
        # Step 4: Decrease in coverage due to acute salinity stress
        #
        #######################################################

        acute_sal_stress  = SpeciesModel.params.acuteSal[loc]

        if acute_sal_stress:
            # Step 4.1: Eliminate coverage of flotant marsh. It follows the same process as step 2 with thin mat becoming open water and thick becoming bareground_flt
            deadThinFloating    = spCoverList['ELBA2_Flt']
            deadThickFloating   = spCoverList['PAHE2_Flt']
            spCoverList['ELBA2_Flt'] -= spCoverList['ELBA2_Flt']
            spCoverList['PAHE2_Flt'] -= spCoverList['PAHE2_Flt']
            spCoverList['WATER']    += deadThinFloating
            spCoverList['DEAD_Flt'] += deadThinFloating  #passed to Morph, where it is set as water 1 m deep
            spCoverList['BAREGRND_Flt']    += deadThickFloating

            # Step 4.2: Eliminate coverage of fresh attached marsh coverage.

            for spName, spModel in filter(lambda kv: kv[1].modelType == 'EmergentWetlandModel' and kv[1].habitat == 'Fresh', iter(spModelList.items())):
                spCoverList['BAREGRND_NEW'] += spCoverList[spName]
                spCoverList[spName] -= spCoverList[spName]#wipes out the fresh attached marsh

        #######################################################
        # Step 5: Check sum
        #
        #######################################################

        checkSum = 0.0
       # for spName in filter( lambda spName : spName != 'DEAD_Flt' or spName != 'FIBS', iter(spModelList.keys())):

        for spName, spModel in itertools.filterfalse(lambda kv: kv[1].modelType == 'NullModel', iter(spModelList.items())):
            checkSum += spCoverList[spName]

        tol = 0.005
        if not( (1.0 - tol) < checkSum and checkSum < (1.0 + tol) ):
            print(('PatchMod: Error: checkSum is out of range. checkSum = ' + str(checkSum) + ' should be 1.0'))
      #      if numpy.isnan(checkSum):  
      #      print(loc)
       #     print(spCoverList)
        #   print()



         #######################################################
        # Step 6: Calculate the FFIBS score
        #
        #######################################################

     # FFIBS is Forested, Fresh, Intermediate, Brackish, Saline and does not include the flotant
        ffibsValues = 0.0
        ffibsCover = 0.0
        for spName, spModel in iter(spModelList.items()):
            if spModel.ffibsScore>-9999:
          #  if ~numpy.isnan(spModel.ffibsScore):
                ffibsValues += (spModel.ffibsScore*spCoverList[spName])
                ffibsCover += spCoverList[spName]
        if ffibsValues:
            spCoverList['FFIBS'] = ffibsValues/ffibsCover
        else:
            spCoverList['FFIBS']=-9999

               # fibsList.append(spModel.fibsScore*spCoverList[spName])
                #fibsCover.append(spCoverList[spName])

            #fibsCover.append(spCoverList[spName]*(spModel.fibsScore/spModel.fibsScore)) #rather goofy way to exclude coverages that don't have a FIBS score, but it works
        #if sum(~numpy.isnan(fibsCover)) > 0:
          #  spCoverList['FIBS'] = numpy.nansum(fibsList)/numpy.nansum(fibsCover)
          
         #######################################################
        # Step 7: Calculate the percent of total land for 
        # bottomland hardwood forest, swamp forest, fresh marsh, inter. marsh, brackish marsh, saline marsh
        #
        #######################################################
        BI = 0.0
        for spName, spModel in filter(lambda kv: kv[1].modelType == 'BarrierIslandModel', iter(spModelList.items())):
                BI += spCoverList[spName]     
                
        total_vegetated_land = ffibsCover + BI ####BARRIER ISLAND SPECIES ARE COUNTED WITH SALINE MARSH 
        
        if total_vegetated_land > 0:
            BF = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'BottomlandHardwoodForestModel', iter(spModelList.items())):
                BF += spCoverList[spName]        
            
            SF = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'SwampForestModel', iter(spModelList.items())):
                SF += spCoverList[spName]
    
            FM = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'EmergentWetlandModel' and kv[1].habitat == 'Fresh' , iter(spModelList.items())):
                FM += spCoverList[spName]        
            
            IM = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'EmergentWetlandModel' and kv[1].habitat == 'Intermediate', iter(spModelList.items())):
                IM += spCoverList[spName]
       
            BM = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'EmergentWetlandModel' and kv[1].habitat == 'Brackish', iter(spModelList.items())):
                BM += spCoverList[spName]        
            
            SM = 0.0
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'EmergentWetlandModel' and kv[1].habitat == 'Saline', iter(spModelList.items())):
                SM += spCoverList[spName]
            
            for spName, spModel in filter(lambda kv: kv[1].modelType == 'BarrierIslandModel', iter(spModelList.items())):
                SM += spCoverList[spName] 
             
            spCoverList['pL_BF'] = BF/total_vegetated_land
            spCoverList['pL_SF'] = SF/total_vegetated_land
            numpy.seterr('raise') 
            spCoverList['pL_FM'] = FM/total_vegetated_land
            spCoverList['pL_IM'] = IM/total_vegetated_land
            spCoverList['pL_BM'] = BM/total_vegetated_land
            spCoverList['pL_SM'] = SM/total_vegetated_land

        else:
            spCoverList['pL_BF'] = 0.0
            spCoverList['pL_SF'] = 0.0
            spCoverList['pL_FM'] = 0.0
            spCoverList['pL_IM'] = 0.0
            spCoverList['pL_BM'] = 0.0
            spCoverList['pL_SM'] = 0.0
        
      
        

    #def copy_to_str(self, dataFormat='{}'):
    #    ret = ''
    #    sep = ''
    #    for value in self.itervalues():
    #        ret = sep + dataFormat.format(value)
    #        sep = ', '
    #    return(ret)

    #def copy_to_stream(self, stream = sys.stdout, dataFormat='{}'):
    #    sep=''
    #    for key, value in self.iteritems():
    #        print >> stream, sep + dataFormat.format(value)
    #        sep = ', '

##\class DynamicsModel
##\brief This class coordinates the updating of plant community dynamics
##\role{Ecology}
##\detail This class handles updating the spatial distribution of species
# in response to environmental conditions and plant dispersal ability.
# This class defines the spatial structure of the model using a
# raster grid. The grid structure is inherited from landscape.LandscapePlus.
# Each element of the raster contains an object of type PatchModel.
# Local dynamics, those taking place within each raster cell, are
# computed by objects of class PatchModel.
class DynamicsModel(landscape.LandscapePlus):
    def __init__(self):
        landscape.LandscapePlus.__init__(self)
        self.patchModel = PatchModel()
        self.spModelList = None
        self.locList     = {}
        self.firstyear = 1
    #@profile
    def config(self, params):
        landscape.Landscape.copy(self, params.initCond)
        self.patchModel.config(params)
        self.spModelList = params.spModelList
        for row, col in filter(   lambda rc: params.initCond.has_data_at(rc[0],rc[1]), itertools.product((range(int(params.initCond.nrow))), (range(int(params.initCond.ncol))) )     ):
            patchIndex                 = params.initCond.data[row,col]
            patchInitCond              = params.initCond[(row,col)]
            self.table[patchIndex]     = patchInitCond
            self.locList[patchIndex]   = (row,col)

        patchInitCond = params.initCond.table.iloc[0].to_dict()
        for key in iter(patchInitCond.keys()):
            patchInitCond[key] = 0.0
        self.table  [params.initCond.nodata_value] = patchInitCond
        self.locList[params.initCond.nodata_value] = (0,0)

    def table_to_stream(self, stream=sys.stdout):
        keyNames = list(next(iter(self.table.values())).keys())
        header = 'CELLID'
        for name in keyNames:
            header += ', ' + str(name)
        header += '\n'
        stream.write(header)

        errorMessage = ''
        try:
            for key,value in filter(lambda kv: kv[0]!= self.nodata_value, iter(self.table.items())):
                line = '{:.0f}'.format(key)
                for elt in keyNames:
                    try:
                        line += ', ' + '{:.5f}'.format(float(value[elt]))
                    except TypeError as error:
                        print(('Class type                = ' + str(    value[elt].__class__    ) ))
                        print(('cover                     = ' + str(    value[elt].cover        ) ))
                        print(('modelType                 = ' + str(    value[elt].modelType    ) ))
                        print(('name()                    = ' + str(    value[elt].name()       ) ))
                        print(('Class float_v1_0() fn id  = ' + str( id(value[elt].float_v1_0)  ) ))
                        print(('Class float_v2_0() fn id  = ' + str( id(value[elt].float_v2_0)  ) ))
                        print(('Class __float__()  fn id  = ' + str( id(value[elt].__float__)   ) ))
                        print(('Class floater()    fn id  = ' + str( id(value[elt].floater)     ) ))
                        print(('Class float_v1_0() fn val = ' + str(    value[elt].float_v1_0() ) ))
                        print(('Class float_v2_0() fn val = ' + str(    value[elt].float_v2_0() ) ))
                        print(('Class __float__()  fn val = ' + str(    value[elt].__float__()  ) ))
                        print(('Class floater()    fn val = ' + str(    value[elt].floater()    ) ))
                        raise error

                line += '\n'
                stream.write(line)
        except KeyError as error:
            errorMessage += 'LandscapePlus.table_to_stream(): Error: A species key does not appear to be defined. Additional info follows.\n'
            errorMessage += str(error)


        if len(errorMessage):
            errorMessage += 'LandscapePlus.table_to_stream(): Error: We\'re hosed, time to crash.\n'
            raise RuntimeError(errorMessage)

    def __getitem__(self, item):
        patchIndex = landscape.Landscape.__getitem__(self, item)
        return self.table[patchIndex]

    def update(self):
        for loc in filter( lambda key: key != self.nodata_value, iter(self.table.keys())):
            self.patchModel.update( self.locList[loc], self.table[loc], self.spModelList, self.firstyear)
        self.firstyear = 0.0

##\class DispersalModel
##\brief Computes the dispersal of each species over space
##\role{Ecology}
##\detail This class computes the dispersal kernel for each species.
class DispersalModel(object):
    dynModel = None
    def __init__(self):
        self.patchIndexLandscape = landscape.Landscape()
        self.patchDict           = dict()

    def config(self, params):
        self.patchIndexLandscape.copy(params.initCond)

        #for row in range(0, int(params.initCond.nrow) ):
           #for col in range(0, int(params.initCond.ncol) ):
                # if params.initCond.has_data_at(row,col):

        for row,col in filter( lambda rc: params.initCond.has_data_at(rc[0],rc[1]), itertools.product(   list(range(0,int(params.initCond.nrow))), list(range(0,int(params.initCond.ncol)))   ) ):
               patchIndex                 = params.initCond.data[row,col]
               patchInitCond              = copy.deepcopy(params.initCond[(row,col)])
               neighborList               = [ (nRow,nCol) for nRow,nCol in filter( lambda rc: params.initCond.has_data_at(rc[0],rc[1]), itertools.product(list(range(row-1,row+2)), list(range(col-1,col+2))) ) ] # This is slick as hell. I like python.
               neighborListExt            = [ (nRow,nCol) for nRow,nCol in filter( lambda rc: params.initCond.has_data_at(rc[0],rc[1]), itertools.product(list(range(row-2,row+3)), list(range(col-2,col+3))) ) ]
               newPatch                   = { 'spFreq':patchInitCond, 'neighborList':neighborList, 'neighborListExt': neighborListExt }
               self.patchDict[patchIndex] = newPatch

        patchInitCond = params.initCond.table.iloc[0].to_dict()
        for key in iter(patchInitCond.keys()):
            patchInitCond[key] = 0.0
        #newPatch = {'loc':(-1,-1), 'spFreq':patchInitCond, 'neighborList':[]}
        newPatch = {'spFreq':patchInitCond, 'neighborList':[]}
        self.patchDict[params.initCond.nodata_value] = newPatch

    def summary(self):
        print ('DispersalModel: Msg: Summary info start')
        print((self.patchIndexLandscape.header()))
        print(('len(patchDict) = ' + str(len(self.patchDict))))
        print((next( iter(self.patchDict.values()) )))
        patchAddress = self.patchIndexLandscape.data[241,44]
        print(('DispersalModel: patchAddress = ' + str(patchAddress)))
        print(('DispersalModel: Stuff at patchAddress = ' + str( self.patchDict[patchAddress] )))
        print(('DispersalModel: Stuff at patchAddress[\'spFreq\'] = ' + str( self.patchDict[patchAddress]['spFreq'])))
        print ('DispersalModel: Msg: Summary info end')

    def __getitem__(self, item):
        patchAddress = self.patchIndexLandscape.data[item]
        return self.patchDict[patchAddress]['spFreq']

    def has_data_at(self,row,col):
        return self.patchIndexLandscape.has_data_at(row,col)

    def compute_local_dsp(self, ret):

        #row,col      = ret['loc']
        ptr          = ret['spFreq']
        neighborListMoore = ret['neighborList']
        neighborListExt = ret['neighborListExt']

  #      for spName, spModel in iter(WetlandMorphModel.dynModel.spModelList.items()):

       # dspClass     = spModel.dspClass
       # total = 0.0

        for sp in iter(ptr.keys()):
            ptr[sp] = 0

        for sp in iter(ptr.keys()):

            dspClass = WetlandMorphModel.dynModel.spModelList[sp].dspClass
            neighborList = neighborListMoore

            if dspClass  > 1:
                neighborList = neighborListExt
            elif dspClass == 1:
                neighborList = neighborListMoore
            else:
                neighborList = []

            for neighbor in neighborList:
                patchCoverList = DispersalModel.dynModel[neighbor]
                ptr[sp] += patchCoverList[sp]

            ptr[sp] /= max(len(neighborList),1)

                #total += patchCoverList[sp]

                ##

       # for neighbor in neighborListExt:
        #    patchCoverList = DispersalModel.dynModel[neighbor]
         #   if maddie == 10:
          #      print('I am in the compute_local_dsp and the patch cover list is:')
           #     print(patchCoverList)
            #for sp,spCover in iter(patchCoverList.items()):
             #   ptr[sp] += spCover
              #  total   += spCover

        #for offsetRow in range(-1,2):
        #    for offsetCol in range(-1,2):
        #        neighborRow = row + offsetRow
        #        neighborCol = col + offsetCol
        #        if DispersalModel.dynModel.has_data_at(neighborRow, neighborCol):
        #            patchModel = DispersalModel.dynModel[neighborRow,neighborCol]
        #            for sp,spModel in patchModel.iteritems():
        #                try:
        #                    ptr[sp] += spModel.cover
        #                    total   += spModel.cover
        #                except exceptions.TypeError as error:
        #                    print 'DispersalModel::_compute_local_dsp(): Error: ret         = ' + str(ret)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: row         = ' + str(row)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: col         = ' + str(col)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: ptr         = ' + str(ptr)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: total       = ' + str(total)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: offsetRow   = ' + str(offsetRow)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: offsetCol   = ' + str(offsetCol)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: neighborRow = ' + str(neighborRow)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: neighborCol = ' + str(neighborCol)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: patch       = ' + str(patch)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: sp          = ' + str(sp)
        #                    print 'DispersalModel::_compute_local_dsp(): Error: cover       = ' + str(cover)
        #                    raise error

       # if total:
        #    for sp in iter(ptr.keys()):
         #       ptr[sp] /= total


    def update(self):
        for loc in filter( lambda key: key != self.patchIndexLandscape.nodata_value, iter(self.patchDict.keys())):
            self.compute_local_dsp(self.patchDict[loc])





##\class WetlandMorphModel
##\brief This class reads data from the wetland morph model and updates the vegetation accordingly
##\role{Ecology}
##\detail This class reads data from the wetland morph model and updates the vegetation accordingly
class WetlandMorphModel:
    dynModel = None
    def __init__(self):
        self.landWater = None

    def config(self, params):
        self.landWater = params.landWater

    def update_patch(self, newLand, spCoverList, spModelList):
        # Step 1: Figure out how much land there is currently. For 2023, morph land = vegetated area + bareground + flotant
        
        curLand = 0.0
        for spName, spModel in itertools.filterfalse( lambda kv: kv[0]=='SAV' or kv[0]=='WATER' or kv[0]=='NOTMOD' or kv[1].modelType == 'NullModel', iter(spModelList.items())):
        #for spName, spModel in itertools.filterfalse (lambda kv: kv[0] == 'SAV' or kv[0] == 'WATER' or kv[1].modelType == 'FloatingMarshModel' or kv[1].modelType == 'NullModel', iter(spModelList.items())):
            curLand += spCoverList[spName]

        # Step 1.5: Reset the bareground ages and Dead flotant. New bareground should be zeroed and added to old
        spCoverList['BAREGRND_OLD'] += spCoverList['BAREGRND_NEW']
        spCoverList['BAREGRND_NEW'] = 0.0
        spCoverList['DEAD_Flt'] = 0.0 

        # Step 2: Adjust the types based on the differences between
        # the current land cover and the new land cover provided by the land/water file.
        flotant = 0.0
        for spName, spModel in filter( lambda kv: kv[1].modelType == 'FloatingMarshModel', iter(spModelList.items())):
            flotant += spCoverList[spName]

      #  curLand = curLand+flotant

        # Step 2.1: If the current land and new land are the same, there is nothing to do.
        if curLand == newLand: return

        # Step 2.2: LAND GAIN - If there is less land in the current veg model state than indicated by the
        # land/water file, then increase the bareground and decrease the area
        # of water types. All the new land is classified as new bareground at this point.
        if curLand < newLand: # curWater > newWater

           curWater                 = 1.0 - (curLand+spCoverList['NOTMOD'])
           newWater                 = 1.0 - (newLand+spCoverList['NOTMOD'])
         
           # curWater                 = 1.0 - (curLand+flotant)
           # newWater                 = 1.0 - (newLand+flotant)         
           if curWater == 0:
               print('The current water is zero, but it is trying to reduce it due to land gain. Ignore for now.')
           else:
               deltaLand                = newLand - curLand
               spCoverList['BAREGRND_NEW'] += deltaLand
               scaleWater               = newWater/curWater # scaleWater < 1.0
               spCoverList['WATER']    *= scaleWater
               spCoverList['SAV']      *= scaleWater

         #   if giveoutput == 1:
         #       print('the cur water, newwater, and scaleWater are:')
         #       print(curWater)
         #       print(newWater)
         #       print(scaleWater)
                
          #      checkSum = 0.0
          #      for spName, spModel in itertools.filterfalse(lambda kv: kv[1].modelType == 'NullModel', iter(spModelList.items())):
          #          checkSum += spCoverList[spName]
          #      tol = 0.005
          #      if not( (1.0 - tol) < checkSum and checkSum < (1.0 + tol) ):
          #         print(('MorphMod, land gain: Error: checkSum is out of range. checkSum = ' + str(checkSum) + ' should be 1.0'))
          #         for spName, spModel in itertools.filterfalse(lambda kv: kv[1].modelType == 'NullModel', iter(spModelList.items())):
          #             print(spName)
          #             print(spCoverList[spName]) 



        # Step 2.3: LAND LOSS - If there is more land in the current veg model state than indicated by the
        # land/water file, then decrease the area of all land types and increase the amount of water.
        # The bareground is reduced first, and if that is less than the total change, the vegetated land is reduced
        else: # curLand > newLand
            deltaLand             = curLand - newLand
            flotant = 0.0
            for spName, spModel in filter( lambda kv: kv[1].modelType == 'FloatingMarshModel', iter(spModelList.items())):
               flotant += spCoverList[spName]
               
            if ((spCoverList['WATER']+spCoverList['NOTMOD'])==1):
                print('The whole cell is water and not mod, but it is trying to reduce land. Ignore for now.') 
                    
            elif spCoverList['BAREGRND_OLD'] >= deltaLand:
                    spCoverList['BAREGRND_OLD'] -= deltaLand
                    spCoverList['WATER'] += deltaLand
            else:
                if (curLand - spCoverList['BAREGRND_OLD'] - flotant)<= 0.0:
                    print('The current land minus BG old and flotant are 0, but it is trying to reduce the land. Ignore for now.')
                else:
                    scaleLand = (newLand-flotant)/(curLand - spCoverList['BAREGRND_OLD'] - flotant)# scaleLand < 1.0
                    #scaleLand = (newLand-spCoverList['NOTMOD'])/(curLand - spCoverList['BAREGRND_OLD'] - spCoverList['NOTMOD'])# scaleLand < 1.0
                    spCoverList['BAREGRND_OLD'] = 0.0
                    # because BG has been knocked out, it's ok to exclude all NullModel_Coverage types from this loop: only vegetated land should be reduced
                    for spName, spModel in itertools.filterfalse( lambda kv: kv[0]=='SAV' or kv[1].modelType == 'NullModel' or kv[1].modelType == 'NullModel_Coverage' or kv[1].modelType == 'FloatingMarshModel', iter(spModelList.items())):
                        spCoverList[spName] *= scaleLand
                    spCoverList['WATER'] += deltaLand   
                     
            #elif (curLand - spCoverList['BAREGRND_OLD'] - flotant) <= 0:

        # Done
        return

    def update(self):
        for loc in filter( lambda key: key != WetlandMorphModel.dynModel.nodata_value, iter(WetlandMorphModel.dynModel.table.keys())):
            pos = WetlandMorphModel.dynModel.locList[loc]
       
            newLand = self.landWater[pos]/100.0
            if newLand < 0: continue
           # if loc == 154037:
           #     giveoutput = 1
           # else:
           #     giveoutput = 0
            self.update_patch(newLand, WetlandMorphModel.dynModel.table[loc], WetlandMorphModel.dynModel.spModelList)


    def act(self):
        self.update()


##\class ModelUpdateEvent
##\brief An event class to update the dynamics model and the dispersal model.
##\role{Machinery}
##\detail This class will likely be removed in an upcoming version of the model.
#
class ModelUpdateEvent(event.Event):
    def __init__(self, time, name, model):
        event.Event.__init__(self, time, name)
        self.model = model

    def act(self):
        print((self.name))
        self.model.update()


class CloseStreamEvent(event.Event):
    def __init__(self, time, name, stream):
        event.Event.__init__(self, time, name)
        self.stream = stream

    def act(self):
        self.stream.close();

##\class StopWatch
##\brief A class to measure elapsed wall clock time.
##\role{Machinery}
##\detail This class is used to measure how much time
# (wall clock time) is used for different aspects of the model.
# This provides some rough profiling information so we can estimate
# how long future runs of the model might take.
class StopWatch:
    def __init__(self):
        self.running  = 0
        self._start   = 0
        self._end     = 0

    def start(self):
        self.running = 1
        self._start = time.time()
        self._stop  = 0

    def stop(self):
        self._end   = time.time()
        self.running = 0
        print(('StopWatch: Msg: Delta t = ' + str(self._end - self._start)))

    def act(self):
        if self.running:
            self.stop()
        else:
            self.start()

    def __str__(self):
        ret  = str(id(self))     + ' '
        ret += str(self.running) + ' '
        ret += str(self._start)  + ' '
        ret += str(self._end)    + ' '
        ret += str(self._end - self._start)
        return ret

#class StopWatchEvent(event.Event):
    #def __init__(self, time, name='StopWatchEvent', stopwatch=None):
        event.Event.__init__(self, time, name)
        self.stopwatch = stopwatch
#
    #def act(self):
        #if self.stopwatch.running:
            #self.stopwatch.stop()
        #else:
            #self.stopwatch.start()


# This is the top level class for the LAVegMod. To get the model going you need
# to do four things.
#  1) Include the model module in the Python file that will be used
#  to coordinate the running of the various models.
#  The include statement should look something like:
#  import model
#  Note that this assumes that the LAVegMod code is in the same
#  directory as the toplevel script.
#
#  2) Instantiate an object of class Model
#  This will create a copy of the model and all of its subcomponents.
#  Instantiating model should look something like:
#
#  laVegMod = model.Model()
#
#
#  3) Call Model.config(self, argv) for the instantated object
#  Model.config takes a single (non-self) argument that contains the
#  name of the configuration file for the model. The call to Model.config()
#  should look something like:
#
#  laVegMod.config( <configFilename> )
#
#  where <configFilename> is a Python string containing the full path
#  and filename of the configuration file for the LAVegMod.
#
#  The call to Model.config() will bring the model up to a state where it
#  is ready to run. It is possible that more or more errors may occur while
#  configuring the model. You can capture these and processes them yourself
#  or you can just let them go and they will halt the code. To capture the
#  errors from the LAVegMod you need code that looks something like:
#
#  try:
#      laVegMod.config( <configFilename> )
#  except exception.RunTimeError as error:
#      <do something with the error>
#
#  4) Call Model.step() repeatedly for the instantated object.
#  Each time you call model.step() the model will advance by one year.
#
#  The call to Model.step() should look something like:
#
#  laVegMod.step()
#
#  Again, this code may throw exceptions if something goes wrong. You can
#  capture these errors using code that looks like this:
#
#  try:
#      laVegMod.step()
#  except exception.RunTimeError as error:
#      <do something with the error>
#
#  Note that the LAVegMod only throws, (or should only throw)
#  exception.RunTimeError()
#
#  The thrown error contains one or error messages stored in a string that describe what
#  went wrong. You can print these with something like:
#
#  print error
#
#  If the model has thrown an exception, then something is broken and
#  there is no way to fix it during run time. You or I, or someone, will
#  have to chase down the error and fix it. Often the error is going to be
#  in the input files. So I would check those first.
#
#  I have tried to make the error checking and reporting as comprehensive
#  and as detailed as possible.
#
#

class Model(object):
    def __init__(self):
        self.currentTime           = event.Time(0,0)
        self.params                = Params()
        self.eventQueue            = event.EventQueue()
        self.dynModel              = DynamicsModel()
        self.dspModel              = DispersalModel()
        self.welModel              = WetlandMorphModel()
        self.plantingModel         = plantingmodel.PlantingModel()

        SpeciesModel.params        = self.params    # I'm not sure I like this.
        SpeciesModel.dspModel      = self.dspModel  #
        DispersalModel.dynModel    = self.dynModel  #
        WetlandMorphModel.dynModel = self.dynModel
        plantingmodel.PlantingModel.dynModel     = self.dynModel

    def config(self, argv):
        print ('Model: Msg: This is LAVegMod V3 (model_v3.py)')

        print ('Model: Msg: Reading configuration information')
        self.params.config(argv)

        print ('Model: Msg: Configuring the dynamics model')
        self.dynModel.config(self.params)

        print ('Model: Msg: Configuring the dispersal model')
        self.dspModel.config(self.params)
        # self.dspModel.summary()

        print ('Model: Msg: Configuring the wetland morph process')
        self.welModel.config(self.params)

        print ('Model: Msg: Configuring the plantings model')
        self.plantingModel.config(self.params)

        print ('Model: Msg: Configuring the update queue')
        self.eventQueue.clear()
        perYearSW = StopWatch()
        totalSW   = StopWatch()

        self.eventQueue.add_event(              event.GenericEvent(event.Time(self.params.startYear,    0) , name='StopWatch', callable=totalSW) )
    #   self.eventQueue.add_event(              event.MsgEvent(event.Time(self.params.startYear,  100),                                                                       stream=self.params.outputStrm['OutputFile'] , msg='# Year = ' + str(self.params.startYear))    )
    #    self.eventQueue.add_event(landscape.WriteASCIIGridPlus(event.Time(self.params.startYear,  200), name='WriteASCIIGrid: Msg: Writing model output',                     stream=self.params.outputStrm['OutputFile'],              landscape=self.dynModel))
        for year in range(self.params.startYear+1, self.params.endYear+1):
            self.eventQueue.add_event(          event.GenericEvent(event.Time(year,  100), name='StopWatch', callable=perYearSW ) )
            self.eventQueue.add_event(         event.MsgEvent(event.Time(year,  200), msg='Year = ' + str(year))    )
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  300), name='ReadASCIIGrid: Msg: Reading wave amp data',                     stream=self.params.inputStrm['WaveAmplitudeFile'],        landscape=self.params.waveAmp))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  400), name='ReadASCIIGrid: Msg: Reading mean salinity data',                stream=self.params.inputStrm['MeanSalinityFile'],         landscape=self.params.meanSal))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  500), name='ReadASCIIGrid: Msg: Reading summer mean water depth data',      stream=self.params.inputStrm['SummerMeanWaterDepthFile'], landscape=self.params.smDepth))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  600), name='ReadASCIIGrid: Msg: Reading summer mean salinity data',         stream=self.params.inputStrm['SummerMeanSalinityFile'],   landscape=self.params.smSal))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  700), name='ReadASCIIGrid: Msg: Reading summer mean temperature data',      stream=self.params.inputStrm['SummerMeanTempFile'],       landscape=self.params.smTemp))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  800), name='ReadASCIIGrid: Msg: Reading water height above ground data',    stream=self.params.inputStrm['HeightAboveWaterFile'],     landscape=self.params.heightAboveWater))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  850), name='ReadASCIIGrid: Msg: Reading bi water height above ground data', stream=self.params.inputStrm['BarrierIslandHeightAboveWaterFile'],     landscape=self.params.biHeightAboveWater))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,  900), name='ReadASCIIGrid: Msg: Reading acute salinity stress data',       stream=self.params.inputStrm['AcuteSalinityStressFile'],          landscape=self.params.acuteSal))
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year, 1000), name='ReadASCIIGrid: Msg: Reading tree establishment data',           stream=self.params.inputStrm['TreeEstCondFile'],          landscape=self.params.treeEstCond))
            self.eventQueue.add_event(       ModelUpdateEvent(event.Time(year, 1100), name='ModelUpdateEvent: Msg: Updating veg dynamics',                  model=self.dynModel  )  )
            self.eventQueue.add_event(       ModelUpdateEvent(event.Time(year, 1200), name='ModelUpdateEvent: Msg: Updating dispersal model',               model=self.dspModel  )  )
           # self.eventQueue.add_event(         event.MsgEvent(event.Time(year, 1300),                                                                       stream=self.params.outputStrm['OutputFile'] , msg='# Year = ' + str(year))    )
           # self.eventQueue.add_event(landscape.WriteASCIIGridPlus(event.Time(year, 1400), name='WriteASCIIGrid: Msg: Writing model output',                stream=self.params.outputStrm['OutputFile'],              landscape=self.dynModel))
            self.eventQueue.add_event(         event.GenericEvent(event.Time(year, 1600), name='StopWatch', callable=perYearSW ) )
            self.eventQueue.add_event(       event.PauseEvent(event.Time(year, 2000)))

        for year in self.params.wetlandMorphYears:
            self.eventQueue.add_event(landscape.ReadASCIIGrid(event.Time(year,450), name='ReadASCIIGrid: Msg: Reading wetland morph data',                  stream=self.params.inputStrm['WetlandMorphLandWaterFile'], landscape=self.params.landWater ))
            self.eventQueue.add_event(     event.GenericEvent(event.Time(year,1050), name='ModelUpdateEvent: Msg: Adding the effects from the wetland morph data', callable=self.welModel  ))

        for year, submodel in iter(self.plantingModel.eventDict.items()):
            self.eventQueue.add_event( ModelUpdateEvent(event.Time(year, 1075), name='ModelUpdateEvent: Msg: Processing plantings for year ' + str(year),  model=submodel ))

        for yearKey,yearStream in filter( lambda kv: re.match(r'^__Single_',kv[0]) != None, iter(self.params.outputStrm.items()) ):
            year = int( yearKey.replace('__Single_','') )
            self.eventQueue.add_event(landscape.WriteASCIIGridPlus(event.Time(year, 1500), name='WriteASCIIGrid: Msg: Writing model output for year  ' + str(year), stream=yearStream, landscape=self.dynModel))
            self.eventQueue.add_event(CloseStreamEvent            (event.Time(year, 1501), name='WriteASCIIGrid: Msg: Closing stream for output year ' + str(year), stream=yearStream ))

        self.eventQueue.add_event(             event.GenericEvent(event.Time(self.params.endYear,3000), name='StopWatch', callable=totalSW) )

    def step(self):
        self.eventQueue.run(while_condition=(lambda arg: arg.name != 'PauseEvent'))
        return 0

    def run(self):
        try:
            self.config(sys.argv)
        except RuntimeError as error:
            print (error)
            return 1
        #except:
        #    print 'Model: Error: Caught an unknown error : ' , sys.exc_info()[0]
        #    return 1

        self.eventQueue.run()

        self.params.done()

        return 0;
