#!/usr/bin/env python

##\file
# This file defined the base Function class.
# This class is used to provide a piece-wise
# linear function for use in the model. The
# primary use of this these classes is presenting
# the information from the species establishment
# and senescenes tables.
#

# Standard Python Modules
from builtins import Exception as exceptions
import itertools
import math
import numpy
import scipy
import scipy.interpolate as interp
import sys
import xlrd

# Third Party Modules
import pandas

class Function(object):
    xValue = list()
    yValue = list()

    def __init__(self, xValue, yValue):
        if len(xValue) != len(yValue):
            msg  = 'Function: Error: Number of elements in xValue an yValue must match'
            msg += ' len(xValue) = ' + str(len(xValue)) + ' len(yValue) = ' + str(len(yValue))
            raise RuntimeError(msg)

        self.xValue = xValue
        self.yValue = yValue

    def toIndex(self, x, vec):
        if ( x < vec[0] or vec[-1] < x ):
            msg  = 'Function: Error: x out of range\n'
            msg += 'x = ' + str(x)             + '\n'
            msg += 'vec[0]  = ' + str(vec[0])  + '\n'
            msg += 'vec[-1] = ' + str(vec[-1]) + '\n'
            raise RuntimeError(msg)

        b = 0                          # b = the "bottom" index
        t = len(vec)-1                 # t = the "top" index
        m = int(scipy.ceil((b+t)/2.0)) # m = the middle index

        if ( x == vec[t] ):
            b     = t-1
            scale = 1.0
            return b,scale

        while ( vec[b+1] <= x ):
            if ( vec[b] <= x and x < vec[m] ):
                t = m
            else:
                b = m
            m = int(scipy.ceil((b+t)/2.0))

        scale = (x - vec[b])/(vec[b+1] - vec[b])

        return b,scale

    def __getitem__(self, x):
        b,scale = self.toIndex(x, self.xValue)
        return self.yValue[b]*(1.0 - scale) + self.yValue[b+1]*scale


class ReadVegTable1D(object):
    def __init__(self):
        pass

    def read(self, filename, modelType, species):
        try:
            data = pandas.read_excel(filename, modelType, index_col=None )
            xValue = list( data.ix[:,0] )
            yValue = list( data[species] )
            return {'elvValue':xValue, 'rate':yValue}
        except IOError as error:
            errorMessage  = 'ReadVegTable1D: Error: Could not open file for reading : ' + filename + '\n'
            errorMessage += 'ReadVegTable1D: Error: Addition error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage)
        except xlrd.biffh.XLRDError as error:
            errorMessage  = 'ReadVegTable1D: Error: Could not find sheet named ' + modelType + ' in file ' + filename + '\n'
            errorMessage += 'ReadVegTable1D: Error: Additional error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage)
        except KeyError as error:
            errorMessage  = 'ReadVegTable1D: Error: Could not find species named ' + species + ' in sheet ' + modelType + ' in file ' + filename +'\n'
            errorMessage += 'ReadVegTable1D: Error: Additional error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage)

class Function2DFast(object):

    def __init__(self, xValue, yValue, data):
        # xValue = wave amplitudes
        # yValue = mean salinity
        # data   = senescence and growth rates
        if ( len(yValue) != numpy.size(data, 0) ):
            msg  = 'Function2DFast: Error: number of elements in yValue must match number of rows in data'
            msg += ' len(yValue) = ' + str(len(yValue)) + ' data.shape = ' + str(data.shape)
            raise RuntimeError(msg)

        if ( len(xValue) != numpy.size(data, 1) ):
            msg  = 'Function2DFast: Error: number of elements in xValue must match number of columns in data'
            msg += ' len(xValue) = ' + str(len(Value)) + ' size(data) = ' + str(data.shape)
            raise RuntimeError(msg)

        self.xValue = xValue
        self.yValue = yValue
        self.data   = data
        self.verbose = False

    def toIndex(self, x, vec):
        if ( x < vec[0] or vec[-1] < x ):
            msg  = 'Function2DFast: Error: x out of range\n'
            msg += 'x = ' + str(x)             + '\n'
            msg += 'vec[0]  = ' + str(vec[0])  + '\n'
            msg += 'vec[-1] = ' + str(vec[-1]) + '\n'
            raise RuntimeError(msg)

        b = 0                          # b = the "bottom" index
        t = len(vec)-1                 # t = the "top" index
        m = int(scipy.ceil((b+t)/2.0)) # m = the middle index

        if ( x == vec[t] ):
            b     = t-1
            scale = 1.0
            return b,scale

        while ( vec[b+1] <= x ):
            if ( vec[b] <= x and x < vec[m] ):
                t = m
            else:
                b = m
            m = int(scipy.ceil((b+t)/2.0))

        scale = (x - vec[b])/(vec[b+1] - vec[b])

        return b,scale

    def __getitem__(self, item):
        i, i_scale = self.toIndex(item[0], self.xValue)
        j, j_scale = self.toIndex(item[1], self.yValue)

        x0 = self.data[j  ,i] * (1.0 - i_scale) + self.data[j  ,i+1]*i_scale
        x1 = self.data[j+1,i] * (1.0 - i_scale) + self.data[j+1,i+1]*i_scale

        ret = x0*(1.0 - j_scale) + x1*j_scale

        return(ret)

class Function2DFaster(Function2DFast):
    def special_round(self,x):
        xprime = math.floor(x)
        if xprime == x:
            return xprime-1
        return xprime

    def compute_mult(self, valueList):
        base      = 10
        mult      = 1

        for p in range(5):
            mult     = base**p
            multList = set( [ int(value * mult) for value in valueList ] )
            if len(multList) == len(valueList):
                break

        return(mult)

    def compute_index_table(self, xList, mult):
        minX       = min(xList)
        maxX       = max(xList)
        minXprime  = int(minX*mult)
        maxXprime  = int(maxX*mult)

        f                 = Function( xList, range(0,len(xList))   )
        xPrimeList        = range(minXprime, maxXprime+1, 1)
        xPrimeToIndexList = numpy.zeros( maxXprime + 1, dtype=numpy.int32 )

        for xPrime in xPrimeList[1:]:
            xPrimeToIndexList[xPrime] = self.special_round( f[float(xPrime)/mult] )

        return( (xPrimeToIndexList, xList) )

    def __init__(self, xList, yList, data):
        multX = self.compute_mult(xList)
        multY = self.compute_mult(yList)

        self.mult   = max([multX, multY])
        self.xValue = self.compute_index_table(xList, self.mult)
        self.yValue = self.compute_index_table(yList, self.mult)

        self.data = data

    def toIndex(self, x, vec):
        xPrimeToIndexList = vec[0]
        xList             = vec[1]

        xPrime = int(x*self.mult)
        i      = xPrimeToIndexList[xPrime]
        i      = i + (xList[i+1] < x)
        xlow   = xList[i]
        xhigh  = xList[i+1]
        scale  = (x - xlow)/(xhigh - xlow)

        return i,scale


class Function2D(object):

    def __init__(self, xValue, yValue, data = 0):
        # xValue = wave amplitudes
        # yValue = mean salinity
        # data   = senescence and growth rates
        if len(yValue) != numpy.size(data, 0):
            msg  = 'Function2D: Error: number of elements in yValue must match number of rows in data'
            msg += ' len(yValue) = ' + str(len(yValue)) + ' data.shape = ' + str(data.shape)
            raise RuntimeError(msg)

        if len(xValue) != numpy.size(data, 1):
            msg  = 'Function2D: Error: number of elements in xValue must match number of columns in data'
            msg += ' len(xValue) = ' + str(len(Value)) + ' size(data) = ' + str(data.shape)
            raise RuntimeError(msg)

        self.point = numpy.array( [ [x,y] for y in yValue for x in xValue ] )
        self.value = numpy.array( [elt for row in data for elt in row] )
        self.minX  = min( xValue )
        self.maxX  = max( xValue )
        self.minY  = min( yValue )
        self.maxY  = max( yValue )

    def __getitem__(self, item):
        if ( item[0] < self.minX or self.maxX < item[0] ):
            raise RuntimeError('Function2D: Error: x coordinate out of range. ' + str(item[0]) + ' not in ' + '[ ' + str(self.minX) + ', ' + str(self.maxX) + ' ]'    )

        if ( item[1] < self.minY or self.maxY < item[1] ):
            raise RuntimeError('Function2D: Error: y coordinate out of range. ' + str(item[1]) + ' not in ' + '[ ' + str(self.minY) + ', ' + str(self.maxY) + ' ]'    )

        return interp.griddata(self.point, self.value, ([item[0]],[item[1]] ), method='linear')[0]

class ReadVegTable2D(object):
    def __init__(self):
        pass

    def read(self, filename, species):
        try:
            data       = pandas.read_excel(filename, species, index_col=None, parse_cols=range(1,23), skiprows=1)
            waValue    = data.columns.values.tolist()[1:len(data.columns)]
            spcode     = data.columns.values.tolist()[0]
            salValue   = data[spcode].tolist()
            rate       = numpy.asarray(data)[:,1:22]
            return {'spcode':spcode, 'waValue':waValue, 'salValue':salValue, 'rate':rate }
        except IOError as error:
            errorMessage  = 'ReadVegTable2D: Error: Could not open file for reading : ' + filename + '\n'
            errorMessage += 'ReadVegTable2D: Error: Addition error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage)
        except xlrd.biffh.XLRDError as error:
            errorMessage  = 'ReadVegTable2D: Error: Could not find sheet named ' + species +  ' in file : ' + filename + '\n'
            errorMessage += 'ReadVegTable2D: Error: Addition error info : ' + str(error) + '\n'
            raise RuntimeError(errorMessage)


def check():
    reader = ReadVegTable2D()
    data   = reader.read('./S11_G01/LaVegMod2_Establishment_Tables_JMV.xlsx','AVGE')
    f      = Function2DFaster(data['waValue'], data['salValue'], data['rate'])

    salList = numpy.arange(0.0, 45.0, 0.1)
    waList  = numpy.arange(0.0, 0.8, 0.01)

    grid    = numpy.zeros( shape=(len(salList), len(waList) ))

    for i,sal in itertools.izip(range(len(salList)), salList):
        for j,wa in itertools.izip(range(len(waList)), waList):
            grid[i,j] = f[wa,sal]

    numpy.savetxt('temp.csv',grid,delimiter=',')

if __name__ == "__main__":
    exit( check() )
