#!/usr/bin/env python

# STD Python modules
import copy
import exceptions
import itertools
import numpy
import re
import StringIO
import sys

# Third party modules
import pandas

# Model modules
import event

class NoDataValueException(exceptions.RuntimeError):
    def __init__(self):
        exceptions.RuntimeError.__init__(self)

class Landscape(object):
    def __init__(self, nrow=0, ncol=0, yllcorner=0, xllcorner=0, cellsize=0, nodata_value=-9999):
        self.nrow         = copy.copy(nrow)
        self.ncol         = copy.copy(ncol)
        self.yllcorner    = copy.copy(yllcorner)
        self.xllcorner    = copy.copy(xllcorner)
        self.cellsize     = copy.copy(cellsize)
        self.nodata_value = copy.copy(nodata_value)
        self.data         = numpy.zeros( shape=(nrow, ncol) )

    def __getitem__(self, item):
        # row = item[0]
        # col = item[1]
        return self.data[item[0],item[1]];

    def __setitem__(self, key, value):
        self.data[key[0], key[1]] = value

    def has_data_at(self, row, col):
        return 0 <= row and row < self.nrow and 0 <= col and col < self.ncol and self.data[row,col] != self.nodata_value

    def resize(self, nrow=0, ncol=0, yllcorner=0, xllcorner=0, cellsize=0, nodata_value=-9999):
        self.nrow         = copy.copy(nrow)
        self.ncol         = copy.copy(ncol)
        self.yllcorner    = copy.copy(yllcorner)
        self.xllcorner    = copy.copy(xllcorner)
        self.cellsize     = copy.copy(cellsize)
        self.nodata_value = copy.copy(nodata_value)
        self.data         = numpy.zeros( shape=(nrow, ncol ));

    def copy(self,orig):
        self.nrow         = copy.copy(orig.nrow);
        self.ncol         = copy.copy(orig.ncol);
        self.yllcorner    = copy.copy(orig.yllcorner);
        self.xllcorner    = copy.copy(orig.xllcorner);
        self.cellsize     = copy.copy(orig.cellsize);
        self.nodata_value = copy.copy(orig.nodata_value);
        self.data         = copy.deepcopy(orig.data);

    def header_to_stream(self, stream = sys.stdout):
        stream.write('nrows '        + '{:.0f}'.format(self.nrow)  + '\n')
        stream.write('ncols '        + '{:.0f}'.format(self.ncol)  + '\n')
        stream.write('yllcorner '    + str(self.yllcorner)         + '\n')
        stream.write('xllcorner '    + str(self.xllcorner)         + '\n')
        stream.write('cellsize '     + '{:.0f}'.format(self.cellsize) + '\n')
        stream.write('nodata_value ' + str(self.nodata_value)      + '\n')

    def data_to_stream(self, stream = sys.stdout, dataFormat = '{}'):
        for row,col in itertools.product(range(0,int(self.nrow)), range(0,int(self.ncol)) ):
                stream.write(dataFormat.format(self.data[row,col]))
                stream.write( [' ','\n'][col+1 == self.ncol])

    def copy_to_stream(self, stream = sys.stdout, dataFormat='{}'):
        self.header_to_stream(stream)
        self.data_to_stream  (stream, dataFormat)


    def __str__(self):
        ret = StringIO.StringIO()
        self.copy_to_stream(ret)
        return ret

class LandscapePlus(Landscape):
    def __init__(self, nrow=0, ncol=0, yllcorner=0, xllcorner=0, cellsize=0, nodata_value=-9999):
        Landscape.__init__(self, nrow, ncol, yllcorner, xllcorner, cellsize, nodata_value)
        self.table = dict()

    def __getitem__(self, item):
        address = Landscape.__getitem__(self, item)
        if address == self.nodata_value:
            raise NoDataValueException()

        try:
            ret = self.table[ self.table['CELLID'] == address].to_dict('record')[0]
            del ret['CELLID']
            return(ret)
        except exceptions.KeyError as error:
            msg = 'LandscapePlus: Error: CELLID column not defined in initial conditions. This is odd\n'
            msg += str(error)
            raise exceptions.RuntimeError(msg)

    def resize(self, nrow=0, ncol=0, yllcorner=0, xllcorner=0, cellsize=0, nodata_value=-9999):
        Landscape.resize(self, nrow, ncol, yllcorner, xllcorner, cellsize, nodata_value)

    def copy(self,orig):
        Landscape.copy(self, orig)
        self.table = copy.deepcopy(orig.table)

    def __str__(self):
        return Landscape.__str__(self)

    def table_to_stream(self, stream=sys.stdout):
        keyNames = self.table.itervalues().next().keys()

        header = 'CELLID'
        for name in keyNames:
            header += ', ' + str(name)
        header += '\n'
        stream.write(header)

        errorMessage = ''
        try:
            for key,value in itertools.ifilter(lambda (k,v): k!= self.nodata_value, self.table.iteritems()):
                line = '{:.0f}'.format(key)
                for elt in keyNames:
                    try:
                        line += ', ' + '{:.5f}'.format(float(value[elt]))
                    except exceptions.TypeError as error:
                        print('Class type                = ' + str(    value[elt].__class__    ) )
                        print('cover                     = ' + str(    value[elt].cover        ) )
                        print('modelType                 = ' + str(    value[elt].modelType    ) )
                        print('name()                    = ' + str(    value[elt].name()       ) )
                        print('Class float_v1_0() fn id  = ' + str( id(value[elt].float_v1_0)  ) )
                        print('Class float_v2_0() fn id  = ' + str( id(value[elt].float_v2_0)  ) )
                        print('Class __float__()  fn id  = ' + str( id(value[elt].__float__)   ) )
                        print('Class floater()    fn id  = ' + str( id(value[elt].floater)     ) )
                        print('Class float_v1_0() fn val = ' + str(    value[elt].float_v1_0() ) )
                        print('Class float_v2_0() fn val = ' + str(    value[elt].float_v2_0() ) )
                        print('Class __float__()  fn val = ' + str(    value[elt].__float__()  ) )
                        print('Class floater()    fn val = ' + str(    value[elt].floater()    ) )
                        raise error

                line += '\n'
                stream.write(line)
        except exceptions.KeyError as error:
            errorMessage += 'LandscapePlus.table_to_stream(): Error: A species key does not appear to be defined. Additional info follows.\n'
            errorMessage += str(error)


        if len(errorMessage):
            errorMessage += 'LandscapePlus.table_to_stream(): Error: We\'re hosed, time to crash.\n'
            raise exceptions.RuntimeError(errorMessage)

    def copy_to_stream(self, stream = sys.stdout, dataFormat='{}'):
        Landscape.copy_to_stream(self, stream, dataFormat='{:.0f}')
        self.table_to_stream(stream)

    def extract_layer(self, layerName):
        try:
            ret = Landscape()
            ret.copy(self)
            for row,col in itertools.ifilter( lambda (r,c): self.has_data_at(r,c), itertools.product(   range(int(ret.nrow)), range(int(ret.ncol))    ) ):
                ret[row,col] = self[row,col][layerName]
        except exceptions.KeyError as error:
            errorMessage = 'LandscapePlus.extract_layer(): Error: Requested layer name does not exist: ' + str(error) + '\n'
            raise exceptions.RuntimeError(errorMessage)

        return ret

class ReadASCIIGrid(event.Event):
    def __init__(self, time=event.Time(0,0), name='ReadASCIIGrid', stream=None, landscape=None):
        event.Event.__init__(self, time, name);
        self.stream    = stream
        self.landscape = landscape

    def read(self, stream, landscape ):
        errorMessage = ''
        symbols      = ['nrows', 'ncols', 'yllcorner', 'xllcorner', 'cellsize', 'nodata_value']
        metadata     = dict()

        for line in stream:
            line = re.sub(r'#.*', '', line)
            line = re.sub(r'\t', ' ', line)
            line = re.sub(r'^ ', '',  line)
            line = re.sub(r' $', '',  line)
            line = re.sub(r'  ', ' ', line)

            if re.match(r'^ *$', line) :
                continue

            key, value      = line.split(' ')
            key             = key.lower()

            metadata[key]   = float(value)

            try:
                symbols.remove(key)
            except exceptions.ValueError as error:
                errorMessage += 'ReadASCIIGrid: Error: Got an unknown or repeated key from metadata: ' + str(key) +'\n'
                errorMessage += str(error) + '\n'

            if len(symbols) == 0:
                break

            if len(metadata) > 6:
                errorMessage += 'ReadASCIIGrid: Error: Expected six lines of metadata but got more\n'

        if len(errorMessage) != 0:
            raise exceptions.RuntimeError(errorMessage);

        landscape.resize(metadata['nrows'], metadata['ncols'], metadata['yllcorner'], metadata['xllcorner'], metadata['cellsize'], metadata['nodata_value'] )

        count   = metadata['nrows']
        strStrm = StringIO.StringIO()
        for line in stream:
            strStrm.write(line)
            count -= 1
            if not(count): break

        strStrm.seek(0,0)
        landscape.data = numpy.loadtxt(strStrm, delimiter=' ')

        #row = 0
        #for line in stream:
        #    line = re.sub(r'\t', ' ',  line)
        #    line = re.sub(r'^ ',  '',  line)
        #    line = re.sub(r' $',  '',  line)
        #    line = re.sub(r'  ', ' ',  line)

        #    rowData = line.split(' ')
        #    landscape.data[row,] = [ float(value) for value in rowData ]

        #    row += 1
        #    if row == landscape.nrow:
        #        break

    def act(self):
        print self.name
        self.read(self.stream, self.landscape)

class ReadASCIIGridPlus(event.Event):
    def __init__(self, time = 0, priority = 100, stream = None, landscape=None):
        event.Event.__init__(self, time, priority)
        self.stream    = stream
        self.landscape = landscape

    def read(self, strm, landscape):
        reader          = ReadASCIIGrid()
        reader.read(strm, landscape)

        strStrm = StringIO.StringIO()
        for line in strm:
            strStrm.write(line)

        strStrm.seek(0)
        landscape.table = pandas.io.parsers.read_csv(strStrm, skipinitialspace=True)

        #print('ReadASCIIGridPlus.read(): Msg: repackaging the table data')
        #errorMessage = ''
        #try:
        #    for r in range(len(table)):
        #        cellID = table.ix[r]['CELLID']
        #        landscape.table[cellID] = table.ix[r].to_dict()
        #        del landscape.table[cellID]['CELLID']
        #except exceptions.KeyError as error:
        #    errorMessage += 'ReadASCIIGridPlus.read(): Error: CELLID is not defined.\n'
        #    errorMessage += str(error)
        #    raise exceptions.RuntimeError(errorMessage)


    def act(self):
        self.read(self.stream, self.landscape)

class WriteASCIIGrid(event.Event):
    def __init__(self, time = event.Time(0, 0), name='WriteAsciiGrid', stream=None, landscape=None, dataFormat='{}'):
        event.Event.__init__(self, time, name)
        self.stream     = stream
        self.landscape  = landscape
        self.dataFormat = dataFormat

    def write(self, stream, landscape, dataFormat='{}'):
        landscape.copy_to_stream(stream, dataFormat)

    def act(self):
        self.write(self.stream, self.landscape, self.dataFormat)

class WriteASCIIGridPlus(WriteASCIIGrid):
    def __init__(self, time=(0,0), name = 'WriteAsciiGridPlus', stream = None, landscape = None, dataFormat='{:.0f}'):
        WriteASCIIGrid.__init__(self, time, name, stream, landscape, dataFormat)
