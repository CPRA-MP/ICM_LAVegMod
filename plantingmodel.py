
from builtins import Exception

import itertools
import re

#from osgeo import ogr


class LAVegModNoneValueError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class LAVegModNoMatchError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PlantingModel:
    dynModel = None

    def __init__(self):
        self.eventDict    = dict()
        self.plantingList = list()

    def __del__(self):
        #self.outSource.Destroy()
        pass

    def setup_shapefile(self, params):
        print ('PlantingModel: Msg: Setting up intersection shapefile')

        outFilename   = './Intersection/intersections.shp'
        self.driver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(outFilename):
            self.driver.DeleteDataSource(outFilename)

        source = next(iter(params.plantingDict.values()))

        self.outSource  = self.driver.CreateDataSource(outFilename)
        outSrs          = source.GetLayer().GetSpatialRef()
        self.outLayer   = self.outSource.CreateLayer('intersections', outSrs, ogr.wkbPolygon)
        outFieldName    = ogr.FieldDefn("ID", ogr.OFTInteger)
        outFieldName.SetWidth(10)
        self.outLayer.CreateField(outFieldName)

    def get_bbox(self, geom):
        line   = geom.GetBoundary().ExportToWkt()
        line   = re.sub(r'(MULTI){0,1}LINESTRING ', '' , line)
        line   = re.sub('\('         , '' , line)
        line   = re.sub('\)'         , '' , line)
        line   = re.sub(','          , ' ', line)
        line   = [ float(x) for x in re.split(' ', line)]
        lat    = line[1::2]
        lon    = line[0::2]
        return max(lat), min(lat), max(lon), min(lon)

    def get_location_list(self, geom, map):
        ret                      = list()

        #print 'PlantingModel: Msg: Getting bounding box'
        north, south, east, west = self.get_bbox(geom)

        #print 'PlantingModel: Msg: Computing start/end row/col'
        northing = map.yllcorner + map.nrow * map.cellsize
        colStart = int( max(0,        math.floor( (west  - map.xllcorner)/map.cellsize ) ) )
        colEnd   = int( min(map.ncol, math.ceil ( (east  - map.xllcorner)/map.cellsize ) ) )
        rowStart = int( max(0,        math.floor( (northing - north)/map.cellsize ) ) )
        rowEnd   = int( min(map.nrow, math.ceil ( (northing - south)/map.cellsize ) ) )

        #print 'PlantingModel: Msg: Finding intersections'
        count = 0
        for row,col in itertools.product( list(range(rowStart,rowEnd)), list(range(colStart,colEnd)) ):
            boxN = northing - row * map.cellsize
            boxS = northing - (row+1) * map.cellsize
            boxE = map.xllcorner + (col+1) * map.cellsize
            boxW = map.xllcorner + col     * map.cellsize

            cell = ogr.Geometry( ogr.wkbLinearRing )
            cell.AddPoint(boxW, boxS)
            cell.AddPoint(boxE, boxS)
            cell.AddPoint(boxE, boxN)
            cell.AddPoint(boxW, boxN)
            cell.AddPoint(boxW, boxS)
            cellPoly = ogr.Geometry(ogr.wkbPolygon)
            cellPoly.AddGeometry(cell)

            if not geom.Intersect(cellPoly):
                continue

            inter = geom.Intersection(cellPoly)

            # If the total area to be adjusted is less that 10x10 meter square,
            # then don't bother with this cell.
            if inter.Area() * 500 * 500 < 100:
                continue

            frac = inter.Area()/cellPoly.Area()
            ret.append((row,col,frac))


            #outFeature = ogr.Feature(self.outLayer.GetLayerDefn())
            #outFeature.SetField("ID", count)
            #count += 1
            #outFeature.SetGeometry(cellPoly)
            #self.outLayer.CreateFeature(outFeature)
            #outFeature = None

            # Note that this is the appropriate way to cleanup the memory allocated by the obj.Geometry() calls.
            # When the reference is reset using = None, python will clean up the memory allocated for these objects.
            inter    = None
            cellPoly = None
            cell     = None


        return ret

    def get_field(self, feature, label, fullPattern, parsePattern ):

        try:
            line = feature.GetField(label) # If the label does not exist, this function throws ValueError
        except ValueError as error:
            raise ValueError('PlantingModel: Warning: Label is not present in shapefile. label = ' + label)

        if line == None: # The field exists, but does not have a value
            raise LAVegModNoneValueError('PlantingModel: Warning: No value associated with field named ' + label)

        ret  = re.match(fullPattern, line)

        if ret == None:
            errMsg  = 'PlantingModel: Warning: Field value does not match expected form. Details follow\n'
            errMsg += 'PlantingModel: Warning: Field code = ' + label     + '\n'
            errMsg += 'PlantingModel: Warning: Value      = ' + line      + '\n'
            errMsg += 'PlantingModel: Warning: regex      = ' + fullPattern
            raise LAVegModNoMatchError( errMsg )

        return re.finditer(parsePattern, line)

    def config_with_planting_list(self, allYearPlantingList ):
        print ('PlantingModel: Msg: Building submodels for each planting year.')
        for planting in allYearPlantingList:
            year = planting[0]
            if year not in self.eventDict:
                print(('Adding planting submodel for year ' + str(year)))
                self.eventDict[year] = PlantingModel()
            self.eventDict[year].plantingList.append( planting )

    def config_with_txt(self, params):
        if params.plantingsStrm == None:
            return

        print ('PlantingModel: Msg: Configuring planting model from text file.')
        allYearPlantingList = list()
        parsePattern = '([a-zA-Z0-9\-\.]+)'
        plantingStrm = params.plantingsStrm

        for line in plantingStrm:
            parseIter = re.finditer(parsePattern, line)
            eltList = [ elt.group(0) for elt in parseIter ]
            year = int( float(eltList[0]) )
            row  = int( float(eltList[1]) )
            col  = int( float(eltList[2]) )
            frac = float( eltList[3] )

            spCoverDict = dict()
            for i in range(4,len(eltList),2):
                spCoverDict[ eltList[i] ] = float( eltList[i+1])

            allYearPlantingList.append([year, row, col, frac, spCoverDict])

        self.config_with_planting_list(allYearPlantingList)

    def config_with_ogr(self, params):

        #self.setup_shapefile(params)

        allYearPlantingList = list()
        for key, dataSource in params.plantingDict.items():
            print(('PlantingModel: Msg: Reading configuration for ' + key))

            layer = dataSource.GetLayer()
            for feature in layer:
                try:
                    objectID     = feature.GetField('OBJECTID')

                    print(('PlantingModel: Msg: Working on object id ' + str(objectID)))
                    value        = self.get_field(feature, 'Date_Pla_3', r'^ *YR +([0-9]+) *$',              r'([0-9]+)')
                    year         = float( value.next().group(0) )

                    value        = self.get_field(feature, 'Veg_Type_M', r'^ *([A-Z0-9]+ *,{0,1} *)+$',      r'([A-Z0-9]+)' )
                    spList       = [ m.group(0) for m in value]

                    value        = self.get_field(feature, 'Frac_Spe_1', r'^ *([0-9]*\.{0,1}[0-9]+ *,{0,1} *)+$', r'([0-9]*\.{0,1}[0-9]+)' )
                    coverList    = [ float(m.group(0)) for m in value]
                except ValueError as error:
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: The reason follows.')
                    print((str(error)))
                    continue
                except LAVegModNoneValueError as error:
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: The reason follows.')
                    print((str(error)))
                    continue
                except LAVegModNoMatchError as error:
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: The reason follows.')
                    print((str(error)))
                    continue

                if year < params.startYear or params.endYear < year:
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: Planting year is out of range of model simulation period')
                    print(('PlantingModel: Warning: Planting year is ' + str(year)))
                    continue

                if len(spList) != len(coverList):
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: Number of species in Veg_Type_M does not match number of cover values in Frac_Spe_1')
                    continue

                #print 'PlantingModel: Msg: Building planting dictionary from field information'
                spCoverDict = dict()
                total       = 0.0
                for sp,cover in zip(spList, coverList ):
                    spCoverDict[sp] = cover
                    total          += cover

                if total > 1.0:
                    print(('PlantingModel: Warning: Planting information will not be used for object id = ' + str(objectID)))
                    print ('PlantingModel: Warning: Total fraction of area to be planted is larger than 1.0')
                    print(('PlantingModel: Warning: total = ' + str(total)))
                    continue

                #print 'PlantingModel: Msg: Getting geometry from feature'
                geom       = feature.GetGeometryRef()

                #print 'PlantingModel: Msg: Building location list'
                locList    = self.get_location_list(geom, params.initCond)

                #print 'PlantingModel: Msg: Adding information'
                for elt in locList:
                    row        = elt[0]
                    col        = elt[1]
                    frac       = elt[2]
                    allYearPlantingList.append([year, row, col, frac, spCoverDict])

        self.config_with_planting_list(allYearPlantingList)

    def config(self, params):
        self.config_with_txt(params)

    def update(self):
        for planting in self.plantingList:
            row            = planting[1]
            col            = planting[2]
            frac           = planting[3]
            plantingDict   = planting[4]

            spCoverDict    = PlantingModel.dynModel[(row,col)]

            totalWater = spCoverDict['SAV'] + spCoverDict['WATER']
            if totalWater == 1.0:
                continue

            totalLand  = 1.0 - totalWater
            frac       = min( frac, totalLand )
            scale      = (totalLand - frac)/totalLand

            for sp in itertools.filterfalse( lambda k : k == 'SAV' or k == 'WATER', iter(spCoverDict.keys()) ):
                spCoverDict[sp] *= scale

            for sp,cover in plantingDict.items():
                spCoverDict[sp] += cover * totalLand * ( 1.0 - scale )

            total = 0.0
            for cover in spCoverDict.values():
                total += cover

            if (total - 1.0) > 0.01:
                errorMsg  = 'PlantingModel: Error: Total cover in cell > 1.0. Details follow.' + '\n'
                errorMsg += 'PlantingModel: Error: row,col      = ' + str(row) + ', ' + str(col) + '\n'
                errorMsg += 'PlantingModel: Error: frac         = ' + str(frac)         + '\n'
                errorMsg += 'PlantingModel: Error: total        = ' + str(total)         + '\n'
                errorMsg += 'PlantingModel: Error: totalWater   = ' + str(totalWater)   + '\n'
                errorMsg += 'PlantingModel: Error: totalLand    = ' + str(totalLand)    + '\n'
                errorMsg += 'PlantingModel: Error: scale        = ' + str(scale)        + '\n'
                errorMsg += 'PlantingModel: Error: spCoverDict  = ' + str(spCoverDict)  + '\n'
                errorMsg += 'PlantingModel: Error: plantingDict = ' + str(plantingDict) + '\n'
                raise RuntimeError(errorMsg)

            if total > 1.0:
                for sp in spCoverDict.keys():
                    spCoverDict[sp] /= total
            else:
                spCoverDict['BAREGRND'] += 1.0 - total
