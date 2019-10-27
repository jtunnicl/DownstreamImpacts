import arcpy, os, traceback, sys, numpy
from arcpy import env
from arcpy.sa import *
env.overwriteOutput = True

class Toolbox(object):
    def __init__(self):
        self.label = "Trace Downstream Impacts"
        self.alias = ""
        # List of tool classes associated with this toolbox
        self.tools = [DownstreamImpacts]

class DownstreamImpacts(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Trace Downstream Impacts from Polygon"
        self.description = "Identify the likely point of deposition within a drainage network, downstream from a site of mass wasting (polygon). \n\r The Search Radius constrains how far from the polygon centroid the algorithm will search; this should be large enough to span a few stream orders, but small enough as to not impact performance - 1000-2000 cells (e.g. 2-4 km for a 2m DEM) should ensure reasonable performance. Other raster items should all have the same extent and resolution, having been generated via the Hydrology toolset, in Spatial Analyst Tools. The resultant point file will record the flow accumulation, elevation and stream order at the site of deposition."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Parameter definitions"""

        searchRadius = arcpy.Parameter(
            displayName="Width of Search Window (Number of Cells)",
            name="search_radius",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        # Set Default Value
        searchRadius.value = 200

        depoGradient = arcpy.Parameter(
            displayName="Channel Slope Gradient below which deposition occurs",
            name="depo_gradient",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        # Set Default Value
        depoGradient.value = 0.15

        slumpPolys = arcpy.Parameter(
            displayName="Disturbance Polygon Shapes",
            name="disturbance_polys",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        slumpPolys.filter.list = ["Polygon"]

        slumpRaster = arcpy.Parameter(
            displayName="Disturbance Polygon Raster",
            name="disturbance_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        elev = arcpy.Parameter(
            displayName="Digital Elevation Model",
            name="dem_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        fdir = arcpy.Parameter(
            displayName="Flow Direction Raster",
            name="fdir_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        facc = arcpy.Parameter(
            displayName="Flow Accumulation Raster",
            name="facc_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        order = arcpy.Parameter(
            displayName="Strahler Stream Order Raster",
            name="strahler_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        newfile = arcpy.Parameter(
            displayName="Output Point File",
            name="output_points",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Output")

        parameters = [searchRadius, depoGradient, slumpPolys, slumpRaster, elev, fdir, facc, order, newfile]
        return parameters

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        # Get the parameter 'values'
        searchRadius = int(parameters[0].valueAsText)
        depoGradient = parameters[1]
        slumpPolys = parameters[2].valueAsText
        polyPath = (arcpy.Describe(slumpPolys).path + '\\' + arcpy.Describe(slumpPolys).name)
        slumpRaster = arcpy.Raster(parameters[3].valueAsText)
        elev = arcpy.Raster(parameters[4].valueAsText)
        fdir = arcpy.Raster(parameters[5].valueAsText)
        facc = arcpy.Raster(parameters[6].valueAsText)
        order = arcpy.Raster(parameters[7].valueAsText)
        newfile = parameters[8].valueAsText

        # Describe data
        dx = elev.meanCellWidth
        width = int(elev.extent.width/dx)
        originX = elev.extent.lowerLeft.X
        originY = elev.extent.lowerLeft.Y
        height = int(elev.extent.height/dx)
        spatialRef = arcpy.Describe(elev).spatialReference
        NoData = arcpy.Describe(elev).noDataValue

        # Output files
        arcpy.management.CreateFeatureclass(os.path.dirname(newfile),os.path.basename(newfile),"POINT","","DISABLED","DISABLED", spatialRef)
        arcpy.AddField_management(newfile, "AccumArea", field_type="DOUBLE")
        arcpy.AddField_management(newfile, "DistFrmSrc", field_type="DOUBLE")
        arcpy.AddField_management(newfile, "StrmOrder", field_type="DOUBLE")

# Helper arrays for determining distances and routing, based on flow direction raster
# Some of this technique borrowed from FelixIP
# https://gis.stackexchange.com/questions/136715/getting-cell-value-along-flow-direction-using-arcpy

        fDirs=(1,2,4,8,16,32,64,128)
        fOrtho = (1, 4, 16, 64)
        fDiag = (2,8,32,138)
        dCol=(1,  1,  0, -1, -1,-1, 0,1)
        dRow=(0,  1,  1,  1,  0, -1, -1,-1)

# loop through the polys
        with arcpy.da.SearchCursor(polyPath, ["OID@", "SHAPE@"]) as cursor:
            for row in cursor:
                offsetX = int((row[1].centroid.X - originX)/dx)   # Get the centroid point mapped onto the raster grid
                CentrX = originX + offsetX * dx
                offsetY = int((row[1].centroid.Y - originY)/dx)
                CentrY = originY + offsetY * dx
                poly_index = int(row[0])
                idx = 0

                ## Arrays for collecting all coordinates and values along the trace route
                steps = int(searchRadius * 0.4)   # Less than half the search radius, otherwise it can run off the domain
                streamOrder = numpy.zeros(steps, dtype=numpy.float32)
                runningDistance = numpy.zeros(steps, dtype=numpy.float32)
                bedElev = numpy.zeros(steps, dtype=numpy.float32)
                upstrArea = numpy.zeros(steps, dtype=numpy.float32)
                nRow = numpy.zeros(steps)    
                nCol = numpy.zeros(steps)

                # Window sections of the raster and carry out tracing
                try:
                    lowerLeft = arcpy.Point(CentrX - int((searchRadius/2)*dx), CentrY - int((searchRadius/2)*dx))
                    # Turn rasters into Numpy Arrays
                    elevArray = arcpy.RasterToNumPyArray(elev,lowerLeft,searchRadius,searchRadius,NoData)
                    dirArray = arcpy.RasterToNumPyArray(fdir,lowerLeft,searchRadius,searchRadius,NoData)
                    faccArray = arcpy.RasterToNumPyArray(facc,lowerLeft,searchRadius,searchRadius,NoData)
                    slumpArray = arcpy.RasterToNumPyArray(slumpRaster,lowerLeft,searchRadius,searchRadius,NoData)
                    orderArray = arcpy.RasterToNumPyArray(order,lowerLeft,searchRadius,searchRadius,NoData)
                    dataRows = elevArray.shape[0] - 1
                    dataCols = elevArray.shape[1] -1
                    vals = faccArray[numpy.where(slumpArray == poly_index)]
                    ans = numpy.where(slumpArray == poly_index) #       two tuples with R, C coordinates within slump
                    maxVal = numpy.amax(faccArray[ans])
                    maxValCoordinate = numpy.where(faccArray[ans] == maxVal) [0][0] # Highest facc value within the slump polygon
                    ## Start point for tracing
                    nR = ans[0][maxValCoordinate]
                    nC = ans[1][maxValCoordinate]
                            # Note this can be problematic with quite large polygons, 
                            # where search origin is possibly much closer to bounds
                    nRow[0] = nR                 # initiate array of coordinates
                    nCol[0] = nC
                    if orderArray[nR,nC] == NoData: streamOrder[0] = 0    # Treat anything less than Order 1 as Order 'zero'
                    else: streamOrder[0] = orderArray[nR,nC]

                    runningDistance[0] = 0       # Initialise the first element of these arrays at the starting point
                    bedElev[0] = elevArray[nR, nC]
                    upstrArea[0] = faccArray[nR, nC]
                    for b in range(1,steps-1):
                        direction=dirArray[nR, nC]
                        if direction<1:break
                        i=fDirs.index(direction)
                        dX=dCol[i]; nC+=dX
                        if nCol[b]<0 or nCol[b]==dataCols: break
                        dY=dRow[i]; nR+=dY
                        if nRow[b]<0 or nRow[b]==dataRows: break
                        nRow[b] = nR
                        nCol[b] = nC
                        if orderArray[nR,nC] == NoData: streamOrder[b] = 0
                        else: streamOrder[b] = orderArray[nR,nC]
                        runningDistance[b] = runningDistance[b-1] + numpy.any(fOrtho == direction)  * dx + numpy.any(fDiag == direction) * 1.4142 * dx
                        bedElev[b] = elevArray[nR, nC]
                        upstrArea[b] = faccArray[nR, nC]

                    bump_coord = numpy.where(numpy.diff(streamOrder)>0)[0]   # point(s) where stream order changes; can be none, or more than one
                    # Case where no junction shows up; give up and use a point shortly downstream of origin
                    if bump_coord.size == 0:
                        idx = 5   # Deposition a short distance (5 cells) downslope of the polygon
                        Yval = nRow[idx]
                        Xval = nCol[idx]
                        vertex = arcpy.Point( lowerLeft.X + ( Xval * dx ) + ( dx/2 ), lowerLeft.Y + ( ( dataRows - Yval ) * dx ) + ( dx/2 ) )

                    # Case where there is only one junction in the whole path; deposition occurs at this transition
                    elif bump_coord.size == 1:
                        idx = bump_coord[0] + 1   # +1 ensures point on path is placed beyond initial order, into the next
                        Yval = nRow[idx]
                        Xval = nCol[idx]
                        vertex = arcpy.Point( lowerLeft.X + ( Xval * dx ) + ( dx/2 ), lowerLeft.Y + ( ( dataRows - Yval ) * dx ) + ( dx/2 ) )

                    # Case with more than one junction; deposition may occur multiple links downstream, if first is steep enough
                    else:
                        depo = 0
                        #for c in range(0, bump_coord.size):
                        while depo < bump_coord.size:          # 4 order transitions is the upper plausible limit
                            # Estimate the slope between origin and the first junction
                            idx = bump_coord[depo] + 1
                            slope = (bedElev[0]- bedElev[idx]) / runningDistance[idx]
                            if slope < depoGradient:   # slope is shallow enough for deposition
                                Yval = nRow[idx]
                                Xval = nCol[idx]
                                break
                            else:
                                depo = depo + 1
                            idx = bump_coord[depo-1] + 1  # Give up, deposit on previous transition. All links are too steep for deposition. 
                            Yval = nRow[idx]
                            Xval = nCol[idx]
                        vertex = arcpy.Point( lowerLeft.X + ( Xval * dx ) + ( dx/2 ), lowerLeft.Y + ( ( dataRows - Yval ) * dx ) + ( dx/2 ) )

                except:  # No evident path downslope, just place the point at the disturbance centroid, fields are zeros
                    vertex = arcpy.Point(CentrX, CentrY)
                    pass

                with arcpy.da.InsertCursor(newfile, ("SHAPE@", "AccumArea", "DistFrmSrc", "StrmOrder")) as cursor:
                    cursor.insertRow((vertex, upstrArea[idx], runningDistance[idx], streamOrder[idx]))

                arcpy.AddMessage('Polygon ID %i' %poly_index + ', vertex x=%f' %vertex.X + ', y=%f' %vertex.Y )

        del cursor

        return
