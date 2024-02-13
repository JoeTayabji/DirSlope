import arcpy, numpy, math
from arcpy.sa import *

arcpy.env.workspace=r'E:\DIRECTIONAL_SLOPE\DaraData'
arcpy.env.overwriteOutput=True
inFile='FishHatchLine_UTM13.shp'
inDEM='FishHatchDEM_1m.tif'

sr=arcpy.Describe(inDEM).spatialReference

#Split Line at Vertices (only available in Advanced license)
arcpy.management.SplitLine(inFile, "LineSplit.shp")

#Add a field to hold the line segment bearing
arcpy.management.AddField("LineSplit.shp", "Bear", "DOUBLE")

#Calculate the line segment bearing into the Bear field
##LINE_BEARING—An attribute will be added to store the start-to-end bearing
##of each line feature. Values range from 0 to 360,
##with 0 meaning north, 90 east, 180 south, 270 west, and so on.
arcpy.management.CalculateGeometryAttributes("LineSplit.shp", [["Bear","Line_Bearing"]])

if arcpy.Exists("LineBearing.tif"):
    arcpy.Delete_management("LineBearing.tif")

print ("Calculated bearing, now making rasters")

#Feature to Raster, set output cell size to DEM, snap raster, OCS to DEM in env
with arcpy.EnvManager(outputCoordinateSystem=sr, cellSize=inDEM, snapRaster=inDEM):
    arcpy.conversion.FeatureToRaster("LineSplit.shp", "Bear", "LineBearing.tif", inDEM)


fdir="LineBearing.tif"

#This next section of code is adapted from Sterling Quinn reference from Matt's method
#   Convert rasters to numpy arrays
fdirArray = arcpy.RasterToNumPyArray(fdir, "", "", "", -9999)
elevArray = arcpy.RasterToNumPyArray(inDEM, "", "", "", -9999)

nRows,nCols=fdirArray.shape
cellsTotal=nCols*nRows
##print(cellsTotal)
d=arcpy.Describe(fdir)
sr = d.spatialReference
origin=d.extent.lowerLeft
cSize=arcpy.Raster(fdir).meanCellHeight

#   Constants to help find cell neighbor
fDirs=(1,2,4,8,16,32,64,128)
rookDirs = (1,4,16,64)
diagDirs = (2,8,32,128)
dCol=(1,  1,  0, -1, -1,-1, 0,1)
dRow=(0,  1,  1,  1,  0, -1, -1,-1)


def forceFlow(num):
    if num >=337.5 or (num >=0 and num<22.5):
        flow=64
    elif num>=22.5 and num<67.5:
        flow=128
    elif num>=67.5 and num<112.5:
        flow=1
    elif num>=112.5 and num<157.5:
        flow=2
    elif num>=157.5 and num<202.5:
        flow=4
    elif num>=202.5 and num<247.5:
        flow=8
    elif num>=247.5 and num<292.5:
        flow=16
    elif num>=292.5 and num<337.5:
        flow=32
    return flow

#   Loop through raster cells
arcpy.SetProgressor("step", "", 0, cellsTotal)
blankArray = numpy.empty((nRows,nCols))
for nRow in range (nRows):
    for nCol in range (nCols):
        arcpy.SetProgressorPosition()

        # Get value of DEM pixel
        elevPixel=elevArray[nRow, nCol]
        if elevPixel == -9999: continue

        # Get value of corresponding flow direction pixel
        fdirPixel=fdirArray[nRow,nCol]
        if fdirPixel == -9999: continue
        ##else: i = fDirs.index(fdirPixel)
        else:
            ##print(fdirPixel)
            fdirPixel=forceFlow(fdirPixel)
            i = fDirs.index(fdirPixel)

        # Get the column of the cell to compare
        nC = nCol + dCol[i]
        if nC < 0 or nC == nCols: continue

        # Get the row of the cell to compare
        nR = nRow + dRow[i]
        if nR < 0 or nR == nRows: continue

        # Get the elevation pixel value of the cell to compare
        dirElevPixel = elevArray[nR, nC]

        # Calculate the difference between elevation pixels
        elevPixelDiff = elevPixel - dirElevPixel

        maxSlope = -9999

        # If cell is in diagonal direction, divide the difference by the orthogonal size
        if fdirPixel in diagDirs:
            maxSlope = ((elevPixelDiff * 1.0) / (math.sqrt(2 * (cSize * cSize))))
        # If cell is in rook direction, divide the difference by the cell size
        elif fdirPixel in rookDirs:
            maxSlope = ((elevPixelDiff * 1.0) / cSize)
        else:
            maxSlope = -9999

        # Write the output cell value to the output array
        blankArray[nRow,nCol] = maxSlope * 100


# Convert the output array to a raster
myRaster = arcpy.NumPyArrayToRaster(blankArray,origin,cSize,cSize)
arcpy.DefineProjection_management(myRaster, sr)
myRaster.save("OutDirSlope.tif")
del fdirArray,elevArray,blankArray

print ("Completed dslope raster")

#Sample the directional slope raster values at the line splits
#This creates points at the centers of the line segments
Sample("OutDirSlope.tif", "LineSplit.shp", "DirSlopePoints.shp",
       generate_feature_class = "FEATURE_CLASS")

#Perform a spatial join to get the dslope values back to the line
#segments
arcpy.analysis.SpatialJoin("LineSplit.shp", "DirSlopePoints.shp",
                           "DirSlopeLineOutput.shp")
print ("All done!")
