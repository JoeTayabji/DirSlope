import arcpy, numpy, math, os
from arcpy.sa import *
from math import *
from shapely.geometry import MultiLineString, Point
from shapely import wkt
from shapely.ops import split
from osgeo import gdal, ogr, osr
import geopandas as gpd

arcpy.env.workspace=r'E:\DIRECTIONAL_SLOPE\WorkingScripts\SampleData'
arcpy.env.overwriteOutput=True
inFile='FishHatchLine_UTM13.shp'
inDEM='FishHatchDEM_1m.tif'
sr=arcpy.Describe(inDEM).spatialReference

lineChoice = input("Enter 1 to sample evenly along your input line. \nEnter 2 to split your line at its vertices. \n")

if lineChoice=="1":
    sampDist=float(input("Enter the number in meters to sample your line by. \n"))

    # Read line file, read fields in ogr
    line_in_ds = ogr.Open(os.path.join(arcpy.env.workspace,inFile), 0)
    gpdFile = gpd.read_file(os.path.join(arcpy.env.workspace,inFile))
    line_in_layer = line_in_ds.GetLayer()
    line_in_spref = line_in_layer.GetSpatialRef()
    line_in_layer_def = line_in_layer.GetLayerDefn() # schema of input line from the .dbf
    num_fields = line_in_layer_def.GetFieldCount()
    field_name = [line_in_layer_def.GetFieldDefn(i).GetName() for i in range(num_fields)]

    # create segmented line (output)
    line_out_path=os.path.join(arcpy.env.workspace, "LineSplit.shp")
    line_out_ds = ogr.GetDriverByName("ESRI Shapefile").CreateDataSource(line_out_path)
    line_out_layer = line_out_ds.CreateLayer("Split Line", geom_type = ogr.wkbLineString, srs = line_in_spref)

    # get input field names and copy them to the output line
    for i in range(num_fields):
        field_defn = line_in_layer_def.GetFieldDefn(i)
        new_field = ogr.FieldDefn(field_defn.GetName(), field_defn.GetType())
        line_out_layer.CreateField(new_field)

    # iterate through lines in input and split at sampling distance
    for p in line_in_layer:
        geom = p.GetGeometryRef().ExportToWkt()
        shapely_line = wkt.loads(geom)
        # Point at distance 0
        first_point = Point(list(shapely_line.coords)[0])
        previous_x = first_point.xy[0][0]
        previous_y = first_point.xy[1][0]

        # Go through the input line in increments equal to sampDist
        line_length = shapely_line.length
        print("The input line length is {:.2f} m.".format(line_length))
        distance = sampDist
        while distance < line_length:
            # get point coordinates for current location on the input line
            point = shapely_line.interpolate(distance)
            x, y = point.xy[0][0], point.xy[1][0]
            # add line to layer
            line_out_geom = ogr.Geometry(ogr.wkbLineString)
            line_out_geom.AddPoint(previous_x, previous_y)
            line_out_geom.AddPoint(x, y)
            line_out_feat = ogr.Feature(line_out_layer.GetLayerDefn())
            line_out_feat.SetGeometry(line_out_geom)
            
            # Copy attributes from the input line to the output line
            for fn in field_name:
                line_out_feat.SetField(fn, p.GetField(fn))
                
            line_out_layer.CreateFeature(line_out_feat)

            # set up for next iteration
            distance += sampDist
            previous_x = x
            previous_y = y
        # do the above steps one last time for the last point on the line
        last_point = Point(list(shapely_line.coords)[-1])
        x, y = last_point.xy[0][0], last_point.xy[1][0]
        line_out_geom = ogr.Geometry(ogr.wkbLineString)
        line_out_geom.AddPoint(previous_x, previous_y)
        line_out_geom.AddPoint(x, y)
        line_out_feat = ogr.Feature(line_out_layer.GetLayerDefn())
        line_out_feat.SetGeometry(line_out_geom)
        for fn in field_name:
            line_out_feat.SetField(fn, p.GetField(fn))
        line_out_layer.CreateFeature(line_out_feat)
        print("Output line with {0} segments created.".format(line_out_layer.GetFeatureCount()))   

    # reset the output shapefile so it has data
    del line_out_ds, line_out_layer
    outFileName = inFile.split(".")[0]+"_DS_Samp"+str(sampDist)+"m.shp"

elif lineChoice=="2":
    #Split Line at Vertices (only available in Advanced license)
    arcpy.management.SplitLine(inFile, "LineSplit.shp")
    outFileName = inFile.split(".")[0]+"_DS_SV.shp"

else:
    print ("You typed an incorrect choice. Stopping script. Try again.")
    exit()

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
with arcpy.EnvManager(outputCoordinateSystem=sr, cellSize=inDEM, snapRaster=inDEM, mask=inDEM,
                      extent=inDEM):
    arcpy.conversion.FeatureToRaster("LineSplit.shp", "Bear", "LineBearing.tif", inDEM)


fdir="LineBearing.tif"

##print(fdir.rowCount)

#This next section of code is adapted from Sterling Quinn reference from Matt's method
#   Convert rasters to numpy arrays
fdirArray = arcpy.RasterToNumPyArray(fdir, "", "", "", -9999)
elevArray = arcpy.RasterToNumPyArray(inDEM, "", "", "", -9999)


#Getting number of rows and columns
#Getting total cells, origin (lower left), cell size
print(elevArray.shape)
print(fdirArray.shape)
nRows,nCols=elevArray.shape
cellsTotal=nCols*nRows
##print(cellsTotal)
d=arcpy.Describe(fdir)
sr = d.spatialReference
origin=d.extent.lowerLeft
cSize=arcpy.Raster(fdir).meanCellHeight



#   Loop through raster cells
arcpy.SetProgressor("step", "", 0, cellsTotal)
#Creating a blank array with the same rows/columns
blankArray = numpy.empty((nRows,nCols))
for row in range (nRows):
    for col in range (nCols):
        arcpy.SetProgressorPosition()

        # Get value of DEM pixel
        elevPixel=elevArray[row, col]
        if elevPixel == -9999: continue

        # Get value of corresponding flow direction pixel
        fdirPixel=fdirArray[row, col]
        if fdirPixel == -9999: continue

        if row+1>=nRows or row-1<0: continue
        if col+1>=nCols or col-1<0: continue

        ##TopLeft = a
        a=elevArray[row-1, col-1]
        ##TopCenter = b
        b=elevArray[row-1, col]
        ##TopLeft = c
        c=elevArray[row-1, col+1]
        ##Left = d
        d=elevArray[row, col-1]
        ##THE PIXEL = e
        e=elevArray[row, col]
        ##Right=f
        f=elevArray[row, col+1]
        ##BottomLeft=g
        g=elevArray[row+1, col-1]
        ##Bottom=h
        h=elevArray[row+1, col]
        ##BottomRight=i
        i=elevArray[row+1, col+1]

        #Right col - left col
        dzdx = ((c+(2*f)+i)-(a+(2*d)+g))/(8*cSize)

        #Top row - bottom row
        dzdy = ((c+(2*b)+a)-(i+(2*h)+g))/(8*cSize)

        #Bottom row - top row
        #dzdy = ((i+(2*h)+g)-(c+(2*b)+a))/(8*cSize)


        direction_rad = fdirPixel*pi / 180.0
        tan_alpha = ((dzdx*sin(direction_rad)) + (dzdy*cos(direction_rad)))

        # output in degrees - downward slopes: positive, upward slopes: negative
        slope = numpy.arctan(tan_alpha)*(-180.0)/pi
        

        # Write the output cell value to the output array
        blankArray[row,col] = slope


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
                           outFileName)
print ("All done!")
