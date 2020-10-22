import os
import numpy as np
import arcpy

#############################################################################################
# RB_to_terrain.py converts the RB output table (.csv) to terrain (.asc)
#   Written by Anzy Lee, Ph.D., Utah State University
#   Date: 10/21/2020
#############################################################################################
# Input variables: RB_name, asc_name, cell_size, execute
RB_name = 'S1'    # Name of the river builder case
asc_name = 'S1'    # Name of the ascii terrain
cell_size = '1'                 # Cell size of ascii terrain
execute = np.array([1,          # 1 if you want to execute "Table to point",
                    1,          # 1 if you want to execute "Create TIN",
                    1,          # 1 if you want to execute "TIN to Raster",
                    1])         # 1 if you want to execute "Raster to asc"
#############################################################################################
# Workspace setting
arcpy.env.workspace = "./"+RB_name+"/"+RB_name
sr = arcpy.SpatialReference(3857, 115700) #  WGS_1984_web_mercator, WGS 1984
#sr = arcpy.SpatialReference(4759, 115700) # WGS 1984, WGS 1984
arcpy.CheckOutExtension("3D")

if execute[0] == 1:
    # 1 Table to point
    in_Table = arcpy.env.workspace+"/SRVtopo.csv"
    output_point = RB_name+'_xyz.shp'
    x_coords = "X"
    y_coords = "Y"
    z_coords = "Z"

    # Make the XY event layer...
    arcpy.management.XYTableToPoint(in_Table, output_point,
                                    x_coords, y_coords, z_coords,#arcpy.SpatialReference(4759, 115700))
                                    sr)
    print(arcpy.GetCount_management(output_point))

    # Points should be adjusted to create a "RUNWAY'
    print('Points should be adjusted to create a RUNWAY')
    os.system("pause")

if execute[1] == 1:
    # 2 Create TIN
    in_point = RB_name+'_xyz.shp'
    output_TIN = RB_name+'_TIN'

    arcpy.ddd.CreateTin(output_TIN, sr, in_point+" Z masspoints")

if execute[2] == 1:
    # 3 TIN to Raster
    in_TIN = RB_name+'_TIN'
    out_tif = RB_name+'.tif'
    # Set variables for TIN to Raster
    dataType = "FLOAT"  # Default
    method = "LINEAR"  # Default
    sampling = "CELLSIZE " + cell_size
    zfactor = "1"

    arcpy.ddd.TinRaster(in_TIN, out_tif, dataType,
                    method, sampling, zfactor)

if execute[3] == 1:
    # 4 Raster to ascii
    in_tif = RB_name+'.tif'
    out_ascii = asc_name + '.asc'

    arcpy.RasterToASCII_conversion(in_tif, out_ascii)
#############################################################################################