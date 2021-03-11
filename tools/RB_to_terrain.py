#############################################################################################
# RB_to_terrain.py converts the RB output table (.csv) to terrain (.asc)
#   Written by Anzy Lee, Postdoctoral Scholar, Utah State University
#   Date: 10/21/2020
#############################################################################################

import os
import numpy as np
import arcpy
import pandas as pd

#############################################################################################
# Input variables: RB_path, asc_name, RB_unit, asc_unit, cell_size, execute
case_name = 'S1_noSlope'         # The name of case
NAME = 'S1'                      # Name of the folder

RB_unit = 'meter'                # Unit of the river archetype
asc_unit = 'meter'              # Unit of the ascii terrain
cell_size = '1'                 # Cell size of ascii terrain in asc_unit

execute = np.array([1,          # 1 if you want to execute "Table to point",
                    1,          # 1 if you want to execute "Create TIN",
                    1,          # 1 if you want to execute "TIN to Raster",
                    1])         # 1 if you want to execute "Raster to asc"

#############################################################################################
# Workspace setting
RB_path = os.path.abspath("..\samples\\"+ NAME + "\\" + case_name)  # path to SRVtopo directory
arcpy.env.workspace = RB_path
sr = arcpy.SpatialReference(3857, 115700) #  WGS_1984_web_mercator, WGS 1984
#sr = arcpy.SpatialReference(4759, 115700) # WGS 1984, WGS 1984
arcpy.CheckOutExtension("3D")

if RB_unit == 'foot':
    RB_conv = 3.28084
else:
    RB_conv = 1
if asc_unit == 'foot':
    asc_conv = 3.28084
else:
    asc_conv = 1
conv_factor = asc_conv/RB_conv

#############################################################################################
if execute[0] == 1:
    # 0 Unit conversion
    print('1. Converting units')
    in_Table = arcpy.env.workspace + "\\SRVtopo.csv"
    out_Table = arcpy.env.workspace + "\\SRVtopo_xyz.csv"
    df = pd.read_csv(in_Table)
    offset = 100 # to prevent minus values
    df.X = df.X*conv_factor +offset
    df.Y = df.Y*conv_factor +offset
    df.Z = df.Z*conv_factor +offset
    df.to_csv(out_Table)

    # 1 Table to point
    in_Table = arcpy.env.workspace+"/SRVtopo_xyz.csv"
    output_point = case_name+'_xyz.shp'
    x_coords = "X"
    y_coords = "Y"
    z_coords = "Z"

    # Make the XY event layer...
    print("2. Running Table to point")
    arcpy.management.XYTableToPoint(in_Table, output_point,
                                    x_coords, y_coords, z_coords,
                                    sr)
    print("# of points = " + str(arcpy.GetCount_management(output_point)))

    # Points can be adjusted to create a "RUNWAY'
    # print('Points should be adjusted to create a RUNWAY')
    # os.system("pause")

#############################################################################################
if execute[1] == 1:
    # 2 Create TIN
    in_point = case_name+'_xyz.shp'
    output_TIN = case_name+'_TIN'

    print("3. Running Create TIN")
    arcpy.ddd.CreateTin(output_TIN, sr, in_point+" Z masspoints")

#############################################################################################
if execute[2] == 1:
    # 3 TIN to Raster
    in_TIN = case_name+'_TIN'
    out_tif = case_name+'.tif'
    # Set variables for TIN to Raster
    dataType = "FLOAT"  # Default
    method = "LINEAR"  # Default
    sampling = "CELLSIZE " + cell_size
    zfactor = "1"

    print("4. Running TIN Raster")
    arcpy.ddd.TinRaster(in_TIN, out_tif, dataType,
                    method, sampling, zfactor)
#############################################################################################

if execute[3] == 1:
    # 4 Raster to ascii
    in_tif = case_name+'.tif'
    out_ascii = case_name + '.asc'

    print("5. Running Raster To ASCII")
    arcpy.RasterToASCII_conversion(in_tif, out_ascii)
#############################################################################################
