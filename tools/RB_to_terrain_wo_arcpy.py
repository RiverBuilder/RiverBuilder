#############################################################################################
# RB_to_terrain.py converts the RB output table (.csv) to terrain (.asc)
#   Written by Anzy Lee, Postdoctoral Scholar, Utah State University
#   Date: 10/21/2020
#############################################################################################

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import rasterio
import geopandas as gpd
import shapely
import pathlib
import math

#############################################################################################
# Input variables: RB_path, asc_name, RB_unit, asc_unit, cell_size, execute
#NAME = 'sfe_209'                               # Name of the case
case_name = 'sfe_322_s2'                        # The name of terrain
directory = "../samples/sfe_322/sfe_322_s2"     # Folder Name

RB_unit = 'meter'               # Unit of the river archetype
asc_unit = 'meter'              # Unit of the ascii terrain
cell_size = 1.0                 # Cell size of ascii terrain in asc_unit
crs = 'EPSG:3857'               # Set CRS, e.g., 'EPSG:3857' (WSG84 Pseudo-Mercator),
                                # 'EPSG:4326' (WGS84), or 'EPSG:2056' (CH1903+/LV95)

NAME = 'sfe_'+case_name.split('_')[1]

execute = np.array([1,          # 1 if you want to execute "Table to point",
                    1,          # 1 if you want to execute "Create TIN",
                    1,          # 1 if you want to execute "TIN to Raster",
                    1])         # 1 if you want to execute "Raster to asc"

#############################################################################################
# Workspace setting
RB_path = pathlib.Path(directory).absolute().resolve() / case_name  # path to SRVtopo directory

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
    in_Table = RB_path / "SRVtopo.csv"
    out_Table = RB_path / "SRVtopo_xyz.csv"
    df = pd.read_csv(in_Table)
    offset = 100 # to prevent minus values
    df.X = df.X*conv_factor +offset
    df.Y = df.Y*conv_factor +offset
    df.Z = df.Z*conv_factor +offset
    # df.to_csv(out_Table,index=False)
    df.to_csv(out_Table)

    # 1 Table to point
    output_point = RB_path / (case_name +'_xyz.shp')

    # Make the XY event layer...
    print("2. Running Table to point")
    print("# of points = " + str(len(df)))

    # Create a GeoDataFrame for the points
    geometry = [shapely.geometry.Point(xy) for xy in zip(df.X.values, df.Y.values)]
    gdf = gpd.GeoDataFrame({'Z': df['Z'].to_numpy(),'LABEL': df['Label'].to_numpy()}, geometry=geometry, crs=crs)

    # Save the GeoDataFrame to a shapefile
    gdf.to_file(output_point)

    # Points can be adjusted to create a "RUNWAY'
    # print('Points should be adjusted to create a RUNWAY')
    # os.system("pause")

#############################################################################################
if execute[1] == 1:
    # 2 Create TIN
    output_TIN = RB_path / (case_name + '_tin_triangles.shp')

    print("3. Running Create TIN")
    TIN_delaunay = scipy.spatial.Delaunay(df[['X','Y']].to_numpy())

    # # Create a GeoDataFrame for the TIN triangles
    # triangles = []

    # for ii,simplex in enumerate(TIN_delaunay.simplices):
    #     print(f"Progress: {ii/len(TIN_delaunay.simplices)*100.0:.2}/100%", end='\r')

    #     pts = df[['X','Y','Z']].to_numpy()[simplex]
    #     poly = shapely.geometry.Polygon([(pts[0][0], pts[0][1]), (pts[1][0], pts[1][1]), (pts[2][0], pts[2][1])])
    #     triangles.append(poly)

    # tin_gdf = gpd.GeoDataFrame(geometry=triangles, crs=crs)
    # tin_gdf.to_file(output_TIN)

#############################################################################################
if execute[2] == 1:
    # 3 TIN to Raster
    out_tif = RB_path / (case_name+'.tif')

    print("4. Running TIN Raster")

    # Calculate the extent of the point dataset
    min_x, max_x = df.X.values.min(), df.X.values.max()
    min_y, max_y = df.Y.values.min(), df.Y.values.max()
    num_x = int(math.ceil((max_x - min_x) / cell_size)) + 1
    num_y = int(math.ceil((max_y - min_y) / cell_size)) + 1
    max_x = min_x + (num_x-1)*cell_size
    max_y = min_y + (num_y-1)*cell_size
    # Create the grid based on the extent and cell size
    grid_x, grid_y = np.mgrid[min_x:max_x:num_x*1j, max_y:min_y:num_y*1j]
    grid_x = grid_x.T
    grid_y = grid_y.T
    grid_x = grid_x + 0.5*cell_size - 10**(-6)*cell_size
    grid_y = grid_y + 0.5*cell_size - 10**(-6)*cell_size
    interpolator = scipy.interpolate.LinearNDInterpolator(TIN_delaunay, df[['Z']].to_numpy())
    grid_z = np.squeeze(interpolator(grid_x,grid_y))

    # Define the GeoTIFF metadata
    transform = rasterio.transform.from_origin(min_x , max_y+cell_size, cell_size, cell_size)
    raster_data = rasterio.open(
                    out_tif,
                    'w',
                    driver='GTiff',
                    height=grid_z.shape[0],
                    width=grid_z.shape[1],
                    count=1,
                    dtype=grid_z.dtype,
                    crs=crs,  # Set CRS to EPSG:4326 (WGS84)
                    transform=transform,
                )

    raster_data.write(grid_z, 1)
    raster_data.close()


#############################################################################################

if execute[3] == 1:
    # 4 Raster to ascii
    out_ascii = RB_path / (case_name+'.asc')
    # delete Ascii if it exists, otherwise rasterio gives an error
    if out_ascii.exists():
        out_ascii.unlink()
    print("5. Running Raster To ASCII")
    raster_data = rasterio.open(
                out_ascii,
                'w',
                driver='AAIGrid',
                height=grid_z.shape[0],
                width=grid_z.shape[1],
                count=1,
                dtype=grid_z.dtype,
                crs=crs,  # Set CRS to EPSG:4326 (WGS84)
                transform=transform,
                force_cellsize=True
                )

    raster_data.write(grid_z, 1)
    raster_data.close()
#############################################################################################
