import sys
import pathlib
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import euclidean
import rasterio
import rasterio.features
import geopandas as gpd
import shapely
from shapely.affinity import scale
from shapely.ops import unary_union, polygonize, split
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString, GeometryCollection, shape, mapping


"""
This script convert RiverBuilder output to a suitable format for 
generating a computational mesh with BASEmesh.



Dependencies:
- numpy
- pandas
- scipy
- geopandas
- pandas
- shapely

Author: Matthias Bürgler,
        ETH Zurich,
        Laboratory of Hydraulics, Hydrology and Glaciology (VAW)
Date:   23.07.2024
"""

def detect_outliers_ransac(group, window_size=7, outlier_threshold=10):
    outliers = np.zeros(len(group), dtype=bool)
    X = group[['X']].values
    y = group['Y'].values
    
    for i in range(len(group) - window_size + 1):
        window_X = X[i:i + window_size]
        window_y = y[i:i + window_size]
        
        # Fit RANSAC regressor
        ransac = RANSACRegressor()
        ransac.fit(window_X, window_y)
        inlier_mask = ransac.inlier_mask_
        outlier_mask = np.logical_not(inlier_mask)
        
        # Calculate the residuals
        predictions = ransac.predict(window_X)
        residuals = np.abs(window_y - predictions)
        
        # Detect outliers
        window_outliers = residuals > outlier_threshold
        outliers[i:i + window_size] |= window_outliers
    
    return outliers

def calculate_min_distances(group):
    distances = []
    indices = group.index.tolist()

    for i in range(len(group)):
        if i == 0:
            # Only one neighbor (next)
            next_point = group.iloc[i + 1][['X', 'Y']]
            distance_next = scipy.spatial.distance.euclidean(group.iloc[i][['X', 'Y']], next_point)
            distances.append(distance_next)
        elif i == len(group) - 1:
            # Only one neighbor (previous)
            prev_point = group.iloc[i - 1][['X', 'Y']]
            distance_prev = scipy.spatial.distance.euclidean(group.iloc[i][['X', 'Y']], prev_point)
            distances.append(distance_prev)
        else:
            # Both neighbors
            prev_point = group.iloc[i - 1][['X', 'Y']]
            next_point = group.iloc[i + 1][['X', 'Y']]
            distance_prev = scipy.spatial.distance.euclidean(group.iloc[i][['X', 'Y']], prev_point)
            distance_next = scipy.spatial.distance.euclidean(group.iloc[i][['X', 'Y']], next_point)
            distances.append(min(distance_prev, distance_next))

    return distances

def determine_outliers(df,factor=3):
    # Calculate distances and add to DataFrame
    df['Min_Distance'] = np.nan
    for name, group in df.groupby('ID'):
        min_distances = calculate_min_distances(group)
        df.loc[group.index, 'Min_Distance'] = min_distances
    
    # Calculate Z-scores
    df['Z_Score'] = df['Min_Distance'].transform(lambda x: scipy.stats.zscore(x, nan_policy='omit'))
    
    # Determine outliers based on Z-score
    df['Outlier'] = np.abs(df['Z_Score']) > 3

    # Return DataFrame with all points and a DataFrame with outliers removed
    df_filtered = df[~df['Outlier']]
    return df_filtered

def points_to_lines(gdf_points):
    # Group by 'ID' and sort points if needed (e.g., by an additional attribute like 'timestamp' or order)
    lines = []
    
    for name, group in gdf_points.groupby('ID'):
        # Ensure points are ordered correctly (if necessary, e.g., by an attribute)
        group = group.sort_index()
        # Create a LineString from the ordered points
        line = shapely.geometry.LineString(group.geometry.tolist())
        # Append the LineString to the list
        lines.append({'ID': name, 'geometry': line})
    
    # Create a GeoDataFrame from the list of lines
    gdf_lines = gpd.GeoDataFrame(lines, geometry='geometry', crs=gdf_points.crs)
    return gdf_lines

def simplify_lines(gdf_lines, tolerance):
    # Apply the simplify method to each LineString with the given tolerance
    gdf_lines['geometry'] = gdf_lines['geometry'].apply(lambda geom: geom.simplify(tolerance))
    return gdf_lines

def lines_to_polygon(gdf_boundary):
    # Ensure the boundary GeoDataFrame has the same CRS
    if not gdf_boundary.crs:
        raise ValueError("GeoDataFrame must have a CRS")

    # Convert lines to a single polygon (assuming they form a closed boundary)
    boundary_union = unary_union(gdf_boundary.geometry)
    
    # Ensure the union is a polygon
    if not isinstance(boundary_union, shapely.geometry.Polygon):
        if boundary_union.geom_type == 'MultiLineString':
            boundary_union = boundary_union.buffer(0)  # Buffer with 0 to close the polygon
        elif boundary_union.geom_type == 'MultiPolygon':
            boundary_union = boundary_union.convex_hull  # Use convex hull to create a single polygon

    if not isinstance(boundary_union, shapely.geometry.Polygon):
        raise ValueError("The boundary lines do not form a valid closed polygon")

    return boundary_union

def clip_lines_with_boundary(gdf_lines, gdf_boundary):
    # Ensure both layers have the same CRS
    if gdf_lines.crs != gdf_boundary.crs:
        gdf_lines = gdf_lines.to_crs(gdf_boundary.crs)

    # Clip lines with the boundary polygon
    clipped_lines = []
    
    for _, line_row in gdf_lines.iterrows():
        line_geom = line_row['geometry']
        
        # Perform the clipping with the boundary polygon
        clipped_geom = line_geom.intersection(gdf_boundary.iloc[0].geometry)
        
        if clipped_geom.is_empty:
            continue
        elif clipped_geom.geom_type == 'LineString':
            clipped_lines.append({'ID': line_row['ID'], 'geometry': clipped_geom})
        elif clipped_geom.geom_type == 'MultiLineString':
            # Keep the longest part
            longest_part = max(clipped_geom.geoms, key=lambda part: part.length)
            clipped_lines.append({'ID': line_row['ID'], 'geometry': longest_part})
    
    # Create a GeoDataFrame from the clipped lines
    gdf_clipped_lines = gpd.GeoDataFrame(clipped_lines, geometry='geometry', crs=gdf_lines.crs)
    
    return gdf_clipped_lines
def create_lines_from_points(gdf_lines):
    # Extract start and end points
    start_points = []
    end_points = []

    for _, row in gdf_lines.iterrows():
        line_geom = row['geometry']
        start_points.append({'ID': row['ID'], 'geometry': shapely.geometry.Point(line_geom.coords[0])})
        end_points.append({'ID': row['ID'], 'geometry': shapely.geometry.Point(line_geom.coords[-1])})

    # Create GeoDataFrames for points
    gdf_start_points = gpd.GeoDataFrame(start_points, geometry='geometry', crs=gdf_lines.crs)
    gdf_end_points = gpd.GeoDataFrame(end_points, geometry='geometry', crs=gdf_lines.crs)

    # Sort IDs based on the specified order
    def sort_key(id):
        if id.startswith('LB'):
            return (0, -int(id[2:]))  # LBXX descending order
        elif 'CL' in id:
            return (1, 0)  # CL comes after LB
        elif id.startswith('RB'):
            return (2, int(id[2:]))  # RBXX ascending order

    sorted_start_points = gdf_start_points.sort_values(by='ID', key=lambda ids: [sort_key(id) for id in ids])
    sorted_end_points = gdf_end_points.sort_values(by='ID', key=lambda ids: [sort_key(id) for id in ids])
    
    # Create LineString from sorted points
    start_line = shapely.geometry.LineString([point for point in sorted_start_points.geometry])
    end_line = shapely.geometry.LineString([point for point in sorted_end_points.geometry])

    # Create GeoDataFrame for new lines
    gdf_new_lines = gpd.GeoDataFrame({
        'ID': ['inflow', 'outflow'],
        'geometry': [start_line, end_line]
    }, geometry='geometry', crs=gdf_lines.crs)

    return gdf_new_lines

def append_lines(existing_gdf, new_lines):
    # Ensure both GeoDataFrames have the same CRS
    if existing_gdf.crs != new_lines.crs:
        new_lines = new_lines.to_crs(existing_gdf.crs)

    # Concatenate the existing GeoDataFrame with the new_lines
    combined_gdf = pd.concat([existing_gdf, new_lines], ignore_index=True)

    return combined_gdf

def remove_valley_breaklines(gdf_lines):
    # Filter out lines with 'LV' or 'RV' in the ID
    gdf_lines_filtered = gdf_lines[~gdf_lines['ID'].str.contains('LV|RV')]
    return gdf_lines_filtered


def create_bounding_polygon(gdf_lines):
    """Create a bounding polygon from multiple lines."""
    # Combine all lines into a single MultiLineString
    all_lines = unary_union(gdf_lines.geometry)
    
    # Ensure that the result is a MultiLineString
    if all_lines.geom_type == 'LineString':
        all_lines = MultiLineString([all_lines])
    elif all_lines.geom_type == 'MultiLineString':
        all_lines = MultiLineString(all_lines.geoms)
    
    # Convert MultiLineString to a polygon using polygonize
    polygons = list(polygonize(all_lines))
    
    # Combine all resulting polygons into a single polygon if needed
    if len(polygons) > 1:
        bounding_polygon = unary_union(polygons)
    elif len(polygons) == 1:
        bounding_polygon = polygons[0]
    else:
        bounding_polygon = None
    
    return bounding_polygon

def save_lines_to_csv_with_ids(line_gdf, output_csv):
    """Save lines with their IDs and points to a CSV file with 'X' and 'Y' columns."""
    # Ensure the DataFrame contains 'ID' and 'geometry'
    if 'ID' not in line_gdf.columns or 'geometry' not in line_gdf.columns:
        raise ValueError("GeoDataFrame must contain 'ID' and 'geometry' columns.")
    
    # Initialize a list to store rows for the CSV
    rows = []
    
    # Iterate over each line feature
    for _, row in line_gdf.iterrows():
        line_id = row['ID']
        line = row['geometry']
        
        # Extract points from each LineString or MultiLineString
        if line.geom_type == 'LineString':
            for coord in line.coords:
                rows.append([line_id, coord[0], coord[1]])
        elif line.geom_type == 'MultiLineString':
            for segment in line:
                for coord in segment.coords:
                    rows.append([line_id, coord[0], coord[1]])
    
    # Create a DataFrame from rows
    points_df = pd.DataFrame(rows, columns=['ID', 'X', 'Y'])
    
    # Save to CSV
    points_df.to_csv(output_csv, index=False)



def main():
    """
        Main function.
    """
    #- User input ----------------------#
    case_name = 'sfe_322_s2'
    directory = "../samples/sfe_322/sfe_322_s2"

    # Set CRS, e.g., 'EPSG:3857' (WSG84 Pseudo-Mercator),
    # 'EPSG:4326' (WGS84), or 'EPSG:2056' (CH1903+/LV95)
    crs = 'EPSG:3857'

    offset = 100

    # Set tolerance value
    simplification_tolerance = 0.5

    #-----------------------------------#

    # Resolve path
    case_path = pathlib.Path(directory).absolute().resolve() / case_name  # path to SRVtopo directory

    # Read in the elevation data 
    dem_file = case_path / (case_name+'.tif')
    with rasterio.open(dem_file) as src:
        dem = src.read(1)  # Read the first band
        transform = src.transform  # Get the affine transformation
        no_data_value = src.nodata  # Get the no-data value if it exists

    # Create a mask for valid (non-NaN and non-no-data) values
    valid_mask = np.isfinite(dem)
    if no_data_value is not None:
        valid_mask &= (dem != no_data_value)

    shapes = rasterio.features.shapes(valid_mask.astype(np.uint8), transform=transform)

    # Create a list of shapely polygons
    polygons = [shapely.geometry.shape(geom) for geom, value in shapes if value == 1]

    # Dissolve all polygons into a single shape
    dissolved_polygon = shapely.ops.unary_union(polygons)

    # Ensure the dissolved shape is a Polygon (buffer to close if necessary)
    if isinstance(dissolved_polygon, LineString):
        dissolved_polygon = dissolved_polygon.buffer(0)
    elif isinstance(dissolved_polygon, MultiLineString):
        dissolved_polygon = dissolved_polygon.convex_hull
    elif isinstance(dissolved_polygon, MultiPolygon):
        dissolved_polygon = dissolved_polygon.convex_hull

    if not isinstance(dissolved_polygon, Polygon):
        raise ValueError("The extracted boundary is not a valid polygon")

    # Create a GeoDataFrame from the polygon
    valid_dem_region_shp = gpd.GeoDataFrame({'ID': ['boundary'], 'geometry': [dissolved_polygon]}, crs=crs)

    # Save the outline
    valid_dem_region_shp_file = case_path / 'valid_dem_region_outline.shp'
    valid_dem_region_shp.to_file(valid_dem_region_shp_file)

    breakline_points_file = case_path / "SRVlevels_xy_labeled.csv"
    breakline_points = pd.read_csv(breakline_points_file)

    # Pre-processing of the breaklines

    # Add offset of 100 to match with DEM
    breakline_points[['X','Y']] = breakline_points[['X','Y']] + offset

    # Save Breaklines as shapefile
    breaklines_points_shp_file = case_path / (case_name +'SRVlevels_xy_labeled.shp')
    geometry = [shapely.geometry.Point(xy) for xy in zip(breakline_points['X'].to_numpy(), breakline_points['Y'].to_numpy())]
    breaklines_points_shp = gpd.GeoDataFrame({'ID': breakline_points['ID'].to_numpy()}, geometry=geometry, crs=crs)
    breaklines_points_shp.to_file(breaklines_points_shp_file)


    # Determine outliers and get the DataFrames
    breakline_points_cleaned = determine_outliers(breakline_points,factor=3)

    # Save Breaklines as shapefile
    breaklines_points_cleaned_shp_file = case_path / (case_name +'SRVlevels_xy_labeled_cleaned.shp')
    geometry = [shapely.geometry.Point(xy) for xy in zip(breakline_points_cleaned['X'].to_numpy(), breakline_points_cleaned['Y'].to_numpy())]
    breaklines_points_cleaned_shp = gpd.GeoDataFrame({'ID': breakline_points_cleaned['ID'].to_numpy()}, geometry=geometry, crs=crs)
    breaklines_points_cleaned_shp.to_file(breaklines_points_cleaned_shp_file)

    # Convert points to lines
    breaklines_lines_cleaned_shp = points_to_lines(breaklines_points_cleaned_shp)
    breaklines_lines_cleaned_shp_file = case_path / (case_name +'SRVlevels_xy_labeled_cleaned_lines.shp')
    breaklines_lines_cleaned_shp.to_file(breaklines_lines_cleaned_shp_file)

    # Simplify the lines
    breaklines_lines_cleaned_simplified_shp = simplify_lines(breaklines_lines_cleaned_shp, simplification_tolerance)
    breaklines_lines_cleaned_simplified_shp_file = case_path / (case_name +'SRVlevels_xy_labeled_cleaned_simplified_lines.shp')
    breaklines_lines_cleaned_simplified_shp.to_file(breaklines_lines_cleaned_simplified_shp_file)

    # Clip the lines with the boundary
    breaklines_lines_cleaned_simplified_clipped_shp = clip_lines_with_boundary(breaklines_lines_cleaned_simplified_shp, valid_dem_region_shp)
    breaklines_lines_cleaned_simplified_clipped_shp_file = case_path / (case_name +'SRVlevels_xy_labeled_cleaned_simplified_clipped_lines.shp')
    breaklines_lines_cleaned_simplified_clipped_shp.to_file(breaklines_lines_cleaned_simplified_clipped_shp_file)

    # Remove valley breaklines
    breaklines_lines_cleaned_simplified_clipped_channel_shp = remove_valley_breaklines(breaklines_lines_cleaned_simplified_clipped_shp)
    breaklines_lines_cleaned_simplified_clipped_channel_shp_file = case_path / (case_name +'SRVlevels_xy_labeled_cleaned_simplified_clipped_channel_lines.shp')
    breaklines_lines_cleaned_simplified_clipped_channel_shp.to_file(breaklines_lines_cleaned_simplified_clipped_channel_shp_file)

    # get inflow and outlow lines
    nodestrings = create_lines_from_points(breaklines_lines_cleaned_simplified_clipped_channel_shp)
    nodestrings_file = case_path / (case_name +'nodestrings.shp')
    nodestrings.to_file(nodestrings_file)
    nodestrings_file = case_path / ('nodestrings.csv')
    save_lines_to_csv_with_ids(nodestrings, nodestrings_file)

    breaklines_lines_cleaned_simplified_clipped_channel_closed_shp = append_lines(breaklines_lines_cleaned_simplified_clipped_channel_shp, nodestrings)
    breaklines_lines_cleaned_simplified_clipped_channel_closed_shp_file = case_path / (case_name +'breaklines.shp')
    breaklines_lines_cleaned_simplified_clipped_channel_closed_shp.to_file(breaklines_lines_cleaned_simplified_clipped_channel_closed_shp_file)
    breaklines_lines_cleaned_simplified_clipped_channel_closed_csv_file = case_path / ('breaklines.csv')
    save_lines_to_csv_with_ids(breaklines_lines_cleaned_simplified_clipped_channel_closed_shp, breaklines_lines_cleaned_simplified_clipped_channel_closed_csv_file)


    # bounding_polygon = create_bounding_polygon(breaklines_lines_cleaned_simplified_clipped_channel_closed_shp)
    # gdf_boundary_polygon = gpd.GeoDataFrame([{'geometry': bounding_polygon}], crs=breaklines_lines_cleaned_simplified_clipped_channel_closed_shp.crs)
    # file = case_path / (case_name +'bounding_polygon.shp')
    # gdf_boundary_polygon.to_file(file)


if __name__ == "__main__":
    main()

