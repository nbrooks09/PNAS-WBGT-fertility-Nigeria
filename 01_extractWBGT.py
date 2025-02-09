##################################################################################
#
#       Count Days
#       By Cascade Tuholske, cascade dot tuholske 1 at montana dot edu 
#
#       This script opens finds the daily maximum wet bulb globe temperatures (WBGTmax)
#       for a set of household clusters from DHS surveys.
#
#       The clusters are buffered by 5-km to accommodate DHS jittering and area-averaged. 
#       
#       To replicate this, you must request access to the DHS data, including the GOS data
#
#       The WBGTmax data can be found here: 
#       https://data.chc.ucsb.edu/people/cascade/UHE-daily/wbgtmax/
#
#       Note: Be sure to check NAN value and split strings on rasters before running; 
#       change as needed. 
#
#################################################################################

# Dependencies 
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point 
import glob
from multiprocessing import Pool
import time
from rasterstats import zonal_stats, gen_zonal_stats
import rasterio
import matplotlib.pyplot as plt
import time
import multiprocessing as mp 

# Functions
def km_to_d(km):
    """
    Convert a distance in kilometers to an approximate distance in decimal degrees.

    This function provides a rough conversion based on the approximation that 
    1 degree of latitude is approximately 111 km. 

    Reference: http://wiki.gis.com/wiki/index.php/Decimal_degrees

    Args:
        km (float): Distance in kilometers.

    Returns:
        float: Approximate distance in decimal degrees.
    
    """
    d = 0.01 * km / 1.11
    
    return d

def buffer_points(csv_fn, lat_col, long_col, dist):
    """
    Create a buffer around survey cluster points.

    This function reads a CSV file containing latitude and longitude coordinates,
    converts them into a GeoDataFrame, and applies a buffer of approximately the 
    given distance in kilometers (converted to decimal degrees).

    Reference: http://wiki.gis.com/wiki/index.php/Decimal_degrees

    Args:
        csv_fn (str): Path to the CSV file containing survey points.
        lat_col (str): Name of the latitude column in the CSV file.
        long_col (str): Name of the longitude column in the CSV file.
        dist (float): Buffer distance in kilometers.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame with buffered geometries.
    
    Note:
        Ensure that the latitude and longitude column names match the actual 
        column names in the CSV file before running the function.
    """
    
    # open data frame
    df = pd.read_csv(csv_fn)
    
    # convert to gdf + buffer
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[long_col], df[lat_col]), crs = 'EPSG:4326')
    gdf['geometry'] = gdf['geometry'].buffer(km_to_d(dist)) # ~5km buffer
    
    return gdf

# Use raster stats zone to assign a mode to the points
def run_zone(polys_in, rst_fn, stats_type, col_stats): #col_merge,
    """
    Calculate zonal statistics for polygons using a raster and return the results as a GeoDataFrame.

    This function applies zonal statistics to a set of input polygons using a given raster
    and returns statistical values as a GeoDataFrame. The results can be renamed for clarity.

    Args:
        polys_in (geopandas.GeoDataFrame): Input polygons for which zonal statistics will be calculated.
        rst_fn (str): Path to the raster file used for zonal statistics.
        stats_type (str): Type of statistic to compute (e.g., "mode", "mean", "sum").
        col_stats (str): New name for the statistics column in the output GeoDataFrame.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame containing the calculated zonal statistics.

    Note:
        - The function does not set a specific `nodata` value for `rasterstats`, and it defaults to -9999.
          Modify this if needed.
        - Ensure that `stats_type` matches the available statistics in `zonal_stats`.
        - The `col_merge` argument was removed but can be reinstated if needed for merging results.
    """
    

    # Run Zonal Stats, July23 changed 'nan' to -9999 
    zs_feats = zonal_stats(polys_in, rst_fn, stats=stats_type, geojson_out=True, nodata = -9999) 
    
    # Turn into GeoDataFrame 
    zgdf = gpd.GeoDataFrame.from_features(zs_feats, crs=polys_in.crs)

    # Rename columns 
    zgdf = zgdf.rename(columns={stats_type: col_stats})

    #polys_out = polys_out.merge(zgdf[[col_merge, col_stats]], on = col_merge, how = 'inner')
    
    return zgdf #polys_out

def point_extract_parallel(start):
    
    """
    Extract daily WBGTmax (Wet Bulb Globe Temperature maximum) for a list of points buffered by 5 km over a given year.

    This function runs in parallel, processing data from multiple raster files for a specified year, and extracts 
    mean values for each buffered point. The results are saved as a CSV file.

    Args:
        start (tuple): A tuple containing:
            - year (int): Year for which to extract data.
            - cluster_fn (str): File path to the CSV file with cluster points (latitude and longitude).
            - rst_path (str): Directory path containing raster files (.tif) for the specified year.
            - out_path (str): Directory path to save the output CSV file.

    Returns:
        None: Outputs a CSV file with extracted data.
    
    Note:
        - Ensure that the 'latnum' and 'longnum' columns exist in the cluster file.
        - The function assumes raster filenames include dates in the format 'max.YYYYMMDD.tif'.
        - Buffers each point by 5 km before extracting mean values from rasters.
        - Designed for use with multiprocessing for parallel execution.
    """
    
    # check process, year
    print(mp.current_process(), start[0])
    
    # args
    year = str(start[0])
    data = start[1]
    cluster_fn = start[2]
    rst_path = start[3]
    out_path = start[4]
    print(year) 
    
    # Open clusters, buffer clusters, turn into GDF, drop un-needed column
    buff = buffer_points(cluster_fn, lat_col = 'latnum', long_col = 'longnum', dist = 5)
    
    # Get rst_fns
    rst_fns = sorted(glob.glob((rst_path + year + '/*.tif')))
    print(rst_fns[0])
    
    
    # Extract data mean for each cluster
    for i, rst_fn in enumerate(rst_fns):
    
        # Get Date
        print(rst_fn)
        date = rst_fn.split('max.')[1].split('.tif')[0] # updated for Tmax
        
        print(date)
   
        stats = 'mean'
        
        if i == 0: # first date
            gdf_out = run_zone(polys_in = buff, rst_fn = rst_fn, stats_type = 'mean', col_stats = date)

        if i > 0:
            gdf_tmp = run_zone(polys_in = buff, rst_fn = rst_fn, stats_type = 'mean', col_stats = date)

            gdf_out[date] = gdf_tmp[date]
  
    # drop col, covert to df and write csv
    fn_out = os.path.join(out_path, 'ng_dhs_gps_'+data+'_'+year+'.csv') 
    df_out = pd.DataFrame(gdf_out.drop(columns='geometry'))
    df_out.to_csv(fn_out, index = False)
    
def parallel_loop(function, start_list, cpu_num):
    
    """
    Executes a given function in parallel using multiple CPU cores.

    This function leverages Python's multiprocessing capabilities to run the 
    specified function concurrently across a given number of CPU cores. The 
    function is applied to each element of the provided `start_list`.

    Args:
        function (callable): The function to execute in parallel. It should be 
            able to accept elements from the `start_list` as arguments.
        start_list (list): A list of arguments to pass to the `function` in parallel.
        cpu_num (int): The number of CPU cores to utilize for parallel processing.

    Note:
        - Ensure that the `function` is defined at the top level of the module, 
          as `multiprocessing` requires functions to be pickleable.
        - The `start_list` should contain all the required inputs for the function 
          to operate correctly.
    """
    
    start = time.time()
    pool = Pool(processes = cpu_num)
    pool.map(function, start_list)
    pool.close()

    end = time.time()
    print(end-start)

# Main 
if __name__ == "__main__":
    
    # Make Lists 
    year_list = list(range(1983,2016+1)) # arg 0
    
    data = 'wbgtmax' # arg 1
    data_list = [data] * len(year_list)
    
    fn = os.path.join('../data/raw/dhs/ng_dhs_gps.csv')  # arg 2
    cluster_fn_list = [fn] * len(year_list)
    
    path = os.path.join('PATH/TO/FILES_IN')  # arg 3
    
    path = os.path.join('PATH/TO/FILES_OUT') # arg 4
    out_path_list = [path] * len(year_list)
    
    # Zip them into a list to pass to parallel func
    args_list = list(zip(year_list, data_list, cluster_fn_list, rst_path_list, out_path_list))
                     

     
    # Run it
    parallel_loop(function = point_extract_parallel, start_list = args_list, cpu_num = 16) # number of CPUs available 