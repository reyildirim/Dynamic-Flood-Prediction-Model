import arcpy
import os
import arcpy.mp as mp
from arcpy import env
from arcpy.sa import *
import time
import datetime

# Record the start time
start_time = time.time()

env.workspace = r"D:\tez_ver_10\tez_ver_10.gdb"

# Open the current ArcGIS Pro project
proj = mp.ArcGISProject("CURRENT")

# Get all maps within the project
maps = proj.listMaps()

# Path to the folder containing WRF NetCDF files
netcdf_folder_path = r"D:\tez_ver_10\WRF_NC"

# Path to the polygon layer for clipping buffer zones
clip_polygon_path = r"D:\tez_ver_10\tez_ver_10.gdb\Samsun_Buffer5k"

# Path to the study area boundary polygon
study_area_mask_path = r"D:\tez_ver_10\tez_ver_10.gdb\Samsun"

def clear_table_of_contents():
    """Cleans up temporary and intermediate layers from the Table of Contents."""
    for map_obj in maps:
        layers = map_obj.listLayers()

        # Remove layers starting with specific prefixes to keep the workspace clean
        for layer in layers:
            if (layer.name.startswith("f_") or layer.name.startswith("n_") 
                or layer.name.startswith("s_") or layer.name.startswith("p_") 
                or layer.name.startswith("normalized_") 
                or layer.name.startswith("out_")):
                map_obj.removeLayer(layer)       

    print("Table of Contents cleared.")


### STEP 1: Process NetCDF files and export to Geodatabase (GDB) ###

# Filter NetCDF files in the directory
netcdf_file_list = [f for f in os.listdir(netcdf_folder_path) if f.endswith(".nc")]

# Loop through NetCDF files to extract meteorological variables
for netcdf_file in netcdf_file_list:
    netcdf_file_path = os.path.join(netcdf_folder_path, netcdf_file)
    
    # Extract file metadata and format the output name
    base_name = os.path.splitext(netcdf_file)[0]
    short_name = base_name[-13:].replace("-", "_")
    export_gdb_name = "p_{}".format(short_name)
    
    # Create a new Point Feature Class in the GDB
    output_feature_class = arcpy.CreateFeatureclass_management(env.workspace, export_gdb_name, "POINT")
    print(f"Created {export_gdb_name} in Geodatabase.")
    
    temp_layer_name = "f_{}".format(short_name)
            
    # Create NetCDF feature layer (Pressure, Specific Humidity, Cumulative Rain, Temperature)
    arcpy.md.MakeNetCDFFeatureLayer(netcdf_file_path, "PSFC;Q2;RAINNC;T2", "XLONG", "XLAT", temp_layer_name, "west_east;south_north", '', '', None, "BY_VALUE")
    print(f"Added NetCDF layers: {temp_layer_name}")
    
    # Save the feature layer as permanent features in GDB
    arcpy.CopyFeatures_management(temp_layer_name, output_feature_class)
    print(f"Exported {export_gdb_name} to GDB.")
    
    # Perform clipping for the specific study buffer
    clipped_feature_class = arcpy.analysis.Clip(export_gdb_name, clip_polygon_path, f"s_{export_gdb_name}")
    print(f"Clipped {export_gdb_name} features.")
    
    clear_table_of_contents()
    print("Processing step for current file completed.")
    
print("STEP 1: All WRF files processed and stored.")
clear_table_of_contents()
    
      
### STEP 2: Calculate specific Humidity field for all WRF data ###

prefix_s = "s_"  
field_humidity = "nem" 
field_type = "DOUBLE"

# Humidity calculation formula derived from WRF variables
humidity_expression = "!Q2! / ((379.90516/!PSFC!)*math.exp(17.2693882*(!T2!-273.16)/(!T2!-35.86)))"

# List all relevant feature classes
feature_class_list = arcpy.ListFeatureClasses(feature_dataset=None, wild_card=f"{prefix_s}*")

for fc in feature_class_list:
    try:
        fields = [f.name for f in arcpy.ListFields(fc)]

        if field_humidity not in fields:
            arcpy.AddField_management(fc, field_humidity, field_type)
            print(f"Added {field_humidity} field to {fc}.")
        else:
            print(f"{field_humidity} field already exists in {fc}.")
        
        # Run calculation
        arcpy.CalculateField_management(fc, field_humidity, humidity_expression, "PYTHON")
        print(f"Calculated humidity for {fc}.")
        
    except Exception as e:
        print(f"Calculation Error in Step 2: {e}")
        
print("STEP 2: Humidity calculations finished.")


### STEP 3: Initialize Precipitation field ###

field_precip = "yagis"

feature_class_list = arcpy.ListFeatureClasses(feature_dataset=None, wild_card=f"{prefix_s}*")

for fc in feature_class_list:
    try:
        fields = [f.name for f in arcpy.ListFields(fc)]

        if field_precip not in fields:
            arcpy.AddField_management(fc, field_precip, field_type)
            arcpy.CalculateField_management(fc, field_precip, 0, "PYTHON")
            print(f"Initialized precipitation field for {fc}.")
       
    except Exception as e:
        print(f"Error in Step 3: {e}")
        
print("STEP 3: Precipitation fields ready.")


### STEP 4: Calculate Hourly Precipitation Intensity (Temporal Difference) ###

target_fcs = arcpy.ListFeatureClasses(feature_dataset=None, wild_card=f"{prefix_s}*")

field_cum_rain = "RAINNC"

# Iterative calculation of hourly rainfall by subtracting cumulative values of previous time step
for i, target_fc in enumerate(target_fcs):
    if i > 0:
        source_fc = target_fcs[i - 1]
        
        with arcpy.da.UpdateCursor(target_fc, [field_cum_rain, "yagis"]) as update_cursor:
            with arcpy.da.SearchCursor(source_fc, [field_cum_rain]) as search_cursor:
                for row_target, row_source in zip(update_cursor, search_cursor):
                    # Hourly Rain = Current Cumulative - Previous Cumulative
                    intensity = (row_target[0] - row_source[0])
                    row_target[1] = intensity
                    update_cursor.updateRow(row_target)

    print(f"Calculated hourly rainfall intensity for {target_fc}.")

print("STEP 4: Temporal intensity analysis completed.")


### STEP 5: Spatial Interpolation (IDW) for Rainfall and Humidity ###

# Define output cell size (100m resolution)
arcpy.env.cellSize = 0.00117238464 

# Interpolate Rainfall
field_to_analyze = "yagis"
feature_class_list = arcpy.ListFeatureClasses(feature_dataset=None, wild_card=f"{prefix_s}*")

for fc in feature_class_list:
    try:
        output_raster_name = "yagis_" + fc.replace(" ", "_")
        arcpy.sa.Idw(fc, field_to_analyze).save(f"{env.workspace}/{output_raster_name}")
        print(f"Generated Rainfall IDW for {fc}.")
    except Exception as e:
        print(f"IDW Error: {e}")
        
# Interpolate Humidity
field_to_analyze = "nem"
for fc in feature_class_list:
    try:
        output_raster_name = "nem_" + fc.replace(" ", "_")
        arcpy.sa.Idw(fc, field_to_analyze).save(f"{env.workspace}/{output_raster_name}")
        print(f"Generated Humidity IDW for {fc}.")
    except Exception as e:
        print(f"IDW Error: {e}")

print("STEP 5: Spatial interpolation outputs generated.")


### STEP 6 & 7: Masking/Clipping Rasters to Study Area ###

def clip_rasters(pattern, prefix_out):
    rasters = arcpy.ListRasters(wild_card=f"{pattern}*")
    for r in rasters:
        out_name = f"{prefix_out}{r}"
        out_path = os.path.join(env.workspace, out_name)
        clipped = ExtractByMask(r, study_area_mask_path)
        clipped.save(out_path)
        print(f"Masked raster saved: {out_name}")

clip_rasters("nem_", "s_")
print("STEP 6: Humidity masking completed.")

clip_rasters("yagis_", "s_")
print("STEP 7: Precipitation masking completed.")


### STEP 8: Raster Normalization (Scale 0-1) ###

# Normalize Precipitation (Based on a fixed threshold of 96.7 mm)
rain_rasters = arcpy.ListRasters(wild_card="s_yagis_*")
for r in rain_rasters:
    in_raster = Raster(r)
    norm_name = "norm_{}".format(r)
    normalized = (in_raster - in_raster.minimum) / (96.7)
    normalized.save(norm_name)

# Normalize Humidity (Min-Max Scaling)
humidity_rasters = arcpy.ListRasters(wild_card="s_nem_*")
for r in humidity_rasters:
    in_raster = Raster(r)
    norm_name = "norm_{}".format(r)
    normalized = (in_raster - in_raster.minimum) / (in_raster.maximum - in_raster.minimum)
    normalized.save(norm_name)

print("STEP 8: Raster normalization to 0-1 scale finished.")


### STEP 9: Multi-Criteria Weighted Sum (Dynamic Flood Risk Calculation) ###

# AHP-derived weights for dynamic and static criteria
w_rain = 0.27487649
w_slope = 0.264867894
w_stream = 0.117766046
w_landuse = 0.110783599
w_humidity = 0.066558585
w_soil = 0.066233803
w_geology = 0.038041881
w_elevation = 0.035812175
w_aspect = 0.025059528

# Static normalized criteria layers
lyr_elevation = "norm_s_yukseklik"
lyr_slope = "norm_s_egim"
lyr_aspect = "norm_s_baki"
lyr_stream = "norm_s_akarsu"
lyr_landuse = "norm_s_ak"
lyr_geology = "norm_s_jeoloji"
lyr_soil = "norm_s_toprak"

# Spatial Reference
utm_proj = arcpy.SpatialReference(32636) # WGS 84 / UTM zone 36N

def run_dynamic_weighted_sum():
    rain_list = arcpy.ListRasters(wild_card="norm_s_yagis_s*")
    hum_list = arcpy.ListRasters(wild_card="norm_s_nem_s*") 
        
    for rain_r in rain_list:
        time_suffix = rain_r[-13:]
        matched_hum = next((h for h in hum_list if h.endswith(time_suffix)), None)

        if matched_hum:
            output_risk_map = "SUM_" + time_suffix.replace(" ", "_")
            
            # Combine static and dynamic criteria using Weighted Sum
            ws_table = WSTable([[rain_r, "VALUE", w_rain],
                                [matched_hum, "VALUE", w_humidity],
                                [lyr_elevation, "VALUE", w_elevation],
                                [lyr_slope, "VALUE", w_slope],
                                [lyr_aspect, "VALUE", w_aspect],
                                [lyr_stream, "VALUE", w_stream],
                                [lyr_landuse, "VALUE", w_landuse],
                                [lyr_geology, "VALUE", w_geology],
                                [lyr_soil, "VALUE", w_soil]])

            risk_output = WeightedSum(ws_table)
            risk_output.save(output_risk_map)
            print(f"Flood Risk Map generated for timestamp: {time_suffix}")

run_dynamic_weighted_sum()
print("STEP 9: Dynamic flood risk mapping complete.")
clear_table_of_contents()


### STEP 10, 11 & 12: Mosaic Dataset and Temporal Management ###

mosaic_name = "FloodRisk_TimeDataset"
arcpy.CreateMosaicDataset_management(env.workspace, mosaic_name, utm_proj)
arcpy.management.AddRastersToMosaicDataset(mosaic_name, "Raster Dataset", arcpy.ListRasters("SUM_*"))

# Add Time metadata for temporal visualization
arcpy.AddField_management(mosaic_name, "Timestamp", "DATE")
base_date = datetime.datetime(2023, 7, 10, 0, 0, 0)

with arcpy.da.UpdateCursor(mosaic_name, ["OBJECTID", "Timestamp"]) as cursor:
    for row in cursor:
        # Increment time by hour based on ObjectID
        row[1] = base_date + datetime.timedelta(hours=row[0] - 1)
        cursor.updateRow(row)

# Optimize for viewing
arcpy.management.BuildOverviews(mosaic_name)

end_time = time.time()
print(f"Workflow completed in {(end_time - start_time)/60:.2f} minutes.")
print("Dynamic Flood Risk Simulation Finished Successfully.")