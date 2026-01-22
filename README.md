# Dynamic Flood Prediction Model
Python-based automated workflow for dynamic flood risk and inundation modeling using WRF meteorological data, MCDM (AHP), and GIS spatial analysis. Developed as part of Ph.D. research at Ondokuz Mayıs University.


## Overview
This repository contains a Python-based automated geospatial workflow developed during my Ph.D. research at **Ondokuz Mayıs University**. The model is designed to process high-resolution meteorological data (WRF) and perform dynamic flood risk mapping using Multi-Criteria Decision Making (MCDM) techniques.

The project demonstrates advanced proficiency in **ArcPy**, automated data processing, and temporal GIS analysis—skills directly applicable to large-scale environmental monitoring and river mobility studies.



## Key Features
* **Automated NetCDF Processing:** Batch processing of WRF (Weather Research and Forecasting) outputs to extract atmospheric variables such as pressure, temperature, and precipitation.
* **Temporal Intensity Analysis:** Algorithmic calculation of hourly rainfall intensity from cumulative data.
* **Multi-Criteria Decision Analysis (MCDA):** Implementation of the **Analytic Hierarchy Process (AHP)** to weigh static (elevation, slope, land use, geology) and dynamic (hourly rainfall, humidity) criteria.
* **Spatial Interpolation:** Automated Inverse Distance Weighting (IDW) to generate high-resolution (100m) continuous surfaces from point data.
* **Temporal Mosaic Integration:** Creation of time-enabled Mosaic Datasets for dynamic visualization and simulation of flood risk over time.

## Methodology
The model follows a structured pipeline:
1. **Data Ingestion:** Reading raw `.nc` files and exporting formatted feature classes to a Geodatabase.
2. **Dynamic Parameter Calculation:** Using Python math libraries to derive specific humidity and hourly rainfall.
3. **Normalization:** Scaling all input rasters to a 0-1 range for multi-criteria compatibility.
4. **Weighted Sum Integration:** Combining 9 different geographical and meteorological layers based on calculated influence weights.
5. **Post-Processing:** Generating temporal overviews for seamless data visualization in GIS environments.

## Technical Stack
* **Language:** Python 3.x
* **Core Library:** ArcPy (ArcGIS Pro)
* **Methods:** AHP, IDW, Spatial Masking, Normalization.

## Research Context
This tool was a core component of my doctoral thesis: *"Development of a Semi-Dynamic Flood and Inundation Prediction Model Based on Geographic Information Systems"*. It highlights my ability to bridge the gap between atmospheric science and geospatial engineering.

## Author
**Dr. Rıdvan Ertuğrul YILDIRIM** 
Geomatics Engineer 
Email: ridvan.yildirim@omu.edu.tr
