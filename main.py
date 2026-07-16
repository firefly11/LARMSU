# -*- coding: utf-8 -*-
import os
import sys
from site import addsitedir
from sys import executable
from os import path
from xml.etree.ElementTree import ElementTree

interpreter = executable
sitepkg = path.dirname(interpreter) + "\\site-packages"
print(sitepkg)
addsitedir(sitepkg)

import numpy as np
import tkinter as tk
from tkinter import filedialog, scrolledtext, ttk
import tkMessageBox
import arcpy
from arcpy.sa import *
from osgeo import gdal
import osgeo.osr as osr
from osgeo.gdalconst import *
import core, core_dinf
import math
from threading import Thread


# Define the numeric validation function
def validate_number_input(input_value):
    if input_value == "":
        return True
    return input_value.isdigit()


# print the console text to log textbox
class RedirectText:
    def __init__(self, textbox):
        self.textbox = textbox

    def write(self, message):
        self.textbox.insert(tk.END, message)
        self.textbox.see(tk.END)

    def flush(self):
        pass


# Function to check if all widgets are filled
def validate_inputs():
    # Check if the required fields are filled
    if not shapefile_path_var.get():
        tkMessageBox.showwarning("Input Error", "Shapefile path cannot be empty!")
        return False
    if not tiff_path_var.get():
        tkMessageBox.showwarning("Input Error", "TIFF path cannot be empty!")
        return False
    if not param_demRes.get():
        tkMessageBox.showwarning("Input Error", "Parameters-demRes cannot be empty!")
        return False
    if not param_travelAg.get():
        tkMessageBox.showwarning("Input Error", "Parameters-travelAg cannot be empty!")
        return False
    if not param_maxV.get():
        tkMessageBox.showwarning("Input Error", "Parameters-maxV cannot be empty!")
        return False
    if not param_diffAg.get():
        tkMessageBox.showwarning("Input Error", "Parameters-diffAg cannot be empty!")
        return False
    if not param_mcTimes.get():
        tkMessageBox.showwarning("Input Error", "Parameters-mcTimes cannot be empty!")
        return False
    if not param_toCF.get():
        tkMessageBox.showwarning("Input Error", "Parameters-toCF cannot be empty!")
        return False
    if not option_var.get():
        tkMessageBox.showwarning("Input Error", "Please select an option!")
        return False
    if not outfile_dir_var.get():
        tkMessageBox.showwarning("Input Error", "Output directory cannot be empty!")
        return False
    if not param_proName.get():
        tkMessageBox.showwarning("Input Error", "proName cannot be empty!")
        return False

    return True


def select_shapefile():
    filepath = filedialog.askopenfilename(
        filetypes=[("Shapefiles", "*.shp")],
        title="Select a Shapefile"
    )
    if filepath:
        shapefile_path_var.set(filepath)


# Function for selecting a TIFF file
def select_tiff():
    filepath = filedialog.askopenfilename(
        filetypes=[("TIFF files", "*.tif")],
        title="Select a TIFF File"
    )
    if filepath:
        tiff_path_var.set(filepath)


# Outfile directory picker label and button
def select_outdir():
    outdir = filedialog.askdirectory(
        title="Select Output Directory"
    )
    if outdir:
        outfile_dir_var.set(outdir)


# Function to execute the main logic
def execute(gridNew, aspect, flow_direction, slope, aspectZonal_ali):
    shapefile = shapefile_path_var.get()
    tiff_file = tiff_path_var.get()
    selected_option = option_var.get()
    outfile_path = outfile_dir_var.get()
    if not validate_inputs():
        return

    if selected_option == 'Dinf':
        try:
            core_dinf.SUsM_dinf(param_proName.get(), shapefile, tiff_file, int(param_demRes.get()),
                                math.tan(math.radians(float(param_travelAg.get()))),
                                int(param_maxV.get()), int(param_diffAg.get()), int(param_mcTimes.get()),
                                outfile_path, gridNew, aspect, flow_direction,
                                slope, aspectZonal_ali)
        except Exception as e:
            print(e.message)
    else:
        try:
            core.SUsM(param_proName.get(), shapefile, tiff_file, int(param_demRes.get()),
                      math.tan(math.radians(float(param_travelAg.get()))),
                      int(param_maxV.get()), int(param_diffAg.get()), int(param_mcTimes.get()),
                      outfile_path, gridNew, aspect, slope, aspectZonal_ali)
        except Exception as e:
            print(e.message)
    print("Done.")

    btn_execute["state"] = tk.NORMAL
    btn_execute["text"] = "Execute"
    root.update()


def mappingSUs(pro_name, DEMPath, dx, flow_acc, outFilePath):
    outFilePath = outFilePath + '\\' + pro_name
    filePath = os.path.dirname(DEMPath)
    DEMName = os.path.basename(DEMPath)

    # Process data path setting
    if not os.path.exists(outFilePath):
        os.makedirs(outFilePath)
    slopePath = outFilePath + '\\' + 'slope.tif'
    aspectPath = outFilePath + '\\' + 'aspect.tif'
    aspectZonal = outFilePath + '\\' + 'aspect_zonal.tif'
    gridInvalidPath = outFilePath + '\\' + 'grid_invalid.tif'
    alignment_path = outFilePath + '\\' + 'alignment_zonal.tif'

    # arcpy setting
    arcpy.env.workspace = filePath
    arcpy.env.overwriteOutput = True
    # check SRS
    sr = arcpy.Describe(DEMName).spatialReference
    print("Spatial Reference System:" + sr.name)
    # check out any necessary licenses
    print ("Spatial Analyst Extension Available")
    print(arcpy.CheckOutExtension("spatial"))

    # register gdal
    gdal.AllRegister()

    # *************
    # Process.1 generate SUs
    # *************
    DEM = arcpy.Raster(DEMName)
    # fill Sinks
    fill = Fill(DEM)
    # generate FlowDirection
    flowdir = FlowDirection(fill, "NORMAL")
    # generate FlowAccumulation
    flowacc = FlowAccumulation(flowdir, "", "FLOAT")
    # extract streams
    streams = SetNull(flowacc, 1, "VALUE <= " + str(flow_acc))  # flowacc <=90 -> null, 90+ -> 1
    # streamlink
    streamlink = StreamLink(streams, flowdir)
    # watershed
    watershedrs = Watershed(flowdir, streamlink, "VALUE")
    arcpy.RasterToPolygon_conversion(watershedrs, "watershed", "NO_SIMPLIFY", "VALUE")  # watershed polygon saved

    # repeat the above steps with reversed DEM
    rDEM = DEM.maximum - DEM
    rfill = Fill(rDEM)
    rflowdir = FlowDirection(rfill, "NORMAL")
    rflowacc = FlowAccumulation(rflowdir, "", "FLOAT")
    rstreamrs = SetNull(rflowacc, 1, "VALUE <= " + str(flow_acc))  # flowacc <=90 -> null, 90+ -> 1
    rstreamlink = StreamLink(rstreamrs, rflowdir)
    rwatershedrs = Watershed(rflowdir, rstreamlink, "VALUE")
    arcpy.RasterToPolygon_conversion(rwatershedrs, "rwatershed", "NO_SIMPLIFY", "VALUE")  # watershed poly
    arcpy.Union_analysis(["watershed.shp", "rwatershed.shp"], "slopeunits")
    print("Process.1 done, Check 'slopeunits' in selected DEM Workspace.")
    root.update()
    # *************
    # Process.2 calculate slope, aspect and flow direction.
    # *************
    data = gdal.Open(DEMPath, GA_ReadOnly)
    if data is None:
        print('Cannot open this file:' + DEMPath)
        sys.exit(1)
    Ref = osr.SpatialReference()
    Ref.ImportFromWkt(data.GetProjectionRef())
    transform = data.GetGeoTransform()

    # project
    gridNew = arcpy.RasterToNumPyArray(fill, nodata_to_value=0)
    gridNew = gridNew.astype(np.float)
    dx1, dy1, dx2, dy2 = core_dinf.calcFiniteSlopes(gridNew, dx)
    slope, aspect = core_dinf.CacSlopAsp(dx1, dy1, dx2, dy2)
    grid_invalid = np.full((gridNew.shape[0], gridNew.shape[1]), -999)
    flow_direction = core_dinf.calculate_flow_direction(gridNew)
    flow_direction = (360 - flow_direction) % 360
    # output
    driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(slopePath):
        os.remove(slopePath)
    if os.path.exists(aspectPath):
        os.remove(aspectPath)
    if os.path.exists(gridInvalidPath):
        os.remove(gridInvalidPath)

    ds1 = driver.Create(slopePath, slope.shape[1], slope.shape[0], 1, GDT_Float32)
    ds1.SetProjection(Ref.ExportToWkt())
    ds1.SetGeoTransform(transform)
    band = ds1.GetRasterBand(1)
    band.WriteArray(slope, 0, 0)

    ds2 = driver.Create(aspectPath, aspect.shape[1], aspect.shape[0], 1, GDT_Float32)
    ds2.SetProjection(Ref.ExportToWkt())
    ds2.SetGeoTransform(transform)
    band = ds2.GetRasterBand(1)
    band.WriteArray(aspect, 0, 0)

    ds3 = driver.Create(gridInvalidPath, grid_invalid.shape[1], grid_invalid.shape[0], 1, GDT_Float32)
    ds3.SetProjection(Ref.ExportToWkt())
    ds3.SetGeoTransform(transform)
    band = ds3.GetRasterBand(1)
    band.WriteArray(grid_invalid, 0, 0)

    del ds1
    del ds2
    del ds3
    print("Process.2 done, Check 'aspect.tif, slope.tif' in Output Workspace.")
    root.update()
    # *************
    # Process.3 ZonalStatistics
    # *************
    arcpy.gp.ZonalStatistics_sa("slopeunits.shp", "FID", aspectPath, aspectZonal, "MEAN", "NODATA")

    print("Process.3 done, Check 'aspect_zonal' in Output Workspace.")
    root.update()
    core_dinf.alignment_data(gridInvalidPath, aspectZonal, alignment_path)
    aspectZonal_ali = arcpy.Raster(alignment_path)
    aspectZonal_ali = arcpy.RasterToNumPyArray(aspectZonal_ali, nodata_to_value=0)
    aspectZonal_ali = aspectZonal_ali.astype(np.float)
    return gridNew, aspect, flow_direction, slope, aspectZonal_ali


def start_task():
    print("Program running...")
    btn_execute["state"] = tk.DISABLED
    btn_execute["text"] = "modelling"
    root.update()

    tiff_file = tiff_path_var.get()
    outfile_path = outfile_dir_var.get()
    gridNew, aspect, flow_direction, slope, aspectZonal_ali = mappingSUs(param_proName.get(), tiff_file,
                                                                         int(param_demRes.get()),
                                                                         int(param_toCF.get()), outfile_path)

    Thread(target=execute(gridNew, aspect, flow_direction, slope, aspectZonal_ali)).start()


# Create main window
root = tk.Tk()
root.title("LARMSU Version 1.0")

# Register validation functions
validate_cmd = root.register(validate_number_input)

# Variables
shapefile_path_var = tk.StringVar()
tiff_path_var = tk.StringVar()
outfile_dir_var = tk.StringVar()
param_demRes = tk.StringVar()
param_travelAg = tk.StringVar()
param_maxV = tk.StringVar()
param_diffAg = tk.StringVar()
param_mcTimes = tk.StringVar()
param_toCF = tk.StringVar()
param_proName = tk.StringVar()
option_var = tk.StringVar(value="Dinf")

# Layout
frame = ttk.Frame(root, padding=10)
frame.grid(row=0, column=0, sticky="NSEW")

# Shapefile Picker
ttk.Label(frame, text="Select Source Area: ").grid(row=0, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=shapefile_path_var, width=50).grid(row=0, column=1, sticky="W", padx=2, pady=2)
ttk.Button(frame, text="Browse", command=select_shapefile).grid(row=0, column=2, sticky="W", padx=2, pady=2)

# TIFF Picker
ttk.Label(frame, text="Select DEM: ").grid(row=1, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=tiff_path_var, width=50).grid(row=1, column=1, sticky="W", padx=2, pady=2)
ttk.Button(frame, text="Browse", command=select_tiff).grid(row=1, column=2, sticky="W", padx=2, pady=2)

# Parameter Input Boxes group
# Dem resolution
ttk.Label(frame, text="Dem resolution (m): ").grid(row=2, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_demRes, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=2, column=1, sticky="W", padx=2, pady=2)
# Travel angle
ttk.Label(frame, text="Travel angle (°): ").grid(row=3, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_travelAg, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=3, column=1, sticky="W", padx=2, pady=2)
# Max velocity
ttk.Label(frame, text="Max velocity (m/s): ").grid(row=4, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_maxV, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=4, column=1, sticky="W", padx=2, pady=2)
# Diffusion angle
ttk.Label(frame, text="Diffusion angle (°): ").grid(row=5, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_diffAg, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=5, column=1, sticky="W", padx=2, pady=2)
# MCMC times
ttk.Label(frame, text="MCMC times: ").grid(row=6, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_mcTimes, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=6, column=1, sticky="W", padx=2, pady=2)
# ToCF
ttk.Label(frame, text="ToCF: ").grid(row=7, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=param_toCF, width=30, validate="key",  # 在键入时触发验证
          validatecommand=(validate_cmd, "%P")).grid(row=7, column=1, sticky="W", padx=2, pady=2)

# Propagation direction (Radio Buttons)
ttk.Label(frame, text="Direction method:").grid(row=8, column=0, sticky="W", padx=2, pady=2)
options = ["Dinf", "aspect"]
for i, option in enumerate(options):
    ttk.Radiobutton(frame, text=option, variable=option_var, value=option).grid(row=8, column=1 + i, sticky="W", padx=2,
                                                                                pady=2)

# Outfile directory picker label and button
ttk.Label(frame, text="Select Output Directory:").grid(row=9, column=0, sticky="W", padx=2, pady=2)
ttk.Entry(frame, textvariable=outfile_dir_var, width=50).grid(row=9, column=1, sticky="W", padx=2, pady=2)
ttk.Button(frame, text="Browse", command=select_outdir).grid(row=9, column=2, sticky="W", padx=2, pady=2)

# Log Textbox
log_textbox = scrolledtext.ScrolledText(frame, width=70, height=15, wrap=tk.WORD)
log_textbox.grid(row=10, column=0, columnspan=3, pady=10)

# ProjectName
ttk.Label(frame, text="ProjectName: ").grid(row=11, column=0, sticky="W")
ttk.Entry(frame, textvariable=param_proName, width=30).grid(row=11, column=1, sticky="W")

# Execute Button
btn_execute = ttk.Button(frame, text="Execute", command=start_task)
btn_execute.grid(row=11, column=2, columnspan=3, pady=10)

sys.stdout = RedirectText(log_textbox)

# Main loop
root.mainloop()
