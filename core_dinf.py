# -*- coding:utf-8 -*-

import os
import sys
import math
import numpy as np
import pandas as pd
import mcmc
import arcpy
from osgeo import gdal
from osgeo import ogr
import osgeo.osr as osr
from osgeo.gdalconst import *


# shp to xy
def shp_to_xy(shpPath, tifPath, xyPath):
    outPath = xyPath
    dataset = gdal.Open(tifPath, gdal.GA_ReadOnly)

    geo_transform = dataset.GetGeoTransform()
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize

    shp = ogr.Open(shpPath, 0)
    m_layer = shp.GetLayerByIndex(0)
    target_ds = gdal.GetDriverByName('GTiff').Create(outPath, xsize=cols, ysize=rows, bands=1,
                                                     eType=gdal.GDT_Byte)
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(dataset.GetProjection())

    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.FlushCache()
    gdal.RasterizeLayer(target_ds, [1], m_layer)
    del dataset
    del target_ds


def assignBCs(elevGrid):
    ny, nx = elevGrid.shape
    Zbc = np.zeros((ny + 2, nx + 2))
    Zbc[1:-1, 1:-1] = elevGrid

    Zbc[0, 1:-1] = elevGrid[0, :]
    Zbc[-1, 1:-1] = elevGrid[-1, :]
    Zbc[1:-1, 0] = elevGrid[:, 0]
    Zbc[1:-1, -1] = elevGrid[:, -1]

    Zbc[0, 0] = elevGrid[0, 0]
    Zbc[0, -1] = elevGrid[0, -1]
    Zbc[-1, 0] = elevGrid[-1, 0]
    Zbc[-1, -1] = elevGrid[-1, 0]

    return Zbc


# Calculate the dx,dy
def calcFiniteSlopes(elevGrid, dx):
    sizex = dx
    sizey = dx
    npgrid = assignBCs(elevGrid)

    nx, ny = npgrid.shape
    s_dx = np.zeros((nx, ny))
    s_dy = np.zeros((nx, ny))
    a_dx = np.zeros((nx, ny))
    a_dy = np.zeros((nx, ny))
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            s_dx[i, j] = ((npgrid[i - 1, j + 1] + 2 * npgrid[i, j + 1] + npgrid[i + 1, j + 1]) - (
                    npgrid[i - 1, j - 1] + 2 * npgrid[i, j - 1] + npgrid[i + 1, j - 1])) / (8 * sizex)
            s_dy[i, j] = ((npgrid[i + 1, j - 1] + 2 * npgrid[i + 1, j] + npgrid[i + 1, j + 1]) - (
                    npgrid[i - 1, j - 1] + 2 * npgrid[i - 1, j] + npgrid[i - 1, j + 1])) / (8 * sizey)

    a_dx = s_dx * sizex
    a_dy = s_dy * sizey
    s_dx = s_dx[1:-1, 1:-1]
    s_dy = s_dy[1:-1, 1:-1]
    a_dx = a_dx[1:-1, 1:-1]
    a_dy = a_dy[1:-1, 1:-1]

    return s_dx, s_dy, a_dx, a_dy


# Calculate slope and aspect using dx dy.
def CacSlopAsp(s_dx, s_dy, a_dx, a_dy):
    # slope convert to (°)
    slope = (np.arctan(np.sqrt(s_dx * s_dx + s_dy * s_dy))) * 180 / math.pi

    # aspect
    a = np.zeros((a_dy.shape[0], a_dy.shape[1]))
    for i in range(0, a_dx.shape[0]):
        for j in range(0, a_dx.shape[1]):
            a[i, j] = math.atan2(a_dy[i, j], -a_dx[i, j]) * 180 / math.pi
    aspect = a
    x, y = a.shape[0], a.shape[1]
    for m in range(0, x):
        for n in range(0, y):
            if a[m, n] < 0:
                aspect[m, n] = 90 - a[m, n]
            elif a[m, n] > 90:
                aspect[m, n] = 360.0 - a[m, n] + 90.0
            else:
                aspect[m, n] = 90.0 - a[m, n]

    return slope, aspect


# Energy change in downhill movement (straight line)
def downhill1(w, slope, u, maxV, dx):
    slope = slope / 180 * np.pi
    w1 = np.tan(slope) - u
    maxW = (maxV * maxV) / (2 * 9.8 * dx)
    wt = w1 + w
    if wt > maxW:
        wt = maxW

    if wt > 0:
        return wt
    else:
        return 0


# Energy change in downhill movement (diagonal line)
def downhill2(w, slope, u, maxV, dx):
    slope = slope / 180 * np.pi
    w1 = np.tan(slope) - u * np.sqrt(2)
    wt = w1 + w
    maxW = (maxV * maxV) / (2 * 9.8 * dx)
    if wt > maxW:
        wt = maxW

    if wt > 0:
        return wt
    else:
        return 0


# Energy change in uphill movement (straight line)
def uphill1(w, slope, u, maxV, dx):
    slope = slope / 180 * np.pi
    w1 = (-np.tan(slope)) + (-u)
    wt = w1 + w
    maxW = (maxV * maxV) / (2 * 9.8 * dx)
    if wt > maxW:
        wt = maxW
    if wt > 0:
        return wt
    else:
        return 0


# Energy change in uphill movement (diagonal line)
def uphill2(w, slope, u, maxV, dx):
    slope = slope / 180 * np.pi
    w1 = (-np.tan(slope)) + (- u * np.sqrt(2))
    wt = w1 + w
    maxW = (maxV * maxV) / (2 * 9.8 * dx)
    if wt > maxW:
        wt = maxW
    if wt > 0:
        return wt
    else:
        return 0


def findOutline(input):
    index = np.where(input != 0)
    x = index[0]
    list_x = np.split(x, np.nonzero(np.diff(x))[0] + 1)
    length = 0
    xxx = []
    result = np.zeros(input.shape)
    for list in list_x:
        xx = [0, 0, 0]
        y = index[1][length:(length + list.shape[0]), ]
        xx[0] = list[0]
        xx[1] = np.min(y)
        xx[2] = np.max(y)
        xxx.append(xx)
        length = list.shape[0] + length
    for i in xxx:
        for j in range(i[1], i[2] + 1):
            result[i[0], j] = 1
    return result


def outlinerReplaceByBoxplot(input):
    list1 = list(np.array(input).flatten())
    array1 = np.array(list1)
    array1 = array1[np.where(array1 != 0)]
    dt = pd.Series(array1)
    Q1 = dt.quantile(q=0.25)
    Q3 = dt.quantile(q=0.75)
    IQR = Q3 - Q1
    UL = Q3 + 1.5 * IQR
    replace_vaule = input[input < UL].max()
    input[input > UL] = replace_vaule
    return input


def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range


# Data Alignment
def alignment_data(null_raster, input_raster, outputpath):
    merge_files = [null_raster, input_raster]
    merge_options = gdal.WarpOptions(format='GTiff', resampleAlg='max')
    dsMer = gdal.Warp(outputpath, merge_files,  options=merge_options)
    dsMer = None


# D∞
def calculate_flow_direction(dem):
    rows, cols = dem.shape
    flow_direction = np.zeros((rows, cols), dtype=np.float32)

    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            dzdx = (dem[i + 1, j - 1] + 2 * dem[i + 1, j] + dem[i + 1, j + 1] - dem[i - 1, j - 1] - 2 * dem[i - 1, j] -
                    dem[i - 1, j + 1]) / 8
            dzdy = (dem[i - 1, j + 1] + 2 * dem[i, j + 1] + dem[i + 1, j + 1] - dem[i - 1, j - 1] - 2 * dem[i, j - 1] -
                    dem[i + 1, j - 1]) / 8
            slope = np.sqrt(dzdx ** 2 + dzdy ** 2)
            if slope == 0:
                flow_direction[i, j] = -1
            else:
                flow_direction[i, j] = math.atan2(dzdy, dzdx) * 180 / np.pi
                if flow_direction[i, j] < 0:
                    flow_direction[i, j] += 360

    return flow_direction


def regionGrowing(im_array, aspect, flowDir, slope, aspectZonal, shppath, xyPath, demFilename,
                  dx, u, maxV, diff_angle, times):
    shp_to_xy(shppath, demFilename, xyPath)
    data = gdal.Open(xyPath, GA_ReadOnly)
    orignArray = data.ReadAsArray().astype(np.int)
    [m, n] = im_array.shape
    a = orignArray / 255

    aspect_a = a * aspect
    aspect_a[aspect_a < 0] = 0
    exist = (a > 0)
    mean_value = aspect_a.sum() / exist.sum()

    f = a
    s = np.zeros(a.shape)
    s = s + a
    # Iterate over the source area
    for c in range(1, m - 1):
        for d in range(1, n - 1):
            if a[c, d] == 1:
                for i in range(times):  # The Monte Carlo cycle

                    w = np.zeros((m, n))
                    w[c, d] = np.tan(slope[c, d] / 180 * np.pi) - u

                    if w[c, d] <= 0 or c == m - 1 or d == n - 1 or c == 0 or d == 0:
                        break
                    if slope[c, d] != 0:
                        maxW = (maxV * maxV) / (2 * 9.8 * dx)
                        if w[c, d] > maxW:
                            w[c, d] = maxW
                    ii = c
                    jj = d
                    isC_flag = -1
                    is_flat = -999
                    last_aspect = -9999
                    while True:
                        if aspectZonal[ii, jj] < 0:
                            break
                        if is_flat == -999:
                            aspect1 = flowDir[ii, jj]
                        else:
                            aspect1 = is_flat
                        # isC = 0
                        if aspect1 < 0:
                            aspect1 = mean_value
                        # Downhill movements
                        if isC_flag == -1:
                            x, y = mcmc.mcmc_single(diff_angle, aspect1)
                            # Case 1: Propagation within the same slope unit
                            if aspectZonal[ii, jj] == aspectZonal[ii + x, jj + y]:
                                if x != 0 and y != 0:
                                    isC = downhill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                else:
                                    isC = downhill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                            else:
                                # Case 2: Enter a new slope unit, and the slope difference between
                                # the new unit and the previous unit is less than or equal to 90
                                if abs(aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) < 90 or abs(
                                        aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) > 180:
                                    if x != 0 and y != 0:
                                        isC = downhill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                    else:
                                        isC = downhill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                if 90 <= abs(aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) <= 180:
                                    # Case 3: Enter the reverse slope unit
                                    isC_flag = flowDir[ii, jj] if 90 <= abs(
                                        flowDir[ii, jj] - aspectZonal[ii + x, jj + y]) <= 180 else aspectZonal[
                                        ii, jj]

                                    if x != 0 and y != 0:
                                        isC = uphill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                    else:
                                        isC = uphill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                        # Uphill movements
                        else:
                            aspect_cp = aspect1 if 90 <= abs(aspect1 - isC_flag) <= 180 else aspectZonal[
                                ii, jj]
                            x, y, isC_flag = mcmc.mcmc_countorslope(diff_angle, aspect_cp, isC_flag, slope[ii, jj], dx,
                                                                    u, w[ii, jj])
                            if x == 0 and y == 0:
                                isC = 0
                            else:
                                # Case 1: Propagation within the same slope unit
                                if aspectZonal[ii, jj] == aspectZonal[ii + x, jj + y]:
                                    if x != 0 and y != 0:
                                        isC = uphill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                    else:
                                        isC = uphill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                else:
                                    # Case 2: Enter a new slope unit, and the slope difference between
                                    # the new unit and the previous unit is less than or equal to 90
                                    if abs(aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) < 90 or abs(
                                            aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) > 180:
                                        if x != 0 and y != 0:
                                            isC = uphill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                        else:
                                            isC = uphill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                    if 90 <= abs(aspectZonal[ii, jj] - aspectZonal[ii + x, jj + y]) <= 180:
                                        # Case 3: Enter the reverse slope unit
                                        isC_flag = -1
                                        if x != 0 and y != 0:
                                            isC = downhill2(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                                        else:
                                            isC = downhill1(w[ii, jj], slope[ii + x, jj + y], u, maxV, dx)
                        if 90 <= abs(last_aspect - aspect1) <= 180:
                            break
                        last_aspect = aspect1
                        w[ii, jj] = 0
                        f[ii + x, jj + y] = 1
                        s[ii + x, jj + y] += 1
                        if isC != 0:
                            w[ii + x, jj + y] = isC
                        else:
                            break
                        # Fixed the flat slope error
                        if slope[ii + x, jj + y] == 0.0:
                            is_flat = aspect[ii, jj]
                        else:
                            is_flat = -999
                        ii = ii + x
                        jj = jj + y
                        if ii == m - 1 or jj == n - 1 or ii == 0 or jj == 0 or aspectZonal[ii, jj] < 0:
                            break
    data = im_array * f
    return data, f, s


def SUsM_dinf(pro_name, sourceShpPath, DEMPath, dx, u, maxV, diff_angle, times, outFilePath,
              gridNew, aspect, flow_direction, slope, aspectZonal_ali):
    outFilePath = outFilePath + '\\' + pro_name
    filePath = os.path.dirname(DEMPath)
    DEMName = os.path.basename(DEMPath)

    # Process data path setting
    if not os.path.exists(outFilePath):
        os.makedirs(outFilePath)

    xyPath = outFilePath + '\\' + 'xy.tif'
    resultPath = outFilePath + '\\' + 'result.tif'
    resultfPath = outFilePath + '\\' + 'resultf.tif'
    shpClassPath = outFilePath + '\\' + 'resultClass.tif'
    tgt = outFilePath + '\\' + 'resultShp.shp'

    # arcpy setting
    arcpy.env.workspace = filePath
    arcpy.env.overwriteOutput = True

    # register gdal
    gdal.AllRegister()
    data = gdal.Open(DEMPath, GA_ReadOnly)
    if data is None:
        print('Cannot open this file:' + DEMPath)
        sys.exit(1)
    Ref = osr.SpatialReference()
    Ref.ImportFromWkt(data.GetProjectionRef())
    transform = data.GetGeoTransform()

    # *************
    # Process.4 regionGrowing
    # *************

    # cell propagation
    driver = gdal.GetDriverByName('GTiff')

    result, f, s = regionGrowing(gridNew, aspect, flow_direction, slope, aspectZonal_ali, sourceShpPath, xyPath,
                                 DEMPath, dx, u, maxV, diff_angle, times)
    ds4 = driver.Create(resultPath, result.shape[1], result.shape[0], 1, GDT_Float32)
    ds4.SetProjection(Ref.ExportToWkt())
    ds4.SetGeoTransform(transform)
    band = ds4.GetRasterBand(1)
    band.WriteArray(result, 0, 0)

    outline = findOutline(f)

    s = outlinerReplaceByBoxplot(s)
    s = normalization(s)
    ds5 = driver.Create(resultfPath, s.shape[1], s.shape[0], 1, GDT_Float32)
    ds5.SetProjection(Ref.ExportToWkt())
    ds5.SetGeoTransform(transform)
    band = ds5.GetRasterBand(1)
    band.WriteArray(s, 0, 0)

    ds6 = driver.Create(shpClassPath, outline.shape[1], outline.shape[0], 1, GDT_Float32)
    ds6.SetProjection(Ref.ExportToWkt())
    ds6.SetGeoTransform(transform)
    band = ds6.GetRasterBand(1)
    band.WriteArray(outline, 0, 0)

    # Raster to Vector
    mask = band
    driver = ogr.GetDriverByName("ESRI Shapefile")
    shp = driver.CreateDataSource(tgt)
    tgtLayer = "extract"
    srs = osr.SpatialReference()
    srs.ImportFromWkt(data.GetProjectionRef())
    layer = shp.CreateLayer(tgtLayer, srs=srs)
    fd = ogr.FieldDefn("DN", ogr.OFTInteger)
    layer.CreateField(fd)
    dst_field = 0
    extract = gdal.Polygonize(band, mask, layer, dst_field, [], None)
    data = None
    projData = None
    del ds4
    del ds5
    del ds6
    print("Process.4 done, Check All results in Current Workspace.")
