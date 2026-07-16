# coding=utf-8

from random import *
from math import *
import numpy as np


def getSpreadPosition(aspect):
    i = 0
    j = 0
    if aspect <= 45:
        if abs(aspect - 0) < abs(aspect - 45):
            i = -1
            r[2] = r[2] + 1
        else:
            i = -1
            j = 1
            r[5] = r[5] + 1
    if (aspect > 45) and (aspect <= 90):
        if abs(aspect - 45) < abs(aspect - 90):
            i = -1
            j = 1
            r[5] = r[5] + 1
        else:
            j = 1
            r[8] = r[8] + 1
    if (aspect > 90) and (aspect <= 135):
        if abs(aspect - 90) < abs(aspect - 135):
            j = 1
            r[8] = r[8] + 1
        else:
            i = 1
            j = 1
            r[11] = r[11] + 1
    if (aspect > 135) and (aspect <= 180):
        if abs(aspect - 135) < abs(aspect - 180):
            i = 1
            j = 1
            r[11] = r[11] + 1
        else:
            i = 1
            r[14] = r[14] + 1
    if (aspect > 180) and (aspect <= 225):
        if abs(aspect - 180) < abs(aspect - 225):
            i = 1
            r[14] = r[14] + 1
        else:
            i = 1
            j = -1
            r[17] = r[17] + 1
    if (aspect > 225) and (aspect <= 270):
        if abs(aspect - 225) < abs(aspect - 270):
            i = 1
            j = -1
            r[17] = r[17] + 1
        else:
            j = -1
            r[20] = r[20] + 1
    if (aspect > 270) and (aspect <= 315):
        if abs(aspect - 270) < abs(aspect - 315):
            j = -1
            r[20] = r[20] + 1
        else:
            i = -1
            j = -1
            r[23] = r[23] + 1
    if aspect > 315:
        if abs(aspect - 315) < abs(aspect - 360):
            i = -1
            j = -1
            r[23] = r[23] + 1
        else:
            i = -1
            r[2] = r[2] + 1
    return i, j


def uniform_two(a1, a2, b1, b2):
    delta_a = a2 - a1
    delta_b = b2 - b1
    if uniform(0, 1) < delta_a / (delta_a + delta_b):
        return uniform(a1, a2)
    else:
        return uniform(b1, b2)


r = [-1, 0, 0, -1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, -1, 0, 0, -1, 0, -1, -1, 0]  # 用于存储mc结果


def change_name():
    global r
    r = [-1, 0, 0, -1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, -1, 0, 0, -1, 0, -1, -1, 0]


def mcmc_angle(times, angle, aspect):
    result = []
    if aspect <= (360 - angle) and aspect >= angle:
        for i in range(times):
            aspectRand = uniform(aspect - angle, aspect + angle)
            x, y = getSpreadPosition(aspectRand)
    elif aspect > (360 - angle):
        for i in range(times):
            aspectRand = uniform_two(aspect - angle, 360, 0, aspect + angle - 360)
            x, y = getSpreadPosition(aspectRand)
    elif aspect < angle:
        for i in range(times):
            aspectRand = uniform_two(aspect - angle + 360, 360, 0, aspect + angle)
            x, y = getSpreadPosition(aspectRand)
    for i in range(1, 9):
        if r[3 * i - 1] != 0:
            a = [r[3 * i - 3], r[3 * i - 2], r[3 * i - 1] / times]
            result.append(a)
    change_name()
    return result


def mcmc_single(angle, aspect):
    x = 0
    y = 0
    if aspect <= (360 - angle) and aspect >= angle:
        aspectRand = uniform(aspect - angle, aspect + angle)
        x, y = getSpreadPosition(aspectRand)
    elif aspect > (360 - angle):
        aspectRand = uniform_two(aspect - angle, 360, 0, aspect + angle - 360)
        x, y = getSpreadPosition(aspectRand)
    elif aspect < angle:
        aspectRand = uniform_two(aspect - angle + 360, 360, 0, aspect + angle)
        x, y = getSpreadPosition(aspectRand)
    return x, y


def mcmc_countorslope(angle, aspect, last_aspect, slope, dx, u, w):
    aspect_d = 0

    c_aspect = aspect + 180 if aspect + 180 <= 360 else aspect + 180 - 360
    beta = abs(last_aspect - c_aspect)
    beta = min(beta, 360 - beta)
    beta_rd = beta * np.pi / 180
    slope = slope * np.pi / 180
    wx = 2 * 9.8 * dx * (w * np.sin(beta_rd) * np.sin(beta_rd) - u * np.cos(slope) * np.sin(beta_rd))
    wy = 2 * 9.8 * dx * (w * np.cos(beta_rd) * np.cos(beta_rd) - np.sin(slope) - u * np.cos(slope) * np.cos(beta_rd))

    if wx > 0 and wy > 0:
        vdx = sqrt(wx)
        vdy = sqrt(wy)
        new_beta = atan(vdx / vdy)
        new_beta = new_beta / np.pi * 180
        diff_beta = abs(new_beta - beta)

        if aspect > last_aspect:
            if abs(last_aspect - aspect) > 180:
                aspect_d = last_aspect - diff_beta + 360 if last_aspect - diff_beta < 0 else last_aspect - diff_beta
            elif abs(last_aspect - aspect) < 180:
                aspect_d = last_aspect + diff_beta - 360 if last_aspect + diff_beta > 360 else last_aspect + diff_beta
            else:
                aspect_d = last_aspect
        if aspect < last_aspect:
            if abs(last_aspect - aspect) > 180:
                aspect_d = last_aspect + diff_beta - 360 if last_aspect + diff_beta > 360 else last_aspect + diff_beta
            elif abs(last_aspect - aspect) < 180:
                aspect_d = last_aspect - diff_beta + 360 if last_aspect - diff_beta < 0 else last_aspect - diff_beta
            else:
                aspect_d = last_aspect
        x, y = mcmc_single(angle, aspect_d)
    else:
        x = 0
        y = 0
    return x, y, aspect_d
