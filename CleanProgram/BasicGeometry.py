# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 11:06:34 2018

@author: GelbJ
"""


######################################################################
## Import des librairies
######################################################################

from math import degrees, atan2,sqrt, pi, radians, sin, cos
from numba import jit


######################################################################
## Fonctions utilitaires
######################################################################

@jit
def FlyDist(A,B) :
    return sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)

@jit
def PointCoords(Point, d, theta):
    x0 = Point[0]
    y0 = Point[1]
    theta_rad = pi/2.0 - radians(theta)
    return x0 + d*cos(theta_rad), y0 + d*sin(theta_rad)

@jit
def NorthAngle(A, B):
    angle = degrees(atan2(B[1] - A[1], B[0] - A[0]))
    bearing2 = (90 - angle) % 360
    return bearing2