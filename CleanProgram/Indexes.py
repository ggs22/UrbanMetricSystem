# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:00:54 2018

@author: GelbJ
"""

import sys,os
sys.path.append("E:/Python/___JBasics")
from JQgis import JVectorLayer as JV
import JGeom, JDecorators
from JPlot import JPlot
from numba import jit
from JHtml import Balise
import itertools
import multiprocessing as mp
from copy import deepcopy

from math import pi, sqrt

import scipy
from scipy.optimize import minimize
import numpy as np

from BasicGeometry import FlyDist, PointCoords, NorthAngle

from FindBorder import FindVectorBorder,CatmullRomChain



###################
## Index of sprawl
###################

#def Power(Distance) :
#    """
#    Fonction estimant la force totale selon une distance et Beta
#    Beta = 6 : rayon de 10km
#    Beta = 11 : rayon de 20km
#    Beta = 21 : rayon de 40km
#    """
#    Dist = Distance/1000.0
#    Repulsion = 1.0 / (float(6.0) + Dist/2.0)
#    Attraction = 1.0/Dist
#    return Attraction - Repulsion


def CircleArea(Radius) :
    return pi*Radius**2

def Radius(Area) :
    return sqrt(Area/pi)

@JDecorators.memoized
def MaxSprawlArea(Beta,Dist,NbEdges=120,RealPoly=False) :
    """
    Permet de determiner l'aire maximum que peut prendre une ville
    """
    Radius = Dist*1000
    if RealPoly == False :
        Center = (0,0)
    else : 
        Centro = RealPoly.Centroid()
        Center = (Centro.GetX(),Centro.GetY())
    CenterPt = JGeom.OgrPoint(Center)
    ForceField=[]
    Angle = 0
    StartPoint = PointCoords(Center,(Radius),Angle)
    Plot1 = JPlot([CenterPt,JGeom.OgrPoint(StartPoint)],"Center","Points",{"Color":(0,0,1),"Marker":"o","Size":25},1)
#    ## Generation d'un cercle de NbEdges rayons avec pour chacun un poid de 1000 au bout
#    for e in range(NbEdges) :
#        Pt2 = PointCoords(Center,Radius,Angle)
#        ForceField.append((Pt2,1))
#        Angle+=360.0/NbEdges
## Utilisation du champ de force pour trouver la ligne de demarcation
    ForceField = [(Center,100000)]
    Border = FindVectorBorder(Center,StartPoint,NbEdges,5,None,Beta,None,ClockWise=True,Plot1=Plot1,ForceField=ForceField)
    #Plot1.Draw()
    RoundedBorder = CatmullRomChain(Border)
    BorderPts = [JGeom.OgrPoint(Pt) for Pt in RoundedBorder]
    print(np.mean([Pt.Distance(CenterPt) for Pt in BorderPts]))
    Poly = JGeom.PolyFromPoints(BorderPts)
    if RealPoly != False :
        Plot1.AddLayer([JGeom.OgrPoint(Pt[0]) for Pt in ForceField],"Center","Points",{"Color":(0,1,0),"Marker":"o","Size":15},3)
        Plot1.AddLayer([Poly],"BackPoly","Polygone",{"BackGroundColor":(0,0,0,0),"BorderColor":(0,0,0),"LineWidth":3},1)
        Plot1.AddLayer([RealPoly],"RealPoly","Polygone",{"BackGroundColor":(0,0,0,0),"BorderColor":(1,0,0),"LineWidth":3},1)
        Plot1.AddLayer([CenterPt,JGeom.OgrPoint(StartPoint)],"Center","Points",{"Color":(0,0,1),"Marker":"o","Size":25},1)
        Plot1.Draw()
    return Poly.Area()


#Poly=MaxSprawlArea(6.0,10)



#def Ei(RealBorder,Beta) :
#    return (RealBorder.Area()/CircleArea(Beta*1000))-1
#
#def Emax(Beta,Dist) :
#    return (MaxSprawlArea(Beta,Dist)/CircleArea(Beta*1000))-1
#
#def SprawlIndex(Polygon,ThisBeta,Dist) :
#    return Ei(Polygon,ThisBeta) / Emax(ThisBeta,Dist)
#    

def SprawlIndex(Polygon,ThisBeta,Dist,ThisRealPoly=None) : 
    ThisArea = Polygon.Area()
    TheoricArea = MaxSprawlArea(ThisBeta,Dist,120,ThisRealPoly)
    return (ThisArea/TheoricArea)-1

###################
## Main Parameters
###################
    
Root = "G:/Projets/Luc_Normand"
Sortie = Root+"/IndexResults/Table.html"

Parameters =[
        {
        "City" : "Montreal",
        "Year" : "2016",
         "InputShp" : Root+"/Montreal Datas/Results/2016/CleanResults/Beta_6_10KM/MergedBorder.shp",
         "Beta" : 6,
         "Dist" : 10},
                
#        {
#        "City" : "Montreal",
#        "Year" : "2016",
#         "InputShp" : Root+"/Montreal Datas/Results/2016/CleanResults/Beta_11_20KM/MergedBorders.shp",
#         "Beta" : 11,
#         "Dist" : 20},
#                
#        {
#        "City" : "Montreal",
#        "Year" : "2016",
#         "InputShp" : Root+"/Montreal Datas/Results/2016/CleanResults/Beta_21_40KM/MergedBorders.shp",
#         "Beta" : 21,
#         "Dist" : 40},
                
        ]

###################
## Executing
###################

Tableau = Balise("table")
Header = Balise("tr")
Header.innerHTML = '<td>City</td><td>Year</td><td>Beta</td><td>Area</td><td>Index</td>'
Tableau.Children.append(Header)

for Params in Parameters : 
    ## recuperation du polygone en entree
    Row = Balise("tr")
    Layer = JV.JFastLayer(Params["InputShp"])
    Layer.Initialize(ID="OID",GeomIndex = False)
    Feat = Layer.NiceFeat(Layer[0])
    Border = Feat["Geom"]
    Poly =  JGeom.PolyFromPoints(JGeom.PointsFromLine(Border))
    Index = SprawlIndex(Poly,Params["Beta"],Params["Dist"],ThisRealPoly = Poly)
    print(Index)
    Row.innerHTML = '<td>'+Params["City"]+'</td><td>'+Params["Year"]+'</td><td>'+str(Params["Beta"])+'</td><td>'+str(round(Poly.Area()/1000000.0,2))+'</td><td>'+str(round(Index,2))+'</td>'
    Tableau.Children.append(Row)
    
Tableau.Write(Sortie)

