# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 15:38:32 2018

@author: gelbj
"""
import sys
sys.path.append("H:/Python/___JBasics")
from JQgis import JVectorLayer as JV
import JGeom
import numpy as np
import json


def GetIndicators(PopFile,MaskFile,ResultFolder,Sens) :
    CircleFile = ResultFolder+"/Circle.shp"
    BorderFile = ResultFolder+"/Borders.shp"
    #ouverture du layer de population
    PopLayer = JV.JFastLayer(PopFile)
    PopLayer.Initialize(ID="OID",GeomIndex=True)
    #ouverture du layer de mask (terre)
    Mask = JV.JFastLayer(MaskFile)
    Mask.Initialize(ID="OID",GeomIndex=True)
    Terre = Mask.Geoms[0]
    #ouverture du layer de bordure
    BorderLayer = JV.JFastLayer(BorderFile)
    BorderLayer.Initialize(ID="OID",GeomIndex=True)
    Border = BorderLayer.GetFirsFeatureLike("Sens == '"+Sens+"'")["Geom"]
    Ring = JGeom.PolyFromPoints(JGeom.PointsFromLine(Border))
    #ouverture du layer cercle
    CircleLayer = JV.JFastLayer(CircleFile)
    CircleLayer.Initialize(ID="OID",GeomIndex=True)
    Circle = CircleLayer.Geoms[0]
    PolyCircle = JGeom.PolyFromPoints(JGeom.PointsFromLine(Circle))
    #recuperation de tous les points dans la bordure
    InterPop = PopLayer.SpatialFilter(Ring)
    TotPop = np.sum(InterPop.AttrTable.GetVector("Population"))
    #recuperation de la superficie
    InterGeom = Terre.Intersection(Ring)
    if InterGeom is None :
        InterGeom = Ring.Intersection(Terre)
        if InterGeom is None :
            print("Are they intersecting ? : "+str(Terre.Intersect(Ring)))
            raise ValueError("The intersection seem to be empty")
    
    InterCircle = Terre.Intersection(PolyCircle)
    if InterCircle is None :
        InterCircle = PolyCircle.Intersection(Terre)
        if InterCircle is None :
            print("Are they intersecting ? : "+str(Terre.Intersect(PolyCircle)))
            raise ValueError("The intersection of the circle seem to be empty")
            
    Area = InterGeom.Area()/1000000.0
    CircleArea = InterCircle.Area()/1000000.0
    #calcul de la densite
    Density = float(TotPop)/Area
    return {"Pop":TotPop,"Area":Area,"Density":Density,"CircleArea":CircleArea}

##################################################################################################################
#####
# Calcul pour Montreal 2006
#####
    
Root = "G:/Projets/Luc_Normand/Montreal Datas/"
Population = Root+"2006/Dot_density_100.shp"
MaskFile = Root+"/Mask.shp"

Mtl2006_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_6","Moyen")
Mtl2006_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_11","Moyen")
Mtl2006_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_21","Moyen")



#####
# Calcul pour Montreal 2016
#####
Population = Root+"2016/PopulationPoints.shp"

Mtl2016_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_6","Moyen")
Mtl2016_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_11","Moyen")
Mtl2016_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_21","Moyen")


##################################################################################################################
#####
# Calcul pour Quebec 2006
#####
Root = "G:/Projets/Luc_Normand/Quebec Datas/"
Population = Root+"2006/PoinPopulationtDensity_100.shp"
MaskFile = Root+"/MASK_terre.shp"

QC2006_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_6","Moyen")
QC2006_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_11","Moyen")
QC2006_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_21","Moyen")


#####
# Calcul pour Quebec 2016
#####
Population = Root+"2016/PointPop_100.shp"

QC2016_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_6","Moyen")
QC2016_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_11","Moyen")
QC2016_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_21","Moyen")


##################################################################################################################
#####
# Calcul pour Toronto 2006
#####
Root = "G:/Projets/Luc_Normand/Toronto Datas/"
Population = Root+"2006/Density_dot_pop100.shp"
MaskFile = Root+"/MASK_terre.shp"

Trt2006_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_6","Moyen")
Trt2006_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_11","Moyen")
Trt2006_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2006/IterativeMethod/Beta_21","AntiHoraire")


#####
# Calcul pour Toronto 2016
#####
Population = Root+"2016/PointPop_100.shp"

Trt2016_Index_Beta6 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_6","Moyen")
Trt2016_Index_Beta11 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_11","Moyen")
Trt2016_Index_Beta21 = GetIndicators(Population,MaskFile,Root+"Results/2016/IterativeMethod/Beta_21","AntiHoraire")


Datas = {'MTL_2006_B6':Mtl2006_Index_Beta6,
         'MTL_2006_B11' : Mtl2006_Index_Beta11,
         'MTL_2006_B21' : Mtl2006_Index_Beta21,
         
         'MTL_2016_B6':Mtl2016_Index_Beta6,
         'MTL_2016_B11' : Mtl2016_Index_Beta11,
         'MTL_2016_B21' : Mtl2016_Index_Beta21,
         
         'QC_2006_B6':QC2006_Index_Beta6,
         'QC_2006_B11' : QC2006_Index_Beta11,
         'QC_2006_B21' : QC2006_Index_Beta21,
         
         'QC_2016_B6':QC2016_Index_Beta6,
         'QC_2016_B11' : QC2016_Index_Beta11,
         'QC_2016_B21' : QC2016_Index_Beta21,
         
         'Trt_2006_B6':Trt2006_Index_Beta6,
         'Trt_2006_B11' : Trt2006_Index_Beta11,
         'Trt_2006_B21' : Trt2006_Index_Beta21,
         
         'Trt_2016_B6':Trt2016_Index_Beta6,
         'Trt_2016_B11' : Trt2016_Index_Beta11,
         'Trt_2016_B21' : Trt2016_Index_Beta21,
         }

for key,Dico in Datas.items() :
    for key2,Value in Dico.items() : 
        Datas[key][key2] = float(Value)

with open("G:/Projets/Luc_Normand/Index.txt", 'w') as outfile:
    json.dump(Datas, outfile)