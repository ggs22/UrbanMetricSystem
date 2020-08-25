# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:13:01 2018

@author: GelbJ
"""
import sys
sys.path.append("H:/Python/___JBasics")
from JQgis import JVectorLayer as JV
import JGeom
import numpy as np


#Source = "G:/Projets/Luc_Normand/Toronto Datas/2006/TES_2006_32188.shp"
Source = "G:/Projets/Luc_Normand/Toronto Datas/2016/TES_2016_32188.shp"
Distance = 200
CountField = "TOTEMPCRT"
#Sortie = "G:/Projets/Luc_Normand/Toronto Datas/2006/TES_2006_joined100.shp"
Sortie = "G:/Projets/Luc_Normand/Toronto Datas/2016/TES_2016_joined200.shp"

Layer = JV.JFastLayer(Source)
Layer.Initialize(ID="OID",GeomIndex = True)


OutPutLayer = JV.JFastLayer("")
OutPutLayer.MakeItEmpty({"SpatialRef":Layer.SpatialRef,"GeomType":"POINT"})
OutPutLayer.AttrTable.AddField("Population","int32",0)

TotPts = Layer.FeatureCount

########################################
# Regroupement de tous les points a moins de 100 les uns des autres
########################################
OID = 0
while Layer.FeatureCount>1 : 
    #recuperation de la 1ere observation
    Feat = Layer.GetRow(0)
    #recuperation de toutes les observations a X m
    Neighbours = Layer.SpatialFilter(Feat["Geom"].Buffer(Distance))
    #recuperation de la somme
    Tot = np.sum(Neighbours.AttrTable.GetVector(CountField))
    #recuperation du point central
    X = np.mean([F["Geom"].GetX() for F in Neighbours.Iterate(True)])
    Y = np.mean([F["Geom"].GetY() for F in Neighbours.Iterate(True)])
    Pt = JGeom.OgrPoint((X,Y))
    #enregistrement de cette nouvelle Feature
    Feature = {"OID":OID,"Population":Tot}
    OutPutLayer.AppendFeat(Feature,Pt)
    OID+=1
    #suppression de ces enregistrements du layer
    Layer = Layer.SpatialFilter(Feat["Geom"].Buffer(Distance),Reverse=True)
    print("Remaining Points : "+str(round(float(Layer.FeatureCount)/float(TotPts)*100,2)))
    
OutPutLayer.Save(Sortie)