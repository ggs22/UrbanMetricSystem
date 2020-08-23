# -*- coding: utf-8 -*-
"""
Created on Wed May 02 18:53:35 2018

@author: gelbj
"""

###############################################################################
## Import des fonctions et packages
###############################################################################

from ForceField import CalculateForceAtPoints
from FindBorderV2 import FindVectorBorder, MeanBorder ,CatmullRomChain, FindCenter, OptimalPosition
import sys,os
sys.path.append("G:/Python/___JBasics")
from JQgis import JVectorLayer as JV
import JGeom
from JPlot import JPlot
import numpy as np
import json



def CreateRing(OuterPath,SensOuter,InnerPath,SensInner,MaskPath) : 
    ## ouverture des layers
    print("_opening some layers")
    LayerOuter = JV.JFastLayer(OuterPath)
    LayerOuter.Initialize(ID="OID",GeomIndex=False)
    LayerInner = JV.JFastLayer(InnerPath)
    LayerInner.Initialize(ID="OID",GeomIndex=False)
    MaskLayer = JV.JFastLayer(MaskPath)
    MaskLayer.Initialize(ID="OID",GeomIndex=False)
    ##recuperation des geometries
    print("_extracting geometries")
    OuterRing = LayerOuter.GetFirsFeatureLike("Sens == '"+SensOuter+"'")["Geom"]
    OuterRingGeom = JGeom.PolyFromPoints(JGeom.PointsFromLine(OuterRing))
    InnerRing = LayerInner.GetFirsFeatureLike("Sens == '"+SensInner+"'")["Geom"]
    InnerRingGeom = JGeom.PolyFromPoints(JGeom.PointsFromLine(InnerRing))
    MaskGeom = MaskLayer.Geoms[0]
    ## creation du contours
    print("_calculating difference")
    Ring = OuterRingGeom.Difference(InnerRingGeom)
    print("_calculating Intersection")
    OkRing = Ring.Intersection(MaskGeom)
    print(OkRing.Area())
    return OkRing
    
    

###############################################################################
## Parametres principaux
###############################################################################


################
### Config Montreal 2006
################
#
#Root = "G:/Projets/Luc_Normand/Montreal Datas/Results/2006/IterativeMethod"
#Population = "G:/Projets/Luc_Normand/Montreal Datas/2006/Dot_density_100.shp"
#Mask = "G:/Projets/Luc_Normand/Montreal Datas/MASK.shp"
#OutPut = Root+"/DensityRings"


###############
## Config Montreal 2016
###############

Root = "J:/Projets/Luc_Normand/Montreal Datas/Results/2016/IterativeMethod_B"
Population = "J:/Projets/Luc_Normand/Montreal Datas/2016/PopulationPoints.shp"
Mask = "J:/Projets/Luc_Normand/Montreal Datas/MASK.shp"
OutPut = Root+"/DensityRings"

Parameters = [
        {"OuterRing" : "24", "InnerRing":"23", "SensOuter" : "Moyen","SensInner" : "Moyen" },
        {"OuterRing" : "23", "InnerRing":"22", "SensOuter" : "Moyen","SensInner" : "Moyen" },
        {"OuterRing" : "22", "InnerRing":"21", "SensOuter" : "Moyen","SensInner" : "Moyen" },
        {"OuterRing" : "21", "InnerRing":"20", "SensOuter" : "Moyen","SensInner" : "Moyen" },
        ]


###############################################################################
## Execution de l'algorithme
###############################################################################

LayerPop = JV.JFastLayer(Population)
LayerPop.Initialize(ID = "OID",GeomIndex = True)

Results={}

i=0

Rings={}

for Params in Parameters :
    print("____Iterating on the params "+str(i))
    i+=1
    #calcul de l'anneau
    OuterPath = Root+"/Beta_"+Params["OuterRing"]+"/Borders.shp"
    InnerPath = Root+"/Beta_"+Params["InnerRing"]+"/Borders.shp"
    Ring = CreateRing(OuterPath,Params["SensOuter"],InnerPath,Params["SensInner"],Mask)
    #recuperation de la population
    Filtered = LayerPop.SpatialFilter(Ring)
    Pop = np.sum(Filtered.AttrTable.GetVector("Population"))
    #recuperation de l'aire et de la density
    Area = Ring.Area()/1000000.0
    Density = Pop/Area
    #enregistrement des valeurs
    Code = Params["OuterRing"]+"_"+Params["InnerRing"]
    Results[Code]  = {"Population" : int(Pop), "Area":Area, "Density":Density}
    Rings[Code] = Ring


with open(OutPut+".js", 'w') as outfile:
    json.dump(Results, outfile)
    
#enregistrement des geometries
MaskLayer = JV.JFastLayer(Mask)
MaskLayer.Initialize(ID="OID",GeomIndex=False)

SavingLayer = JV.JFastLayer("")
SavingLayer.MakeItEmpty({"SpatialRef":MaskLayer.SpatialRef,"GeomType":"POLYGON"})
SavingLayer.AttrTable.AddField("Code","|S25","")

i=0
for Code,Geom in Rings.items() :
    Feat = {"Code":Code,"OID":i}
    i+=1
    SavingLayer.AppendFeat(Feat,Geom)
    
SavingLayer.Save(OutPut+".shp")