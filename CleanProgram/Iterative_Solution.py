# -*- coding: utf-8 -*-
"""
Created on Sat Apr 07 22:31:06 2018

@author: gelbj
"""

###############################################################################
## Import des fonctions et packages
###############################################################################

from ForceField import CalculateForceAtPoints
from FindBorderV2 import FindVectorBorder, MeanBorder ,CatmullRomChain, FindCenter, OptimalPosition
import sys,os
# sys.path.append("H:/Python/___JBasics")
sys.path.append("C:/repos/UrbanMetricSystem")
from JQgis import JVectorLayer as JV
from JQgis import JGeom
import JPlot

###############################################################################
##definitio des parametres generaux
###############################################################################
Distances = [44,42,30,25,22,20,15,12,10]


#TODO: un-hardcode these parameters
#####
# Montreal Config 2016
#####
StartCenter = (297684.293078,5043750.56996)
StartBorder = (331334.803927,5063512.03874)

Root = "E:/Projets/Luc_Normand/Montreal Datas"
SmallGrid = Root+'/Grid_200m_zoom10k.shp'
LargeGrid = Root+"/Grid_1000m.shp"

Population = Root+"/2016/Dot_density_100Tot.shp"
PonderField = "Population"

OutPut = Root+"/Results/2016/IterativeMethod_C"


######
## Montreal Config 2006
######
#StartCenter = (297684.293078,5043750.56996)
#StartBorder = (331334.803927,5063512.03874)
#
#Root = "E:/Projets/Luc_Normand/Montreal Datas"
#SmallGrid = Root+'/Grid_200m_zoom10k.shp'
#LargeGrid = Root+"/Grid_1000m.shp"
#
#Population = Root+"/2006/Dot_density_100Tot.shp"
#PonderField = "Population"
#
#OutPut = Root+"/Results/2006/IterativeMethod_C"


######
## Toronto Config
######
#StartCenter = (-169878.298979,4858132.98613)
#StartBorder = (-171357.456644,4874392.39626)
#
#Root = "G:/Projets/Luc_Normand/Toronto Datas"
#SmallGrid = Root+"/Grid_400m.shp"
#LargeGrid = Root+"/Grid_1000m.shp"
#
#Population = Root+"/2006/Dot_density_100Tot.shp"
#PonderField = "Population"
#
#OutPut = Root+"/Results/2006/IterativeMethod"

######
## Toronto Config 2016
######
#StartCenter = (-169878.298979,4858132.98613)
#StartBorder = (-171357.456644,4874392.39626)
#
#Root = "G:/Projets/Luc_Normand/Toronto Datas"
#SmallGrid = Root+"/Grid_400m.shp"
#LargeGrid = Root+"/Grid_1000m.shp"
#
#Population = Root+"/2016/Dot_density_100Tot.shp"
#PonderField = "Population"
#
#OutPut = Root+"/Results/2016/IterativeMethod"

######
## Quebec Config 2006
######
#StartCenter = (476948.723606,5189021.62289)
#StartBorder = (466118.146029,5228341.32844)
#
#Root = "G:/Projets/Luc_Normand/Quebec Datas"
#SmallGrid = Root+"/Grid_1000m.shp"
#LargeGrid = Root+"/Grid_1000m.shp"
#
#Population = Root+"/2006/PoinPopulationtDensity_100.shp"
#PonderField = "Population"
#
#OutPut = Root+"/Results/2006/IterativeMethod"

######
## Quebec Config 2016
######
#StartCenter = (476948.723606,5189021.62289)
#StartBorder = (466118.146029,5228341.32844)
#
#Root = "G:/Projets/Luc_Normand/Quebec Datas"
#SmallGrid = Root+"/Grid_1000m.shp"
#LargeGrid = Root+"/Grid_1000m.shp"
#
#Population = Root+"/2016/PointPop_100.shp"
#PonderField = "Population"
#
#OutPut = Root+"/Results/2016/IterativeMethod"


Parameters =[
        {"Distance" : 48, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 46, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 44, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 41, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 42, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 39.86, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        #{"Distance" : 39.30, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 38, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 41, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 40.5, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 40, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":True,"CalculateCircile":True},   # Beta = 21
#        {"Distance" : 39, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 38, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 37, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 36, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 36, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 34, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 33, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 32, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 31, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 30, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 29, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 28, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 27, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 26, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 25, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 24, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 23, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 22, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 21, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 20, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":True,"CalculateCircile":True},   # Beta = 11
#        {"Distance" : 19, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 18, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 17, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 16, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 15, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 14, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 13, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
        {"Distance" : 12, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False},
#        {"Distance" : 11, "SaveIt":True,"Grid":LargeGrid,"CalculateVectors":False,"CalculateCircile":False}, 
        {"Distance" : 10, "SaveIt":True,"Grid":SmallGrid,"CalculateVectors":True,"CalculateCircile":True},  # Beta = 16
        ]

###############################################################################
##Algorithme principal
###############################################################################
i=0
ActualCenter = StartCenter

## Generation du ForceField
LayerSector = JV.JFastLayer(Population)
LayerSector.Initialize(ID="OID",GeomIndex = True)
ForceField = [((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
PointSelection = LayerSector

## mise en place de l'arriere plan pour le plot
Background = Root+"/MASK.shp"
LayerB = JV.JFastLayer(Background)
LayerB.Initialize(ID="OID",GeomIndex=False)
BackPoly = LayerB.Geoms.values()
Plot1 = JPlot(BackPoly,"BackPoly","Polygone",{"BackGroundColor":(0,0,0,0),"BorderColor":(0,0,0),"LineWidth":1},1)

if __name__ == "__main__" :
    for Params in Parameters :
        print("Iterating on configuration : "+str(i))
        print(Params)
        i+=1
        Beta = Params["Distance"]/2.0+1
        
        ## Etape1 : trouver le centre
        print("___finding the real center")
        RealCenter = FindCenter(ActualCenter,"",Beta,"",Tolerance = 5, MaxIter = 100, ForceField=ForceField)
        
        ## Etape 2 : trouver les bordures
        print("___finding the borders")
        Pts1 = FindVectorBorder(RealCenter,"",StartBorder,Beta,"", Tolerance = 5, NbEdge =180,ClockWise=True,ForceField=ForceField,Plot=Plot1)
        Pts2 = FindVectorBorder(RealCenter,"",StartBorder,Beta,"", Tolerance = 5, NbEdge = 180,ClockWise=False,ForceField=ForceField,Plot=Plot1)
        Pts3 = MeanBorder(Pts1,Pts2,RealCenter)
        B1 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts1)])
        B2 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts2)])
        B3 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts3)])
        
        ## enregistrement des elements si demande
        print("___Savings...")
        if Params["SaveIt"] :
            Folder = OutPut+"/Beta_"+str(Beta)
            if os.path.isdir(Folder)==False :
                os.makedirs(Folder)
            # enregistrement des lignes
            LayerBorders = JV.JFastLayer("")
            LayerBorders.MakeItEmpty({"SpatialRef":LayerB.SpatialRef,"GeomType":"LINESTRING"})
            LayerBorders.AttrTable.AddField("Sens","|S25","")
            F1 = {"OID":1,"Sens":"Horaire"}
            F2 = {"OID":2,"Sens":"AntiHoraire"}
            F3 = {"OID":3,"Sens":"Moyen"}
            for Feat,Geom in zip([F1,F2,F3],[B1,B2,B3]) :
                LayerBorders.AppendFeat(Feat,Geom)
            LayerBorders.Save(Folder+"/Borders.shp")
            #enregistrement des centres
            LayerCenter = JV.JFastLayer("")
            LayerCenter.MakeItEmpty({"SpatialRef":LayerB.SpatialRef,"GeomType":"POINT"})
            LayerCenter.AttrTable.AddField("Type","|S10","None")
            F1 = {"OID":1,"Type":"ProposedCenter"}
            F3 = {"OID":2,"Type":"FoundCenter"}
            LayerCenter.AppendFeat(F1,JGeom.OgrPoint(ActualCenter))
            LayerCenter.AppendFeat(F3,JGeom.OgrPoint(RealCenter))
            LayerCenter.Save(Folder+"/Centers.shp")
            if Params["CalculateCircile"] :
                #calcul du cercle
                CircleBorder = FindVectorBorder(RealCenter,PointSelection,StartBorder,Beta,PonderField, Tolerance = 10, NbEdge = 180,ClockWise=True,ForceField=[(RealCenter,1000)],Plot=Plot1)
                B5 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(CircleBorder)])
                LayerCircle = JV.JFastLayer("")
                LayerCircle.MakeItEmpty({"SpatialRef":LayerB.SpatialRef,"GeomType":"LINESTRING"})
                LayerCircle.AttrTable.AddField("Sens","|S25","")
                F1 = {"OID":1,"Sens":"Horaire"}
                LayerCircle.AppendFeat(F1,B5)
                LayerCircle.Save(Folder+"/Circle.shp")
            if Params["CalculateVectors"] :
                #calcul et enregistrement des vecteurs
                print("______Saving grids")
                LayerGrid = JV.JFastLayer(Params["Grid"])
                LayerGrid.Initialize(ID="OID",GeomIndex=False)
                print("Generaing the vector map (may require some time)")
                VectorGrid = CalculateForceAtPoints(LayerGrid,PointSelection,PonderField=PonderField,Beta = Beta,Cores = 6)
                VectorGrid.Save(Folder+"/VectorDatas.shp")
            
        
        ## Etape 3 : preparer le forcefield et les autres elements necessaires a la prochain iteration
        Poly = JGeom.PolyFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts3)])
        PointSelection = LayerSector.SpatialFilter(Poly)
        ActualCenter = RealCenter
        ForceField = [((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in PointSelection.Iterate(True)]
        StartBorder = Pts1[-4]
        
        if Params["Distance"] == 39.30 :
            break
    
    
    