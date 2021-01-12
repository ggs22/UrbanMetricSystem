# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 11:21:39 2018

@author: GelbJ
"""

###############################################################################
# import des fonctions et librairies
###############################################################################

from ForceField import CalculateForceAtPoints
from FindBorderV2 import FindVectorBorder, MeanBorder ,CatmullRomChain, FindCenter, OptimalPosition
import sys
sys.path.append("H:/Python/___JBasics")
from JQgis import JVectorLayer as JV
import JGeom
import JPlot

###############################################################################
# Definition des parametres principaux
###############################################################################

# Background = "E:/Projets/Luc_Normand/Montreal Datas/2006/MASK.shp"
Background = "../Data/Montreal/MASK.shp"
LayerB = JV.JFastLayer(Background)
LayerB.Initialize(ID="OID",GeomIndex=False)
BackPoly = LayerB.Geoms.values()
Plot1 = JPlot(BackPoly,"BackPoly","Polygone",{"BackGroundColor":(0,0,0,0),"BorderColor":(0,0,0),"LineWidth":1},1)

Parameters =[
        {"LayerSector":"E:/Projets/Luc_Normand/Montreal Datas/2006/Dot_density_100.shp",
         "LayerGrid":"E:/Projets/Luc_Normand/Montreal Datas/Grid_200m_zoom10k.shp",
         "Beta":6.0,
         "OutPut":"E:/Projets/Luc_Normand/Montreal Datas/Results/2006/Beta_6_10KM",
         "PonderField":"Population",
         "GenerateMapVector":False,
         "GenerateBorders" : True,
         "CenterStatic" : False,
         "GenerateCircle" : True,
         "CenterPoint" : (298269.300934,5043268.31072),
         "StartPoint" : (290523.932712,5051741.96683),
         "NbEdges" : 180,
         "Tolerance" : 5,
         "Cores" : 4
         },
         
#         {"LayerSector":"J:/Projets/Luc_Normand/Montreal Datas/2006/Dot_density_100.shp",
#         "LayerGrid":"J:/Projets/Luc_Normand/Montreal Datas/Grid_1000m.shp",
#         "Beta":11.0,
#         "OutPut":"J:/Projets/Luc_Normand/Montreal Datas/Results/2006/Beta_11_20KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (296967.309429,5043846.5265),
#         "StartPoint" : (283876.083443,5029455.80382),
#         "NbEdges" : 120,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
          
#         {"LayerSector":"J:/Projets/Luc_Normand/Montreal Datas/2006/Dot_density_100.shp",
#         "LayerGrid":"J:/Projets/Luc_Normand/Montreal Datas/Grid_1000m.shp",
#         "Beta":21.0,
#         "OutPut":"J:/Projets/Luc_Normand/Montreal Datas/Results/2006/Beta_21_40KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (296695.978194,5044053.10126),
#         "StartPoint" : (291420.747641,5085021.77491),
#         "NbEdges" : 120,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
          
#          {"LayerSector":"J:/Projets/Luc_Normand/Montreal Datas/2016/PopulationPoints.shp",
#         "LayerGrid":"J:/Projets/Luc_Normand/Montreal Datas/Grid_200m_zoom10k.shp",
#         "Beta":6.0,
#         "OutPut":"J:/Projets/Luc_Normand/Montreal Datas/Results/2016/Beta_6_10KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "GenerateCircle":True,
#         "CenterPoint" : (298362.29177,5043261.16236),
#         "StartPoint" : (290479.726928,5051355.54276),
#         "NbEdges" : 120,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
#
#         {"LayerSector":"E:/Projets/Luc_Normand/Montreal Datas/2016/PopulationPoints.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Montreal Datas/Grid_1000m.shp",
#         "Beta":11.0,
#         "OutPut":"E:/Projets/Luc_Normand/Montreal Datas/Results/2016/Beta_11_20KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (297147.682467,5043844.8735),
#         "StartPoint" : (298508.947897,5063619.04501),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
#
#         {"LayerSector":"E:/Projets/Luc_Normand/Montreal Datas/2016/PopulationPoints.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Montreal Datas/Grid_1000m.shp",
#         "Beta":21.0,
#         "OutPut":"E:/Projets/Luc_Normand/Montreal Datas/Results/2016/Beta_21_40KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (296497.274525,5044146.47316),
#         "StartPoint" : (291477.352803,5085307.65728),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
#        
#          {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2006/PoinPopulationtDensity_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":6.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2006/Beta_6_10KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (477517.987996,5189245.22045),
#         "StartPoint" : (474346.87052,5199381.47096),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
#
#         {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2006/PoinPopulationtDensity_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":21.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2006/Beta_21_40KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (476948.723606,5189021.62289),
#         "StartPoint" : (466118.146029,5228341.32844),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },

           
#         {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2006/PoinPopulationtDensity_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":11.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2006/Beta_11_20KM",
#         "PonderField":"Population",
#         "GenerateMapVector":True,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (477066.447275,5189139.34656),
#         "StartPoint" : (454698.950106,5191022.92527),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },

#        {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2016/PointPop_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":6.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2016/Beta_6_10KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterPoint" : (477468.692896,5188865.50158),
#         "StartPoint" : (473431.130927,5199602.75429),
#         "CenterStatic" : False,
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
#         
#         {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2016/PointPop_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":11.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2016/Beta_11_20KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "CenterPoint" : (477585.039903,5189197.31732),
#         "StartPoint" : (471435.734881,5209226.48225),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },
          
#          {"LayerSector":"E:/Projets/Luc_Normand/Quebec Datas/2016/PointPop_100.shp",
#         "LayerGrid":"E:/Projets/Luc_Normand/Quebec Datas/2006/Grid_1000m.shp",
#         "Beta":21.0,
#         "OutPut":"E:/Projets/Luc_Normand/Quebec Datas/Results/2016/Beta_21_40KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterStatic" : False,
#         "GenerateCircle" : True,
#         "CenterPoint" : (476595.552598,5189139.34656),
#         "StartPoint" : (441631.622813,5214214.48812),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 4
#         },

#          {"LayerSector":"H:/Projets/Luc_Normand/Toronto Datas/2006/PointPop_200.shp",
#         "LayerGrid":"H:/Projets/Luc_Normand/Toronto Datas/Grid_1000Km.shp",
#         "Beta":6.0,
#         "OutPut":"H:/Projets/Luc_Normand/Toronto Datas/Results/2006/Beta_6_10KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterPoint" : (-176972.27496,4859841.31781),
#         "StartPoint" : (-171067.489501,4847293.64871),
#         "NbEdges" : 120,
#         "Tolerance" : 5,
#         "Cores" : 8
#         },
           
#         {"LayerSector":"H:/Projets/Luc_Normand/Toronto Datas/2006/PointPop_200.shp",
#         "LayerGrid":"H:/Projets/Luc_Normand/Toronto Datas/Grid_1000Km.shp",
#         "Beta":11.0,
#         "OutPut":"H:/Projets/Luc_Normand/Toronto Datas/Results/2006/Beta_11_20KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterPoint" : (-170751.161708,4855939.9417),
#         "StartPoint" : (-157887.164814,4845501.12455),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 8
#         },
          
#          {"LayerSector":"H:/Projets/Luc_Normand/Toronto Datas/2006/PointPop_200.shp",
#         "LayerGrid":"H:/Projets/Luc_Normand/Toronto Datas/Grid_1000Km.shp",
#         "Beta":21.0,
#         "OutPut":"H:/Projets/Luc_Normand/Toronto Datas/Results/2006/Beta_21_40KM",
#         "PonderField":"Population",
#         "GenerateMapVector":False,
#         "GenerateBorders" : True,
#         "CenterPoint" : (-171383.817293,4856361.71209),
#         "StartPoint" : (-215775.150839,4877450.23159),
#         "NbEdges" : 180,
#         "Tolerance" : 10,
#         "Cores" : 8
#         },
        ]

###############################################################################
# Execution
###############################################################################

if __name__ == '__main__':
    i=0
    for Params in Parameters :
        i+=1
        print("______Iterating on configuration : "+str(i))
        #chargement des layers
        print("loading layers")
        LayerGrid = JV.JFastLayer(Params["LayerGrid"])
        LayerSector = JV.JFastLayer(Params["LayerSector"])
        LayerGrid.Initialize(ID="OID",GeomIndex=False)
        LayerSector.Initialize(ID="OID",GeomIndex=False)
        
        #generation de la carte de champs de force
        if Params["GenerateMapVector"] :
            print("Computing the vector Map")
            VectorGrid = CalculateForceAtPoints(LayerGrid,LayerSector,Params["PonderField"],Params["Beta"], Params["Cores"])
            VectorGrid.Save(Params["OutPut"]+"/VectorDatas.shp")
        if Params["GenerateBorders"] :
            #Ajustement du centre
            print("Looking for real center")
            if Params["CenterStatic"] == False :
                RealCenter = FindCenter(Params["CenterPoint"],LayerSector,Params["Beta"],Params["PonderField"],Params["Tolerance"], MaxIter = 300)
            else : 
                RealCenter = Params["CenterPoint"]
            #generation du shp avec le centre
            LayerCenter = JV.JFastLayer("")
            LayerCenter.MakeItEmpty({"SpatialRef":LayerGrid.SpatialRef,"GeomType":"POINT"})
            LayerCenter.AttrTable.AddField("Type","|S10","None")
            F1 = {"OID":1,"Type":"ProposedCenter"}
            F2 = {"OID":2,"Type":"Start"}
            F3 = {"OID":3,"Type":"FoundCenter"}
            LayerCenter.AppendFeat(F1,JGeom.OgrPoint(Params["CenterPoint"]))
            LayerCenter.AppendFeat(F2,JGeom.OgrPoint(Params["StartPoint"]))
            LayerCenter.AppendFeat(F3,JGeom.OgrPoint(RealCenter))
            LayerCenter.Save(Params["OutPut"]+"/Centers.shp")
            print("Generating the border")
            Plot1.AddLayer([JGeom.OgrPoint(RealCenter),JGeom.OgrPoint(Params["StartPoint"])],"StartinPoints","Points",{"Color":(1,0,0),"Marker":"o","Size":15},5)
            Plot1.Draw()
            #Pts1 = FindVectorBorder(RealCenter,Params["StartPoint"],Params["NbEdges"],Params["Tolerance"],LayerSector,Params["Beta"],Params["PonderField"],ClockWise=True,Plot = Plot1)
            #Pts2 = FindVectorBorder(RealCenter,Params["StartPoint"],Params["NbEdges"],Params["Tolerance"],LayerSector,Params["Beta"],Params["PonderField"],ClockWise=False,Plot1 = Plot1)
            

            Pts1 = FindVectorBorder(RealCenter,LayerSector,Params["StartPoint"],Params["Beta"],Params["PonderField"], Tolerance = Params["Tolerance"], NbEdge = Params["NbEdges"],ClockWise=True,ForceField=None,Plot=Plot1)
            Pts2 = FindVectorBorder(RealCenter,LayerSector,Params["StartPoint"],Params["Beta"],Params["PonderField"], Tolerance = Params["Tolerance"], NbEdge = Params["NbEdges"],ClockWise=False,ForceField=None,Plot=Plot1)
            
            Pts3 = MeanBorder(Pts1,Pts2,RealCenter)
            #B1 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts1)])
            #B2 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts2)])
            #B3 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(Pts3)])
            
            
            B1 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in Pts1])
            B2 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in Pts2])
            B3 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in Pts3])
            
            LayerBorders = JV.JFastLayer("")
            LayerBorders.MakeItEmpty({"SpatialRef":LayerGrid.SpatialRef,"GeomType":"LINESTRING"})
            LayerBorders.AttrTable.AddField("Sens","|S25","")
            F1 = {"OID":1,"Sens":"Horaire"}
            F2 = {"OID":2,"Sens":"AntiHoraire"}
            F3 = {"OID":3,"Sens":"Moyen"}
            for Feat,Geom in zip([F1,F2,F3],[B1,B2,B3]) :
                LayerBorders.AppendFeat(Feat,Geom)
            LayerBorders.Save(Params["OutPut"]+"/Borders.shp")
            
            if Params["GenerateCircle"] :
                CircleBorder = FindVectorBorder(RealCenter,LayerSector,Params["StartPoint"],Params["Beta"],Params["PonderField"], Tolerance = Params["Tolerance"], NbEdge = Params["NbEdges"],ClockWise=True,ForceField=[(RealCenter,1000)],Plot=Plot1)
                B3 = JGeom.LineFromPoints([JGeom.OgrPoint(Pt) for Pt in CatmullRomChain(CircleBorder)])
                LayerCircle = JV.JFastLayer("")
                LayerCircle.MakeItEmpty({"SpatialRef":LayerGrid.SpatialRef,"GeomType":"LINESTRING"})
                LayerCircle.AttrTable.AddField("Sens","|S25","")
                F1 = {"OID":1,"Sens":"Horaire"}
                LayerCircle.AppendFeat(F1,B3)
                LayerCircle.Save(Params["OutPut"]+"/Circle.shp")
            
    print("Finished !!")
    print("\a")