# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 11:05:28 2018

@author: GelbJ
"""


######################################################################
## Import des librairies
######################################################################

import sys,os
sys.path.append("H:/Python/___JBasics")
from JQgis import JVectorLayer as JV
from numba import jit
import itertools
import JGeom
import multiprocessing as mp
from copy import deepcopy
import numpy as np

from BasicGeometry import FlyDist, PointCoords, NorthAngle

######################################################################
## Fonctionsd e calcul des champs de force
######################################################################

@jit
def Attraction(Dist) : 
    return 1.0/(1+Dist/1000.0)

@jit
def Repulsion(Beta,Dist) : 
    return 1.0 / (Beta + (Dist/1000.0)/2.0)

@jit    
def Power(Distance,Weight,Beta) :
    """
    Fonction estimant la force totale selon une distance et Beta
    Beta = 6 : rayon de 10km
    Beta = 11 : rayon de 20km
    Beta = 21 : rayon de 40km
    """
#    Dist = Distance/1000.0
#    Repulsion = 1.0 / (float(Beta) + Dist/2.0)
#    Attraction = 1.0/Dist
    return round((Attraction(Distance) - Repulsion(Beta,Distance)) * Weight,15)


@jit
def MSum(List) :
    """
    version amelioree de la somme
    """
    Arr = np.array(List)
    Values=[]
    i=0
    for e in List[0] :
        Values.append(np.sum(Arr[:,i]))
        i+=1
    return Values

@jit
def Sum(List) :
    """
    version amelioree de la somme
    """
    Tot = 0
    for element in List :
        Tot+=element
    return Tot
    
#def SplitSeq(iterable, size):
#    """
#    fonction pour splitter une liste en X (Size) parts
#    """
#    it = iter(iterable)
#    item = list(itertools.islice(it, size))
#    while item:
#        yield item
#        item = list(itertools.islice(it, size))
    

def SplitSeq(seq, num) :
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out
        
#def VectorForce(A,ForceField,Beta) :
#    """
#    Fonction permettant de calculer la vecteur resultante en un point considerant 
#    un ensemble de point d'attraction (ForceField)
#    """
#    Coords=[]
#    for Dest in ForceField :
#        B = Dest[0]
#        Dist = FlyDist(A,B)
#        if Dist ==0 :
#            Dist = 1
#        Angle = NorthAngle(A,B)        
#        PtPower = Power(Dist,Dest[1],Beta)
#        ProjPoint = PointCoords((0,0), PtPower, Angle)
#        Coords.append(ProjPoint)
#    X = Sum(Coord[0] for Coord in Coords)
#    Y = Sum(Coord[1] for Coord in Coords)
#    Angle = NorthAngle((0,0),(X,Y))
#    PtPower = FlyDist((0,0),(X,Y))
#    return (Angle,PtPower,(A[0]+X,A[1]+Y))


def VectorForce(A,ForceField,Beta) :
    """
    Fonction permettant de calculer la vecteur resultante en un point considerant 
    un ensemble de point d'attraction (ForceField)
    """
    def ApplyForce(Dest) : 
        B = Dest[0]
        Dist = FlyDist(A,B)
        Angle = NorthAngle(A,B)        
        PtPower = Power(Dist,Dest[1],Beta)
        ProjPoint = PointCoords((0,0), PtPower, Angle)
        return (ProjPoint)
    
    Coords=map(ApplyForce,ForceField)
    X,Y = MSum(Coords)
    Angle = NorthAngle((0,0),(X,Y))
    PtPower = FlyDist((0,0),(X,Y))
    return (Angle,PtPower,(A[0]+X,A[1]+Y))


#def __CalculateVector(Params) :
#    """
#    sous fonction appelee dans le programme en mode multiprocessing
#    """
#    Origins = Params[0]
#    Destinations = Params[1]
#    Queue = Params[2]
#    PonderField = Params[3]
#    Beta = Params[4]
#    for Start in Origins :
#        A = (Start["Geom"].GetX(),Start["Geom"].GetY())
#        Coords =[]
#        for Dest in Destinations :
#            B = (Dest["Geom"].GetX(),Dest["Geom"].GetY())
#            Dist = FlyDist(A,B)
#            Angle = NorthAngle(A,B)            
#            PtPower = Power(Dist,Dest[PonderField],Beta)
#            ProjPoint = PointCoords((0,0), PtPower, Angle)
#            Coords.append(ProjPoint)
#        X = Sum(Coord[0] for Coord in Coords)
#        Y = Sum(Coord[1] for Coord in Coords)
#        Angle = NorthAngle((0,0),(X,Y))
#        PtPower = FlyDist((0,0),(X,Y))
#        Queue.put((Start["OID"],Angle,PtPower))
    
def __CalculateVector(Params) :
    """
    sous fonction appelee dans le programme en mode multiprocessing
    """
    Origins = Params[0]
    Destinations = Params[1]
    Queue = Params[2]
    PonderField = Params[3]
    Beta = Params[4]
    ForceField = [((Pt["Geom"].GetX(),Pt["Geom"].GetY()),Pt[PonderField]) for Pt in Destinations]
    for Start in Origins :
        A = (Start["Geom"].GetX(),Start["Geom"].GetY())
        Angle,PtPower,Junk = VectorForce(A,ForceField,Beta)
        Queue.put((Start["OID"],Angle,PtPower))
        
        
def CalculateForceAtPoints(LayerGrid,LayerSector,PonderField,Beta,Cores=1) :
    """
    Fonction pour calculer sur un ensemble de point (LayerGrid) les vecteurs resultantes
    en fonction d'un autre layer comprenant les points d'attraction (LayerSector)
    """
    Origins = [Feat for Feat in LayerGrid.Iterate(True)]
    Destinations = [Feat for Feat in LayerSector.Iterate(True)]
    LayerGrid.AttrTable.AddField("Force","float32",0)
    LayerGrid.AttrTable.AddField("Angle","float32",0)
    LayerGrid.AttrTable.AddField("Attracted","int32",0)
    Startings = list(SplitSeq(Origins,Cores))    
    m = mp.Manager()
    Resultats = m.Queue()
    Pool = mp.Pool()
    print("Poll prepared")
    AllParams = list([[Starts,Destinations,Resultats,PonderField,Beta] for Starts in Startings])
    Pool.map(__CalculateVector,AllParams)
    print("Pooling is finished")
    #sortie des resultats
    while not Resultats.empty():
        Pt = Resultats.get()
        LayerGrid.AttrTable.SetValue(Pt[0],"Angle",Pt[1])
        LayerGrid.AttrTable.SetValue(Pt[0],"Force",Pt[2])
    Pool.terminate()
    Pool.join()
    return LayerGrid


if __name__=="__main__" : 
    P0 = (0,0)
    P1 = (0,10000)
    P2 = (10000,0)
    P3 = (0,-10000)
    P4 = (-10000,0)
    Beta = 6
    Power1 = Power(FlyDist(P0,P1),1,Beta)
    Power2 = Power(FlyDist(P0,P2),1,Beta)
    Power3 = Power(FlyDist(P0,P3),1,Beta)
    Power4 = Power(FlyDist(P0,P4),1,Beta)
    print((Power1,Power2,Power3,Power4))
    ##cas 1 (toutes les forces sont a 0)
    ForceField=[(P1,1),(P2,1),(P3,1),(P4,1)]
    Resultante = VectorForce(P0,ForceField,Beta)
    ##cas 2 toutes les forces doivent s'annuler
    ForceField=[(P1,1),(P2,1),(P3,1),(P4,1)]
    Resultante2 = VectorForce(P0,ForceField,11)
    ##cas 3 le vecteur doit pointer vers le bas
    ForceField=[(P2,1),(P3,1),(P4,1)]
    Resultante3=  VectorForce(P0,ForceField,11)