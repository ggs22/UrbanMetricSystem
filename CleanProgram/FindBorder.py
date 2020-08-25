# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 11:23:31 2018

@author: GelbJ
"""

######################################################################
## Import des librairies
######################################################################
import sys
from BasicGeometry import FlyDist, PointCoords, NorthAngle
import numpy as np
from ForceField import VectorForce
sys.path.append("G:/Python/___JBasics")
#from JPlot import JPlot
import JGeom
#from JQgis import JVectorLayer as JV

######################################################################
## Import des librairies
######################################################################

def ChaikinAlgo(Points,Iter) :
    """
    Fonction permettant d'appliquer l'aglorithme de simplification du contours
    d'un polygone
    """
    for u in range(Iter) :
        NewPoints=[]
        for e in range(len(Points)-1) :
            Pi = np.array(Points[e])
            Pi2 = np.array(Points[e+1])
            Qi = 0.75*Pi + 0.25*Pi2
            Ri = 0.25*Pi + 0.75*Pi2
            NewPoints.append(Qi)
            NewPoints.append(Ri)
        Points = NewPoints
        Points.append(Points[0])
    return NewPoints


def CatmullRomSpline(P0, P1, P2, P3, nPoints=100):
    # source : https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
    """
    P0, P1, P2, and P3 should be (x,y) point pairs that define the Catmull-Rom spline.
    nPoints is the number of points to include in this curve segment.
    """
    # Convert the points to numpy so that we can do array multiplication
    P0, P1, P2, P3 = map(np.array, [P0, P1, P2, P3])

    # Calculate t0 to t4
    alpha = 0.5
    def tj(ti, Pi, Pj):
        xi, yi = Pi
        xj, yj = Pj
        return ( ( (xj-xi)**2 + (yj-yi)**2 )**0.5 )**alpha + ti

    t0 = 0
    t1 = tj(t0, P0, P1)
    t2 = tj(t1, P1, P2)
    t3 = tj(t2, P2, P3)

    # Only calculate points between P1 and P2
    t = np.linspace(t1,t2,nPoints)

    # Reshape so that we can multiply by the points P0 to P3
    # and get a point for each value of t.
    t = t.reshape(len(t),1)
    A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
    A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
    A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
    B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
    B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

    C    = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
    return C


def CatmullRomChain(P):
    # source : https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
    """
    Calculate Catmull Rom for a chain of points and return the combined curve.
    """
    P.append(P[0])
    P.append(P[1])
    P.append(P[2])
    sz = len(P)
    # The curve C will contain an array of (x,y) points.
    C = []
    for i in range(sz-3):
        c = CatmullRomSpline(P[i], P[i+1], P[i+2], P[i+3])
        C.extend(c)

    return C

def Direction(CenterPoint,Pt,Resultante) :
    """
    Permet d'indiquer si la resultante calculee indique un rapprochement ou un
    eloignement du centre
    """
    ##
    #Old methode
    ##
    D1 = FlyDist(CenterPoint,Pt)
    D2 = FlyDist(CenterPoint,Resultante[2])
    
    if D1>D2 :
        ## le prochain point s'eloigne du centre
        return "Outside"
    elif D1<D2 :
        ## le prochain point s'approche du centre
        return "Inside"
    elif D1==D2 :
        ## le prochain point est a egal distance du centre
        return "Equilibrium"
    
def OptimalPosition(CenterPoint,StartPoint,Tolerance,ForceField,Beta,MaxIter = 500,Plot1=None) :
    #on regarde ce que donne le point de reference (point probable)
    Resultante = VectorForce(StartPoint,ForceField,Beta)
    Dir = Direction(CenterPoint,StartPoint,Resultante)
    ActualPoint = StartPoint
    Angle = NorthAngle(CenterPoint,StartPoint)
    Jump = Tolerance * 10
    #algorithme permettant de trouver la position optimale
    i=0
    Points=[]
    while i<MaxIter :
        i+=1
        if Dir == "Outside" :
            ## dans ce cas il faut se rapprocher
            Pt2 = PointCoords(CenterPoint, FlyDist(CenterPoint,ActualPoint)-Jump, Angle)
        elif Dir == "Inside" :
            ## dans ce cas, il faut s'eloigner
            Pt2 = PointCoords(CenterPoint, FlyDist(CenterPoint,ActualPoint)+Jump, Angle)
        elif Dir == "Equilibrium" :
            ## on a trouve l'equilibre ! stop !
            print("Fuck, equilibrium !!")
            if i>1 :
                return Pt2
            else : 
                return ActualPoint
        Resultante2 = VectorForce(Pt2,ForceField,Beta)
        Dir2 = Direction(CenterPoint,Pt2,Resultante2)
        
        if Dir2==Dir :
            ## on va dans la meme direction, il faut recommencer
            Dir = Dir2
            ActualPoint = Pt2
            Points.append(Pt2)
        else :
            ## on a change de direction ! on se rapproche du but :
            if Jump>Tolerance :
                Jump = Jump/2.0
                Dir = Dir2
            else :
                ## on a atteint le pallier demande, le travail est fini !!
                if Dir == "Inside" :
                    ## on cherchait a se rapprocher
                    MidDist =  FlyDist(CenterPoint,Pt2)+(FlyDist(ActualPoint,Pt2)/2.0)
                else :
                    ## on cherchait a s'eloigner
                    MidDist =  FlyDist(CenterPoint,Pt2)-(FlyDist(ActualPoint,Pt2)/2.0)
                MidPt = PointCoords(CenterPoint, MidDist, Angle)
                return MidPt
    print("Sorry I haven't found a point that match the tolerance ...")
    Plot1.AddLayer([JGeom.OgrPoint(Pt) for Pt in Points],"ErrorPts"+str(len(Plot1.Layers.keys())),"Points",{"Color":(0,0,1),"Marker":"o","Size":5},1)
    Plot1.AddLayer([JGeom.OgrPoint(CenterPoint),JGeom.OgrPoint(StartPoint)],"StartPts","Points",{"Color":(1,0,0),"Marker":"o","Size":15},2)
    Plot1.Draw()
    raise ValueError("Algorithm didn't converge")
    

def FindVectorBorder(CenterPoint,StartPoint,NbEdges,Tolerance,LayerSector,Beta,PonderField,ClockWise=True,Plot1=None,ForceField=None) :
    if ForceField is None :
        ForceField =[((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
    ### Setting des premiers parametres
    ProbaPoint = StartPoint
    StartAngle = NorthAngle(CenterPoint,StartPoint)
    StepAngle = 360.0/NbEdges
    OkPoints=[]
    Missed=0
    ##lancement de l'algorithme
    for e in range(NbEdges) :
        print("avancement : "+str(round(e/float(NbEdges)*100,2))+"%")
        #algorithme permettant de trouver la position optimale
        OptimalPt = OptimalPosition(CenterPoint,ProbaPoint,Tolerance,ForceField,Beta,MaxIter = 500, Plot1=Plot1)
        OkPoints.append(OptimalPt)
        Dist = FlyDist(CenterPoint,OptimalPt)
        #ajustement des variables pour la prochaine boucle
        if ClockWise :
            ProbaPoint = PointCoords(CenterPoint,Dist,StartAngle+(e+1)*StepAngle)
        else :
            ProbaPoint = PointCoords(CenterPoint,Dist,StartAngle-(e+1)*StepAngle)
    return OkPoints

def MeanBorder(B1,B2,Center) :
    B2_5 = [B2.pop(0)]
    B2.reverse()
    B2_5.extend(B2)
    Border=[]
    for Pt1,Pt2 in zip(B1,B2_5) :
        D1 = FlyDist(Center,Pt1)
        D2 = FlyDist(Center,Pt2)
        D3 = FlyDist(Pt1,Pt2)
        Min = min([D1,D2])
        Angle1 = NorthAngle(Center,Pt1)
        Angle2 = NorthAngle(Center,Pt2)
        Border.append(PointCoords(Center,Min+(D3/2.0),Angle1))
    return Border


def FindCenter(Start,LayerSector,Beta,PonderField,Tolerance = 5, MaxIter = 100) :
    """
    Fonction pour trouver de facon iterative le centre
    
    Cas 1 : signifie que l'on doit continuer a boucler
    """
    ForceField =[((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
    Jump = Tolerance *100
    PrecPt = None
    ActualPt = Start
    i=0
    while i<MaxIter :
        i+=1
        Resultante = VectorForce(ActualPt,ForceField,Beta)
        NewPt = PointCoords(ActualPt, Jump, Resultante[0]+180)
        if not (PrecPt is None):
            D1 = FlyDist(NewPt,ActualPt)
            D2 = FlyDist(NewPt,PrecPt)
            if D1<D2 :
                Case = 1
            else : 
                Case = 2
        else :
            Case = 1
        #si l'on est dans le cas 1, on prepare la nouvelle boucle
        if Case ==1 :
            PrecPt = ActualPt
            ActualPt = NewPt
        #si l'on est dans le cas 2, on doit verifier que Jump est descendu en dessous de Tolerance
        else :
            if Jump>Tolerance : 
                Jump = Jump/2.0
            else :
                FinalPt = PointCoords(ActualPt, Jump/2.0, Resultante[0])
                return FinalPt
            
        