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
sys.path.append("H:/Python/___JBasics")
#from JPlot import JPlot
# import JGeom
from JQgis import JGeom
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


######################################################################
## Definition des fonctions
######################################################################
    

def Orientation(Center,A,B) : 
    Dist1 = FlyDist(Center,A)
    Dist2 = FlyDist(Center,B)
    if Dist1 == Dist2 :
        return "Equilibrium"
    elif Dist1 < Dist2 : 
        return "Outside"
    elif Dist1 > Dist2 : 
        return "Inside"


def OptimalPosition(Center,StartPoint,ForceField,Beta,Tolerance = 10, Plot=None,MaxIter = 500) :
    """
    Fonction permettant de trouver par iteration le point tombant a la frontiere
    a partir d'un point estime
    """
    ######
    # definition des variables initiales
    ######
    MaxJump = 1000
    Jump = 100
    ActualPoint = StartPoint
    MainAngle = NorthAngle(Center,StartPoint)
    
    ######
    # definition de la situation initiale
    ######
    Angle,Power,Junk = VectorForce(StartPoint,ForceField,Beta)
    Pt2 = PointCoords(StartPoint,5,Angle)
    Direction = Orientation(Center,StartPoint,Pt2)
    
    ######
    # Debut des iterations
    ######
    Pts = []
    i=0
    while i < MaxIter :
        i+=1
        if Direction == "Outside" :
            #dans ce cas, il faut que le prochain point se rapproche
            ProjPt = PointCoords(Center,FlyDist(Center,ActualPoint)-Jump,MainAngle)
        elif Direction == "Inside" :
            #dans ce cas, il faut que le prochain point s'eloigne du centre
            ProjPt = PointCoords(Center,FlyDist(Center,ActualPoint)+Jump,MainAngle)
        elif Direction == "Equilibrium" : 
            #oua, point d'equilibre, on s'arrete
            return ActualPoint
        
        #####
        #Ici, observation de la situation actuelle
        #####
        Angle,Power,Junk = VectorForce(ProjPt,ForceField,Beta)
        Pt2 = PointCoords(ProjPt,5,Angle)
        NewDirection = Orientation(Center,ProjPt,Pt2)
        Pts.append(ProjPt)
        
        if NewDirection == Direction :
            #on continue dans le meme sens, il faut donc refaire un tour de boucle
            Direction = NewDirection
            ActualPoint = ProjPt
        elif NewDirection != Direction and Jump>Tolerance :
            #on a croise le centre d'equilibre, on va donc devoir essayer de s'en approcher
            Direction = NewDirection
            ActualPoint = ProjPt
            Jump = Jump/2.0
        elif NewDirection != Direction and Jump<Tolerance :
            #on a atteind le point d'equilibre attendu !
            D1 = FlyDist(Center,ActualPoint)
            D2 = FlyDist(Center,ProjPt)
            D3 = FlyDist(ActualPoint,ProjPt)
            if D1>D2 : 
                return PointCoords(Center,D2+D3/2.0,MainAngle)
            else :
                return PointCoords(Center,D1+D3/2.0,MainAngle)
        
    print("Error, impossible to converge")
    Plot.AddLayer([JGeom.OgrPoint(Center),JGeom.OgrPoint(StartPoint)],"IterPOINT","Points",{"Color":(0,1,0),"Marker":"o","Size":20},15)
    Plot.AddLayer([JGeom.OgrPoint(Pt) for Pt in Pts],"IterPOINT","Points",{"Color":(1,0,0),"Marker":"o","Size":10},8)
    Plot.Draw()
    raise ValueError("Can't converge after "+str(MaxIter)+" iterations")        
    
    
    

def FindVectorBorder(Center,LayerSector,StartPoint,Beta,PonderField, Tolerance = 10,NbEdge = 180,ClockWise=True,Plot=None, ForceField=None) :
    ######
    # Creation du forcefield si il n'est pas defini
    ######
    if ForceField is None :
        ForceField =[((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
        
    #####
    #Defintion des variables initiles
    #####
    StartAngle = NorthAngle(Center,StartPoint)
    StepAngle = 360.0/NbEdge
    i=0
    ProbaPoint = StartPoint
    Pts=[]
    while i<NbEdge :
        if i%10 == 0 :
            print("Avancement : "+str(round(float(i)/float(NbEdge)*100,2)))
        #trouver la meilleure position pour le point probabale
        OkPt = OptimalPosition(Center,ProbaPoint,ForceField,Beta,Tolerance,MaxIter=500,Plot = Plot)
        Pts.append(OkPt)
        i+=1
        #definition du prochain angle de vise
        if ClockWise :
            ThisAngle = StartAngle+(i*StepAngle)
        else : 
            ThisAngle = StartAngle-(i*StepAngle)
        ProbaPoint = PointCoords(Center,FlyDist(Center,OkPt),ThisAngle)
    return Pts
    
    
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


#def FindCenter(Start,LayerSector,Beta,PonderField,Tolerance = 5, MaxIter = 100, ForceField=None) :
#    """
#    Fonction pour trouver de facon iterative le centre
#    
#    Cas 1 : signifie que l'on doit continuer a boucler
#    """
#    if ForceField is None :
#        ForceField =[((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
#    Jump = Tolerance *100
#    PrecPt = None
#    ActualPt = Start
#    i=0
#    while i<MaxIter :
#        i+=1
#        Resultante = VectorForce(ActualPt,ForceField,Beta)
#        NewPt = PointCoords(ActualPt, Jump, Resultante[0])
#        if not (PrecPt is None):
#            D1 = FlyDist(NewPt,ActualPt)
#            D2 = FlyDist(NewPt,PrecPt)
#            if D1>D2 :
#                Case = 1
#            else : 
#                Case = 2
#        else :
#            Case = 1
#        #si l'on est dans le cas 1, on prepare la nouvelle boucle
#        if Case ==1 :
#            PrecPt = ActualPt
#            ActualPt = NewPt
#        #si l'on est dans le cas 2, on doit verifier que Jump est descendu en dessous de Tolerance
#        else :
#            if Jump>Tolerance : 
#                Jump = Jump/2.0
#            else :
#                FinalPt = PointCoords(ActualPt, Jump/2.0, Resultante[0])
#                return FinalPt
    

def FindCenter(Start,LayerSector,Beta,PonderField,Tolerance = 5, MaxIter = 500, ForceField=None) :
    """
    Fonction pour trouver de facon iterative le centre
    
    Cas 1 : signifie que l'on doit continuer a boucler
    """
    if ForceField is None :
        ForceField =[((Feat["Geom"].GetX(),Feat["Geom"].GetY()),Feat[PonderField]) for Feat in LayerSector.Iterate(True)]
    Jump = Tolerance *100
    ActualPt = Start
    i=0
    PrecPt = None    
    while i<MaxIter :
        Resultante = VectorForce(ActualPt,ForceField,Beta)
        ProjPt = PointCoords(ActualPt,Jump,Resultante[0])
        # si on est a la premiere iteration il suffit de passer a la prochaine
        if PrecPt is None :
            Case = 1
        # sinon, on va s'interesser aux distances entre les points
        else :
            D1 = FlyDist(PrecPt,ActualPt)
            D2 = FlyDist(PrecPt,ProjPt)
            if D2 > D1 :
                #si le nouveau point est plus eloigne du precedent, il suffit de continuer
                Case = 1
            elif D1>D2 :
                #si le nouveau point est plus proche du precedent, on divise le Jump par 2
                if Jump > Tolerance :
                    Case = 2
                    Jump = Jump/2.0
                else :
                    #on a fini
                    return ProjPt
                    
        if Case == 1 :
            #dans le cas 1 : il faut continuer a appliquer la methode
            PrecPt = ActualPt
            ActualPt = ProjPt
            
        elif Case ==2 :
            #dans le cas 2 on doit recommencer de la on on est
            PrecPt = None
            ActualPt = ProjPt
            
    raise ValueError("Fing center couldn't converge")
            