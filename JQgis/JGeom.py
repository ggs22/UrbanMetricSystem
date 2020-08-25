#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      jeremy
#
# Created:     27/05/2016
# Copyright:   (c) jeremy 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from osgeo import ogr
import shapely

from osgeo import ogr,osr
import shapely
import itertools
import math
from shapely import wkt,ops,geometry,affinity
from collections import defaultdict
import matplotlib.pyplot as plt
import copy
import numpy as np
from scipy.stats import norm, chi2
from scipy import stats
from sklearn import linear_model

from math import degrees, atan2,sqrt, pi, radians, sin, cos

import matplotlib 
from matplotlib.patches import Ellipse,Polygon

###################################################
###################################################
#_______________________________fonctions utilitaires

def comp(v1, v2):
    if v1[1]<v2[1]:
        return -1
    elif v1[1]>v2[1]:
        return 1
    else:
        return 0
        
def DrawLine(line,color="r",width=2) :
    if line.GetGeometryName()=="MULTILINESTRING" :
        Subss =[LineSequence(line.GetGeometryRef(g)) for g in range(line.GetGeometryCount())]
        for Subs in Subss :
            for element in Subs :
                Extremites = GetExtremites(element)
                plt.plot([Extremites[0].GetX(),Extremites[1].GetX()],[Extremites[0].GetY(),Extremites[1].GetY()],color=color,mew=width)
    Subs = LineSequence(line)
    for element in Subs :
        Extremites = GetExtremites(element)
        plt.plot([Extremites[0].GetX(),Extremites[1].GetX()],[Extremites[0].GetY(),Extremites[1].GetY()],color=color,mew=width)

def DrawPoints(points,color="ro",width=2) :
    X = [pt.GetX() for pt in points]
    Y = [pt.GetY() for pt in points]
    plt.plot(X,Y,color,mew=width)

def DrawPolygone(Geom,BackGroundColor=(1,0,0),BorderColor=(0.5,0,0),LineWidth=2) :
    if "Z" in Geom.GetGeometryName() :
        Geom.FlattenTo2D()
    if Geom.GetGeometryName()=="MULTIPOLYGON" :
        for e in range(Geom.GetGeometryCount()) :
            G = Geom.GetGeometryRef(e)
            Bounds = G.GetBoundary()
            ring_mixed = [Bounds.GetPoint_2D(e) for i in range(Bounds.GetPointCount())]
            ring_patch = Polygon(ring_mixed,facecolor=BackGroundColor,edgecolor=BorderColor,linewidth=LineWidth)
            plt.gca().add_patch(ring_patch)
    else :
        Bounds = Geom.GetBoundary()
        ring_mixed = [Bounds.GetPoint_2D(e) for e in range(Bounds.GetPointCount())]     
        ring_patch = Polygon(ring_mixed,facecolor=BackGroundColor,edgecolor=BorderColor,linewidth=LineWidth)
        plt.gca().add_patch(ring_patch)
    

    
def DrawFigure(Geom,Color) :
    if "LINESTRING" in Geom.GetGeometryName() :
        DrawLine(Geom,Color)
    elif "POLYGON" in Geom.GetGeometryName() :
        DrawPolygone(Geom,Color)
    elif "POINT" in Geom.GetGeometryName() :
        DrawPoints([Geom],Color)
        
def SetView(Box) :
    plt.axis(Box)
    
def Clear() :
    plt.clf()

def RoundCoords(Geom,Decimals=6) :
    if Geom.GetGeometryName()=="MULTILINESTRING" :
        Lines = [Geom.GetGeometryRef(e) for e in range(Geom.GetGeometryCount())]
        RoundedLines =[RoundCoords(Line) for Line in Lines]
        NewMultiLine = ogr.Geometry(ogr.wkbMultiLineString)
        for Line in RoundedLines :
            NewMultiLine.AddGeometry(Line)
        return NewMultiLine
    elif Geom.GetGeometryName()=="LINESTRING" :
        Pts = [Geom.GetPoint(e) for e in range(Geom.GetPointCount())]
        RoundedPts = [(round(Coords[0],Decimals),round(Coords[1],Decimals)) for Coords in Pts]
        NewLine = ogr.Geometry(ogr.wkbLineString)
        for Coords in RoundedPts :
            NewLine.AddPoint(Coords[0],Coords[1])
        return NewLine
        
###################################################
###################################################
#_______________________________fonctions de geoemtrie classique
def Pytagore(A,B) :
    X2 = (abs(A[0])-abs(B[0]))**2
    Y2 = (abs(A[1])-abs(B[1]))**2
    return math.sqrt(X2+Y2)

def FlyDist(A,B) :
    return sqrt((A[0]-B[0])**2 + (A[1]-B[1])**2)

def PointCoords(Point, d, theta):
    x0 = Point[0]
    y0 = Point[1]
    theta_rad = pi/2.0 - radians(theta)
    return x0 + d*cos(theta_rad), y0 + d*sin(theta_rad)

def NorthAngle_T(A, B):
    angle = degrees(atan2(B[1] - A[1], B[0] - A[0]))
    bearing2 = (90 - angle) % 360
    return bearing2
    


###################################################
###################################################
#_______________________________fonctions de convertion
def ToShapely(Geom) :
    WKT= Geom.ExportToWkt()
    try :
        ShGeom = shapely.wkt.loads(WKT)
        return ShGeom
    except :
        print(WKT)

def ToOgr(Geom) :
    WKT = shapely.wkt.dumps(Geom)
    OgrGeom = ogr.CreateGeometryFromWkt(WKT)
    return OgrGeom
      
###################################################
###################################################
#_______________________________fonctions informative
def GetExtent(Geoms) :
    MinX = None
    MaxX= None
    MinY = None
    MaxY = None
    for G in Geoms :
        Env=G.GetEnvelope()
        if MinX==None or MinX>Env[0] :
            MinX = Env[0]
        if MaxX==None or MaxX<Env[1] :
            MaxX = Env[1]
        if MinY==None or MinY>Env[2] :
            MinY=Env[2]
        if MaxY==None or MaxY<Env[3] :
            MaxY=Env[3]
    return (MinX,MaxX,MinY,MaxY)


###################################################
###################################################
#_______________________________fonctions de projection
def ReProj(Geom,Source,Target):
    try :
        S1 = osr.SpatialReference()
        S1.ImportFromWkt(Source.ExportToWkt())
    except AttributeError :
        S1 = osr.SpatialReference()
        Source.ImportFromEPSG(int(Source))
    try :
        S2 = osr.SpatialReference()
        S2.ImportFromWkt(Target.ExportToWkt())
    except AttributeError :
        S2 = osr.SpatialReference()
        S2.ImportFromEPSG(int(Target))
    transform = osr.CoordinateTransformation(S1, S2)
    Geom.Transform(transform)
    return Geom

def ReScale(Geom,Factor,Center=None) :
    Shape = ToShapely(Geom)
    if Center==None :
        NewShape = shapely.affinity.scale(Shape, xfact=Factor, yfact=Factor, zfact=Factor, origin='center')
    else :
        NewShape = shapely.affinity.scale(Shape, xfact=Factor, yfact=Factor, zfact=Factor, origin=(Center.GetX(),Center.GetY()))
    return ToOgr(NewShape)

###################################################
###################################################
#_______________________________fonctions de Validation

def Intersect(G1,G2) :
    return G1.Intersect(G2)
    
def CentroidInside(G1,G2) :
    C1 = G1.Centroid()
    return C1.Intersect(G2)
    
def ShareSegment(G1,G2) :
    Inter = G1.Intersection(G2)
    if Inter == None :
        return False
    else :
        Name = Inter.GetGeometryName()
        if "LINE" in Name or "POLYGON" in Name :
            return True
        else :
            return False
        
        
def GetSimpleGeom(Geom) : 
    if "MULTI" not in Geom.GetGeometryName() : 
        if Geom.GetGeometryCount() == 0 :
            return [Geom]
        else : 
            return [Geom.GetGeometryRef(0)]
    else : 
        Herited=[]
        for i in range(Geom.GetGeometryCount()) :
            G = Geom.GetGeometryRef(i)
            if G.GetGeometryCount() == 0 :
                Herited.append(G)
            else : 
                Herited.append(G.GetGeometryRef(0))
        return Herited



###################################################
###################################################
#_______________________________fonctions de generation

def EllipseGeom(Points) :
    Xs=[]
    Ys=[]
    for element in Points :
        Xs.append(element.GetX())
        Ys.append(element.GetY())
    X = np.mean(Xs)
    Y = np.mean(Ys)
    #formule obtenues sur 
    #http://resources.esri.com/help/9.3/arcgisdesktop/com/gp_toolref/spatial_statistics_tools/how_directional_distribution_colon_standard_deviational_ellipse_spatial_statistics_works.htm
    #voir la formule de Wang http://www.portailsig.org/content/ellipsewangpy    
    #calcul des SDE
    SDEX = math.sqrt(np.sum([(Xe-X)**2 for Xe in Xs])/float(len(Xs)))
    SDEY = math.sqrt(np.sum([(Ye-Y)**2 for Ye in Ys])/float(len(Ys)))
    #calcul de l'angle
    print("Xmean : "+str(X))
    print("Ymean : "+str(Y))
    A = np.sum([(Xe-X)**2 for Xe in Xs])-np.sum([(Ye-Y)**2 for Ye in Ys])
    C = np.sum(([(Xe-X)*(Ye-Y) for Xe,Ye in zip(Xs,Ys)]))
    B = math.sqrt(A**2 +4*(C**2))
    print("A : "+str(A))
    print("B : "+str(B))
    print("C : "+str(C))
    Angle = math.atan((A+B)/(2*C))
    print("Angle : "+str(Angle))
    #calcul des deviations
    DevX = math.sqrt(2*(np.sum([((Xe-X)*math.cos(Angle)-(Ye-Y)*math.sin(Angle))**2 for Xe,Ye in zip(Xs,Ys)])/float(len(Xs))))
    DevY = math.sqrt(2*(np.sum([((Xe-X)*math.sin(Angle)-(Ye-Y)*math.cos(Angle))**2 for Xe,Ye in zip(Xs,Ys)])/float(len(Xs))))
    
    El = Ellipse(xy=(X,Y), width=DevX*3, height=DevY*3 ,angle=math.degrees(Angle))
    NewPoints = []
    for coord in El.get_verts() :
        Pt = ogr.Geometry(ogr.wkbPoint)
        Pt.AddPoint(coord[0],coord[1])
        NewPoints.append(Pt)

    return NewPoints
    
def EllipseWang(Points,nstd,**kwargs) :
    #creation de la matrice de covariance des points
    Coords=[]
    for element in Points :
        Coords.append([element.GetX(),element.GetY()])
    Coords = np.array(Coords)
    pos = Coords.mean(axis=0)
    cov = np.cov(Coords, rowvar=False)
    #Calcul des parametres
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * nstd * np.sqrt(vals)
    El = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)
    NewPoints = []
    for coord in El.get_verts() :
        Pt = ogr.Geometry(ogr.wkbPoint)
        Pt.AddPoint(coord[0],coord[1])
        NewPoints.append(Pt)
    return NewPoints
    
def EllipseWang2(Points,q=None, nsig=None, **kwargs) :
    """
    q : float, optional, Confidence level, should be in (0, 1)
    nsig : int, optional
        Confidence level in unit of standard deviations. 
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.
    """
    #creation de la matrice de covariance des points
    #http://stackoverflow.com/questions/12301071/multidimensional-confidence-intervals
    Coords=[]
    for element in Points :
        Coords.append([element.GetX(),element.GetY()])
    Coords = np.array(Coords)
    pos = Coords.mean(axis=0)
    cov = np.cov(Coords, rowvar=False)
    
    if q is not None:
        q = np.asarray(q)
    elif nsig is not None:
        q = 2 * norm.cdf(nsig) - 1
    else:
        raise ValueError('One of `q` and `nsig` should be specified.')
    r2 = chi2.ppf(q, 2)
    
    val, vec = np.linalg.eigh(cov)
    width, height = 2 * np.sqrt(val[:, None] * r2)
    rotation = np.degrees(np.arctan2(*vec[::-1, 0]))
    
    El = Ellipse(xy=pos, width=width, height=height, angle=rotation, **kwargs)
    NewPoints = []
    for coord in El.get_verts() :
        Pt = ogr.Geometry(ogr.wkbPoint)
        Pt.AddPoint(coord[0],coord[1])
        NewPoints.append(Pt)
    return NewPoints,OgrPoint((pos[0],pos[1]))


def DivideSquare(Geom,Number) :
    Envelope = Geom.GetEnvelope()
    L = int(math.sqrt(Number))
    if L != math.sqrt(Number) :
        print("La division doit contenir un nombre dont la racine carree est un entier")
    width = abs(abs(Envelope[0])-abs(Envelope[1])) / L
    Height = abs(abs(Envelope[2])-abs(Envelope[3])) / L
    Geoms=[]
    for i in range(L) :
        for j in range(L) :
            geom = Extent_geom((Envelope[0]+width*j,Envelope[0]+width*(j+1),Envelope[2]+i*Height,Envelope[2]+(i+1)*Height))
            Geoms.append(geom)
    return Geoms

def Extent_geom(geom,AsOgr=True) :
    """
    GetEnveloppe ==> (minX,minY,maxX,maxY)
    """
    try :
        geom[0]
        env = geom
    except TypeError :
        env = geom.GetEnvelope()
    TXT = "POLYGON(("+str(env[0])+" "+str(env[2])+","+str(env[1])+" "+str(env[2])+","+str(env[1])+" "+str(env[3])+","+str(env[0])+" "+str(env[3])+","+str(env[0])+" "+str(env[2])+"))"
#    new_geom = ogr.Geometry(ogr.wkbPolygon)
#    ring = ogr.Geometry(ogr.wkbLinearRing)
#    ring.AddPoint(env[0],env[2])
#    ring.AddPoint(env[1],env[2])
#    ring.AddPoint(env[1],env[3])
#    ring.AddPoint(env[0],env[3])
#    ring.AddPoint(env[0],env[2])
#    new_geom.AddGeometry(ring)
    if AsOgr :
        new_geom = ogr.CreateGeometryFromWkt(TXT)
    else :
        new_geom = shapely.wkt.loads(TXT)
    return new_geom
    
    
def OgrPoint(coords,AsOgr=True) :
    TXT = "POINT("
    for c in coords :
        TXT+=str(c)+" "
    TXT = TXT[:-1]+")"
    if AsOgr :
        Pt = ogr.CreateGeometryFromWkt(TXT)
    else :
        Pt = shapely.wkt.loads(TXT)
#    Pt = ogr.Geometry(ogr.wkbPoint)
#    if len(coords)==2:
#        Pt.AddPoint(coords[0],coords[1])
#    else :
#        Pt.AddPoint(coords[0],coords[1],coords[2])
        
    return Pt
    

def Translate3D(Geom,vector) :
    Shape= ToShapely(Geom)
    Translated = shapely.affinity.translate(Shape, xoff=vector[0], yoff=vector[1], zoff=vector[1])
    return ToOgr(Translated)
    

###################################################
###################################################
#_______________________________fonctions d'agregation
def CascadedUnion(Geoms) :
    WKT = [Geom.ExportToWkt() for Geom in Geoms]
    ShGeoms = [shapely.wkt.loads(e) for e in WKT]
    ShUnion = shapely.ops.cascaded_union(ShGeoms)
    Union=shapely.wkt.dumps(ShUnion)
    return ogr.CreateGeometryFromWkt(Union)

def MergeLines(Geoms) :
    if len(Geoms)>1 :
        ShGeoms=[ToShapely(Geom) for Geom in Geoms]
        ShLine = shapely.ops.linemerge(ShGeoms)
        Union = shapely.wkt.dumps(ShLine)
        NewLine=ogr.CreateGeometryFromWkt(Union)
        return NewLine
    else :
        return Geoms[0]
        
def AppendLine(L1,L2) :
    Extremites1 = GetExtremites(L1)
    buffera = Extremites1[0].Buffer(0.01)
    bufferb = Extremites1[1].Buffer(0.01)
    Extremites2 = GetExtremites(L2)
    #cas de figure ou les lignes ne se touche pas a une extremites
    if LineTouches(L1,L2)==False :
        print("ces deux lignes ne se touchent pas")
        raise ValueError
    #cas de figure ou la premiere ligne doit etre inversee
    if Extremites2[0].Intersect(buffera) or Extremites2[1].Intersect(buffera)  :
        L1 = GeomReverse(L1)
    #verifions si la deuxieme ligne doit etre inversee
    if Extremites2[1].Intersect(buffera) or Extremites2[1].Intersect(bufferb) :
        L2 = GeomReverse(L2)
    #combinaison des points
    AllPoints = []
    for coords in L1.GetPoints() :
        AllPoints.append(OgrPoint(coords))
    for coords in L2.GetPoints() :
        AllPoints.append(OgrPoint(coords))
    return LineFromPoints(AllPoints)
        

def MergeLines2(Lines) :
    L1 = Lines.pop(0)
    i = 0
    Continue=True
    while Continue==True :
        for element in Lines :
            if LineTouches(element,L1) :
                i+=1
                Lines.pop(Lines.index(element))
                L1 = AppendLine(L1,element)
        if i==0 :
            Continue=False
        else :
            i=0
    return L1,Lines

###################################################
###################################################
#_______________________________petites fonction sur geometrie

def GeomReverse(geom,Shapes=False) :
    if Shapes ==True :
        Coords = list(geom.coords)
        Coords.reverse()
        NewGeom = shapely.geometry.LineString(Coords)
        return NewGeom
    else :
        NpPt = geom.GetPointCount()
        NewGeom = ogr.Geometry(ogr.wkbLineString)
        for e in range(NpPt) :
            coords=geom.GetPoint(NpPt-1-e)
            NewGeom.AddPoint(coords[0],coords[1])
        return NewGeom


def TestSens (Ligne,Depart) :
    Dep = GetExtremites(Ligne)[0]
    if Dep.Intersects(Depart.Buffer(0.001))==False :
        return GeomReverse(Ligne)
    else :
        return Ligne

def LineFromPoints(Points) :
    Line = ogr.Geometry(ogr.wkbLineString)
    for element in Points :
        if abs(element.GetX())<1.1e+307 :
            Line.AddPoint(element.GetX(),element.GetY())
    return Line
    
def PolyFromPoints(Points,D3=False) :
    Poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for element in Points :
        if D3 :
            ring.AddPoint(element.GetX(),element.GetY(),element.GetZ())
        else :
            ring.AddPoint(element.GetX(),element.GetY())
    if D3 :
        ring.AddPoint(Points[0].GetX(),Points[0].GetY(),Points[0].GetZ())
    else :
        ring.AddPoint(Points[0].GetX(),Points[0].GetY())
    Poly.AddGeometry(ring)
    return Poly
    
def PointsFromLine(line,D3=False) :
    if line.GetGeometryName()=="LINESTRING" :
        Points=[]
        for i in range(line.GetPointCount()) :
            Pt = ogr.Geometry(ogr.wkbPoint)
            coords=line.GetPoint(i)
            if D3 :
                Pt.AddPoint(coords[0],coords[1],coords[2])
            else :
                Pt.AddPoint(coords[0],coords[1])
            Points.append(Pt)
    elif line.GetGeometryName()=="MULTILINESTRING" :
        Points=[]
        for e in range(line.GetGeometryCount()) :
            element=line.GetGeometryRef(e) 
            for i in range(element.GetPointCount()) :
                Pt = ogr.Geometry(ogr.wkbPoint)
                coords=element.GetPoint(i)
                if D3 :
                    Pt.AddPoint(coords[0],coords[1],coords[2])
                else :
                    Pt.AddPoint(coords[0],coords[1])
                Points.append(Pt)
    return Points
            

def LineSequence(line) :
    """
    Deconstruction d'une ligne en ses sous sous arretes
    """
    Count = line.GetPointCount()
    Lines=[]
    for e in range(Count-1) :
        co1 = line.GetPoint(e)
        co2 = line.GetPoint(e+1)
        NewLine = ogr.Geometry(ogr.wkbLineString)
        NewLine.AddPoint(co1[0],co1[1])
        NewLine.AddPoint(co2[0],co2[1])
        Lines.append(NewLine)
    return Lines

#fonction convertissant une ligne en un vecteur absolu

def To_vector(geom) :
    Dep,Arriv = GetExtremites(geom)
    vector=[Arriv.GetX()-Dep.GetX(),Arriv.GetY()-Dep.GetY()]
    return vector

#fonction permettant de projeter perpendiculairement un point sur une ligne
def Projection(point,line) :
    shape_point = ToShapely(point)
    shape_line = ToShapely(line)
    new_pt=shape_line.interpolate(shape_line.project(shape_point))
    result = ogr.Geometry(ogr.wkbPoint)
    result.AddPoint(new_pt.x, new_pt.y)
    return result


#Interpolation d'un point le long d'une Ligne
def Interpolate(line,distance,Start=None) :
    if Start!=None :
        Extremites = GetExtremites(line)
        if Start.Distance(Extremites[0])>Start.Distance(Extremites[1]) :
            line = GeomReverse(line)
    ShLine = shapely.wkt.loads(line.ExportToWkt())
    pt = ShLine.interpolate(distance)
    PT = ogr.CreateGeometryFromWkt(shapely.wkt.dumps(pt))
    return PT

def Prolongate(line,distance) :
    A = line.GetPoint(line.GetPointCount()-1)
    B = line.GetPoint(line.GetPointCount()-2)
    lenAB = Pytagore(A,B)
    Cx = B[0] + (B[0] - A[0]) / lenAB * (distance-lenAB)
    Cy = B[1] + (B[1] - A[1]) / lenAB * (distance-lenAB)
    line.AddPoint(Cx,Cy)
    return line


def SnappOnLineProp(points,line) :
    """Projection d'un ensemble de point sur une ligne par methode 
    des distances proportionnelle"""
    BasicLine = LineFromPoints(points)
    if BasicLine.Length()==0 :
        DrawPoints(points,"ro")
        DrawLine(line)
    #recuperation des points de departs et d'arrivee
    StartPoint=points[0]
    EndPoint=points[-1]
    ProjStartPoint = Projection(StartPoint,line)
    ProjEndPoint = Projection(EndPoint,line)
    DrawPoints([ProjStartPoint,ProjEndPoint],"go",width=20)
    #definition du rapport de distance
    ProjLine = PartOfLine(ProjStartPoint,ProjEndPoint,line)
    Rapport = ProjLine.Length()/BasicLine.Length()
    TotDist=0
    NewPoints=[ProjStartPoint]
    #snapping des points
    for e in range(len(points)-1) :
        if e>0 :
            pt1 = points[e]
            pt2 = points[e+1]
            TotDist += pt1.Distance(pt2) * Rapport
            NewPoints.append(Interpolate(ProjLine,TotDist))
    NewPoints.append(ProjEndPoint)
    return NewPoints
        
def SnappOnLine(points,line) :
    Proj=[Projection(element,line) for element in points]
    return Proj
    
#fonction pour decouper une geoemtrie
def SplitPolygone(Geom,RowNumber,LineNumber) :
    Extent = Geom.GetEnvelope()
    StartX=Extent[0]
    StartY = Extent[2]
    #decoupage de l'etendue
    Width = Extent[1]-Extent[0]
    Height = Extent[3]-Extent[2]
    PartWidth = Width/float(RowNumber)
    PartHeight = Height/float(LineNumber)
    #creation des etendue de decoupage
    Extents=[]
    TotX=0
    for i in range(RowNumber) :
        TotY=0
        for j in range (LineNumber) :
            Ext = Extent_geom((StartX+(i*PartWidth),StartX+((i+1)*PartWidth),StartY+(j*PartHeight),StartY+((j+1)*PartHeight)))
            Inter =  Ext.Intersection(Geom)    
            if Inter!=None and Inter.Area()>0 :
                Extents.append(Inter)
    return Extents


#fonction supprimant les points inutiles d'une ligne
def SimplifyLine(Line,ATolerance=0.001,DTolerance=0.001) :
    AllPoints = PointsFromLine(Line)
    P1 = AllPoints.pop(0)
    LastPoint=AllPoints[-1]
    FinalPoints=[P1]
    while len(AllPoints)>1 :
        Pt1 = FinalPoints[-1]
        Pt = AllPoints.pop(0)
        Pt2 = AllPoints[0]
        if Pt1.Distance(Pt)>DTolerance :
                L1 = LineFromPoints([Pt1,Pt])
                L2 = LineFromPoints([Pt,Pt2])
                v1 = To_vector(L1)
                v2 = To_vector(L2)
                Angle = Calculate_angle(v1,v2)
                if Angle>=ATolerance :
                    FinalPoints.append(Pt)
    FinalPoints.append(LastPoint)
    return LineFromPoints(FinalPoints)


#fonction de regression lineaire sur un ensemble de points OGR
def MeanLine(points) :
    X=[]
    Y=[]
    for element in points :
        X.append(element.GetX())
        Y.append(element.GetY())
    X = np.array(X)
    X = X.reshape(-1,1)
    Y = np.array(Y)
    Y = Y.reshape(-1,1)
    Reg = linear_model.LinearRegression()
    Reg = Reg.fit(X,Y)
    y1 = Reg.predict(np.min(X))[0][0]
    y2 = Reg.predict(np.max(X))[0][0]
    line = ogr.Geometry(ogr.wkbLineString)
    line.AddPoint(min(X)[0],float(y1))
    line.AddPoint(max(X)[0],float(y2))
    return line


#fonction pour calculer la distance orthogonal d'un point vers une ligne
def Distance_orthogonale(ligne,point) :
    sh_ligne = shapely.wkt.loads(ligne.ExportToWkt())
    sh_point = shapely.wkt.loads(point.ExportToWkt())
    projected = sh_ligne.interpolate(sh_ligne.project(sh_point))
    distance = projected.distance(sh_point)
    return distance

#fonction pour calculer un angle entre deux vecteurs
def Calculate_angle(v1,v2) :
    AB = math.sqrt(v1[0]**2+v1[1]**2)
    CD = math.sqrt(v2[0]**2+v2[1]**2)
    if AB*CD !=0 :
        a =abs(v1[0]*v2[0]+v1[1]*v2[1])/(AB*CD)
    else :
        return 0
    if a>1 and a<1.00000000000002 :
        Acos = 0
    elif a<=1 :
        Acos = math.acos(a)
    else :
        raise ValueError("Wrong value for this angle")

    angle = math.degrees(Acos)
    return angle


# fonction pour trouver l'azimuth entre deux point avec comme 0 direction Nord
def NorthAngle(pt1, pt2):
  """Calculate azimuth angle from two points. (Zero is north.)"""
  import math
  # recup coords
  x1, y1, x2, y2 = pt1.GetX(), pt1.GetY(), pt2.GetX(), pt2.GetY()

  dx, dy = (x2 - x1, y2 - y1)
  angle = 90 - math.degrees(math.atan2(dy, dx))
  if angle < 0:
      angle += 360
  return angle

#fonction pour trouver l'angle de brisure entre le centre de trois geometries
def BrisureAngle(Geoms) :
    if len(Geoms)!=3 :
        print ("This treatment need exactly 3 geometry")
    Centres = [g.Centroid() for g in Geoms]
    #trouver la geometrie centrale
    Centre=None
    Liste=[]
    for Id in range(len(Centres)) :
        geom = Centres[Id]
        Dist=0
        for Id2 in range(len(Centres)) :
            geom2 = Centres[Id2]
            if Id!=Id2 :
                Dist+=geom.Distance(geom2)
        Liste.append((geom,Dist))
    Sort = sorted(Liste,key=lambda tup : tup[1])
    Centre=Sort[0][0]
    Extrem=[Sort[1][0],Sort[2][0]]
    #trouver les deux vecteurs
    L1 = ogr.Geometry(ogr.wkbLineString)
    L2 = ogr.Geometry(ogr.wkbLineString)
    L1.AddPoint(Centre.GetX(),Centre.GetY())
    L2.AddPoint(Centre.GetX(),Centre.GetY())
    L1.AddPoint(Extrem[0].GetX(),Extrem[0].GetY())
    L2.AddPoint(Extrem[1].GetX(),Extrem[1].GetY())
    V1 = To_vector(L1)
    V2 = To_vector(L2)
    return Calculate_angle(V1,V2)


###################################################
###################################################
#_______________________________fonctions topologiaues

#verifier si une geometrie a bien ses deux extremites intersectees par d'autres dans une liste
def Check_extremites(geom,geoms,tolerance,check_id=False):
    if check_id==False :
        g1=geom
    else :
        g1=geom[0]
        ID=geom[1]
    extremites=GetExtremites(g1)
    b1=extremites[0].Buffer(tolerance)
    b2=extremites[1].Buffer(tolerance)
    inter_b1=False
    inter_b2=False
    if check_id==True :
        for element in geoms :
            s_ID=element[1]
            if ID!=s_ID :
                s_geom=element[0]
                s_extremites=GetExtremites(s_geom)
                sb1=s_extremites[0].Buffer(tolerance)
                sb2=s_extremites[1].Buffer(tolerance)
                if sb1.Intersect(b1) or sb2.Intersect(b1) :
                    inter_b1=True
                if sb1.Intersect(b2) or sb2.Intersect(b2) :
                    inter_b2=True
        if inter_b1==False or inter_b2==False :
            return False
        else :
            return True
    else :
        for element in geoms :
            s_geom=element[0]
            s_extremites=GetExtremites(s_geom)
            if s_extremites[0].Intersect(b1) or s_extremites[1].Intersect(b1) :
                inter_b1=True
            if s_extremites[0].Intersect(b2) or s_extremites[1].Intersect(b2) :
                inter_b2=True
        if inter_b1==False or inter_b2==False :
            return False
        else :
            return True

#verifier si une geometrie a bien ses deux extremites intersectees dans un layer
def Check_extremites_layer(feat,layer,tolerance) :
    geom = feat.GetGeometryRef()
    ID = feat.GetFID()
    extremites = GetExtremites(geom)
    b1 = extremites[0].Buffer(tolerance)
    b2 = extremites[1].Buffer(tolerance)
    inter_b1=False
    inter_b2=False
    layer.SetSpatialFilter(Extent_geom(b1))
    if layer.GetFeatureCount()>1 :
        inter_b1=True
    layer.SetSpatialFilter(Extent_geom(b2))
    if layer.GetFeatureCount()>1 :
        inter_b2=True
    layer.SetSpatialFilter(None)
    if inter_b1==True and inter_b2==True :
        return True
    else :
        return False
        
#verifier si une geometrie a bien ses extremites intersectee dans un dictionnaire
def CheckExtremityDict(Geom,Geoms) :
    Extremites=GetExtremites(Geom[1])
    E0=False
    E1=False
    for ID,geom in Geoms.items() :
        if ID!=Geom[0] :
            Ex2 = GetExtremites(geom)
            b1 = Ex2[0].Buffer(0.001)
            b2 = Ex2[1].Buffer(0.001)
            if b1.Intersect(Extremites[0]) or b2.Intersect(Extremites[0]) :
                E0=True
            if b1.Intersect(Extremites[1]) or b2.Intersect(Extremites[1]) :
                E1=True
            if E0 and E1 :
                return (True,True)
    return (E0,E1)

#verifier si deux lignes se touchent a leurs extremites
def LineTouches(L1,L2) :
    Extremites1 = GetExtremites(L1)
    buffera = Extremites1[0].Buffer(0.01)
    bufferb = Extremites1[1].Buffer(0.01)
    Extremites2 = GetExtremites(L2)
    #cas de figure ou les lignes ne se touche pas a une extremites
    if Extremites2[0].Intersect(buffera)==False and Extremites2[0].Intersect(bufferb)==False and Extremites2[1].Intersect(buffera)==False and Extremites2[1].Intersect(bufferb)==False :
        return False
    else :
        return True

#fonction permettant de recuperer les extremites d'une geometrie sous forme de point
def GetExtremites(geom) :
    coord1 = geom.GetPoint(0)
    coord2 = geom.GetPoint(geom.GetPointCount()-1)
    p1 = ogr.Geometry(ogr.wkbPoint)
    p2 = ogr.Geometry(ogr.wkbPoint)
    p1.AddPoint(coord1[0],coord1[1])
    p2.AddPoint(coord2[0],coord2[1])
    return (p1,p2)
    
#fonction pour creer un convex hull d'un ensemble de geometries
def ConvexHull(Geoms) :
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for Geom in Geoms :
        geomcol.AddGeometry(Geom)
    convexhull = geomcol.ConvexHull()
    return convexhull

#fonction permettant d'identifier toutes les lignes fusionnables a partir d'une
def LookForFusionnable(Feat,Layer,DiffFields, CheckField,Tolerance) :
    Ids = defaultdict(lambda : None)
    Ids[Feat.GetFID()]=True
    Fusionnable=defaultdict(lambda : None)
    Geom = Feat.GetGeometryRef()
    Extremites = GetExtremites(Geom)
    CheckIt=list(Extremites)
    while len(CheckIt)>0 :
        Ext = CheckIt.pop()
        Layer.SetSpatialFilter(Ext.Buffer(Tolerance))
        if Layer.GetFeatureCount()==2 :
            for voisin in Layer :
                VId = voisin.GetFID()
                if Ids[VId]==None :
                    Same = True
                    Ids[VId]=True
                    for field in CheckField :
                        if field not in DiffFields :
                            if Feat.GetField(field)!=voisin.GetField(field) :
                                Same = False
                    if Same ==True :
                        Fusionnable[VId]=True
                        VGeom = voisin.GetGeometryRef()
                        VExtremites=GetExtremites(VGeom)
                        CheckIt.append(VExtremites[0])
                        CheckIt.append(VExtremites[1])
    return Fusionnable

###################################################
###################################################
#_______________________________Eclatement de ligne


def SplittingOnintersect(geom,intersections) :
    #conversion en shapely
    wkt = geom.ExportToWkt()
    shape_ligne = shapely.wkt.loads(wkt)
    coords = list(shape_ligne.coords)
    shape_intersections = []
    intersections_coord = []
    for feat in intersections :
        element=feat.GetGeometryRef()
        ide=feat.GetFID()
        pt = shapely.wkt.loads(element.ExportToWkt())
        if pt.geometryType()=="MultiPoint" :
            for pts in pt :
                dist = shape_ligne.project(pts)
                shape_intersections.append([pts,dist,ide])
                intersections_coord.append((pts.x,pts.y))
        else :
            dist = shape_ligne.project(pt)
            shape_intersections.append([pt,dist,ide])
            intersections_coord.append((pt.x,pt.y))

    #on cree une nouvelle ligne comprenant tous les points
    passage = 0
    temp_line = ogr.Geometry(ogr.wkbLineString)
    for a in range(len(coords)) :
        pt = shapely.geometry.Point(coords[a])
        distance = shape_ligne.project(pt)
        shape_intersections.append([pt,distance])
    #on ordonne les points et leurs distances
    shape_intersections.sort(cmp=comp)
    for element in shape_intersections :
        pt = element[0]
        temp_line.AddPoint(pt.x,pt.y)
    # separation de la geometrie en sortie
    final_lines=[]
    new_line = ogr.Geometry(ogr.wkbLineString)
    for a in range(temp_line.GetPointCount()) :
        ptt = temp_line.GetPoint(a)
        pt =(ptt[0],ptt[1])
        if pt not in intersections_coord :
            new_line.AddPoint(pt[0],pt[1])
            if a == (temp_line.GetPointCount()-1) :
                final_lines.append([new_line,-1])
        else :
            new_line.AddPoint(pt[0],pt[1])
            index=intersections_coord.index(pt)
            final_lines.append([new_line,shape_intersections[index][2]])
            new_line = ogr.Geometry(ogr.wkbLineString)
            new_line.AddPoint(pt[0],pt[1])
    return final_lines


#eclatement de ligne selon un ensemble de point, nouvelle version
def SplittingOnintersectV2(geom,intersections,Shapes=False) :
    """
    Eclatement d'une Ligne a partir d'un ensemble de point (a reprendre avec la fonction frac)
    """
    #verification si on a pas des multipoints
    intersections2=[]
    for element in intersections :
        if element.GetGeometryName()=="MULTIPOINT" :
            for e in range(element.GetGeometryCount()) :
                intersections2.append(element.GetGeometryRef(e))
        else :
            intersections2.append(element)
    if Shapes==False :
        if geom.GetGeometryName()=="MULTILINESTRING" :
            points = PointsFromLine(geom)
            line = LineFromPoints(points)
            ShLine = ToShapely(line)
        else :
            ShLine = ToShapely(geom)
    else :
        ShLine=copy.deepcopy(geom)
        geom = ogr.CreateGeometryFromWkt(shapely.wkt.dumps(ShLine))
    ListePoints=[]
    FinalLines=[]
    #recuperation des distances pour chaque point au debut de la ligne
    for element in intersections2 :
        if Shapes==False :
            Shelement = shapely.wkt.loads(element.ExportToWkt())
        else :
            Shelement=element
        Distance = ShLine.project(Shelement)
        if Shapes==False :
            ListePoints.append([element,Distance])
        else :
            ListePoints.append([ogr.CreateGeometryFromWkt(shapely.wkt.dumps(Shelement)),Distance])
    #recuperation du depart et de l'arrivee egalement
    Coords = ShLine.coords
    Dep = shapely.geometry.Point(Coords[0])
    Arriv=shapely.geometry.Point(Coords[-1])
    DistanceDep = ShLine.project(Dep)
    DistanceArriv = ShLine.project(Arriv)
    ListePoints.append([ogr.CreateGeometryFromWkt(shapely.wkt.dumps(Dep)),DistanceDep])
    ListePoints.append([ogr.CreateGeometryFromWkt(shapely.wkt.dumps(Arriv)),DistanceArriv])

    #ordonnencement des ces points par rapport a cette distance
    ListePoints.sort(cmp=comp)
    #iteration sur les points deux a deux pour extraire les fragments de ligne
    for i in range(len(ListePoints)) :
        if i<len(ListePoints)-1 :
            pt1 = ListePoints[i][0]
            pt2 = ListePoints[i+1][0]
            Frag = PartOfLine(pt1,pt2,geom)
            FinalLines.append(Frag)
    return FinalLines


def SplitOnIntersectV3(Line, Points):
    line = ToShapely(Line)
    points=[ToShapely(element) for element in Points]
    # First coords of line
    coords = list(line.coords)

    # Keep list coords where to cut (cuts = 1)
    cuts = [0] * len(coords)
    cuts[0] = 1
    cuts[-1] = 1

    # Add the coords from the points
    coords += [list(p.coords)[0] for p in points]    
    cuts += [1] * len(points)        

    # Calculate the distance along the line for each point    
    dists = [line.project(shapely.geometry.Point(p)) for p in coords]    
    # sort the coords/cuts based on the distances    
    # see http://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list    
    coords = [p for (d, p) in sorted(zip(dists, coords))]    
    cuts = [p for (d, p) in sorted(zip(dists, cuts))]          

    # generate the Lines    
    #lines = [LineString([coords[i], coords[i+1]]) for i in range(len(coords)-1)]    
    lines = []        

    for i in range(len(coords)-1):    
        if cuts[i] == 1:    
            # find next element in cuts == 1 starting from index i + 1   
            j = cuts.index(1, i + 1)    
            lines.append(shapely.geometry.LineString(coords[i:j+1]))            
    #get lines as ogr geometry
    Lines =[ToOgr(element) for element in lines]
    return Lines


#decoupage d'une ligne en troncon de ligne egaux
def SplitLineByDist(Line,Dist) :
    if (Line.Length()>=Dist) :
        Nb = int(Line.Length()/Dist)
        Pts=[]
        for e in range(Nb) :
            e=e+1
            dist = Dist*e
            Pts.append(Interpolate(Line,dist))
        Splitted = SplitOnIntersectV3(Line,Pts)
    else :
        Splitted=[Line]
    return Splitted
                


#Extraction d'une partie d'une ligne a partir de deux points
def PartOfLine(pt1,pt2,line) :
    sh_pt1=ToShapely(pt1)
    sh_pt2=ToShapely(pt2)
    sh_line=ToShapely(line)
    name = line.GetGeometryName()
    if name=="MULTILINESTRING" :
        print("shit this is a multilinestring")
        print("tentative de sauvetage 1")
        Rep = Sauvetage2(line)
        if Rep==False :
            print("tentative de sauvetage 2")
            #arrondissement des coordonees
            RoundedLine=RoundCoords(line)
            UnderLines = [RoundedLine.GetGeometryRef(e) for e in range(RoundedLine.GetGeometryCount())]
            Merged=MergeLines(UnderLines)
            if Merged.GetGeometryName()=="LINESTRING" :
                line=Merged
                sh_line=ToShapely(line)
            else :
                print("tentative de sauvetage 3")
                Merged=Sauvetage1(RoundedLine)
                UnderLines = [Merged.GetGeometryRef(e) for e in range(Merged.GetGeometryCount())]
                Merged=MergeLines(UnderLines)
                if Merged.GetGeometryName()=="LINESTRING" :
                    line=Merged
                    sh_line=ToShapely(line)
                else :
                    print ("unable to Save this line")
                    print(line.ExportToWkt())
                    input()
        else :
            line=Rep
            sh_line=ToShapely(line)
        
    d1=sh_line.project(sh_pt1)
    d2=sh_line.project(sh_pt2)
    if d1>d2 :
        premier = d2
        keeped_points=[sh_pt2]
        deuxieme = d1
    else :
        premier =d1
        keeped_points=[sh_pt1]
        deuxieme = d2
    
    coords = sh_line.coords
    for element in coords :
        point = shapely.geometry.Point(element)
        d3 = sh_line.project(point)
        #recuperation des points entre depart et arrivee
        if d3>premier and d3<deuxieme :
            keeped_points.append(point)
    if d1==premier :
        keeped_points.append(sh_pt2)
    else :
        keeped_points.append(sh_pt1)
    geom = ogr.Geometry(ogr.wkbLineString)
    for element in keeped_points :
        geom.AddPoint(element.x,element.y)
    return geom

def Sauvetage1(MultiLineString) :
    Lines = [MultiLineString.GetGeometryRef(e) for e in range(MultiLineString.GetGeometryCount())]
    for e in range(len(Lines)-1) :
        L1 = Lines[e]
        L2=Lines[e+1]
        if L1.Distance(L2)<1 and L1.Intersect(L2)==False :
            AddedLine = ogr.Geometry(ogr.wkbLineString)
            Ex1 = GetExtremites(L1)
            Ex2 = GetExtremites(L2)
            D1 = Ex1[0].Distance(Ex2[0])
            D2 = Ex1[0].Distance(Ex2[1])
            D3 = Ex1[1].Distance(Ex2[0])
            D4 = Ex1[1].Distance(Ex2[1])
            if D1<D2 and D1<D3 and D1<D4 :
                 AddedLine.AddPoint(Ex1[0].GetX(),Ex1[0].GetY())
                 AddedLine.AddPoint(Ex2[0].GetX(),Ex2[0].GetY())
            elif D2<D1 and D2<D3 and D2<D4 :
                 AddedLine.AddPoint(Ex1[0].GetX(),Ex1[0].GetY())
                 AddedLine.AddPoint(Ex2[1].GetX(),Ex2[1].GetY())
            elif D3<D1 and D3<D2 and D3<D4 :
                 AddedLine.AddPoint(Ex1[1].GetX(),Ex1[1].GetY())
                 AddedLine.AddPoint(Ex2[0].GetX(),Ex2[0].GetY())
            elif D4<D1 and D4<D2 and D4<D3 :
                 AddedLine.AddPoint(Ex1[1].GetX(),Ex1[1].GetY())
                 AddedLine.AddPoint(Ex2[1].GetX(),Ex2[1].GetY())
            MultiLineString.AddGeometry(AddedLine)
    return MultiLineString

def Sauvetage2(Geom) :
    Lines = [Geom.GetGeometryRef(e) for e in range(Geom.GetGeometryCount())]
    LinesLen =[(Line,Line.Length()) for Line in Lines]
    Biggest = max(LinesLen,key=lambda x : x[1])
    for line,length in LinesLen :
        if line!=Biggest[0] :
            Coords = [line.GetPoint(e) for e in range(line.GetPointCount())]
            for co in Coords :
                Pt = ogr.CreateGeometryFromWkt("POINT ("+str(co[0])+" "+str(co[1])+")")
                Buffer = Pt.Buffer(0.001)
                if Buffer.Intersect(Biggest[0])==False :
                    #si un des points des sous lignes n'appartient pas a la ligne principale
                    return False
    #sinon on retourne juste la ligne principale
    return Biggest[0]
    

###################################################
###################################################
#_______________________________Proximite


#fonction permettant de trouver l'objet le plus proche d'un autre dans un autre layer
def nearest(feat,layer,distance,interdit=[]) :
    geom = feat.GetGeometryRef()
    BUFFER = geom.Buffer(distance)
    layer.SetSpatialFilter(BUFFER)
    i=1
    while(layer.GetFeatureCount()<len(interdit)+1) :
        i+=1
        BUFFER=geom.Buffer(distance*i)
        layer.SetSpatialFilter(BUFFER)
    feats=[]
    for a in range(layer.GetFeatureCount()) :
        voisin = layer.GetNextFeature()
        voisin_id = voisin.GetFID()
        voisin_geom = voisin.GetGeometryRef()
        dist = voisin_geom.Distance(geom)
        feats.append([dist,voisin_id])

    #verification de qui est le plus proche
    dist_ok =9999999
    id_ok =""
    for element in feats :
        if element[0]<dist_ok and element[1] not in interdit :
            dist_ok=element[0]
            id_ok = element[1]
    return [dist_ok,id_ok]


#fonction permettant de trouver tous les objets dans un rayon x
def Get_in_buffer(geom,layer,x,interdit=[]) :
    ids=[]
    BUFFER = geom.Buffer(x)
    layer.SetSpatialFilter(Extent_geom(BUFFER))
    for e in range(layer.GetFeatureCount()) :
        feat = layer.GetNextFeature()
        if feat.GetFID() not in interdit :
            ids.append(feat.GetFID())
    layer.SetSpatialFilter(None)
    return ids



#Retrouver l'objet le plus proche dans un Layer et retrouner son identifiant
def NearestObject(feat,layer,restrictions={}) :
    """
    A partir d'un objet en entree, on cherche a recuperer l'identifiant de l'objet le plus proche dans un
    autre layer
    restrictions est un dictionnaire qui permet d'exclure certains objets qui ne correspondent pas a certaines conditions :
        champ : 'valeur ... condition'
    """
    search=100
    try :
        geom = feat.GetGeometryRef()
    except :
        geom = feat
    BUFFER = geom.Buffer(search)
    layer.SetSpatialFilter(BUFFER)
    i=1
    while layer.GetFeatureCount()<1 :
        i+=1
        BUFFER = geom.Buffer(search*i)
        layer.SetSpatialFilter(BUFFER)
    #on itere sur les lignes recuperees
    ok_dist = 99999
    ID_ligne = 9999999
    for a in range(layer.GetFeatureCount()) :
        voisin = layer.GetNextFeature()
        #on verifie que les conditions soient respectees
        if restrictions !={} :
            for field,expression in restrictions.items() :
                value = voisin.GetField(field)
                expression = expression.replace('valeur',str(value))
                #si l'expression est respectee on verifie alors la distance
                if eval(expression) :
                    voisin_geom = voisin.GetGeometryRef()
                    dist = voisin_geom.Distance(geom)
                    if dist<ok_dist :
                        ok_dist = dist
                        ID_ligne = voisin.GetFID()
        else :
            voisin_geom = voisin.GetGeometryRef()
            dist = voisin_geom.Distance(geom)
            if dist<ok_dist :
                ok_dist = dist
                ID_ligne = voisin.GetFID()
    return ID_ligne


#obtention de K objets les plus proches d'une feature dans un layer
def K_NearestObject(feat,layer,K,restrictions={}) :
    layer.SetSpatialFilter(None)
    p=100
    try :
        geom = feat.GetGeometryRef()
        checkId=True
    except AttributeError :
        geom=feat
        checkId=False
    layer.SetSpatialFilter(geom.Buffer(p))
    u=1
    while layer.GetFeatureCount()<K+1 :
        u+=1
        layer.SetSpatialFilter(geom.Buffer(p*u))
    Dists=[]
    for Feat in layer :
        if checkId==True :
            if Feat.GetFID()!=feat.GetFID() :
                Dists.append((geom.Distance(Feat.GetGeometryRef()),Feat))
        else :
            Dists.append((geom.Distance(Feat.GetGeometryRef()),Feat))
    return sorted(Dists,key=lambda x :x[0])[0:K]


#obtention de tous les voisins direct et indirect d'une geometrie
def GetNeighbours(feat,layer,Ids=[],K=float("inf")) :
    geom = feat.GetGeometryRef()
    layer.SetSpatialFilter(geom.Buffer(100))
    Ids=[feat.GetFID()]
    News=[]
    for f in layer :
        g = f.GetGeometryRef()
        if g.Intersects(geom)==True and f.GetFID() not in Ids :
            Ids.append(f.GetFID())
            News.append(f.GetFID())
    u=0
    while len(News)>0 and u<K :
        u+=1
        Feats = [layer.GetFeature(e) for e in Ids]
        Geoms = [f.GetGeometryRef() for f in Feats]
        Geom = CascadedUnion(Geoms)
        layer.SetSpatialFilter(Geom.Buffer(0.001))
        News=[]
        for f in layer :
            g = f.GetGeometryRef()
            if g.Intersects(Geom)==True and f.GetFID() not in Ids :
                Ids.append(f.GetFID())
                News.append(f.GetFID())
    return [layer.GetFeature(e) for e in Ids]
    
def GetNeighboursList(geom,Geoms) :
    """
    geom : (ID,Geom)
    Geoms : {ID:Geom}
    """
    Intersected=[]
    for ID,G in Geoms.items() :
        if ID!=geom[0]:
            if G.Intersect(geom[1]):
                Intersected.append((ID,G))
    return Intersected

###############################################################################
## functions for geometries application on roads
###############################################################################

def PerpendicularPoint(Line,Dist,Start=True) :
    Extremites = GetExtremites(Line)
    ax = Extremites[0].GetX()
    ay = Extremites[0].GetY()
    bx = Extremites[1].GetX()
    by = Extremites[1].GetY()
    dx = ax - bx
    dy = ay - by
    
    dist = math.sqrt(dx * dx + dy * dy)
    
    offset = Dist
    
    normX = dx / dist
    normY = dy / dist
    
    xPerp = offset * normX
    yPerp = offset * normY
    
    
    cx = ax + yPerp
    cy = ay - xPerp
    dx = ax - yPerp
    dy = ay + xPerp
    
    ex = bx - yPerp
    ey = by + xPerp
    fx = bx + yPerp
    fy = by - xPerp
    if Start :
        return (OgrPoint((cx,cy)),OgrPoint((dx,dy)))
    else :
        return (OgrPoint((ex,ey)),OgrPoint((fx,fy)))
        

def BissectorPoints(Line1,Line2,Dist) :
    V1 = To_vector(Line1)
    V2 = To_vector(Line2)
    AngleLines=Calculate_angle(V1,V2)
    Extr = GetExtremites(Line2)
    AngleNorth = NorthAngle(Extr[0],Extr[1])
    if Extr[0].GetY()>Extr[1].GetY() :
        HalfAngle = (AngleLines)/2.0
        Sens=True
    else :
        Sens=False
        HalfAngle = (AngleLines-90)/2.0
    
    OkAngle1 = AngleNorth+HalfAngle+90
    OkAngle2 = AngleNorth+HalfAngle-90
    X1 = Extr[0].GetX() + Dist * math.cos(math.radians(OkAngle1))
    Y1 = Extr[0].GetY() + Dist * math.sin(math.radians(OkAngle1))
    X2 = Extr[0].GetX() + Dist * math.cos(math.radians(OkAngle2))
    Y2 = Extr[0].GetY() + Dist * math.sin(math.radians(OkAngle2))
    if Sens :
        return (OgrPoint((X1,Y1)),OgrPoint((X2,Y2)))
    else :
        return (OgrPoint((X2,Y2)),OgrPoint((X1,Y1)))
        
def NearestPoint(Pt,Pts) :
    Dists = [(Pti,Pt.Distance(Pt1)) for Pti in Pts]
    Dists.sort(key=lambda T : T[1])
    return Dists[0][0]

def CreateRoad(Line,Width) :
    """
    Fonction pour creer une forme de route a partir d'une ligne et d'une largeur
    """
    #operation avec shapely
    RefPts = PointsFromLine(Line,True)
    Shape = ToShapely(Line)
    LeftOffseted = Shape.parallel_offset(Width/2.0,"left",join_style=3,mitre_limit=10.0)
    RightOffseted = Shape.parallel_offset(Width/2.0,"right",join_style=3,mitre_limit=10.0)
    RightC=[]
    LeftC=[]
    #tentative de recuperer les geometries 3D
#    for OldCoords,LeftCoords,RightCoords in zip(Shape.coords,LeftOffseted.coords,RightOffseted.coords):
#        RightC.append((RightCoords[0],RightCoords[1],OldCoords[2]))
#        LeftC.append((LeftCoords[0],LeftCoords[1],OldCoords[2]))
    #pour le left
    for Coord in LeftOffseted.coords :
        Pt = NearestPoint(OgrPoint(Coord),RefPts)
        LeftC.append((Coord[0],Coord[1],Pt.GetZ()))
        
    for Coord in RightOffseted.coords :
        Pt = NearestPoint(OgrPoint(Coord),RefPts)
        RightC.append((Coord[0],Coord[1],Pt.GetZ()))
        
    #pour le right
    LeftLine = shapely.geometry.LineString(LeftC)
    RightLine = shapely.geometry.LineString(RightC)
    #retour a ogr    
    Left = ToOgr(LeftLine)
    Right = ToOgr(RightLine)
    #creation de la geometrie
    Pts1 = PointsFromLine(Left,True)
    Pts2 = PointsFromLine(Right,True)
    Pts1.extend(Pts2)
    return (PolyFromPoints(Pts1,True))
    
    
        
    

if __name__=="__main__" :
    
    WKT1 = "LINESTRING(3 3,5 5)"
    WKT2 = "LINESTRING(5 5,7 3)"
    WKT3 = "LINESTRING(3 3,5 5,7 3,8 6)"
    WKT4 = "LINESTRING(7 3,8 6)"
#    Line1 = ogr.CreateGeometryFromWkt(WKT1)
#    Line2 = ogr.CreateGeometryFromWkt(WKT2)
    Line3 = ogr.CreateGeometryFromWkt(WKT3)
#    Line1 = ogr.CreateGeometryFromWkt(WKT2)
#    Line2 = ogr.CreateGeometryFromWkt(WKT4)
#    Pts = BissectorPoints(Line1,Line2,1)
#    DrawPoints(Pts,"bo")

    Poly = CreateRoad(Line3,1)
    DrawPolygone(Poly,(0,0,0,0),(0,1,0))
    
    DrawLine(Line3)
    SetView((0,10,0,10))
    
    
    