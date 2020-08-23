# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 10:26:23 2017

@author: GelbJ
HELP : https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.dtypes.html
"""
#from __future__ import unicode_literals

import os
if 'GDAL_DATA' not in os.environ:
    #os.environ['GDAL_DATA'] = 'C:/Users/gelbj/Anaconda2/lib/site-packages/osgeo/data/gdal' #Ordi Maison
    os.environ['GDAL_DATA'] = 'C:/Users/GelbJ/Documents/New_Anaconda27/Library/share/gdal' #Ordi Labo

import numpy as np
import JGeom
from .QuadTree import Index as QdIndex
from numpy.lib import recfunctions
from osgeo import ogr,osr
import copy
import ast
import sqlite3
from copy import deepcopy
import re
from collections import defaultdict
import chardet

##################################################
##Fonctions utilitaires
##################################################
def multifind(string, value, start = 0, stop = None):
    values = []
    while True:
        found = string.find(value, start, stop)
        if found == -1:
            break
        values.append(found)
        start = found + 1
    return values

    
def ReadStrType(String,ForceString=False) :
    try:
        String = ast.literal_eval(String)
    except SyntaxError as e :
        if ForceString : 
            return String
        else :
            raise SyntaxError(str(e))
    return String

    
def GessNumpyType(String) :
    try :
        Value = float(String)
        if ("." in String) :
            return "float32"
        else :
            return "int32"
    except ValueError :
        return "|S25"
        
def GetSqlLiteType(value) :
    if type(value)==str or type(value)==unicode :
        return "TEXT"
    elif type(value)==int or type(value)==bool:
        return "INTEGER"
    elif type(value)==float :
        return "REAL"

def GetSqlLiteType(TypeDef) :
    if TypeDef.GetTypeName()=="String" :
        return "TEXT"
    if TypeDef.GetTypeName()=="Real" :
        return "REAL"
    if TypeDef.GetTypeName()=="Integer" :
        return "INTEGER"
    if TypeDef.GetTypeName()=="Integer64" :
        return "INTEGER"
    if TypeDef.GetTypeName()=="Date" :
        return "TEXT"
    raise ValueError("This type is not actualy programmed : "+str(TypeDef.GetTypeName()))


def GetNumpyType(TypeDef) :
    if TypeDef.GetTypeName()=="String" :
        Width = TypeDef.GetWidth()
        if Width >0 :
            return '|S'+str(Width)
        else : 
            return '|S600'
    if TypeDef.GetTypeName()=="Real" :
        return 'float32'
    if TypeDef.GetTypeName()=="Integer" :
        return 'int32'
    if TypeDef.GetTypeName()=="Integer64" :
        return 'int64'
    if TypeDef.GetTypeName()=="Date" :
        return 'datetime64[D]'
    raise ValueError("This type is not actualy programmed : "+str(TypeDef.GetTypeName()))
    
def GetOgrType(TEXT) :
    if '|S' in TEXT :
        return ogr.OFTString
    elif "float" in  TEXT :
        return ogr.OFTReal 
    elif "int" in TEXT :
        return ogr.OFTInteger
    elif "datetime64[D]" in TEXT :
        return ogr.OFTDate
    raise ValueError("This type is not actualy programmed : "+TEXT)
    
def GetTextBetween(String,element) :
    Parts=[]
    Part=""
    Open=False
    for e in String :
        if e==element and Open==False :
            Open=True
        elif e==element and Open==True :
            Open=False
            Parts.append(Part)
            Part=""
        if Open==True and e!=element :
            Part+=e
    return Parts


def ListLayerInSdb(File) : 
    Driver=ogr.GetDriverByName('SQLite')
    Source = Driver.Open(File)
    Layers=[]
    for e in range(Source.GetLayerCount()) :
        Layer = Source.GetLayer(e)
        Layers.append(Layer.GetName())
    Source.Destroy()
    return Layers

###############################################################################
##JDatatable from Files
###############################################################################
def JDataTableFromCSV(File,Sep=";",ID="OID",Types=None) : 
    """
    Fonction pour generer une JDataTable from a csv
    """
    CSV = open(File,"r")
    Header = CSV.readline().replace("\n","").split(Sep)
    L1 = CSV.readline().replace("\n","").split(Sep)
    if Types is None :
        Types = [GessNumpyType(D) for D in L1]
    if len(Header) != len(Types) : 
        print("There is an error here, specified types must have the same length as Header")
        print("Header : "+str(Header))
        print("Types : "+str(Types))
        raise ValueError("Error concerning Field Number")
    NeedOID = False
    if ID == "OID" and "OID" not in Header: 
        Header.append("OID")
        Types.append("int32")
        NeedOID = True
    
    Line = CSV.readline()
    Datas=[]
    i=0
    while Line is not None and Line != "" : 
        Data = Line.replace("\n","").split(Sep)
        OkData = [D for D in Data]
        if NeedOID : 
            OkData.append(i)
        i+=1
        Datas.append(tuple(OkData))
        Line = CSV.readline()
        
    print(Header)
    print(Types)
    print(Datas[0:5])
    Table = JDataTable(Header,Types,ID,Datas)
    return Table


###############################################################################
##JDatatable :
# Objet permettant de stocker une table attributaire 
# a l'aide d'un tableau numpy
#
###############################################################################
class JDataTable(object) :

    def __init__(self,Fields,Types,ID,Datas=[]) :
        """
        Fields est une liste de nom de champs
        Types est une liste de type (string) compris par numpy
        Datas est une liste de tuple comprenant les donnees
        ID : nom d'un champ etant l'identifiant unique
        """
        self.Fields = Fields
        self.Types = Types
        self.Defaults={}
        self.ID = ID
        self._Ids ={}
        try :
            self._Array = np.array(Datas,dtype=[(f,t) for f,t in zip(Fields,Types)])
        except : 
            print(Fields)
            print(Types)
            print(Datas)
            L1 = Datas[0]
            
            self._Array = np.array(Datas,dtype=[(f,t) for f,t in zip(Fields,Types)])
        
    def __repr__(self) :
        return str(self._Array)
        
    def FieldDesc(self) :
        TXT=""
        for Field,Type in zip(self.Fields,self.Types) :
            TXT+=Field+" : "+Type+"\n"
        return TXT
        
    def __getitem__(self,ID) :
        return self._Array[self._Ids[ID]]
    
    def GetRow(self,i) : 
        return self._Array[i]
        
    @property
    def FeatureCount(self) :
        return len(self._Array)
		
    def SetValue(self,ID,Field,Value) :
        self._Array[self._Ids[ID]][Field] = Value
        
    def GetModalities (self,Field) :
        return np.unique(self._Array[Field])
        
    def UpdateID(self) :
        """
        permet de mettre a jour le lien en ID et numeros de ligne
        """
        Index = self.Fields.index(self.ID)
        Dico={self._Array[i][Index] : i for i in range(len(self._Array))}
        self._Ids = Dico
        
    def Iterate(self,NiceFeat=False) :
        """
        permet d'iterer sur toutes les entrees
        """
        for i in range (len(self._Array)) :
            if NiceFeat == False :
                yield self._Array[i]
            else : 
                yield {Field:Value for Field,Value in zip(self.Fields,self._Array[i])}
        
    def AddRow(self,Row) :
        """
        Row : liste de valeur a ajouter au tableau
        """
        NewRow = np.array(Row,dtype=self._Array.dtype)
        self._Array = np.append(self._Array,NewRow)
        Index = self.Fields.index(self.ID)
        
    def MakeEmptyTable(self) : 
        Table = JDataTable(self.Fields,self.Types,self.ID)
        return Table
        
    def DeleteRow(self,ID,UpdateID=True) :
        """
        permet de supprimer une entree par son ID
        met a jour en consequence les matching ID et numeros de ligne si demande
        """
        Row = self._Ids[ID]
        del self._Ids[ID]
        self._Array = np.delete(self._Array,Row)
        if UpdateID :
            self.UpdateID()
                
    def DeleteField(self,Field) :
        """
        permet de supprimer un champ
        """
        if Field in self.Fields:
            Index = self.Fields.index(Field)
            del self.Fields[Index]
            del self.Types[Index]
        self._Array = self._Array[self.Fields]
        
    def AddField(self,Field,Type,FillValue) :
        """
        permet de rajouter un champ
        """
        self.Defaults[Field] = FillValue
        self.Fields.append(Field)
        self.Types.append(Type)
        NewCol = np.zeros(len(self._Array),dtype=Type)
        NewCol.fill(FillValue)
        self._Array = recfunctions.append_fields(self._Array,Field,NewCol,usemask=False)
        
    def RenameField(self,OldName,NewName) :
        """
        permet de renommer un champ
        """
        Index = self.Fields.index(OldName)
        self.Fields[Index] = NewName
        self._Array.dtype.names = self.Fields
        
    def CalculateField(self,FieldName,Func,**kwargs) :
        """
        permet de calculer un champ
        def Add(Feat) :
          return Feat["Field1"]+Feat["Field2"]
        """
        for i,Feat in enumerate(self.Iterate()) :
            Value = Func({F:V for F,V in zip(self.Fields,Feat)},**kwargs)
            self._Array[i][FieldName] = Value
            
    def Join(self,SecondTable,Field,JoinType="inner",PostFix1='',PostFix2='jn_') :
        """
        jointype : {'inner', 'outer', 'leftouter'}, optional
        If 'inner', returns the elements common to both r1 and r2.
        If 'outer', returns the common elements as well as the elements of
        r1 not in r2 and the elements of not in r2.
        If 'leftouter', returns the common elements and the elements of r1
        not in r2.
        """
        NewArray = recfunctions.join_by(Field,self._Array,SecondTable,JoinType,PostFix1,PostFix2)
        
    def AttributeFilter(self,Filters) :
        """
        permet de filter par des conditions : "( (Field1 == 18) AND (Field2 == 'titi') ) OR (Field2 == 'toto')"
        NB : chaque condition doit etre specifiee entre parentheses ... (c'est plate je sais)
        """
        Array = self._Array
        for element in self.Fields :
            #Filters=Filters.replace(element,'Array["'+element+'"]')
            Filters=re.sub(r"\b"+element+r"\b",'Array["'+element+'"]',Filters)
#            Places = multifind(Filters," "+element+" ")
#            while len(Places)>0 :
#                Place = Places.pop()
#                Filters = Filters[0:Place]+" Array['"+element+"']"+Filters[Place+1+len(element):len(Filters)]
#                Places = multifind(Filters," "+element+" ")
#                    #Filters[Place+1:(Place+1+len(element))] = "Array["+element+"]"
        Filters = Filters.replace("AND","&").replace("OR","|")
        Filters = 'Array['+Filters+']'
#        Filters = "Array["+Filters+"]"
        NewArray = eval(Filters)
        NewTable = JDataTable(self.Fields,self.Types,self.ID)
        NewTable._Array = NewArray
        NewTable.UpdateID()
        return NewTable
    
    def GetFirsFeatureLike(self,Filters) : 
        """
        Permet de retourner ma premiere feature (ordre d'enregistrement) qui valide les caracteristiques demandee
        "( Field1 == 18 AND Field2 == 'titi' ) OR Field2 == 'toto'"
        """
        Table = self.AttributeFilter(Filters)
        return Table._Array[0]
        
        
    def GetVector(self,Field) :
        return self._Array[Field]
        
        
    def Sort(self,by) :
        """
        permet de trier le tableau de donnees selon une ou plusieurs colonnes
        """
        self._Array = np.sort(self._Array,order=by)
        self.UpdateID()
        
    def RandomSplit(self,Prt) :
        """
        Permet de scinder la table en 2 partie, la premiere correspondant au prt
        indique et la deuxieme au reste de la table
        """
        Nb = int(Prt*self.FeatureCount)
        Old_Array = copy.deepcopy(self._Array)
        np.random.shuffle(Old_Array)
        Array1 = list(Old_Array[:Nb])
        Array2 = list(Old_Array[Nb:(self.FeatureCount-1)])
        print(Array1)
        print(Array2)
        NewTable1 = JDataTable(self.Fields,self.Types,self.ID,Array1)
        NewTable2 = JDataTable(self.Fields,self.Types,self.ID,Array2)
        return(NewTable1,NewTable2)
        
    def SaveAsCsv(self,Path,sep=";") :
        """
        Permet d'exporter la table en CSV
        """
        File = open(Path,"w")
        Header = sep.join(self.Fields)+"\n"
        File.write(Header)
        for Row in self.Iterate() :
            Line = sep.join([str(element) for element in Row])+"\n"
            File.write(Line)
        File.close()
        
            


###############################################################################
##JBigDatatable :
# Objet permettant de stocker une table attributaire 
# a l'aide d'une bd sqlite en memoire vive
#
###############################################################################
class JBigDataTable(object) :

    def __init__(self,Fields,Types,ID,Datas=[]) :
        """
        Fields est une liste de nom de champs
        Types est une liste de type (string) compris par numpy
        Datas est une liste de tuple comprenant les donnees
        ID : nom d'un champ etant l'identifiant unique
        """
        self.Fields = Fields
        self.Types = Types
        self.Defaults={}
        self.ID = ID
        self._Ids ={}
        self.DBase = sqlite3.connect(":memory:")
        #creation de la table dans la BD
        self.Cursor = self.DBase.cursor()
        Request = "CREATE TABLE Data("
        for Field,Type in zip(self.Fields, self.Types) :
            if Field!=ID : 
                Request+=Field+" "+Type+", "
            else :
                Request+=Field+" "+Type+" PRIMARY KEY, "
        Request = Request[0:-2]+");"
        self.Cursor.execute(Request)
        #rajout des donnes dans la table
        for Row in Datas :
            Request = "INSERT INTO Data ("+",".join(self.Fields)+") VALUES ("
            for element in Row :
                if type(element)==str or type(element)==unicode :
                    Request+="'"+str(element)+"', "
                else :
                    Request+=str(element)+", "
            Request = Request[0:-2]+");"
            print(Request)
            self.Cursor.execute(Request)
        #self._Array = np.array(Datas,dtype=[(f,t) for f,t in zip(Fields,Types)])
     
    def FieldDesc(self) :
        TXT=""
        for Field,Type in zip(self.Fields,self.Types) :
            TXT+=Field+" : "+Type+"\n"
        return TXT
        
    def __getitem__(self,Id) :
        if self.Types[self.Fields.index(self.ID)]=="TEXT" :
            Request = "SELECT * FROM Data WHERE "+self.ID+" = '"+Id+"';"
        else :
            Request = "SELECT * FROM Data WHERE "+self.ID+" = "+str(Id)
        self.Cursor.execute(Request)
        Datas =self.Cursor.fetchone()
        return Datas
        
    @property
    def FeatureCount(self) :
        self.Cursor.execute("SELECT Count(*) FROM Data")
        return self.Cursor.fetchone()[0]

		
    def SetValue(self,ID,Field,Value) :
        if type(Value)==str or type(Value)==unicode :
            Value = "'"+Value+"'"
        if self.Types[self.Fields.index(self.ID)]=="TEXT" :
            ID = "'"+ID+"'"
        Request = "UPDATE Data SET "+Field+" = "+str(Value)+" WHERE "+self.ID+"="+str(ID)+";"
        print(Request)        
        self.Cursor.execute(Request)
        
    def GetModalities (self,Field) :
        Request = "SELECT DISTINCT "+Field+" FROM Data;"
        self.Cursor.execute(Request)
        return [element[0] for element in self.Cursor.fetchall()]
        
    def Iterate(self) :
        """
        permet d'iterer sur toutes les entrees
        """
        self.Cursor.execute("SELECT * FROM Data")
        for element in self.Cursor.fetchall() :
            yield element
            
    def AddRow(self,Row) :
        """
        Row : liste de valeur a ajouter au tableau
        """
        Request = "INSERT INTO Data ("+",".join(self.Fields)+") VALUES ("+str(Row)[1:-1]+");"
        self.Cursor.execute(Request)
        
    def DeleteRow(self,ID) :
        """
        permet de supprimer une entree par son ID
        """
        if self.Types[self.Fields.index(self.ID)]=="TEXT" :
            ID = "'"+ID+"'"
        Request = "DELETE FROM Data WHERE "+self.ID+" = "+ID+";"
        self.Cursor.execute(Request)
        
        
    def DeleteField(self,Field) :
        """
        permet de supprimer un champ
        """
        print("Wraning : this operation is not well supported by sqlite3")
        if Field in self.Fields:
            Index = self.Fields.index(Field)
            del self.Fields[Index]
            del self.Types[Index]
        KeepedFields = ",".join(self.Fields)
        Request = "CREATE TABLE Data_backup AS SELECT "+KeepedFields+" FROM Data;"
        print(Request)
        self.Cursor.execute("CREATE TABLE Data_backup AS SELECT "+KeepedFields+" FROM Data;")
        self.Cursor.execute("DROP TABLE Data;")
        self.Cursor.execute("ALTER TABLE Data_backup RENAME TO Data;")      
        
        
    def AddField(self,Field,Type,FillValue) :
        """
        permet de rajouter un champ
        """
        self.Defaults[Field] = FillValue
        self.Fields.append(Field)
        self.Types.append(Type)
        self.Cursor.execute("ALTER TABLE Data ADD COLUMN "+Field+" "+Type+";")
        if type(FillValue)==str :
            FillValue = "'"+FillValue+"'"
        self.Cursor.execute("UPDATE Data SET "+Field+" = "+str(FillValue))
        
    def CalculateField(self,FieldName,Func) :
        """
        permet de calculer un champ
        def Add(Feat) :
          return Feat["Field1"]+Feat["Field2"]
        """
        print("Warning : this function is not optimised and will perform a lot or request")
        for i,Feat in enumerate(self.Iterate()) :
            Feature = {F:V for F,V in zip(self.Fields,Feat)}
            Value = Func(Feature)
            self.SetValue(Feature[self.ID],FieldName,Value)
            
            
    def AttributeFilter(self,Filters) :
        """
        permet de filter par des conditions : "( Field1 == 18 AND Field2 == 'titi' ) OR Field2 == 'toto'"
        """
        Request = "SELECT * FROM Data WHERE "+Filters+";"
        print(Request)
        self.Cursor.execute(Request)
        Datas = list(self.Cursor.fetchall())
        print(Datas)
        NewTable = JBigDataTable(self.Fields,self.Types,self.ID,Datas)
        return NewTable

###############################################################################
##QuadSpatialTree :
# Objet permettant de stocker des geometries dans un index spatial
# Sans dependance C++
###############################################################################
class QuadSpatialTree(object) :
    """
    un objet permettant de ranger des geometries et de les retrouver plus
    vite via un index spatial
    
    Bbox : (xmin,ymin,xmax,ymax)
    
    Ne necessite pas de dependance C++
    Geoms : {"ID":Geom}
    """
    def __init__(self,Geoms) :
        BBox = JGeom.GetExtent(Geoms.values())
        self.Index = QdIndex((BBox[0],BBox[2],BBox[1],BBox[3]))
        self.Geoms = Geoms
        for Id,element in Geoms.items() :
            Env = element.GetEnvelope()
            self.Index.insert(Id,((Env[0],Env[2],Env[1],Env[3])))
            
    def AddGeometry(self,Geom,Id) :
        #python 2.7
#        if self.Geoms.has_key(Id)==False :
#            self.Geoms[Id]=Geom
#            Env = Geom.GetEnvelope()
#            self.Index.insert(Id,(Env[0],Env[2],Env[1],Env[3]))
#        else :
#            raise KeyError("This SpatialTree already know this object")
            
        if Id not in self.Geoms :
            self.Geoms[Id]=Geom
            Env = Geom.GetEnvelope()
            self.Index.insert(Id,(Env[0],Env[2],Env[1],Env[3]))
        else :
            raise KeyError("This SpatialTree already know this object")
    def RemoveGeometry(self,Id) :
        self.Geoms[Id] = None
        
            
    def NearestGeom(self,Geom,MaxDist=float("inf"),Jump=100) :            
        """
        permet de retrouver la geometrie la plus proche de celle demandee
        retourne None si aucune geometrie n'a ete trouvee
        retourne un tuple (ID,Geom,Dist) si une geometrie a ete trouvee        
        """
        Dist = Jump
        Continue = True
        while Continue :
            Buffer = Geom.Buffer(Dist)
            Dist+=Jump
            Results=self.SpatialQuery(Buffer)
            if len(Results)>0 or Dist>MaxDist :
                Continue = False
        if len(Results)==0 :
            return None
        else :
            Geoms =[(ID,geom,geom.Distance(Geom)) for ID,geom in Results.items()]
            Geoms.sort(key=lambda x : x[2])
            return Geoms[0]
        
            
    def SpatialQuery(self,Geom) :        
    #gestion des polygones et points
        if Geom.GetGeometryName()!="LINESTRING" :
            Buffer = Geom.Buffer(0.1)
            Env = Buffer.GetEnvelope()
            Results=list(self.Index.intersect((Env[0],Env[2],Env[1],Env[3])))
            RealResults = {}
            for Id in Results :
                if self.Geoms[Id]!=None :
                    if (Geom.Intersect(self.Geoms[Id])) :
                        RealResults[Id]=self.Geoms[Id]
            return RealResults
        #gestion des linestrings (decoupage pour avoir moins d'intersection)
        elif Geom.GetGeometryName()=="LINESTRING" :
            Parts = JGeom.LineSequence(Geom)
            GeomsDics =[]
            for element in Parts :
                Buff = element.Buffer(0.1)
                Env = Buff.GetEnvelope()
                GeomsDics.append(self.Index.intersect((Env[0],Env[2],Env[1],Env[3])))
            RealResults={}
            for element in GeomsDics :
                for Id in element :
                    if self.Geoms[Id]!=None :
                        if (Geom.Intersect(self.Geoms[Id])==True) :
                            RealResults[Id]=self.Geoms[Id]
            return RealResults


###############################################################################
##JFastLayer
#objet permettant de charger un layer en memoire vive pour le travailler
###############################################################################

class JFastLayer(object) :
    
    def __init__(self,File,LayerName=0) :
        self.Source = File
        self.LayerName = LayerName

        
    def Initialize(self,ID,GeomIndex=True,Params={}) :
        self.ID = ID
        Ext = self.Source.split(".")[-1]
        if Ext=="shp" or Ext in ("sdb","sqlite") :
            Datas,Geoms = self._ReadFromShp(Params)
        elif Ext=="csv" :
            Datas,Geoms = self._ReadFromCsv(Params)
        elif Ext == "geojson" :
            Datas,Geoms = self._ReadFromShp(Params)
        else :
            raise ValueError("Extension de fichier noon geree")
        self._Populate(Datas,Geoms,GeomIndex)
        
    def BuildSpIndex(self) : 
        self.SpatialTree = QuadSpatialTree(self.Geoms)
        
    def MakeItEmpty(self,Params={}, OID = True) :
        self.Source=None
        self.LayerName = None
        self.ID = "OID"
        self.Fields = []
        self.FieldsTypes =[]
        self.SpatialTree = None
        self.Geoms = {}
        for key,value in Params.items() :
            self.__dict__[key]=value
        if "OID" not in self.Fields and OID and self.ID=="OID" :
            self.Fields.append("OID")
            self.FieldsTypes.append("int32")
        self.AttrTable = JDataTable(self.Fields,self.FieldsTypes,Datas=[],ID=self.ID)
        
    def CreateEmptyLayer(self) : 
        NewLayer = JFastLayer("")
        NewLayer.Geoms = {}
        NewLayer.SpatialTree = None
        NewLayer.ID = self.ID
        NewLayer.Fields = self.Fields
        NewLayer.FieldsTypes = self.FieldsTypes
        NewLayer.SpatialRef = self.SpatialRef
        NewLayer.GeomType = self.GeomType
        NewLayer.AttrTable = JDataTable(self.Fields,self.FieldsTypes,Datas=[],ID=self.ID)
        return NewLayer
        
    def AppendFeat(self,Feat,Geom):
        self.AttrTable.AddRow(tuple([Feat[Field] for Field in self.Fields]))
        if self.SpatialTree!=None :
            self.SpatialTree.AddGeometry(Geom,Feat[self.ID])
        else :
            self.Geoms[Feat[self.ID]]=Geom
        
    def RemoveFeat(self,Id,UpdateID=True) :
        self.SpatialTree.RemoveGeometry(Id)
        self.AttrTable.DeleteRow(Id,UpdateID)
        
        
    def BuildContiguity(self,Mode="Queen",Tolerance = 0.1) : 
        """
        Permet de construire une matrice de contiguite
        """
        self.Contiguity = defaultdict(lambda : [])
        for Feat in self.Iterate(True) :
            Buff = Feat["Geom"].Buffer(Tolerance)
            Neighbours = self.SpatialFilter(Buff)
            for Neighbour in Neighbours.Iterate(True) : 
                if Neighbour[self.ID]!=Feat[self.ID] :
                    if Mode == "Queen" : 
                        self.Contiguity[Feat[self.ID]].append(Neighbour[self.ID])
                    elif Mode == "Rook" : 
                        raise ValueError("the Rook mode is not currently implemented (because of tolerance pb)")
                    else : 
                        raise ValueError("the Mode must be either Queen or Rook")
        
        
        
        
        
    def __repr__(self) :
        Txt = """____________________
        Layer Object from : """+self.Source+"""
            FIELDS :
"""+self.FieldDesc()+"""
            Geometry : """+self.GeomType+"""
            Feature Count : """+str(self.FeatureCount)+"""
____________________"""
        return Txt
        
    def FieldDesc(self) :
        return self.AttrTable.FieldDesc()
        
    def Draw(self,color=(1,0,0),border =(0,0,0)) :
        if "LINESTRING" in self.GeomType :
            for Key,Geom in self.Geoms.items() :
                JGeom.DrawLine(Geom,color,border)
        elif "POLYGON" in self.GeomType :
            for Key,Geom in self.Geoms.items() :
                JGeom.DrawPolygone(Geom,color,border)
        elif "POINT" in self.GeomType :
            JGeom.DrawPoints(self.Geoms.values(),color)

    @property        
    def FeatureCount(self) :
        return len(self.AttrTable._Array)
        
    def __getitem__(self,ID) :
        return (self.AttrTable[ID],self.Geoms[ID])
        
    def GetRow(self,Nb) : 
        Index = self.Fields.index(self.ID)
        Row = self.AttrTable._Array[Nb]
        Geom = self.Geoms[Row[Index]]
        return self.NiceFeat((Row,Geom))
        
#    def Iterate(self,Nice=False) :
#        for key in self.Geoms.keys() :
#            if Nice==False :
#                yield(self[key])
#            else :
#                yield(self.NiceFeat(self[key]))
        
#    def Iterate(self,Nice=False) :
#        Index = self.Fields.index(self.ID)
#        for i in range (len(self.AttrTable._Array)) :
#            Row = self.AttrTable._Array[i]
#            if Nice==False :
#                yield(self[Row[Index]])
#            else :
#                yield(self.NiceFeat(self[Row[Index]]))
        
    def Iterate(self,Nice=False ) : 
        Index = self.Fields.index(self.ID)
        for i in range(self.FeatureCount) : 
            Row = self.AttrTable._Array[i]
            if Nice : 
                Geom = self.Geoms[Row[Index]]
                yield self.NiceFeat((Row,Geom))
            else : 
                yield Row
            
            
    def NiceFeat(self,Feat) :
        Dico={}
        for Field,Value in zip(self.Fields,Feat[0]) :
            Dico[Field]=Value
        Dico["Geom"] = Feat[1]
        return Dico
        
    def CalculateFieldWithGeom(self,Field,Func,Verbose=False,**kwargs) :
        """
        Func(Feat,**kwargs)
        Func doit retourner la valeur
        """
        if Verbose :
            i=0
            for key,geom in self.Geoms.items() :
                i+=1
                if i%5 == 0 :
                    print("Avancement : "+str(round(float(i)/float(self.FeatureCount)*100,2)))
                Feat = self.NiceFeat(self[key])
                self.AttrTable._Array[self.AttrTable._Ids[key]][Field] = Func(Feat,**kwargs)
        else :
            for key,geom in self.Geoms.items() :
                Feat = self.NiceFeat(self[key])
                self.AttrTable._Array[self.AttrTable._Ids[key]][Field] = Func(Feat,**kwargs)
            
    def AttributeFilter(self,Filters,GeomIndex=True) : 
        """
        permet de filter par des conditions : "( Field1 == 18 AND Field2 == 'titi' ) OR Field2 == 'toto'"
        """
        Filtered = self.AttrTable.AttributeFilter(Filters)
        Geoms = {ID : self.Geoms[ID] for ID in Filtered._Ids.keys()}
        NewLayer = JFastLayer(File="")
        NewLayer.MakeItEmpty(Params={"Fields":self.Fields,"FieldsTypes":self.FieldsTypes,"SpatialRef":self.SpatialRef,"GeomType":self.GeomType,"ID":self.ID})
        NewLayer.AttrTable = Filtered
        if GeomIndex :
            NewLayer.SpatialTree = QuadSpatialTree(Geoms)
            NewLayer.Geoms = NewLayer.SpatialTree.Geoms
        else :
            NewLayer.Geoms = Geoms
        return NewLayer
    
        
    @property
    def Extent(self) :
        return JGeom.GetExtent(self.Geoms.values())
            
    def _ReadFromCsv(self,Params) :
        """
        recuperation des donnees a partir d'un CSV
        """
        File = open(self.Source,"r")
        SpatialRef = osr.SpatialReference()
        SpatialRef.ImportFromEPSG(Params["EPSG"])
        self.SpatialRef=SpatialRef
        self.Fields = []
        self.FieldsTypes=[]
        if self.ID =="OID" :
            self.Fields.append("OID")
            self.FieldsTypes.append("int32")
        Header = File.readline()
        Fields=[element for element in Header.replace("\n","").split(Params["Sep"])]
        XIndex = Fields.index(Params["X"])
        YIndex = Fields.index(Params["Y"])
        self.Fields = self.Fields+Fields
        Line1 = File.readline()
        self.FieldsTypes=self.FieldsTypes+[GessNumpyType(element) for element in Line1.replace("\n","").split(Params["Sep"])]
        self.GeomType = "POINT"
        Datas=[]
        Geoms={}
        File.seek(0)
        Header = File.readline()
        Line = File.readline()
        i=0
        while Line!="" :
            Elements = Line.replace("\n","").split(Params["Sep"])
            Point = ogr.Geometry(ogr.wkbPoint)
            Point.AddPoint(float(Elements[XIndex]),float(Elements[YIndex]))
            # if Params.has_key("SpatialFilter") : 2.7 syntax
            if "SpatialFilter" in Params :
                if Point.Intersect(Params["SpatialFilter"]) :
                    AddIt=True
                else :
                    AddIt=False
            else :
                AddIt = True
            if AddIt :
                Values = [i]
                Values = Values + [ReadStrType(e) for e in Elements]
                Datas.append(Values)
                Geoms[i] = Point
                i+=1
            Line = File.readline()
        return Datas,Geoms
        
    def _ReadFromShp(self,Params) :
        """
        Recupere les donnees a partir d'un SHP
        """
        ##connexion au fichier
        Shp = ogr.Open(self.Source)
        Layer = Shp.GetLayer(self.LayerName)
        #if Params.has_key("SpatialFilter") : python 2.7 syntax
        if "SpatialFilter" in Params :
            Layer.SetSpatialFilter(Params["SpatialFilter"])
        
        ##recuperation des champs
        Def = Layer.GetLayerDefn()
        if Params.has_key("EPSG") : 
            self.SpatialRef = osr.SpatialReference()
            self.SpatialRef.ImportFromEPSG(Params["EPSG"])
        else :
            self.SpatialRef = Layer.GetSpatialRef().Clone()
        self.Fields=[]
        self.FieldsTypes=[]
        if self.ID =="OID" :
            self.Fields.append("OID")
            self.FieldsTypes.append("int32")
        for e in range(Def.GetFieldCount()) :
            FDef = Def.GetFieldDefn(e)
            if FDef.GetName()!="OID" :
                self.Fields.append(FDef.GetName())
                self.FieldsTypes.append(GetNumpyType(FDef))
        F1 = Layer.GetNextFeature()
        G1 = F1.GetGeometryRef()
        Layer.ResetReading()
        self.GeomType = G1.GetGeometryName()
        ##recuperation des donnees
        Datas =[]
        Geoms ={}
        for Feat in Layer :
            Values=[]
            if Feat.GetGeometryRef()!=None :
                for F,T in zip(self.Fields,self.FieldsTypes) :
                    if F=="OID" :
                        Values.append(Feat.GetFID())
                    elif T == "datetime64[D]" :
                        Values.append(Feat.GetField(F).replace("/","-"))
                    elif "|S" in T : 
                        Value = Feat.GetField(F)
                        if Value is None :
                            Value = 'None'
                        if (Params.has_key("Decode")) :
                            if Params["Decode"]=="NONE" :
                                Values.append(Value.encode("UTF8"))
                            else :
                                try :
                                    V = Value.decode(Params["Decode"])
                                    Values.append(V.encode("utf-8"))
                                except UnicodeDecodeError : 
                                    Codec = chardet.detect(Value)
                                    raise ValueError("The encoding specified seems to be wrong... try this one : "+str(Codec))
                        else :
                            Values.append(str(Value).encode("utf-8"))
                    else :
                        Values.append(Feat.GetField(F))
                    if self.ID == "OID" :
                        Geoms[Feat.GetFID()]=Feat.GetGeometryRef().Clone()
                    else :
                        Geoms[Feat.GetField(self.ID)]=Feat.GetGeometryRef().Clone()
                Datas.append(tuple(Values))
                #Datas.append(tuple([Feat.GetField(F) for F in self.Fields]))
            ##
        return Datas,Geoms
        
    def _Populate(self,Datas,Geoms,GeomIndex) :
        """
        Cree l'index spatial et la table d'attributs
        """
        self.AttrTable = JDataTable(self.Fields,self.FieldsTypes,Datas=Datas,ID=self.ID)
        self.AttrTable.UpdateID()
        if GeomIndex :
            self.SpatialTree = QuadSpatialTree(Geoms)
            self.Geoms = self.SpatialTree.Geoms
        else :
            self.Geoms = Geoms
            
    def UpdateTree(self) :
        self.SpatialTree = QuadSpatialTree(self.Geoms)
        self.Geoms = self.SpatialTree.Geoms
            
    def Save(self,File,GetTxt=False,LayerName=None) :
        if File.split(".")[-1]=="shp":
            self._SaveShp(File)
        elif File.split(".")[-1]=="geojson" :
            self._SaveGeoJson(File,GetTxt)
        elif File.split(".")[-1]=="sdb" :
            self._SaveSdb(File,LayerName)
        else :
            raise ValueError("Actually save only support shapefile format")
            
    def _SaveSdb(self,File,LayerName) : 
        """
        Sauvegarde du layer dans une SBD
        """
        ## si la BD n'existe pas, il faut la creer
        Driver=ogr.GetDriverByName('SQLite')
        if os.path.isfile(File)==False :
            Database = Driver.CreateDataSource(File,options=['SPATIALITE=yes'])

        Source = Driver.Open(File,1)
        
        if self.GeomType == "LINESTRING" :
            geom_type = ogr.wkbLineString25D
        elif self.GeomType == "MULTILINESTRING" : 
            geom_type = ogr.wkbMultiLineString
        if self.GeomType == "POINT" :
            geom_type = ogr.wkbPoint
        elif self.GeomType == "MULTIPOINT" : 
            geom_type = ogr.wkbMultiPoint
        if self.GeomType == "POLYGON" :
            geom_type = ogr.wkbPOLYGON
        elif self.GeomType == "MULTIPOLYGON" : 
            geom_type = ogr.wkbMultiPolygon
        
        print(LayerName)
        #si le layer existe deja, il faut le supprimer
        if Source.GetLayerByName(LayerName) is not None :
            Source.DeleteLayer(LayerName)
        
        SortieLayer = Source.CreateLayer(str(LayerName), srs = self.SpatialRef,geom_type = geom_type ,options=['FORMAT=SPATIALITE'])
        #creation des champs du layer
        for Name,Type in zip(self.Fields,self.FieldsTypes) :
            if ("|S") in Type :
                NewField = ogr.FieldDefn(Name, GetOgrType(Type))
                NewField.width=int(Type.split("S")[1])
            else :
                NewField = ogr.FieldDefn(Name, GetOgrType(Type))
            SortieLayer.CreateField(NewField)
        #remplissage du Shp
        FDefn = SortieLayer.GetLayerDefn()
        for e,row in enumerate(self.Iterate(True)) :
            Feature = ogr.Feature(FDefn)
            for Field,Type in zip(self.Fields,self.FieldsTypes) :
                value = row[Field]
                Feature.SetGeometry(row["Geom"])
                if Type == "float32" or Type == "float64" : 
                    Feature.SetField(Field,float(value))
                elif Type == "int32" or Type == "int64" : 
                    Feature.SetField(Field,int(value))
                elif "|" in Type : 
                    Feature.SetField(Field,str(value))
                else :
                    print("Unknow type : "+Type)
                    Feature.SetField(Field,value)
            SortieLayer.CreateFeature(Feature)
        del SortieLayer
        Source.Destroy()
            
            
    def _SaveShp(self,File) :
        """
        Sauvegarde du layer sous forme d'un SHP
        """
        #creation du shapefile de sortie
        Driver=ogr.GetDriverByName('ESRI Shapefile')
        Sortie = Driver.CreateDataSource(File)
        SortieLayer = Sortie.CreateLayer('layer',self.SpatialRef)
        
        #creation des champs du layer
        for Name,Type in zip(self.Fields,self.FieldsTypes) :
            if ("|S") in Type :
                NewField = ogr.FieldDefn(Name, GetOgrType(Type))
                NewField.width=int(Type.split("S")[1])
            else :
                NewField = ogr.FieldDefn(Name, GetOgrType(Type))
            SortieLayer.CreateField(NewField)
        #remplissage du Shp
        for e,row in enumerate(self.Iterate(True)) :
            Feature = ogr.Feature(SortieLayer.GetLayerDefn())
            Feature.SetFID(int(e))
            for Field,Type in zip(self.Fields,self.FieldsTypes) :
                value = row[Field]
                Feature.SetGeometry(row["Geom"])
                if len(Field)>=10 :
                    Field = Field[0:10]
                if Type == "float32" or Type == "float64" : 
                    Feature.SetField(Field,float(value))
                elif Type == "int32" or Type == "int64" : 
                    Feature.SetField(Field,int(value))
                elif "|" in Type : 
                    Feature.SetField(Field,str(value))
                else :
                    print("Unknow type : "+Type)
                    Feature.SetField(Field,value)
            SortieLayer.CreateFeature(Feature)
        Sortie.Destroy()
        
        
    def _SaveGeoJson(self,File,GetTxt=False) :
        """
        Sauvegarde du layer sous forme d'un geojson
        """
        Text = ""
        ProjWkt = self.SpatialRef.ExportToPrettyWkt().split("\n")[-1]
        try :
            EPSG = GetTextBetween(ProjWkt,'"')[1]
        except IndexError :
            EPSG ="Unknown"
        Text+="""{
"type": "FeatureCollection",
"crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::"""+EPSG+"""" } },
                                                                                                                              
"features": [\n"""
        e=0
        for Feat in self.Iterate(True) :
            Geom = Feat["Geom"]
            e+=1
            ##preparation des attributes
            String="""{"type":"Feature","properties":{"""
            for i,Field in enumerate(self.Fields) :
                if Field!="OID" :
                    value = Feat[Field]
                    ##controle des Nan
                    if str(value)=="nan" or str(value)=="None" :
                        value = -999
                    if "|S" in self.FieldsTypes[i] or "date" in self.FieldsTypes[i]:
                        String+='"'+Field+'"'+' : "'+str(value)+'", '
                    else :
                        String+='"'+Field+'" : '+str(value)+", "
            String = String[:-2]+'},"geometry":'
            ##preparation de la geometrie
            Gtxt = Geom.ExportToJson()
            if e<self.FeatureCount :
                String+=Gtxt+'},\n'
            else :
                String+=Gtxt+'}]\n'
            Text+=String
        
        Text+="}"
        if GetTxt :
            return Text
        else :
            Sortie = open(File,"w")
            Sortie.write(Text)
            Sortie.close()
        
        
    def Reproj(self,EPSG) :
        """
        fonction pour reprojeter les geometries
        """
        if type(EPSG)==int :
            To = osr.SpatialReference()
            res = To.ImportFromEPSG(EPSG)
            if res!=0 :
                raise ValueError("EPSG non connu, impossible d'effectuer la transformation")
                return None
        elif type(EPSG)==str :
            To = osr.SpatialReference()
            To.ImportFromWkt(EPSG)
        else :
            To = EPSG
        From = self.SpatialRef
        transform = osr.CoordinateTransformation(From, To)
        AllGeoms = copy.deepcopy(self.Geoms)
        NewGeoms={}
        for key,Geom in AllGeoms.items() :
            Geom.Transform(transform)
            NewGeoms[key] = Geom
        self.SpatialTree = QuadSpatialTree(NewGeoms)
        self.Geoms = self.SpatialTree.Geoms
        self.SpatialRef = To
        
        
    def SpatialFilter(self,Geom,GeomIndex=True,Reverse=False,Intersections=False) : 
        """
        fonction pour retourner un nouveau layer filtre selon l'intersection
        a une geometrie
        """
        Intersected=self.SpatialTree.SpatialQuery(Geom)
        if len(Intersected)==0 :
            return None
        
        if Reverse :
            SetID = set(Intersected.keys())
            OldIds = set(self.Geoms.keys())
            NewIds = OldIds.difference(SetID)
            Intersected = {ID:self.Geoms[ID] for ID in NewIds}
            
        if Intersections : 
            Inters ={}
            for key,Geom2 in Intersected.items() : 
                Inter = Geom2.Intersection(Geom)
                if Inter is not None : 
                    Inters[key] = Inter
            Intersected = Inters
        
        Datas =[self[ID][0] for ID in Intersected.keys()]
        NewLayer = JFastLayer(File="")
        NewLayer.Fields = self.Fields
        NewLayer.FieldsTypes = self.FieldsTypes
        NewLayer.SpatialRef = self.SpatialRef
        NewLayer.GeomType=self.GeomType
        NewLayer.ID = self.ID
        NewLayer._Populate(Datas,Intersected,GeomIndex)
        return NewLayer
    
    def RoundCoordinates(self,Digits) : 
        for ID,Geom in self.Geoms.items() : 
            NewGeom = JGeom.RoundCoords(Geom,Digits)
            self.Geoms[ID] = NewGeom
        self.UpdateTree()
        
    def Flatten(self) : 
        for ID,Geom in self.Geoms.items() : 
            self.Geoms[ID].FlattenTo2D()
        self.UpdateTree()
        
    def NearestFeature(self,Geom,MaxDist=float("inf"),Jump=100) :
        """
        permet de retrouver la feature la plus proche de la geometrie fournie
        """
        Result = self.SpatialTree.NearestGeom(Geom,MaxDist=MaxDist,Jump=Jump)
        if Result!=None :        
            return self.NiceFeat(self[Result[0]]),Result[2]
        else :
            return None
        
    def GetFirsFeatureLike(self,Filters) : 
        """
        Permet de retourner ma premiere feature (ordre d'enregistrement) qui valide les caracteristiques demandee
        "( Field1 == 18 AND Field2 == 'titi' ) OR Field2 == 'toto'"
        """
        Table = self.AttrTable.AttributeFilter(Filters)
        Feat = Table._Array[0]
        return self.NiceFeat(self[Feat[self.ID]])
        


if __name__ == "__main__" :

    Datas =[(1,"titi"),(2,"toto"),(3,"tutu"),(4,"tata")]
    Fields =["NB","Name"]
    Types =[">i4",'|S10']
    #Types = ["INTEGER","TEXT"]
    Table = JDataTable(Fields,Types,"Name",Datas)
    #Table.AddField("Taille","REAL",0)
    #Table.DeleteField("Taille")
    
    #Table.AddField("Letter","TEXT","--")
    #Func = lambda x : x["Name"][1]
    #Table.CalculateField("Letter",Func)
    
    #NewTable = Table.AttributeFilter("NB >= 2")
    
    #for T in Table.Iterate() :
    #    print(T)

#Table = JDataTable(Fields,Types,"Name",Datas)
#Table.AddRow((10,"toutou"))
#Table.UpdateID()
#
#Table.AddField("Date","|S10","00-00-00")
#Table.RenameField("NB","Age")
#Table2 = Table.AttributeFilter(" (( Age > 2 ) OR ( Name == 'titi' )) AND ( Name != 'toutou' ) ")
#Table["toutou"]["Age"]=25
#Table.AddField("NameMaj","|S10","NONE")
#Table.CalculateField("NameMaj",lambda x : x["Name"].upper())
#
#File = "I:/Python/Temp/training_set.shp"
#Shp = ogr.Open(File)
#Layer = Shp.GetLayer(0)
#LayerDef = Layer.GetLayerDefn()
#FieldDef=LayerDef.GetFieldDefn(1)
#Layer = JFastLayer(File)
#Layer.Initialize(ID="OID",GeomIndex=True)
#
#B1 = Layer[0][1].Buffer(400)
#JGeom.DrawPolygone(B1,(0,0,1))
#Layer2 = Layer.SpatialFilter(B1,True)
#Layer.Draw(color="ro")
#Layer2.Draw(color="go")