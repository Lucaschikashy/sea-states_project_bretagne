#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.4.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/home/H22552/ENPC/SPAR/0-MESH')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Cylinder_1 = geompy.MakeCylinderRH(5, 200)
Box_1 = geompy.MakeBoxDXDYDZ(400, 400, 400)
Translation_1 = geompy.MakeTranslation(Box_1, -200, 0, -10)
Cut_1 = geompy.MakeCutList(Cylinder_1, [Translation_1], True)
Translation_2 = geompy.MakeTranslation(Cut_1, 0, 0, -200)
Group_1 = geompy.CreateGroup(Translation_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_1, [3, 17])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Translation_1, 'Translation_1' )
geompy.addToStudy( Cut_1, 'Cut_1' )
geompy.addToStudy( Translation_2, 'Translation_2' )
geompy.addToStudyInFather( Translation_2, Group_1, 'Group_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Group_1)
Regular_1D = Mesh_1.Segment()
Local_Length_1 = Regular_1D.LocalLength(1,None,1e-07)
NETGEN_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_2D)
Length_From_Edges_1 = NETGEN_2D.LengthFromEdges()
Quadrangle_Preference_1 = NETGEN_2D.SetQuadAllowed()
Local_Length_1.SetLength( 1.3 )
Local_Length_1.SetPrecision( 1e-07 )
isDone = Mesh_1.Compute()
try:
  Mesh_1.ExportDAT( r'/home/H22552/ENPC/SPAR/0-MESH/Spar.dat' )
  pass
except:
  print('ExportDAT() failed. Invalid file name?')


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_2D.GetAlgorithm(), 'NETGEN 2D')
smesh.SetName(Quadrangle_Preference_1, 'Quadrangle Preference_1')
smesh.SetName(Local_Length_1, 'Local Length_1')
smesh.SetName(Length_From_Edges_1, 'Length From Edges_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
