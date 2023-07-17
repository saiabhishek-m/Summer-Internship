#**************************************************************************************************#
# Importing Necessary Modules
#**************************************************************************************************#
import numpy as np
from random import random
import csv

# ABAQUS
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

#**************************************************************************************************#
# Master Function Definition
#**************************************************************************************************#
def Master_Function(le, PR, AR, TR, Vf, r1, l, lc, tv, th, LRVE, BRVE, r2, x1, x2, y1, sheetsize, gridspace, 
                    rhop, Ep, vp, rhom, Em, vm, seedVal, initincval, maxincval, 
                    minincval, maxnumincval, timep, xdisp, ydisp, cpunum, gpunum, fail_p,
                      norm_fail_m, shear_fail_m):
    
    # Defining Module Functions
    def Part_Module(le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1, sheetsize, gridspace):
        mod = mdb.models[modelname]
        mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize, ).setPrimaryObject(option=STANDALONE)
        s = mod.sketches['__profile__']
        s.rectangle(point1=(0.0, 0.0), point2=(LRVE, BRVE))
        session.viewports['Viewport: 1'].view.fitView()
        mod.Part(name=partname, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
        mod.parts[partname].BaseShell(sketch=mod.sketches['__profile__'])
        mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize).unsetPrimaryObject()
        
        # Generating Partitions
        session.viewports['Viewport: 1'].setValues(displayedObject=mod.parts[partname])
        p = mod.parts[partname]
        f = p.faces
        RVEcent = (0.5*LRVE ,0.5*BRVE, 0.0)
        t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=RVEcent)
        # In the above line, the coordinates set as the origin are defined with respect to the origin which was used to when sketching the part, which is
        # the bottom left corner.

        s = mod.ConstrainedSketch(name='__profile__', sheetSize=sheetsize, gridSpacing=gridspace, transform=t)
        s.setPrimaryObject(option=SUPERIMPOSE)
        p = mod.parts[partname]
        p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

        # Note while specifying the necessary curves/shapes for partition sketches, give coordinates defined by the 'RVEcent' (translation of origin)
        s.Line(point1=(x1-RVEcent[0],0-RVEcent[1]),point2=(LRVE-x2-RVEcent[0],0-RVEcent[1]))                  #EF
        s.Line(point1=(LRVE-x2-RVEcent[0],0-RVEcent[1]),point2=(x1+l-RVEcent[0],0.5*le-RVEcent[1]))             #FR
        s.Line(point1=(x1+l-RVEcent[0],0.5*le-RVEcent[1]),point2=(x1+r1*l-RVEcent[0],0.5*lc-RVEcent[1]))      #RU
        s.Line(point1=(x1+r1*l-RVEcent[0],0.5*lc-RVEcent[1]),point2=(x1-RVEcent[0],0.5*le-RVEcent[1]))         #UQ
        s.Line(point1=(x1-RVEcent[0],0.5*le-RVEcent[1]),point2=(x1-RVEcent[0],0-RVEcent[1]))                   #QE

        s.Line(point1=(0-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]),point2=(r1*l-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]))            #MW
        s.Line(point1=(r1*l-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]),point2=(r1*l-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]))             #WZ
        s.Line(point1=(r1*l-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]),point2=(0-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]))              #ZP
        #s.Line(point1=(0-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]),point2=(0-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]))          #PM

        s.Line(point1=(r1*l+tv-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]))            #XN
        #s.Line(point1=(LRVE-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]))    #NO
        s.Line(point1=(LRVE-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]),point2=(r1*l+tv-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]))          #OY 
        s.Line(point1=(r1*l+tv-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]),point2=(r1*l+tv-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]))                   #YX

        s.Line(point1=(x1-RVEcent[0],BRVE-0.5*le-RVEcent[1]),point2=(x1+r1*l-RVEcent[0],BRVE-0.5*lc-RVEcent[1]))                   #TV
        s.Line(point1=(r1*l+x1-RVEcent[0],BRVE-0.5*lc-RVEcent[1]),point2=(x1+l-RVEcent[0],BRVE-0.5*le-RVEcent[1]))            #VS
        s.Line(point1=(x1+l-RVEcent[0],BRVE-0.5*le-RVEcent[1]),point2=(LRVE-x2-RVEcent[0],BRVE-RVEcent[1]))       #SG
        s.Line(point1=(LRVE-x2-RVEcent[0],BRVE-RVEcent[1]),point2=(x1-RVEcent[0],BRVE-RVEcent[1]))                 #GH
        s.Line(point1=(x1-RVEcent[0],BRVE-RVEcent[1]),point2=(x1-RVEcent[0],BRVE-0.5*le-RVEcent[1]))                 #HT

        s.Line(point1=(x1-RVEcent[0],0-RVEcent[1]),point2=(0-RVEcent[0],0-RVEcent[1]))             #EA
        #s.Line(point1=(0-RVEcent[0],0-RVEcent[1]),point2=(0-RVEcent[0],0.5*le+y1-RVEcent[1]))           #AI
        s.Line(point1=(0-RVEcent[0],0.5*le+y1-RVEcent[1]),point2=(x1-RVEcent[0],0.5*le-RVEcent[1]))            #IQ

        s.Line(point1=(x1+l-RVEcent[0],BRVE-0.5*le-RVEcent[1]),point2=(LRVE-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]))        #SK
        #s.Line(point1=(LRVE-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]),point2=(LRVE-RVEcent[0],BRVE-RVEcent[1]))            #KC
        s.Line(point1=(LRVE-RVEcent[0],BRVE-RVEcent[1]),point2=(LRVE-x2-RVEcent[0],BRVE-RVEcent[1]))             #CG

        s.Line(point1=(LRVE-x2-RVEcent[0],0-RVEcent[1]),point2=(LRVE-RVEcent[0],0-RVEcent[1]))          #FB
        #s.Line(point1=(LRVE-RVEcent[0],0-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*le+y1-RVEcent[1]))        #BJ 
        s.Line(point1=(LRVE-RVEcent[0],0.5*le+y1-RVEcent[1]),point2=(x1+l-RVEcent[0],0.5*le-RVEcent[1]))         #JR

        s.Line(point1=(x1-RVEcent[0],BRVE-0.5*le-RVEcent[1]),point2=(0-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]))     #TL
        #s.Line(point1=(0-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]),point2=(0-RVEcent[0],BRVE-RVEcent[1]))          #LD
        s.Line(point1=(0-RVEcent[0],BRVE-RVEcent[1]),point2=(x1-RVEcent[0],BRVE-RVEcent[1]))          #DH

        #s.Line(point1=(0-RVEcent[0],0.5*le+y1-RVEcent[1]),point2=(0-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]))                        #IM
        #s.Line(point1=(LRVE-RVEcent[0],0.5*le+y1-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*(BRVE-lc)-RVEcent[1]))         #JN
        #s.Line(point1=(LRVE-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]),point2=(LRVE-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]))       #KO
        #s.Line(point1=(0-RVEcent[0],BRVE-0.5*le-y1-RVEcent[1]),point2=(0-RVEcent[0],0.5*(BRVE+lc)-RVEcent[1]))         #LP

        s.Line(point1=(r1*l-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]),point2=(x1+r1*l-RVEcent[0],0.5*(BRVE-le)-y1-RVEcent[1]))              #Wa
        s.Line(point1=(x1+r1*l-RVEcent[0],0.5*(BRVE-le)-y1-RVEcent[1]),point2=(r1*l+tv-RVEcent[0],0.5*(BRVE-le)-RVEcent[1]))          #aX
        s.Line(point1=(r1*l-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]),point2=(x1+r1*l-RVEcent[0],0.5*(BRVE+le)+y1-RVEcent[1]))          #Zb
        s.Line(point1=(x1+r1*l-RVEcent[0],0.5*(BRVE+le)+y1-RVEcent[1]),point2=(r1*l+tv-RVEcent[0],0.5*(BRVE+le)-RVEcent[1]))              #bY

        s.Line(point1=(0-RVEcent[0],0-RVEcent[1]),point2=(0-RVEcent[0],BRVE-RVEcent[1]))                  #AD
        s.Line(point1=(LRVE-RVEcent[0],0-RVEcent[1]),point2=(LRVE-RVEcent[0],0-RVEcent[1]))             #BC   


        p = mod.parts[partname]

        # Parametrically defining faces for partition based on above created rectangles. The following are the centroids of those rectangles
        # These coordinates are with respect to the origin used for creating the sketch for the part
        pt1 = (x1+0.5*l,0.25*lc,0)
        f1 = f.findAt(pt1)
        pt2 = (0.5*r1*l,0.5*BRVE,0)                      #1-5 plateletes
        f2 = f.findAt(pt2)                                            #6 - matrix
        pt3 = (LRVE-0.5*r2*l,0.5*BRVE,0)
        f3 = f.findAt(pt3)
        pt4 = (x1+0.5*l,BRVE-0.25*lc,0)
        f4 = f.findAt(pt4) 
        pt5 = (0.5*x1,0.25*le,0)
        f5 = f.findAt(pt5)  
        pt6 = (LRVE-0.5*x2,0.25*le,0)
        f6 = f.findAt(pt6)
        pt7 = (x1+r1*l,0.5*(0.5*(BRVE+lc-le)-y1),0)
        f7 = f.findAt(pt7)
        pt8 = (x1+r1*l,0.5*BRVE,0)
        f8 = f.findAt(pt8)
        pt9 = (x1+r1*l,0.5*(1.5*BRVE+0.5*(le-lc)+y1),0)
        f9 = f.findAt(pt9)
        pt10 = (0.5*x1,BRVE-0.25*le,0)
        f10 = f.findAt(pt10)
        pt11 = (LRVE-0.5*x2,BRVE-0.25*le,0)
        f11 = f.findAt(pt11)


        pickedfaces = (f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, )
        p.PartitionFaceBySketch(faces=pickedfaces, sketch=s)
        s.unsetPrimaryObject()


    def Property_Module(rhop, Ep, vp, rhom, Em, vm, le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1):
        mod = mdb.models[modelname]
        # Platelet Material
        mod.Material(name=matp)
        mod.materials[matp].Density(table=((rhop, ), ))
        mod.materials[matp].Elastic(type=ISOTROPIC, table=((Ep, vp), ))
        # mod.materials[matp].Plastic(table=((240,0), ))

        # Matrix Material
        mod.Material(name=matm)
        mod.materials[matm].Density(table=((rhom, ), ))
        mod.materials[matm].Elastic(type=ISOTROPIC, table=((Em, vm), ))
        # matm_pt = ExcelDataToScriptData('TilAl4V_ELI_Stress_Strain.xlsx', 'Sheet1')         # Usage of the plastic table via 'ExcelDatatoScriptData' fucntion
        # matrix_plastic_table = ((790.0, 0.0), (791.23, 0.00022), (795.05, 0.00085), (798.87, 0.00165), (802.69, 0.00265), (806.51, 0.00391), (810.33, 0.0055), (814.15, 0.0075), (817.97, 0.01003), (821.79, 0.0132), (825.62, 0.01719), (829.44, 0.02221), (833.26, 0.02851), (837.08, 0.03642), (840.9, 0.04635), (844.72, 0.05879), (848.54, 0.07437), (852.36, 0.09387), (856.18, 0.11823))
        # mod.materials[matm].Plastic(table=matrix_plastic_table)

        # Platelet Section
        p = mod.parts[partname]
        psecname = 'Platelet_Section'

        mod.HomogeneousSolidSection(material=matp, name=psecname, thickness=100)

        p1 = (round(x1+0.5*l, 6), round(0.25*lc, 6), 0.0)
        p2 = (round(0.5*r1*l, 6), round(0.5*BRVE, 6), 0.0)
        p3 = (round(LRVE-0.25*r2*l, 6), round(0.5*BRVE, 6), 0.0)
        p4 = (round(x1+0.5*l, 6), round(BRVE-0.25*lc, 6), 0.0)


        p.Set(name='Platelet_Faces', faces=p.faces.findAt((p1, ), (p2, ), (p3, ), (p4, ), ))

        p.SectionAssignment(region=p.sets['Platelet_Faces'], sectionName=psecname, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        # Matrix Section
        p = mod.parts[partname]
        msecname = 'Matrix_Section'

        mod.HomogeneousSolidSection(material=matm, name=msecname, thickness=100)

        m5 = (round(0.5*x1, 6), round(0.25*le, 6), 0)
        m6 = (round(LRVE-0.5*x2, 6), round(0.25*le, 6), 0)
        m7 = (round(x1+r1*l, 6), round(0.5*(0.5*(BRVE+lc-le)-y1), 6), 0)
        m8 = (round(x1+r1*l, 6), round(0.5*BRVE, 6), 0)
        m9 = (round(x1+r1*l, 6), round(0.5*(1.5*BRVE+0.5*(le-lc)+y1), 6), 0)
        m10 = (round(0.5*x1, 6), round(BRVE-0.25*le, 6), 0)
        m11 = (round(LRVE-0.5*x2, 6), round(BRVE-0.25*le, 6), 0)

        p.Set(name='Matrix_Faces', faces=p.faces.findAt((m5, ), (m6, ), (m7, ), (m8, ), (m9, ), (m10, ), (m11, ), ))
        p.Set(name='VI_FACES', faces=p.faces.findAt((m5, ), (m6, ), (m8, ), (m10, ), (m11, ), ))
        p.Set(name='HI_FACES', faces=p.faces.findAt((m7, ), (m9, ), ))
        p.SectionAssignment(region=p.sets['Matrix_Faces'], sectionName=msecname, offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        #mod.Tie(name='Matrix_Section', master=p.sets['VI_FACES'], slave=p.sets['HI_FACES'], positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)

    def Assembly_Module():
        mod = mdb.models[modelname]
        a = mod.rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mod.parts[partname]
        a.Instance(name=instancename, part=p, dependent=OFF)

    def Mesh_Module(seedVal, le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1):
        mod = mdb.models[modelname]
        a = mod.rootAssembly
        c = a.instances[instancename]
        session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
        session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(meshTechnique=ON)

        # Creating Face Sets
        p1 = (round(x1+0.5*l, 6), round(0.25*lc, 6), 0.0)
        p2 = (round(0.5*r1*l, 6), round(0.5*BRVE, 6), 0.0)
        p3 = (round(LRVE-0.25*r2*l, 6), round(0.5*BRVE, 6), 0.0)
        p4 = (round(x1+0.5*l, 6), round(BRVE-0.25*lc, 6), 0.0)

        m5 = (round(0.5*x1, 6), round(0.25*le, 6), 0)
        m6 = (round(LRVE-0.5*x2, 6), round(0.25*le, 6), 0)
        m7 = (round(x1+r1*l, 6), round(0.5*(0.5*(BRVE+lc-le)-y1), 6), 0)
        m8 = (round(x1+r1*l, 6), round(0.5*BRVE, 6), 0)
        m9 = (round(x1+r1*l, 6), round(0.5*(1.5*BRVE+0.5*(le-lc)+y1), 6), 0)
        m10 = (round(0.5*x1, 6), round(BRVE-0.25*le, 6), 0)
        m11 = (round(LRVE-0.5*x2, 6), round(BRVE-0.25*le, 6), 0)

        a.seedPartInstance(regions=(c, ), size=seedVal, deviationFactor=0.0, minSizeFactor=0.1)

        # Setting Element Type and Assigning Mesh Controls
        elemType1 = ElemType(elemCode=CPE4, elemLibrary=STANDARD)
        elemType2 = ElemType(elemCode=CPE3, elemLibrary=STANDARD)

        a.setElementType(elemTypes=(elemType1, elemType2), regions=(c.faces.findAt((p1, ), (p2, ), (p3, ), (p4, ), (m5, ), (m6, ), (m7, ), (m8, ), (m9, ), (m10, ), (m11, ), ), ))
        a.setMeshControls(elemShape=QUAD, regions=c.faces.findAt((p1, ), (p2, ), (p3, ), (p4, ), (m5, ), (m6, ), (m7, ), (m8, ), (m9, ), (m10, ), (m11, ), ), technique=STRUCTURED)
        a.generateMesh(regions=(c, ))

    def Step_Module(initincval, maxincval, minincval, maxnumincval, timep):
        mod = mdb.models[modelname]
        mod.StaticStep(initialInc=initincval, maxInc=maxincval, minInc=minincval, maxNumInc=maxnumincval, name=stepname, previous='Initial', timePeriod=timep)

    def Eqn_Constraints(le, lc, l, tv, th, r1, r2, x1, x2, y1, LRVE, BRVE):
        mod = mdb.models[modelname]
        a = mod.rootAssembly
        p = a.instances[instancename]

        L1 = p.edges.findAt((0.0, round(0.25*le, 6), 0.0))
        L2 = p.edges.findAt((0.0, round(0.5*(0.5*(BRVE+le-lc) + y1), 6), 0.0))
        L3 = p.edges.findAt((0.0, round(0.5*BRVE, 6), 0.0))
        L4 = p.edges.findAt((0.0, round(0.5*(BRVE-0.5*le+0.5*(BRVE+lc)-y1), 6), 0.0))
        L5 = p.edges.findAt((0.0, round(BRVE-0.25*le, 6), 0.0))

        R1 = p.edges.findAt((round(LRVE, 6), round(0.25*le, 6), 0.0))
        R2 = p.edges.findAt((round(LRVE, 6), round(0.5*(0.5*(BRVE+le-lc) + y1), 6), 0.0))
        R3 = p.edges.findAt((round(LRVE, 6), round(0.5*BRVE, 6), 0.0))
        R4 = p.edges.findAt((round(LRVE, 6), round(0.5*(BRVE-0.5*le+0.5*(BRVE+lc)-y1), 6), 0.0))
        R5 = p.edges.findAt((round(LRVE, 6), round(BRVE-0.25*le, 6), 0.0))

        U1 = p.edges.findAt((round(0.5*x1, 6), round(BRVE, 6), 0.0))
        U2 = p.edges.findAt((round(0.5*l+x1, 6), round(BRVE, 6), 0.0))
        U3 = p.edges.findAt((round(LRVE-0.5*x2, 6), round(BRVE, 6), 0.0))

        D1 = p.edges.findAt((round(0.5*x1, 6), 0.0, 0.0))
        D2 = p.edges.findAt((round(x1+0.5*l, 6), 0.0, 0.0))
        D3 = p.edges.findAt((round(LRVE-0.5*x2, 6), 0.0, 0.0))


        q11 = U1.index
        q12 = U2.index
        q13 = U3.index

        q21 = D1.index
        q22 = D2.index
        q23 = D3.index

        q31 = R1.index
        q32 = R2.index
        q33 = R3.index
        q34 = R4.index
        q35 = R5.index

        q41 = L1.index
        q42 = L2.index
        q43 = L3.index
        q44 = L4.index
        q45 = L5.index

        EdUe1 = p.edges[q11:q11+1]
        EdUe2 = p.edges[q12:q12+1]
        EdUe3 = p.edges[q13:q13+1]  


        EdDe1 = p.edges[q21:q21+1]
        EdDe2 = p.edges[q22:q22+1]
        EdDe3 = p.edges[q23:q23+1]  


        EdRe1 = p.edges[q31:q31+1]
        EdRe2 = p.edges[q32:q32+1]
        EdRe3 = p.edges[q33:q33+1]
        EdRe4 = p.edges[q34:q34+1]
        EdRe5 = p.edges[q35:q35+1]


        EdLe1 = p.edges[q41:q41+1]
        EdLe2 = p.edges[q42:q42+1]
        EdLe3 = p.edges[q43:q43+1]
        EdLe4 = p.edges[q44:q44+1]
        EdLe5 = p.edges[q45:q45+1]


        combinedEdgesL = EdLe1 + EdLe2 + EdLe3 + EdLe4 + EdLe5
        a.Set(edges=combinedEdgesL, name='Combined_LeftEdges')
        Leftnodes = a.sets['Combined_LeftEdges'].nodes


        combinedEdgesR = EdRe1 + EdRe2 + EdRe3 + EdRe4 + EdRe5
        a.Set(edges=combinedEdgesR, name='Combined_RightEdges')
        Rightnodes = a.sets['Combined_RightEdges'].nodes


        combinedEdgesU = EdUe1 + EdUe2 + EdUe3 
        a.Set(edges=combinedEdgesU, name='Combined_UpEdges')
        Upnodes = a.sets['Combined_UpEdges'].nodes


        combinedEdgesD = EdDe1 + EdDe2 + EdDe3 
        a.Set(edges=combinedEdgesD, name='Combined_DownEdges')
        Downnodes = a.sets['Combined_DownEdges'].nodes

        Upcoord, Downcoord, Leftcoord, Rightcoord = [], [], [], []

        for node in Upnodes:
            Upcoord = Upcoord + [[node.coordinates[0], node.coordinates[1], node.label]]
        for node in Downnodes:
            Downcoord = Downcoord + [[node.coordinates[0], node.coordinates[1], node.label]]

        for node in Leftnodes:
            Leftcoord = Leftcoord + [[node.coordinates[0], node.coordinates[1], node.label]]
        for node in Rightnodes:
            Rightcoord = Rightcoord + [[node.coordinates[0], node.coordinates[1], node.label]]

        for i in range(0,len(Leftcoord)):
            if (Leftcoord[i][0]<10E-6):
                Leftcoord[i][0] = 0
            else:
                continue


        Upcoord.sort()
        Downcoord.sort()
        Leftcoord.sort()
        Rightcoord.sort()

        NodeTol = seedVal/200

        Num = len(Upcoord)
        for i in range(0, Num):
            if (abs(Upcoord[i][0]-Downcoord[i][0])<NodeTol):
                Nlabel = Upcoord[i][2]
                a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='UpNode-'+str(i))
                Nlabel = Downcoord[i][2]
                a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='DownNode-'+str(i))


        Num = len(Leftcoord)
        for i in range(0, Num):
            if (abs(Leftcoord[i][1]-Rightcoord[i][1])<NodeTol):
                Nlabel = Rightcoord[i][2]
                a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='RightNode-'+str(i))
                Nlabel = Leftcoord[i][2]
                a.Set(nodes=p.nodes[Nlabel-1:Nlabel], name='LeftNode-'+str(i))

        for i in range(1, len(Rightcoord)-1):
            mod.Equation(name='Eqn-LR-X-'+str(i), terms=((-1.0, 'LeftNode-'+str(i), 1), (1.0, 'RightNode-'+str(i), 1), (-1.0, 'RightNode-0', 1)))
        for i in range(1, len(Rightcoord)-1):
            mod.Equation(name='Eqn-LR-Y-'+str(i), terms=((1.0, 'LeftNode-'+str(i), 2), (-1.0, 'RightNode-'+str(i), 2)))

        for i in range(1, len(Upcoord)-1):
            mod.Equation(name='Eqn-UD-Y-'+str(i), terms=((-1.0, 'DownNode-'+str(i), 2), (1.0, 'UpNode-'+str(i), 2), (-1.0, 'UpNode-0', 2)))
        for i in range(1, len(Upcoord)-1):
            mod.Equation(name='Eqn-UD-X-'+str(i), terms=((1.0, 'DownNode-'+str(i), 1), (-1.0, 'UpNode-'+str(i), 1)))


        mod.Equation(name='Eqn-TR-X', terms=((1.0, 'UpNode-'+str(len(Upcoord)-1), 1), (-1.0, 'DownNode-'+str(len(Upcoord)-1), 1)))  # along the X direction
        mod.Equation(name='Eqn-TR-Y', terms=((1.0, 'UpNode-'+str(len(Upcoord)-1), 2), (-1.0, 'UpNode-0', 2)))   # along Y direction


        return [Upcoord, Downcoord, Leftcoord]

    def Boundary_Conditions(xdisp, ydisp, LRVE, BRVE):
        mod = mdb.models[modelname]
        a = mod.rootAssembly
        v = a.instances[instancename].vertices
        ver = v.findAt((0.0, 0.0, 0.0))
        q = ver.index
        FixVer1 = v[q:q+1]
        region = a.Set(vertices=FixVer1, name='Set-Fix-1')
        mod.EncastreBC(name='Fix1', createStepName='Initial', region=region, localCsys=None)

        # ver = v.findAt((0.0, round(BRVE, 6), 0.0))
        # q = ver.index
        # FixVer2 = v[q:q+1]
        # #region = a.Set(vertices=FixVer2, name='Set-Fix-2')
        # mod.EncastreBC(name='Fix-2', createStepName='Initial', region=region, localCsys=None)

        #Fix Upper Left Node along X direction
        ver = v.findAt((0.0, round(BRVE, 6), 0.0))
        q = ver.index
        MoveVer = v[q:q+1]
        region = a.Set(vertices=MoveVer, name='Y') 
        mod.DisplacementBC(name='Y', createStepName=stepname, region=region, u1=0.0, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)

        # Fix Bottom Right Node along Y direction
        ver = v.findAt((round(LRVE, 6), 0.0, 0.0))
        q = ver.index
        MoveVer = v[q:q+1]
        region = a.Set(vertices=MoveVer, name='X')
        mod.DisplacementBC(name='X', createStepName=stepname, region=region, u1=xdisp, u2=0.0, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None) # change the u1 value to 'xdisp' for cases other than validation

    def Job_Module(cpunum, gpunum):
        mdb.Job(name=jobname, model=modelname, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, 
        numCpus=cpunum, numGPUs=gpunum, numDomains=cpunum)
        mdb.jobs[jobname].submit(consistencyChecking=OFF)
        mdb.jobs[jobname].waitForCompletion()

    def Post_Processing(BRVE, fail_p, norm_fail_m, shear_fail_m, Upcoord, Downcoord):
        # Accessing ODB
        odb = session.openOdb(name=odbpath)
        o2 = odb.rootAssembly
        o3 = o2.instances[instancename]

        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT, 'S11'))
        session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(deformationScaling=UNIFORM, uniformScaleFactor=0)

        S11_Plt_val = []
        S11_VI_val = []
        S12_HI_val = []

        S11_plt = odb.steps[stepname].frames[-1].fieldOutputs['S'].getSubset(region=o3.elementSets['PLATELET_FACES'])
        for s in S11_plt.values:
            S11_Plt_val.append(s.data[0])

        S11_VI = odb.steps[stepname].frames[-1].fieldOutputs['S'].getSubset(region=o3.elementSets['VI_FACES'])
        for s in S11_VI.values:
            S11_VI_val.append(s.data[0])

        S12_HI = odb.steps[stepname].frames[-1].fieldOutputs['S'].getSubset(region=o3.elementSets['HI_FACES'])
        for s in S12_HI.values:
            S12_HI_val.append(s.data[2])
        
        def average(list):
            return sum(list)/len(list)
        
        s11_plt_avg = average(S11_Plt_val)
        s11_vi_avg = average(S11_VI_val)
        S12_hi_avg = average(S12_HI_val)

        R_p = s11_plt_avg/fail_p
        R_v = s11_vi_avg/norm_fail_m
        R_h = S12_hi_avg/shear_fail_m

        # print(round(R_p, 6), round(R_v, 6), round(R_h, 6))

        session.Path(name='Middle', type=NODE_LIST, expression=(('RVE', (Downcoord[int(np.ceil(len(Downcoord)*0.5))][2], Upcoord[int(np.ceil(len(Upcoord)*0.5))][2], )), ))
        pth = session.paths['Middle']
        session.XYDataFromPath(name='XYData-Middle', path=pth, includeIntersections=True, projectOntoMesh=False, pathStyle=PATH_POINTS, numIntervals=10, projectionTolerance=0, shape=UNDEFORMED, labelType=TRUE_DISTANCE_Y)

        x0 = session.xyDataObjects['XYData-Middle']
        total=0
        for i in range(0,(len(x0)-1)):
            a1=list(x0[i])
            a2=list(x0[i+1])
            c=0.5*(a2[0]-a1[0])*(a2[1]+a1[1])
            total=c+total
        S11_avg= total/BRVE
        # print(round(S11_avg, 6))

        session.odbs[odbpath].close()
        mdb.close()

        return round(R_p, 6), round(R_v, 6), round(R_h, 6), round(S11_avg, 6)
    
    def Del_Module(Upcoord, Leftcoord):
        mod = mdb.models[modelname]
        a1 = mod.rootAssembly

        # Deleting sets containing individual nodes
        for i in range(len(Upcoord)):
            del a1.sets['UpNode-'+str(i)]
            del a1.sets['DownNode-'+str(i)]

        for i in range(len(Leftcoord)):
            del a1.sets['LeftNode-'+str(i)]
            del a1.sets['RightNode-'+str(i)]
        
        # Deleting all equation constraints
        for i in range(1, len(Upcoord)-1):
            del mod.constraints['Eqn-UD-X-'+str(i)]
            del mod.constraints['Eqn-UD-Y-'+str(i)]

        for i in range(1, len(Leftcoord)-1):
            del mod.constraints['Eqn-LR-X-'+str(i)]
            del mod.constraints['Eqn-LR-Y-'+str(i)]
            
        del mod.constraints['Eqn-TR-X']
        del mod.constraints['Eqn-TR-Y']

    # Calling Module Functions and storing outputs
    Part_Module(le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1, sheetsize, gridspace)
    Property_Module(rhop, Ep, vp, rhom, Em, vm, le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1)
    Assembly_Module()
    Mesh_Module(seedVal, le, l, lc, tv, th, LRVE, BRVE, r1, r2, x1, x2, y1)
    Step_Module(initincval, maxincval, minincval, maxnumincval, timep)
    Eqn_Constraints(le, lc, l, tv, th, r1, r2, x1, x2, y1, LRVE, BRVE)
    Coord = Eqn_Constraints(le, lc, l, tv, th, r1, r2, x1, x2, y1, LRVE, BRVE)
    Upcoord = Coord[0]
    Downcoord = Coord[1]
    Leftcoord = Coord[2]
    Boundary_Conditions(xdisp, ydisp, LRVE, BRVE)
    Job_Module(cpunum, gpunum)
    Del_Module(Upcoord, Leftcoord)
    Post_Processing(BRVE, fail_p, norm_fail_m, shear_fail_m, Upcoord, Downcoord)
    res1 = Post_Processing(BRVE, fail_p, norm_fail_m, shear_fail_m, Upcoord, Downcoord)

    return [le, PR, AR, TR, Vf, r1, res1[0], res1[1], res1[2], res1[3]]

#**************************************************************************************************#
# Defining Random Hypercube Sampling and Parameter Ranges using it
#**************************************************************************************************#
def Random_Latin_Hypercube(num_samples, param_ranges):
    # num_samples = number of samples (an integer value only)
    # param_ranges = a list of parameter ranges, ex: for a two parameter case
    # param_ranges = [(1, 2), (3, 4)]
    
    num_var = len(param_ranges) # number of variables
    data = np.zeros((num_samples, num_var)) # np.zeros((num_rows, num_columns)). output is a matrix

    # Generating random values for each parameter within its range
    for i in range(num_var):    # iterating through columns
        for j in range(num_samples):    # iterating through rows = number of samples = number of stratifications
            data[j, i] = param_ranges[i][0] + (param_ranges[i][1]-param_ranges[i][0])*(j+random())/num_samples

    # shuffling the rows to create a more random Hypercube
    for i in range(num_var):    # iterating through the number of variables, i.e. columns
        order = np.arange(num_samples)  # arranging the the sample number i.e. the row number in an increasing order starting with zero
        np.random.shuffle(order)
        data[:, i] = data[order, i]
    return data

# le_range = (5E-3, 10E-3)
# PR_range = (0.2, 0.9)
# AR_range = (2, 50)
# TR_range = (-2, 2)
# Vf_range = (0.4, 0.8)
# r1_range = (0.1, 0.9)

le_range = (5E-3, 15E-3)
PR_range = (0.4, 0.9)
AR_range = (5, 40)
TR_range = (-2, 2)
Vf_range = (0.4, 0.8)
r1_range = (0.3, 0.8)


param_ranges = [le_range, PR_range, AR_range, TR_range, Vf_range, r1_range]
#num_samples = 8000
num_samples = 10
InpData = Random_Latin_Hypercube(num_samples, param_ranges)
le = []
PR = []
AR = []
TR = []
Vf = []
r1 = []
for i in range(num_samples):
    le.append(InpData[i, 0])
    PR.append(InpData[i, 1])
    AR.append(InpData[i, 2])
    TR.append(InpData[i, 3])
    Vf.append(InpData[i, 4])
    r1.append(InpData[i, 5])

#**************************************************************************************************#
# Defining software specific parameters
#**************************************************************************************************#
modelname = 'Model-1'
partname = 'Part'
sheetsize = 5000
gridspace = 0.5
matp = 'Platelet_Material'
matm = 'Matrix_Material'
instancename = 'RVE'
stepname = 'Step-1'
initincval = 1
maxincval = initincval
minincval = 0.1
timep = 1
maxnumincval = 1
jobname = 'Master'
cpunum = 6
gpunum = 1
odbpath='Master.odb'

#**************************************************************************************************#
# Defining Material Properties | Units: ton, mm, second, N, MPa, N-mm, ton/mm^3
#**************************************************************************************************#
fail_p = 260    # UTS of paltelet in MPa form MatWeb
norm_fail_m = 80    # UTS of matrix in MPa from MatWeb
shear_fail_m = 60   # Shear Strength of matrix in MPa from MatWeb

# Platelet - 99.5% Alumina
Ep = 375000             # Young's Modulus of Platelet in MPa
vp = 0.22               # Poisson's Ratio of Platelet
rhop = 3.89E-9           # Density of Platelet material in tonne/mm3

# Matrix - PMMA
Em = 2500             # Young's Modulus of Matrix in MPa
rhom = 1.185E-9         # Density of Matrix material in tonne/mm3
vm = 0.375              # Poisson's Ratio (0.35 - 0.4) (explicit reference - https://www.mit.edu/~6.777/matprops/pmma.htm)

#**************************************************************************************************#
# Calling Master_Function() for each set of variables
#**************************************************************************************************#
l = []  
b = []
lc = []        # length of platelet
tv = []         # thickness of horizontally aligned matrix layer
th = []         # thickness of vertically aligned matrix layer
LRVE = []       # Length of RVE
BRVE = []       # Breadth of RVE
r2 = []   
x1 = [] 
x2 = []
y1 = []     # Complementary Offset percentage
xdisp = []
ydisp = []
seedVal = []    # Seed value for meshing
Result = [0]*10     # Container for all results

with open('MainRVE_Samples-1.csv', 'w') as file:
    writer = csv.writer(file)
    for i in range(0, len(le)):
        a1=le[i]*AR[i]
        l.append(a1)

        a2 = PR[i]*le[i]
        lc.append(a2)

        b1 = 0.5*(le[i]+lc[i])
        b.append(b1)

        c1 = np.e**TR[i]
        c2 = l[i]+b[i]*np.e**TR[i]
        c3 = l[i]*b[i]*(1-1/Vf[i])
        a3 = (-c2+np.sqrt(c2**2-4*c1*c3))/(2*c1)# Thickness of Horizontally aligned matrix in mm
        th.append(a3)

        a4 = th[i]*np.e**TR[i]    # Thickness of Vertically aligned matrix in mm
        tv.append(a4)

        a5 = l[i]+tv[i]
        LRVE.append(a5)

        a7 = 1-r1[i]
        r2.append(a7)

        a8 = 0.001*LRVE[i]
        xdisp.append(a8)

        a11 = tv[i]*(1-(r1[i]*l[i])/(LRVE[i]-tv[i]))
        x2.append(a11)

        a12 = tv[i]-x2[i]
        x1.append(a12)

        a13 = 0.5*(le[i]-lc[i])*tv[i]/(LRVE[i]-tv[i])
        y1.append(a13)

        a6 = le[i]+ lc[i]+ 2*th[i] + 2*y1[i]
        BRVE.append(a6)

        a9 = 0.001*BRVE[i]
        ydisp.append(a9)

        a10 = BRVE[i]/20
        seedVal.append(a10)

        R = Master_Function(le[i], PR[i], AR[i], TR[i], Vf[i], r1[i], l[i], lc[i], tv[i], th[i], LRVE[i], BRVE[i], r2[i], x1[i], x2[i], y1[i], sheetsize, gridspace, 
                        rhop, Ep, vp, rhom, Em, vm, seedVal[i], initincval, maxincval, 
                        minincval, maxnumincval, timep, xdisp[i], ydisp[i], cpunum, gpunum, fail_p,
                        norm_fail_m, shear_fail_m)

        Result = [R[j] for j in range(10)]
        writer.writerow(Result)
#**************************************************************************************************#
# Exporting Result to CSV file as DataSet for Machine Learning
#**************************************************************************************************#