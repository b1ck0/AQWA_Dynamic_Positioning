# -*- coding: utf-8 -*-
from AqwaServerMgr import *
from math import *

####################################################################################################################################################

def FSTYAW_1(yaw, Rx, Ry, sign): #transofrmation of axis sing = 1 --> FRA to LSA; sign == -1 --> LSA to FRA
    roll = 0.0
    pitch = 0.0
    
    c = cos(yaw)
    s = sin(yaw)
    
    if sign == 1:
        roll = Rx*c + Ry*s
        pitch = -Rx*s + Ry*c
    else:
        roll = Rx*c - Ry*s
        pitch = Rx*s + Ry*c
        
    return roll, pitch 

def sign(x):
    if x > 0:
        return 1.0
    elif x < 0:
        return -1.0
    elif x == 0:
        return 0.0
    else:
        return x

####################################################################################################################################################

class Globs:
    RHO = 1025.0 # [kg/m^3] Water Density
    T = 7.50 # [m] Draft at Midship
    B = 32.0 # [m] Breadth of the Ship
    L = 180.0 # [m] Length of the vessel
    Displ = 22102028.0 # [kg] Vessel Displacement
    Vol = Displ/RHO # [m^3] Vessel Displaced Volume
    Cb = Vol/(L*B*T) # [-] Block Coefficient of the Vessel
    S = L*(1.7*T+Cb*B) # [m^2] Wetted Area of the Vessel
    nu = 0.00001307 # [m^2/s] Kinematic Viscosity of the Water

    StructIdx = 3 # Here comes the index of the tug

    Integrals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    IniPos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	
    Fx = 50.0 * 9806.65 #Metric-Tons to Newtons
    Fy = 50.0 * 9806.65 #Metric-Tons to Newtons

    Kpx = Fx/10.0 
    Kpy = Fy/10.0 

    Kdx = Fx/50.0 
    Kdy = Fy/50.0

    Kix = Fx*0.0
    Kiy = Fy*0.0

    # DIP Coeffs : 27 coeffs
    D   = [
        [-Kdx,0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, -Kdy, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0]
        ]
    
    I   = [
        [-Kix,0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, -Kiy, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0]
        ]
    
    P   = [
        [-Kpx,0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, -Kpy, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0,0.0]
        ]
# End of Global persisting storage

def UF1(Analysis,Mode,Stage,Time,TimeStep,Pos,Vel):

    # We consider only one structure, so Structure index is always going to be '0'
    StructIdx = Globs.StructIdx
    # Check that input file is the one expected without .DAT extension
    Error   = 0
    # We create an empty container for AddMass and Force
    AddMass = BlankAddedMass(Analysis.NOfStruct)
    Force   = BlankForce(Analysis.NOfStruct)

    # Active D.O.F are X,Y,RZ
    DOFs = [0,1]

    # Initialise at t=0
    if (Time == 0.0):
        for dof in range(6):
            Globs.IniPos[dof] = Pos[StructIdx][dof]

    # Integrate Positions from relative position at t=0
    ActivePos = [0.0 ,0.0, 0.0, 0.0, 0.0, 0.0]
    for dof in range(6):
        ActivePos[dof] = Pos[StructIdx][dof] - Globs.IniPos[dof]
        Globs.Integrals[dof] += ActivePos[dof]
    
    # Now calculate response from coeffsi 
    for dof_Force in DOFs:
        for dof_Pos in DOFs:
            Force[StructIdx][dof_Force] += (
                Globs.D[dof_Force][dof_Pos] * Vel[StructIdx][dof_Pos]
                + Globs.I[dof_Force][dof_Pos] * Globs.Integrals[dof_Pos]
                + Globs.P[dof_Force][dof_Pos] * ActivePos[dof_Pos]
                )
    
    #making sure that the tug doesn't deliver more thrust that it really can        
    if abs(Force[StructIdx][0]) - Globs.Fx > 1.0:
         Force[StructIdx][0] = sign(Force[StructIdx][0])*Globs.Fx  

    if abs(Force[StructIdx][1]) - Globs.Fy > 1.0:
         Force[StructIdx][1] = sign(Force[StructIdx][1])*Globs.Fy   
      
    # Now return the results

    return Force,AddMass,Error

########################################################################################################

def UF4(Analysis,Mode,Stage,Time,TimeStep,Pos,Vel):

    Error   = 0
    AddMass = BlankAddedMass(Analysis.NOfStruct)
    Force = BlankForce(Analysis.NOfStruct)
    
    L = Globs.L
    rho = Globs.RHO
    nu = Globs.nu
    S = Globs.S

    VSurge = 0.0
    VSway = 0.0
    
    VSurge, VSway = FSTYAW_1(Pos[0][5],Vel[0][0],Vel[0][1],1)

    if abs(VSurge) > 0.0:
        Rn = abs(VSurge)*L/nu # Reynolds Number for the Ship
        A = (log10(Rn)-2)**2
        FViscous = 0.5*rho*(VSurge**2)*S*(0.075/A) #Viscous Resistance Force of the Ship
        Force[0][0], Force[0][1] = FSTYAW_1(Pos[0][5],-FViscous,0.0,-1)

    return Force,AddMass,Error

########################################################################################################

def UFTotal(Analysis,Mode,Stage,Time,TimeStep,Pos,Vel):
    
    Force1, AddMass1, Error = UF1(Analysis,Mode,Stage,Time,TimeStep,Pos,Vel) #Forces from Dynamic Positioning System
    
    return Force1, AddMass1,Error

Server = AqwaUserForceServer()

Server.Run(UFTotal)