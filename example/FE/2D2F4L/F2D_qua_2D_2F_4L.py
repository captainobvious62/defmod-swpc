#!/usr/bin/python
import numpy as np
import os, sys, netCDF4
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
import pandas as pd
################################################################################
# 2D Model Preprocessing Script
# Intended for Defmod
# 
# Last revision 24 Feb 2017, RLW
#
#plt.style.use('fivethirtyeight')
#plt.style.use('ggplot')
fin = sys.argv[1]
#fin = "2D_2F_4L.exo"
nc = netCDF4.Dataset(fin)

#
#
# 2D Horizontal Plane, injection into a channel
# Slip Weakening Frictional
################################################################################
# Sim Parameters
#===============================================================================
#mpiexec -n 4 defmod -ss 1 -f 2D_2F_4L.inp -pc_type lu -pc_factor_mat_solver_package mumps
# Define model type
explicit = False 
implicit = False 
fault = True 
poro = True
visc = False
dim = 2 
NZ_K_PR = 9                  # K matrix nonzeros per row

# Quasi-Static Time Stepping (implicit)
dt_hr=6                      # dt in hr
t = 3600.*24*30              # Overall sim time, sec
dt = 3600.*dt_hr             # Implicit time step
nviz = 1                     # Snapshot Sample Interval
dsp = 1                      # Output type: 1 absolute, 0 incremental

# Dynamic Time Stepping (explicit)
t_dyn = 3.000                # Explicit t
dt_dyn= 0.00025              # Explicit dt - Timestep length, post dynamic switchover
nviz_dyn=20                  # Explicit Snapshot Sample Interval
t_lim = 3.                   # Quake duration
dsp_hyb=0                    # Explicit output type: 1 absolute, 0 incremental

# Yes/No (1/0) Switches
dsp_str=1                    # Output Stress
bod_frc=1                    # Body forces (gravity)
rsf=0                        # Rate and State Friction
init=1                       # Initial Pressure - KEEP THIS FLAG ON

# For Rate and State Friction (rsf==1)
v_bg=1E-12                   # initial slip rate for rsf, m/s

# Simulation Options
hyb=1                        # Hybrid Solver
nviz_wave=10                 # Waveform Sample Interval
nviz_slip=100                # Slip Rate Sample Interval

# Rayleigh Damping Constants
alpha= 0.
beta = 0.00125
rfrac = 0

################################################################################
# Nodeset / Sideset Info
#===============================================================================

# NS/SS with fault information, C ordering
faultID = np.array([1])
# Non-Conforming Node NS/SS
# Just set something valid and ignore this for now
NCF_node_list = np.array([6,8])
#===============================================================================
# Boundary Conditions

# Fixed
# Fixed NS/SS
# Set at -999 if no fixed value

# Y direction
fixedID_Y = np.array([2]) 
#fixedID_Y = np.array([-999])

# X direction
fixedID_X = np.array([3,5]) 
#fixedID_X = np.array([-999])

# Traction NS/SS
#tractionID = np.array([3, 5])  - 1

# No Flow (no pressure) - Actually think this is set drained condition
fixedID_Poro = np.array([-999])





# If Body Forces are enabled:

if bod_frc == 1:
# traction bc 
# Positive up, east
#            frcx(Pa)        frcy(Pa)     poro    start   end
    trac_val_F = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # Floor
    trac_val_E = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # East 
    trac_val_S = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # Surface
    trac_val_W = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # West

elif bod_frc == 0:
    trac_val_F = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # Floor
    trac_val_E = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # East 
    trac_val_S = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # Surface
    trac_val_W = [ 0.0 ,         0.0,      0.0,    0.0,    t]    # West




# ID numbers for edge BCs
TracBCFloor     =   np.array([2])   # Bottom
TracBCEast      =   np.array([3])   # Right
TracBCSurface   =   np.array([4])   # Top
TracBCWest      =   np.array([5])   # Left



## ID numbers for edge BCs
#BCFloor     =   np.array([4])
#BCWest      =   np.array([3]) 
#BCEast      =   np.array([5]) 
#BCSurface   =   np.array([2])

## Traction, pascals, (x, y)
#trac_val = [50E6, -65E6]

# Values after faults
first_bnd = 1 
last_bnd = 5


################################################################################
## Physical model parameters
################################################################################

# Number of ID Values
numIDs = 10
init = 1
ID_Range = np.linspace(1,numIDs,numIDs)
dp_val = 0 #source/sink 
num_mats = 4                # Set for diagram plotting

#                       MR          FM           WB      EB         WFZ1       EFZ1      WFZ2      EFZ2      WFZ3      EFZ3
perm_k      = np.array([  15.   ,  15.    ,      15.,    5.   ,     15.,       15.,      15.,      15.,      15.,      15.])*9.8692327E-16          # Rock Permeability, m**2 (converted from mD)
mu_fluid    = np.array([0.002   , 0.002   ,   0.002 ,  0.002  ,   0.002,     0.002,    0.002,    0.002,    0.002,    0.002])                        # Fluid Viscosity, Pa * s
vp          = np.array([4200.   , 4200.   ,    4200.,  4200.  ,   4200.,     4200.,    4200.,    4200.,    4200.,    4200.])                        # Vel P, m/s
vs          = np.array([2400.   , 2400.   ,    2400.,  2400.  ,   2400.,     2400.,    2400.,    2400.,    2400.,    2400.])                        # Vel S, m/s
rho_m       = np.array([2700.   , 2700.   ,    2700.,  2700.  ,   2700.,     2700.,    2700.,    2700.,    2700.,    2700.])                        # Matrix Density
rho_f       = np.array([1000.   , 1000.   ,    1000.,  1000.  ,   1000.,     1000.,    1000.,    1000.,    1000.,    1000.])                        # Fluid Density
mu_solid    = np.array([1.9e25  , 1.9e25  ,    1.9e25, 1.9e25 ,  1.9e25,    1.9e25,   1.9e25,   1.9e25,   1.9e25,   1.9e25])                        # Solid Viscosity 
B           = np.array([0.9     , 0.9     ,    0.9  ,  0.9    ,     0.9,       0.9,      0.9,      0.9,      0.9,      0.9])                        # Biot's Constant
r           = np.array([1.0     , 1.0     ,    1.0  ,  1.0    ,     1.0,       1.0,      1.0,      1.0,      1.0,      1.0])                        # Power Law Exponent
cf          = np.array([2.2e9   , 2.2e9   ,    2.2e9,  2.2e9  ,   2.2E9,     2.2E9,    2.2E9,    2.2E9,    2.2E9,    2.2E9])                        # Fluid Compressibility / Fluid Bulk Modulus
K           = perm_k / mu_fluid                                                                                                                     # Fluid Mobility, m**2 / Pa / s
phi         = np.array([0.10    , 0.10    ,    0.10 ,  0.10   ,    0.10,      0.10,     0.10,     0.10,     0.10,     0.10])                        # Matrix Porosity, frac
source      = np.array([dp_val  , dp_val  ,  dp_val , dp_val  ,      0.,        0.,       0.,       0.,       0.,       0.])*6894.7573              # Source/Sink Pressure Delta, Pa
rho         = rho_m*(1-phi) + rho_f*phi                                                                                                                                                 # (converted from psi)

# Calculate Lame Parameters
E_dyn=rho*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
nu_dyn=(vp**2-2*vs**2)/2/(vp**2-vs**2)
E=E_dyn/1.6
nu=nu_dyn

if init==1:
    numblocks = 10
    numvals = 12
    mat = np.zeros([numblocks,numvals])
    mat[:,0] = E
    mat[:,1] = nu
    mat[:,2] = mu_solid
    mat[:,3] = r
    mat[:,4] = rho
    mat[:,5] = K
    mat[:,6] = B
    mat[:,7] = phi
    mat[:,8] = cf
    mat[:,9] = source
    mat[:,10] = E_dyn
    mat[:,11] = nu_dyn

    mat = mat[0:num_mats,:]
else:      
          #         Rock Properties            Fluid Properties
          # ================================ ======================
          # |                              | |                    |
          # E       nu    mu_s    r    rho_s Kf       B    phi   cf
    mat = [[3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],        #
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-11, 0.9, 0.15, 2.2E9],        #
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],           
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9],                                                       
           [3.0E10, 0.25, 1.0E25, 1.0, 3000, 1.0E-15, 0.9, 0.15, 2.2E9]]

#===============================================================================
# Fault Properties (Per Fault)

# Static Friction (Taking SCEC TVP10 values)

S_fc =   1.0#0.85        # Stationary Friction Coefficient
S_fcd =  1.0#0.65        # Dynamic Friction Coefficient
S_dc = 0.05         # Slip Weakening Distance
S_st_init = 0       # Initial Fault Traction
S_frc = np.int(0)   # Frictional (1 = yes) <==== THIS CAN SHUT OFF THE FAULT
S_coh = 0           # Cohesion
S_dcoh = 1          # Cohesion Weakening Slip
# Rate & State (From SCEC TVP102))
tau0_val  = 75E6
sn0_val   = 120E6
a_val     = np.array([.008 ])
b0_val    = np.array([ 0.6 ])
V0_val    = np.array([ 1E-6])
dtau0_val = np.array([ 25E6])
b_val     = np.array([ .012])
L_val     = np.array([ .02 ])

# Fault Permeability and Initial Stress
FaultPerm=0
#FT_INIT_STRESS=[s*2.7E7, -5E7]


#===============================================================================
# Nodal Force / Flux Parameters
# Commence injection after initial time step - after model has stabilized
flux     =  0. #-112.5 250. # Units?? injection positive
flux     = flux * dt_hr/24.      # convert input into vol / day
flux_mid =  0.


      #             x         y         xforce       yforce         flux        tstart       tend
well_a = np.array([[-2.25,      -2.000,    0,           0,          flux    ,     dt*5,      t*0.8],
                   [-2.25,      -2.000,    0,           0,          flux*0  ,     dt*10,     t*0.5],
                   [-2.25,      -2.000,    0,           0,          flux*0  ,    t*0.5,      t*0.8],
                   [-2.25,      -2.000,    0,           0,          flux*0  , 3*t/4,          t]])

well_b = np.array([[-2.25,      -1.800,    0,           0,          flux    ,     dt*5,      t*0.8],
                   [-2.25,      -1.800,    0,           0,          flux*0  ,     dt*10,     t*0.5],
                   [-2.25,      -1.800,    0,           0,          flux*0  ,    t*0.5,      t*0.8],
                   [-2.25,      -1.800,    0,           0,          flux*0  , 3*t/4,          t]])

well_c = np.array([[-2.25,      -2.150,    0,           0,          flux    ,     dt*5,      t*0.8],
                   [-2.25,      -2.150,    0,           0,          flux*0  ,     dt*10,     t*0.5],
                   [-2.25,      -2.150,    0,           0,          flux*0  ,    t*0.5,      t*0.8],
                   [-2.25,      -2.150,    0,           0,          flux*0  , 3*t/4,          t]])



fnodes = np.vstack([well_a,well_b,well_c])

#===============================================================================
# Debug Shots
      #             x         y       xforce       yforce      flux     tstart       tend
shot_a = np.array([[-.5,   -0.680,   -5E7,        -5E7,          0. ,     dt*2.5,    dt*2.5],
                   [-.5,   -0.680,    5E7,        -5E7,          0. ,     dt*2.5,    dt*2.5],
                   [-.5,   -0.680,   -5E7,         5E7,          0. ,     dt*2.5,    dt*2.5],
                   [-.5,   -0.680,    5E7,         5E7,          0. ,     dt*2.5,    dt*2.5]])

shot_b = np.array([[-.5,   -0.4,     -5E7,        -5E7,          0. ,     dt,        dt],
                   [-.5,   -0.4,      5E7,        -5E7,          0. ,     dt,        dt],
                   [-.5,   -0.4,     -5E7,         5E7,          0. ,     dt,        dt],
                   [-.5,   -0.4,      5E7,         5E7,          0. ,     dt,        dt]])

shot_nodes = np.vstack([shot_a,shot_b])


# Flux     x         y         xforce    yforce      flux    tstart   tend
#fnodes = [[0.,      -2.000,    0,        0,          0.      ,     0,    dt*10],
#          [0.,      -2.000,    0,        0,          flux    ,  dt*10,   t*0.5],
#          [0.,      -2.000,    0,        0,          flux*5  , t*0.5,    t*0.8],
#          [0.,      -2.000,    0,        0,          flux*0  , 3*t/4,       t],
#          [5.750,   -0.620,    0,        0,          flux_mid,     0,       t],
#          [5.750,   -0.640,    0,        0,          flux_mid,     0,       t],
#          [5.750,   -0.660,    0,        0,          flux_mid,     0,       t],
#          [5.750,   -0.680,   -5E7,        -5E7,          0.      ,     dt*5,   dt*5],
#          [5.750,   -0.680,    5E7,        -5E7,          0.      ,     dt*5,   dt*5],
#          [5.750,   -0.680,   -5E7,         5E7,          0.      ,     dt*5,   dt*5],
#          [5.750,   -0.680,    5E7,         5E7,          0.      ,     dt*5,   dt*5]]

#===============================================================================
# Observation Points for Recorded Waveforms

# Observation info
ogrid = np.array([[-2.5, 0.],
                  [-2.0, 0.],
                  [-1.5, 0.],
                  [-1.0, 0.],
                  [-0.5, 0.],
                  [-0.25,0.],
                  [ 0.0, 0.],
                  [ 0.25,0.],
                  [ 0.5, 0.],
                  [ 1.0, 0.],
                  [ 1.5, 0.],
                  [ 2.0, 0.],
                  [ 2.5, 0.]])
#ogrid = np.array([])
nobs=len(ogrid)


#===============================================================================
# SCRIPT RUN ###################################################################
#===============================================================================
# Sort mode conflict
if not (explicit or implicit or fault):
    implicit = True
poro = poro and not explicit
visc = visc and not explicit
if explicit: 
#    t = 7.5; dt = 0.001; nviz = 25; dsp=0
    max_slip = np.array([.1, 0.]).reshape(2,1)
    dt_slip = .02
    slip = max_slip/(dt_slip/dt)
    t_rup = .6 
    alpha= 0.; beta = 0.0125; rfrac = 0
    line1 = ["explicit qua 12"]
    line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
    line4 = np.array([alpha, beta, rfrac]).reshape(1,3)
elif fault: 
    # dsp_str=1 output fault stress; bod_frc=1 body force; hyb=1 hybrid; init=1 initial pressure 
    #t = 3600.*24*200; dt = 3600.*24/2; nviz = 1 
    #t_dyn = 0.025; dt_dyn= 0.00025; t_lim = 5.; dsp=1; dsp_hyb=1; dsp_str=1; rsf=0; v_bg=1E-12
    #bod_frc=1; hyb=1; nviz_dyn=120; nviz_wave=10; nviz_slip=100; init=1
    #alpha= 0.; beta = 0.00125; rfrac = 0
    if  poro and visc:
        line1 = ["fault-pv qua 12"]
    elif poro and not visc:
        line1 = ["fault-p qua 12"]
    elif not poro and visc:
        line1 = ["fault-v qua 12"]
    else:
        line1 = ["fault qua 12"]
    line3 = np.array([t, dt, nviz, dsp]).reshape(1,4)
    line4 = np.array([t_dyn, dt_dyn, nviz_dyn, t_lim, dsp_hyb, dsp_str, bod_frc,hyb,rsf,init]).reshape(1,10)
    if hyb == 1:
        if rsf==0:
            line5 = np.array([nviz_wave,nviz_slip]).reshape(1,2)
        else:
            line5 = np.array([nviz_wave,nviz_slip,v_bg]).reshape(1,3)
    line6 = np.array([alpha, beta, rfrac]).reshape(1,3)
else:
    #t = 1500000.; dt = 10000.; nviz = 5
    if poro and visc:
        line1 =["implicit-pv qua 12"]
    elif poro and not visc:
        line1 =["implicit-p qua 12"]
    elif not poro and visc:
        line1 =["implicit-v qua 12"]
    else:
        line1 =["implicit qua 12"]
    line3 = np.array([t, dt, nviz, 1]).reshape(1,4)

# node data
print 'Extracting mesh...'

coord = np.hstack((nc.variables['coordx'][:].\
        reshape(len(nc.variables['coordx']),1),\
        nc.variables['coordy'][:].\
        reshape(len(nc.variables['coordy']),1)))
qd_node = np.empty(shape=[0, 4], dtype=np.uint32)
mat_typ = np.empty(shape = (0,1), dtype=np.uint32)

for i in nc.variables['eb_prop1'][:]:
    cnct = nc.variables["connect"+str(i)][:]
    n_elem = len(cnct)
    cnct = cnct.reshape(n_elem*4)
    cnct = cnct.reshape(n_elem,4) 
    qd_node = np.vstack((qd_node, cnct))
    mat_typ = np.vstack((mat_typ, i*np.ones((len(cnct),1))))
print '%d nodes, %d elements' %(len(coord), len(qd_node))
meshnumnodes = nc.dimensions['num_nodes'].size
meshnumelems = nc.dimensions['num_elem'].size
# add fault nodes
# For each fault nodeset/sideset....
###############################################################################
#===============================================================================
# Fault Processing Section
for fault_set in faultID:
    print 'Forming fault constraints...'
    id_tmp = nc.variables['ss_prop1'][fault_set-1]
    el_flt = nc.variables['elem_ss' + str(id_tmp)]
    sd_flt = nc.variables['side_ss' + str(id_tmp)]
    node_flt = qd_node[el_flt[:]-1,:]
    node_tap = np.empty((0,2),dtype=np.uint32)
    node_flt_p = np.empty((0),dtype=np.uint32)
    crd_flt_p = np.empty((0,2),dtype=float)
    node_flt_n = np.empty((0),dtype=np.uint32)
    crd_flt_n = np.empty((0,2),dtype=float)
    sd_flt_p = np.empty((0),dtype=np.uint32)
    spair_flt = np.empty((0,2),dtype=np.uint32)
    for i in range(len(el_flt)):
        el = el_flt[i]
        sd = sd_flt[i]
        node = node_flt[i,:]
        if   sd ==1: node_on = node[[0,1]];idn=[0,1];node_off=node[[2,3]] 
        elif sd ==2: node_on = node[[1,2]];idn=[1,2];node_off=node[[0,3]] 
        elif sd ==3: node_on = node[[2,3]];idn=[2,3];node_off=node[[0,1]] 
        elif sd ==4: node_on = node[[3,0]];idn=[3,0];node_off=node[[1,2]] 
        # negative side has elements lower right than the fault
        if sum(coord[node_on-1,1]-coord[node_on-1,0]) > sum(coord[node_off-1,1]-coord[node_off-1,0]):
            for j in range(2):
                if node_on[j] in node_tap[:,0]: 
                    node_add = node_tap[node_tap[:,0]==node_on[j],1]
                    qd_node[el-1,idn[j]] = node_add
                else:
                    node_add = len(coord)+len(node_tap)+1
                    qd_node[el-1,idn[j]] = node_add
                    node_tap = np.vstack((node_tap,[node_on[j],node_add]))
                    node_flt_n = np.hstack((node_flt_n, node_add))
                    crd_flt_n = np.vstack((crd_flt_n, coord[node_on[j]-1,:])) 
        else:
            for j in range(2):
                if not (node_on[j] in node_flt_p):
                    node_flt_p = np.hstack((node_flt_p, node_on[j]))
                    crd_flt_p = np.vstack((crd_flt_p, coord[node_on[j]-1,:]))
            sd_flt_p = np.hstack((sd_flt_p,sd))        
            spair_flt = np.vstack((spair_flt,node_on))
    # remove split node at tip
    crd_add=np.empty((0),dtype=bool)
    #crd_add = np.zeros(len(node_tap)).astype(np.bool)
    for i in range(len(node_tap)):

        loc = qd_node==node_tap[i,1]

        loc_flt = node_flt_n==node_tap[i,1]
        
        if sum(sum(loc))<2:
            qd_node[loc] = node_tap[i,0]
            node_flt_n[loc_flt] = node_tap[i,0]
            crd_add=np.hstack((crd_add,False))
        else:
            qd_node[loc] = qd_node[loc]-sum(~crd_add)
            node_flt_n[loc_flt] = node_flt_n[loc_flt]-sum(~crd_add)
            crd_add=np.hstack((crd_add,True))
    coord = np.vstack((coord,crd_flt_n[crd_add,:])) # 258 rows added here


    # Pair fault nodes
    ft_map = np.array(np.array(np.all((crd_flt_p[:,None,:]\
           ==crd_flt_n[None,:,:]),axis=-1).nonzero()).T.tolist())
    node_flt_n = node_flt_n[ft_map[:,1]]
    node_flt_p = node_flt_p[ft_map[:,0]]
    crd_flt_n = crd_flt_n[ft_map[:,1],:]
    crd_flt_p = crd_flt_n[ft_map[:,0],:]
    crd_add = crd_add[ft_map[:,1]]
    node_flt_n = node_flt_n[crd_add]
    node_flt_p = node_flt_p[crd_add]
    crd_flt_n = crd_flt_n[crd_add,:]
    crd_flt_p = crd_flt_p[crd_add,:]

    # Fault's strike and normal vectors
    vec_fn = np.zeros((len(node_flt_p), 2), dtype=float)
    vec_fs = np.zeros((len(node_flt_p), 2), dtype=float)
    for i in range(len(spair_flt)):
        v1 = coord[spair_flt[i,0]-1,:]
        v2 = coord[spair_flt[i,1]-1,:]
        vec1 = v2 - v1
        vec1 = -vec1*np.sign(vec1[1])
        vec2 = np.array([vec1[1], -vec1[0]])
        vec2 = vec2*np.sign(vec2[1])
        row = np.squeeze(np.hstack((np.where(node_flt_p==spair_flt[i,0]),\
              np.where(node_flt_p==spair_flt[i,1]))))
        vec_fs[row,:] += vec1 
        vec_fn[row,:] += vec2
    vec_fs /= (np.ones((2,1))*np.linalg.norm(vec_fs, axis=1)).T
    vec_fn /= (np.ones((2,1))*np.linalg.norm(vec_fn, axis=1)).T
    vecf = np.empty(shape=(0,6))
    xfnode = np.empty(shape=(0,2))
    num_flt_nodes = len(node_flt_p)

    #===============================================================================
    # Fault Permeability and Initial Stress
    #-------------------------------------------------------------------------------

    # Make the fault permeable at certain segment; add initial stress 
    perm    = np.ones((num_flt_nodes,1))
    st_init = np.zeros((num_flt_nodes,2))

    y_min = min(coord[:,1])
    y_max = max(coord[:,1])

    # Loop through all positive nodes in fault
    for node_pos, i in zip(node_flt_p,range(len(node_flt_p))):
        y = coord[node_pos - 1,1]
        s = -np.sign(vec_fn[i,0])
        if y > y_min + (y_max - y_min)*0.6 and y < y_min + (y_max - y_min)*0.75:
            perm[i]=FaultPerm
            #st_init[i,:]=[s*2.7E7, 5E7]
            st_init[i,:]=[0., 0.]
        else:
            perm[i]=FaultPerm
            #st_init[i,:]=[s*2.7E7, 5E7]
            st_init[i,:]=[0., 0.]

    if rsf==1: # Rate-state parameters
        a     =     a_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        b0    =    b0_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        V0    =    V0_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        dtau0 = dtau0_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        b     =     b_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        L     =     L_val[fault_set-1]  * np.ones((num_flt_nodes,1))
        theta_init = L/V0*np.exp((a*np.log(2.*np.sinh(0.7/a))-b0-a*np.log(v_bg/V0))/b)


#===============================================================================
# Locate Non Conformal Nodes
# These refer to nodes that are not exclusively connected by terminating
# bars from elements
#-------------------------------------------------------------------------------


print 'Locating nonconformal nodes...' 
NCF_slave2master = np.zeros((0,3),np.uint32) 
NCF_weight  = np.zeros((0,2),np.float)
hit=[]
pair=np.empty((0,2,),dtype=np.uint32)
for s in NCF_node_list:
    master_nodes  = nc.variables['node_ns'+str(s)][:]
    slave_nodes   = nc.variables['node_ns'+str(s+1)][:]
    master_nodes  = np.setdiff1d(master_nodes,node_flt_n)
    master_nodes  = np.setdiff1d(master_nodes,node_flt_p)
    slave_nodes   = np.setdiff1d(slave_nodes,node_flt_n)
    slave_nodes   = np.setdiff1d(slave_nodes,node_flt_p)
    master_coord  = coord[master_nodes-1,:]
    slave2master  = np.zeros((0,3),np.uint32) 
    weight_master = np.zeros((0,2),np.float)
    slave_nodes   = np.setdiff1d(slave_nodes,np.array(hit))
    for slave in slave_nodes:
        dis2node =  np.linalg.norm(master_coord - np.dot(np.ones((len(master_nodes),1),dtype=np.float),coord[slave-1,:].reshape(1,2)),axis=1)
        master2 = np.argsort(dis2node)[:2]
        if min(dis2node) < 1E-6: # Node on node
            nodetmp = master_nodes[master2[0]]
            hit.append(slave)
            pair=np.vstack((pair,[slave,nodetmp]))
        else:
            onedg=False
            for el, side in zip(nc.variables['elem_ss'+str(s)][:], nc.variables['side_ss'+str(s)][:]):
                el_node = qd_node[el-1,:]
                if   side ==1: e_node = el_node[[0,1]] 
                elif side ==2: e_node = el_node[[1,2]] 
                elif side ==3: e_node = el_node[[2,3]] 
                elif side ==4: e_node = el_node[[3,0]]
                v2j=coord[e_node[0]-1,:]-coord[slave-1,:]
                v2k=coord[e_node[1]-1,:]-coord[slave-1,:]
                onedg=np.linalg.norm(np.cross(v2j,v2k))<1E-6 and np.dot(v2j,v2k)<0 # Node on edge
                if onedg:
                    dis2j=np.linalg.norm(v2j)
                    dis2k=np.linalg.norm(v2k)
                    wj = dis2k/(dis2j+dis2k); wk = dis2j/(dis2j+dis2k)
                    slave2master= np.vstack((slave2master,np.hstack((slave,e_node[0],e_node[1]))))
                    weight_master=np.vstack((weight_master,[wj,wk]))
                    hit.append(slave)
                    break
    NCF_slave2master=np.vstack((NCF_slave2master,slave2master))
    NCF_weight=np.vstack((NCF_weight,weight_master))
print '%d nonconformal nodes located' %(len(NCF_slave2master))


#----------------------Boundary ID-------------------------
# Side/nodesets: 0 fault, 1 lower, 2 right, 3 upper, 4 left
#----------------------------------------------------------
bnd_el = []
# Traction and abs boundaries,id 0 preserved for fault faces
for i in nc.variables['ss_prop1'][first_bnd:last_bnd]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),\
            sides.reshape(len(sides),1))))

#===============================================================================
# fixed nodes

if poro:
    # For 2D poroelastic case, BC format is as follows:
    # [(x flag) (y flag) (poro flag) (nodal pore pressure)]
    # nodal pore pressure only matters if poro flag = 0
    bc_typ = np.ones((len(coord),3), dtype=np.int8)
else:
    bc_typ = np.ones((len(coord),2), dtype=np.int8)

# Y Direction
if fixedID_Y[0] == -999:
    bcy_nodes = []
else:
    for i in fixedID_Y:
        bcy_nodes = nc.variables['node_ns' + str(i)][:] 
    for node in bcy_nodes:
        bc_typ[node - 1, 1] = 0   

# X Direction
if fixedID_X[0] == -999:
    bcx_nodes = []
else:
    for i in fixedID_X:
        bcx_nodes = nc.variables['node_ns' + str(i)][:]
    for node in bcx_nodes:
        bc_typ[node - 1, 0] = 0

# Poro BC - Set only for floor at the moment
if poro:
    if fixedID_Poro[0] == -999:
        bcporo_nodes = []
    else:
        loopcount = 0
        for i in fixedID_Poro:
            bcporo_nodes = nc.variables['node_ns' + str(i)][:]
            for node in bcporo_nodes:
                bc_typ[node - 1, 2] = 0
                #bc_typ[node - 1, 3] = 4000.*9.80*1000.
                

###################
# LM Dirichlet BC (pressure)
# Set fixed pressure at floor, or at least try to
# Assuming hydrostatic gradient, pure water, etc
floor_press = 4000.*9.80*1000.
hydrostatic_nodes = nc.variables['node_ns2']
hydrostatic_press = np.ones(len(hydrostatic_nodes))*floor_press

surface_nodes = nc.variables['node_ns4']
surface_press = np.zeros(len(surface_nodes))


#===============================================================================
# nodal force/flux 
if poro:
    fnode_bc = np.empty((0,6))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((len(coord), 1)),\
                np.array(fnode[:2]).reshape(1,2)), axis=1)
        row = np.argsort(dis)[0]
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[2:]))))
else:
    fnode_bc = np.empty((0,5))
    for fnode in fnodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((len(coord), 1)),\
                np.array(fnode[:2]).reshape(1,2)), axis=1)
        row = np.argsort(dis)[0]
        fnode_bc = np.vstack((fnode_bc, np.hstack((row + 1, fnode[2:]))))

#===============================================================================
# Debug Shots
if len(shot_nodes) > 0:
    shot_info = np.empty((0,6))
    for snode in shot_nodes:
        dis = np.linalg.norm(coord - np.dot(np.ones((len(coord), 1)),\
                np.array(snode[:2]).reshape(1,2)), axis=1)            
        row = np.argsort(dis)[0]
        shot_info = np.vstack((shot_info, np.hstack((row + 1, fnode[2:]))))
        
#===============================================================================
# traction bc 


bnd_el = []
# Traction and abs boundaries,id 0 preserved for fault faces

# Floor traction
floor_bnd_el = []
for j in np.arange(len(TracBCFloor)):
    i = nc.variables['ss_prop1'][TracBCFloor[j]-1]
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    floor_bnd_el.append(np.hstack((els.reshape(len(els),1),\
                sides.reshape(len(sides),1))))
floor_bnd_el = floor_bnd_el[0]                

# West traction
west_bnd_el = []
for j in np.arange(len(TracBCWest)):
    i = nc.variables['ss_prop1'][TracBCWest[j]-1]
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    west_bnd_el.append(np.hstack((els.reshape(len(els),1),\
                sides.reshape(len(sides),1))))
west_bnd_el = west_bnd_el[0]

# East traction
east_bnd_el = []
for j in np.arange(len(TracBCEast)):
    i = nc.variables['ss_prop1'][TracBCEast[j]-1]
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    east_bnd_el.append(np.hstack((els.reshape(len(els),1),\
                sides.reshape(len(sides),1))))
east_bnd_el = east_bnd_el[0]

# Surface traction
surf_bnd_el = []
for j in np.arange(len(TracBCSurface)):
    i = nc.variables['ss_prop1'][TracBCSurface[j]-1]
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    surf_bnd_el.append(np.hstack((els.reshape(len(els),1),\
                sides.reshape(len(sides),1))))
surf_bnd_el = surf_bnd_el[0]

# Vestigal
for i in nc.variables['ss_prop1'][first_bnd-1:last_bnd-1]:
    els = nc.variables['elem_ss' + str(i)][:]
    sides = nc.variables['side_ss' + str(i)][:]
    bnd_el.append(np.hstack((els.reshape(len(els),1),\
            sides.reshape(len(sides),1))))
#trac_el1 = bnd_el[3]
#trac_el2 = bnd_el[2]
#abs_bc1 = bnd_el[0]
#abs_bc2 = bnd_el[1]

if explicit:
    trac_val = [0., 0.]
if poro:
    # Apply predetermined tractions to length of boundaries.
    # Later this will be more specific
    floor_trac_bc = np.zeros(shape=[len(floor_bnd_el), 5])
    west_trac_bc  = np.zeros(shape=[len(west_bnd_el ), 5])    
    east_trac_bc  = np.zeros(shape=[len(east_bnd_el ), 5])        
    surf_trac_bc  = np.zeros(shape=[len(surf_bnd_el ), 5])
    floor_trac_bc[:, 0:] = trac_val_F
    west_trac_bc[:,  0:] = trac_val_W
    east_trac_bc[:,  0:] = trac_val_E
    surf_trac_bc[:,  0:] = trac_val_S
    trac_el = np.vstack( (floor_bnd_el,  west_bnd_el,  east_bnd_el,  surf_bnd_el ))
    trac_bc = np.vstack( (floor_trac_bc, west_trac_bc, east_trac_bc, surf_trac_bc))
else:
    floor_trac_bc = np.zeros(shape=[len(floor_bnd_el), 4])
    west_trac_bc  = np.zeros(shape=[len(west_bnd_el ), 4])    
    east_trac_bc  = np.zeros(shape=[len(east_bnd_el ), 4])        
    surf_trac_bc  = np.zeros(shape=[len(surf_bnd_el ), 4])
    floor_trac_bc[:, 0:] = trac_val_F
    west_trac_bc[:,  0:] = trac_val_W
    east_trac_bc[:,  0:] = trac_val_E
    surf_trac_bc[:,  0:] = trac_val_S
    trac_el = np.vstack( (floor_bnd_el,  west_bnd_el,  east_bnd_el,  surf_bnd_el ))
    trac_bc = np.vstack( (floor_trac_bc, west_trac_bc, east_trac_bc, surf_trac_bc))


#===============================================================================
# absorbing bc 

# I guess for 2D cases we are always assume the free surface condition 
#
# This will need some work as well - for now all of east and west, plus lower 
# are crammed accordingly 

abs_bc_floor = np.hstack((floor_bnd_el, 2*np.ones((len(floor_trac_bc) , 1))))
abs_bc_surf  = np.hstack((surf_bnd_el , 2*np.ones((len(surf_trac_bc), 1))))
abs_bc_west  = np.hstack((west_bnd_el ,   np.ones((len(west_trac_bc) , 1))))
abs_bc_east  = np.hstack((east_bnd_el ,   np.ones((len(east_trac_bc), 1))))

abs_bc = np.vstack((abs_bc_floor, abs_bc_west, abs_bc_east))
#  abs_bc4 upper bound is not absorbing


#===============================================================================
# Constraint Function
#-------------------------------------------------------------------------------

# Total length of constraint function
if poro:
    neqNCF=(dim+1)*len(NCF_slave2master)
    neqFT=dim*num_flt_nodes + sum(perm)
else:
    neqNCF=dim*len(NCF_slave2master)
    neqFT=dim*num_flt_nodes
neq = neqNCF+neqFT 
print '%d NCF and %d fault constraint equations.' %(neqNCF,neqFT)

#===============================================================================
# Eliminate Doubling Nodes
#-------------------------------------------------------------------------------

for i in range(len(pair)):
    slave  = pair[i,0]
    master = pair[i,1]
    qd_node[ qd_node == slave] = master
    NCF_slave2master[NCF_slave2master == slave] = master
    node_flt_n[node_flt_n == slave] = master
    node_flt_p[node_flt_p == slave] = master
    fnode_bc[fnode_bc[:,0] == slave, 0] = master
    coord=np.vstack( (coord[:slave-1,:], coord[slave:,:]) )
    bc_typ=np.vstack ((bc_typ[:slave-1, :], bc_typ[slave:, :]) )
    id_shift = qd_node > slave
    qd_node[id_shift] = qd_node[id_shift] - 1
    id_shift = NCF_slave2master > slave
    NCF_slave2master[id_shift] = NCF_slave2master[id_shift] - 1
    id_shift = node_flt_p > slave
    node_flt_p[id_shift] = node_flt_p[id_shift] - 1
    id_shift = node_flt_n > slave
    node_flt_n[id_shift] = node_flt_n[id_shift] - 1
    id_shift = fnode_bc[:,0] > slave
    fnode_bc[id_shift,0] = fnode_bc[id_shift,0]-1
    id_shift = pair > slave
    pair[id_shift] = pair[id_shift] - 1

#===============================================================================
# Export to file
#-------------------------------------------------------------------------------

# Export to Defmod .inp file
fout = fin.rsplit('.')[0] + '.inp'
print 'Write to ' + fout + '...'
if os.path.isfile(fout): os.remove(fout)

# Consider Dirichlet Constraints
neq = neq + len(hydrostatic_nodes) + len(surface_nodes)


f = open(fout, 'a')
if explicit:
    line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
            len(fnode_bc), len(trac_el), len(abs_bc)]).reshape(1,7)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d %d %d %d %d %d %d')
    np.savetxt(f, line3, delimiter=' ', fmt='%.4f %.4f %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %g')
elif fault:
    line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
            len(fnode_bc), len(trac_el), len(abs_bc), num_flt_nodes, len(ogrid), neqNCF]).reshape(1,10)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d '*10)
    np.savetxt(f, line3, delimiter=' ', fmt='%g %g %d %d')
    np.savetxt(f, line4, delimiter=' ', fmt='%g %g %d %g %d %d %d %d %d %d')
    if hyb == 1 and rsf==0:
        np.savetxt(f, line5, delimiter=' ', fmt='%d %d')
    if hyb == 1 and rsf==1:
        np.savetxt(f, line5, delimiter=' ', fmt='%d %d %g')
    np.savetxt(f, line6, delimiter=' ', fmt='%g %g %g')
else:
#    line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
#            len(fnode_bc), len(trac_el), 0, num_flt_nodes]).reshape(1,8)
    line2 = np.array([len(qd_node), len(coord), len(mat), neq,\
            len(fnode_bc), len(trac_el), 0, nobs]).reshape(1,8)
    np.savetxt(f, line1, fmt='%s')
    np.savetxt(f, line2, delimiter=' ', fmt='%d %d %d %d %d %d %d %d')
    np.savetxt(f, line3, delimiter=' ', fmt='%.1f %.1f %d %d')
np.savetxt(f, np.column_stack((qd_node, mat_typ)), delimiter=' ',\
           fmt='%d %d %d %d %d')
if poro:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g %g %d %d %d')
else:
    np.savetxt(f, np.column_stack((coord, bc_typ)) , delimiter = ' ',\
               fmt='%g %g %d %d')
if init==1:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*12)
else:
    np.savetxt(f, mat, delimiter=' ', fmt = '%g '*9) 

#===============================================================================
# Write NCF constraints
for slave2master, weight in zip(NCF_slave2master,NCF_weight):
    npt = np.count_nonzero(slave2master)
    if npt==3:
        if explicit or not poro:
            vec1  = [[   -1, 0, slave2master[0]], 
                     [weight[0], 0, slave2master[1]],
                     [weight[1], 0, slave2master[2]]]
            vec2  = [[0,    -1, slave2master[0]], 
                     [0, weight[0], slave2master[1]],
                     [0, weight[1], slave2master[2]]]
        elif (fault or implicit) and poro:
            vec1  = [[   -1, 0, 0, slave2master[0]], 
                     [weight[0], 0, 0, slave2master[1]],
                     [weight[1], 0, 0, slave2master[2]]]
            vec2  = [[0,    -1, 0, slave2master[0]], 
                     [0, weight[0], 0, slave2master[1]],
                     [0, weight[1], 0, slave2master[2]]]
            vec3  = [[0, 0,    -1, slave2master[0]], 
                     [0, 0, weight[0], slave2master[1]],
                     [0, 0, weight[1], slave2master[2]]]
    if (explicit or not poro):
        np.savetxt(f, [npt], fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %d') 
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
    elif (fault or implicit) and poro:
        np.savetxt(f, [npt], fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")
        np.savetxt(f, [npt], fmt = '%d')
        np.savetxt(f, vec3, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, [[0.,0.,0.]], delimiter = ' ', fmt = "%1.2E %g %g")

# LM Dirichlet Test
i1 = [1]

if poro:
    if len(hydrostatic_press) > 0:
        for i in np.arange(0,len(hydrostatic_nodes)):
            line1 = [[0, 0, 1, hydrostatic_nodes[i] ]]
            line2 = [[hydrostatic_press[i], 0.0, dt]]
            np.savetxt(f, i1,    fmt = '%d')
            np.savetxt(f, line1, delimiter = ' ', fmt = '%g %g %g %d')
            np.savetxt(f, line2, delimiter = ' ', fmt = "%d %g %d")

    if len(surface_press) > 0:
        for i in np.arange(0,len(surface_nodes)):
            line1 = [[0, 0, 1, surface_nodes[i] ]]
            line2 = [[surface_press[i], 0.0, dt]]
            np.savetxt(f, i1,    fmt = '%d')
            np.savetxt(f, line1, delimiter = ' ', fmt = '%g %g %g %d')
            np.savetxt(f, line2, delimiter = ' ', fmt = "%d %g %d")




#===============================================================================
# fault slip: strike and open
if not explicit:
    slip = np.array([0.0, 0.0]).reshape(2,1) 
    c = 0; dt_slip=0; t_rup=0
n = [2]
ft_neg_nodes_tap = []
j = 0
for node_p, node_n, i in zip(node_flt_p, node_flt_n, range(len(node_flt_p))):
    if explicit or not poro:
        vec1  = [[1,  0, node_p], 
                [-1,  0, node_n]]
        vec2  = [[0,  1, node_p], 
                 [0, -1, node_n]]
    elif (fault or implicit) and poro:
        vec1  = [[1, 0, 0, node_p], 
                [-1, 0, 0, node_n]]
        vec2  = [[0, 1, 0, node_p], 
                 [0,-1, 0, node_n]]
    mat_ft = np.hstack((vec_fs[i,:].reshape(2,1),\
                vec_fn[i,:].reshape(2,1)))
    mat_f = np.matrix.transpose(mat_ft).reshape(1,4)
    val = np.dot(mat_ft,slip)
    y = np.array(coord[node_p - 1,:][1])
    if explicit: 
       c = np.sqrt(1 - (y+2.)**2)
    t_act =  dt_slip+(1-c)*t_rup 
    t_slip = [t_act-dt_slip/2, t_act+dt_slip/2]
    cval1 = np.hstack((c*val[0], t_slip)).reshape(1,3)
    cval2 = np.hstack((c*val[1], t_slip)).reshape(1,3)
    if explicit or not poro:
        np.savetxt(f, n, fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %d') 
        np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, n, fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %d')
        np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
    elif (fault or implicit) and poro:
        np.savetxt(f, n, fmt = '%d') 
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %g %d') 
        np.savetxt(f, cval1, delimiter = ' ', fmt = "%1.2E %g %g") 
        np.savetxt(f, n, fmt = '%d')
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, cval2, delimiter = ' ', fmt = "%1.2E %g %g")
    vecf = np.vstack((vecf,np.hstack(([[node_p, node_n]], mat_f))))
    xfnode = np.vstack((xfnode,coord[node_p-1,:]))
    # permeable fault 
    if poro and perm[j] > 0 and not explicit: 
        vec4 = [[0, 0,  1, node_p], 
                [0, 0, -1, node_n]]
        cval4 =[[0, 0, 0]]
        np.savetxt(f, n, fmt = '%d')
        np.savetxt(f, vec4, delimiter = ' ', fmt = '%g %g %g %d')
        np.savetxt(f, cval4, delimiter = ' ', fmt = "%1.2E %g %g")
    j = j + 1

#===============================================================================
# Write test shot/impulsive forces

    
#===============================================================================    
# Define frictional parameters (static friction)

if fault:
    fc = np.empty((len(vecf),1),dtype=float)
    fcd = np.empty((len(vecf),1),dtype=float)
    dc = np.empty((len(vecf),1),dtype=float)
    frc = np.empty((len(vecf),1),dtype=np.uint32)
    for node_pos, i in zip(node_flt_p,range(len(node_flt_p))):
        x,y = coord[node_pos - 1,:]
        if y > y_min + (y_max - y_min)*0.48 and y < y_min + (y_max - y_min)*0.52:
            fc[i] = S_fc #.7
            fcd[i] = S_fcd #.4
            dc[i] = S_dc #.01
            frc[i] = S_frc # 1
        else:
            fc[i] = S_fc #.7
            fcd[i] = S_fcd #.4
            dc[i] = S_dc #.01
            frc[i] = S_frc # 1
    coh =  np.ones((len(vecf),1),dtype=float)*S_coh       # Cohesion
    dcoh = np.ones((len(vecf),1),dtype=float)*S_dcoh      # Cohesion Weakening Slip
    # Write fault orientation tensor + frictional parameters
    if rsf==1:
        if poro:
            np.savetxt(f, np.hstack((vecf, b0, V0, dtau0, a, b, L, theta_init, perm, st_init, xfnode, frc,coh,dcoh)), delimiter = ' ',\
               fmt = '%d '*2 + '%g '*4 + '%g '*7 + '%d ' + '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
        else:
            np.savetxt(f, np.hstack((vecf, b0, V0, dtau0, a, b, L, theta_init, perm, st_init, xfnode, frc,coh,dcoh)), delimiter = ' ',\
               fmt = '%d '*2 + '%g '*4 + '%g '*7 +         '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
    else:
        if poro:
            np.savetxt(f, np.hstack((vecf, fc, fcd, dc, perm, st_init, xfnode, frc,coh,dcoh)), delimiter = ' ',\
               fmt = '%d '*2 + '%g '*4 + '%g '*3 + '%d ' + '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
        else:
            np.savetxt(f, np.hstack((vecf, fc, fcd, dc, st_init, xfnode, frc,coh,dcoh)), delimiter = ' ', \
               fmt = '%d '*2 + '%g '*4 + '%g '*3 +         '%g '*2 + '%g '*2 + '%d ' + '%g '*2)
# Form rotated constraint matrix for fault model
if fault:
    for node_p, node_n, i in zip(node_flt_p, node_flt_n, range(len(node_flt_p))):
        mat_ft = np.hstack((vec_fs[i,:].reshape(2,1),vec_fn[i,:].reshape(2,1)))
        vec = np.array([1., 0.]).reshape(2,1)
        vec = np.dot(mat_ft, vec).reshape(2,)
        vec1  = [[ vec[0],  vec[1], node_p], 
                 [-vec[0], -vec[1], node_n]]
        vec = np.array([0., 1.]).reshape(2,1)
        vec = np.dot(mat_ft, vec).reshape(2,)
        vec2  = [[ vec[0], vec[1],  node_p], 
                 [-vec[0], -vec[1], node_n]]
        np.savetxt(f, vec1, delimiter = ' ', fmt = '%g %g %d') 
        np.savetxt(f, vec2, delimiter = ' ', fmt = '%g %g %d')

# Point force/source and boundary traction/flux
if poro:
    np.savetxt(f, fnode_bc, delimiter=' ',\
            fmt ='%d %1.2E %1.2E %1.2E %g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %1.2E %g %g')
else:
    np.savetxt(f, fnode_bc, delimiter=' ', fmt ='%d %1.2E %1.2E %g %g')
    np.savetxt(f, np.column_stack((trac_el, trac_bc)), delimiter=' ',\
            fmt ='%d %d %1.2E %1.2E %g %g') 
 
# Observation grid
if len(ogrid) > 0:
    np.savetxt(f, ogrid , delimiter = ' ', fmt='%g '*2)
# Abs boundary
if (explicit or fault): np.savetxt(f, abs_bc, delimiter=' ', fmt='%d %d %d')
f.close(); 
print 'Defmod file ' + fout + ' created'

## h* and CFL condition
#lbd=rho*(vp**2-2*vs**2); mu=rho*vs**2
#if rsf==1:
#    L=np.mean(L); a=np.mean(a); b=np.mean(b)
#    sigma_e=max(map(abs,trac_val))
#    hstar=min(2*mu*L/(b-a)/np.pi/sigma_e)
#    print "Critical RSF distance h*=%0.3f m" %hstar
#Lc=min(dt_dyn*np.sqrt(E/rho))
#print "Critical element length h=%0.3f m" %Lc

################################################################################
##===============================================================================
## Visualize mesh
##-------------------------------------------------------------------------------

## Mesh name
#mesh_name = fin.rsplit('.')[0]

## Element Block Listing
#eb_list = nc.variables['eb_prop1'][:]
#X = nc.variables['coordx']
#Y = nc.variables['coordy']
#xy = np.array([X[:], Y[:]]).T
#connect_list = []
#elem_patches = []
#all_connect = np.empty([0,4])

#for i in eb_list:
#    connect = nc.variables['connect' + str(i)]
#    connect_list.append(nc.variables['connect' + str(i)])
#    all_connect = np.vstack((all_connect,nc.variables['connect'+str(i)][:]))
#    patches = []
#    for plot_coords in xy[connect[:]-1]:
#        plot_quad = Polygon(plot_coords, True)
#        patches.append(plot_quad)
#    elem_patches.append(patches)


#    
#fig, ax = plt.subplots()


##fig = plt.figure()
#colors = 100 * np.random.rand(len(patches))
## Add element block patch collections
#eb1 = PatchCollection(elem_patches[0], cmap=matplotlib.cm.coolwarm, alpha=0.4,label='Basement')
#eb2 = PatchCollection(elem_patches[1], cmap=matplotlib.cm.coolwarm, alpha=0.1,label='Low Permeability Formation')
#eb3 = PatchCollection(elem_patches[2], cmap=matplotlib.cm.coolwarm, alpha=0.3,label='High Permeability Formation')
#eb4 = PatchCollection(elem_patches[3], cmap=matplotlib.cm.coolwarm, alpha=0.4,label='Mudrock')


## Set element block colors
#eb1.set_color((0,0,0))
#eb2.set_color((0,1,0))
#eb3.set_color((0,0,1))
#eb4.set_color((1,1,0))
#eb1_l = plt.Rectangle((0,0),1,1,fc=eb1.get_facecolor()[0])
#eb2_l = plt.Rectangle((0,0),1,1,fc=eb2.get_facecolor()[0])
#eb3_l = plt.Rectangle((0,0),1,1,fc=eb3.get_facecolor()[0])
#eb4_l = plt.Rectangle((0,0),1,1,fc=eb4.get_facecolor()[0])
## Label element blocks

##eb1.set_label('Basement')
##eb2.set_label('Low Permeability Formation')
##eb3.set_label('High Permeability Formation')
##eb4.set_label('Mudrock')
##eb1.set_array(np.array([0.1]))
##eb2.set_array(np.array([0.25]))
##eb3.set_array(np.array([0.5]))
##eb4.set_array(np.array([0.75]))

#ax.add_collection(eb1)
#ax.add_collection(eb2)
#ax.add_collection(eb3)
#ax.add_collection(eb4)

## Draw Faults
#faultID = 1

#plot_flt_nds = nc.variables['node_ns'+ str(faultID)]
#plot_flt_coords = xy[plot_flt_nds[:]-1]
#plot_flt_LHS = plot_flt_coords[np.flatnonzero(plot_flt_coords[:,0] < 0)]
#plot_flt_RHS = plot_flt_coords[np.flatnonzero(plot_flt_coords[:,0] > 0)]

## Highlight faults
#flt_LHS = mpl.lines.Line2D(plot_flt_LHS[:,0],plot_flt_LHS[:,1], 
#                            color='red',label='LHS Fault',linewidth=5,alpha=0.45)
#flt_RHS = mpl.lines.Line2D(plot_flt_RHS[:,0],plot_flt_RHS[:,1], 
#                            color='green',label='RHS Fault',linewidth=5,alpha=0.45)

#ax.add_line(flt_LHS)
#ax.add_line(flt_RHS)


### Indicate well location

#WF_xy = fnodes[0][0:2]

#plot_WF = mpl.patches.Ellipse((WF_xy[:]),0.05,0.05,0.0,color='red')
#ax.add_patch(plot_WF)
## Point out wellfoot
##ax.annotate('$Well Foot$',
##            xy=(WF_xy[0],WF_xy[1]), xycoords='data',
##            xytext=(-3.5, -3.5), textcoords='offset points',
##            bbox=dict(boxstyle="round", fc="0.8"),
##            arrowprops=dict(arrowstyle="->",
##                            patchB=el,
##                            connectionstyle="angle,angleA=90,angleB=0,rad=10"))




## Plot seismic receivers

#ax.scatter(ogrid[:,0],ogrid[:,1],c='black',marker='v',s=32)

##Label BC Traction
#props = dict(boxstyle='square', facecolor='white', alpha=0.6)

## Surface
#plot_surf_trac = 'Trac, X: ' + str(trac_val_S[0]) + ' Pa, Z: ' + str(trac_val_S[1]) + ' Pa'
#plot_surf_fixed = 'NOT Fixed'
#plot_surf_abc = 'NO Absorbing BC'
#plot_surf = plot_surf_trac + '\n' + plot_surf_fixed + '\n' + plot_surf_abc
#ax.text(0.,0.05,plot_surf,ha='center',fontsize=10,
#        bbox=props,fontname='courier',weight='bold')

## Floor
#plot_floor_trac = 'Trac, X: ' + str(trac_val_F[0]) + ' Pa, Z: ' + str(trac_val_F[1]) + ' Pa'
#plot_floor_fixed = 'Fixed'
#plot_floor_abc = 'Absorbing BC'
#plot_floor = plot_floor_trac + '\n' + plot_floor_fixed + '\n' + plot_floor_abc
#ax.text(0.,-4.35,plot_floor,ha='center',fontsize=10,
#        bbox=props,fontname='courier',weight='bold')

## West
#plot_west_trac = 'Trac, X: ' + str(trac_val_W[0]) + ' Pa, Z: ' + str(trac_val_W[1]) + ' Pa'
#plot_west_fixed = 'NOT Fixed'
#plot_west_abc = 'Absorbing BC'
#plot_west = plot_west_trac + '\n' + plot_west_fixed + '\n' + plot_west_abc
#ax.text(-3.2,-2.,plot_west,ha='center',fontsize=10,
#        bbox=props,fontname='courier',rotation=90.,weight='bold')

## East
#plot_east_trac = 'Trac: X: ' + str(trac_val_E[0]) + ' Pa, Z: ' + str(trac_val_E[1]) + ' Pa'
#plot_east_fixed = 'NOT Fixed'
#plot_east_abc = 'Absorbing BC'
#plot_east = plot_east_trac + '\n' + plot_east_fixed + '\n' + plot_east_abc
#ax.text(3.2,-2.,plot_east,ha='center',fontsize=10,
#        bbox=props,fontname='courier',rotation=90.,weight='bold')

## Misc. info box

#if bod_frc == 1:
#    bf_text = 'Sim Parameters' + '\n' + '==============' + '\n' + 'Gravity ON'
#    ax.text(-3,-4.15,bf_text,ha='left',va='top',fontsize=8,
#        bbox=props,fontname='courier',rotation=0.,weight='bold')

#################################################################################
## Finalize Plot
#################################################################################
#xmax = X[:].max() + 0.5
#xmin = X[:].min() - 0.5
#ymax = Y[:].max() + 0.5
#ymin = Y[:].min() - 0.5



## Plot Fault Label



### Rotate Label
##LHS_dx = plot_flt_LHS[-1,0] - plot_flt_LHS[0,0]
##LHS_dy = plot_flt_LHS[-1,1] - plot_flt_LHS[0,1]
##LHS_angle = np.rad2deg(np.arctan(LHS_dx/LHS_dy))

##LHS_label_anc = plot_flt_LHS[-1,:] - 0.1
##RHS_label_anc = plot_flt_RHS[-1,:] + 0.1

##trans_angle = plt.gca().transData.transform_angles(np.array((90.-LHS_angle + 2,)),
##                                                   LHS_label_anc.reshape((1, 2)))[0]



##LHS_label = mpl.text.Text(plot_flt_LHS[0,0], plot_flt_LHS[0,1], 'text rotated correctly', fontsize=16,
##               rotation=trans_angle, rotation_mode='anchor',ha='right')
##               
##ax.add_artist(LHS_label)       
###########################################
#ax.set_xlim([xmin,xmax])
#ax.set_ylim([ymin,ymax])
#ax.set_aspect('equal')




#################################
## Add legend
#legend_items = [eb1_l,eb2_l,eb3_l,eb4_l,flt_LHS,flt_RHS]
#legend_text = ['Basement',
#               'Low Permeability FM',
#               'High Permeability FM',
#               'Mudrock',
#               'LHS Fault',
#               'RHS Fault']
#plt.legend(legend_items,legend_text,fontsize=6,loc=4)
##plt.legend()
## Label figure
#ax.set_xlabel('Easting, km')
#ax.set_ylabel('Depth, km')

## Final figure parameters
#fig.set_figwidth(10)
#fig.set_figheight(10)
#plt.tight_layout()
#plt.savefig(mesh_name + '.pdf',orientation='landscape')
##plt.show()    
#    
#################################################################################
## Plot Material Property Table
##-------------------------------------------------------------------------------
#fig = []; ax = []
##fig, ax = plt.subplots()
#col_label = ['Basement','Low Perm FM', 'Hi Perm FM','Mudrock'  ]
#row_label = ['k, m**2' ,'mu f, Pa*s' , 'Vp, m/s'   , 'Vs, m/s',
#             'rho'     ,'mu s, Pa*s' , 'B'         ,'r'       ,
#             'cf'      , 'Fluid Mobility, m**2 / Pa / s','phi, %',
#             'Source/Sink Delta P, Pa']



#tab_dat = np.array([perm_k,mu_fluid,vp,vs,rho,mu_solid,B,r,cf,K,phi,source])
#table = pd.DataFrame(columns=col_label,index=row_label,data=tab_dat[:,0:4])
#table.to_latex(buf='mat_props.tex',float_format=np.str)

#################################################################################
## Simulation Parameter Table
##-------------------------------------------------------------------------------
#row_label = ['Simulation Type',
#             'Implicit delta T',
#             'Total Sim Time',
#             'Body forces?',
#             '# Faults',]
#            
#             
#             
#             
#             

#################################################################################
## Plot Injection Schedule
##-------------------------------------------------------------------------------






