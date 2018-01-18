import numpy as np
from os import path

def read_par(dir, parfile='flash.par'):
    pars = {}
    with open(path.join(dir,parfile)) as f:
        for line in f.readlines():
            if line.startswith('#') or line.startswith('\n'): continue
            name, value = line.split('=')[:2]
            try:
                pars[name.strip()] = float(value.lstrip(' ').split('#')[0].rstrip())
            except:
                pars[name.strip()] = value.lstrip(' ').split('#')[0].rstrip().strip('"')
    return pars


def calcNozzleCoords(ds, proj_axis):
    # Axial direction of the jet nozzle
    # TODO: Not functional now, since the vectors are the inital runtime parameters, not scalars in FLASH
    zvec = np.array([ds.parameters['nozzlevecx'],\
                     ds.parameters['nozzlevecy'],\
                     ds.parameters['nozzlevecz']])
    # Radial vectors of the jet nozzle
    rxvec = np.cross(zvec, np.array([0.,1E-10,1.]))
    rxvec /= np.sqrt(sum(rxvec**2))
    ryvec = np.cross(zvec, rxvec)

    r = ds.parameters['nozzleradius']+ds.parameters['rfeatherout']
    l = ds.parameters['nozzlehalfl']

    # Number of points in the ring for annotating the nozzle
    # Total number of points will be 3*nPnt
    nPnt = 12
    xx = [np.cos(th) for th in np.linspace(0,2*np.pi, nPnt, endpoint=False)]
    yy = [np.sin(th) for th in np.linspace(0,2*np.pi, nPnt, endpoint=False)]

    # Position vector of the ring points from the center of the ring
    nozzleCorners = np.zeros([3*nPnt,3])
    for i in range(nPnt):
        # Points at the upper edge of the nozzle
        nozzleCorners[i,:] =        r*(xx[i]*rxvec + yy[i]*ryvec) + l*zvec
        # Points at the central cross-section of the nozzle
        nozzleCorners[nPnt+i,:] =   r*(xx[i]*rxvec + yy[i]*ryvec)
        # Points at the lower edge of the nozzle
        nozzleCorners[2*nPnt+i,:] = r*(xx[i]*rxvec + yy[i]*ryvec) - l*zvec

    axisDict = { 'x': 0, 'y': 1, 'z': 2 }
    iaxis = axisDict[proj_axis]
    sizes = nozzleCorners[:,iaxis]/(max(nozzleCorners[:,iaxis])-min(nozzleCorners[:,iaxis]))

    planeDict = { 'x': (1,2), 'y': (2,0), 'z': (0,1) }
    idim = planeDict[proj_axis]

    nozzlePos = [ds.parameters['nozzleposx'],\
                 ds.parameters['nozzleposy'],\
                 ds.parameters['nozzleposz']]

    #nozzleCoords = [([ corner[idim[0]] + nozzlePos[idim[0]],\
    #                   corner[idim[1]] + nozzlePos[idim[1]] ], size)\
    #               for corner, size in zip(nozzleCorners[sizes.argsort()], sorted(sizes)) ]
    nozzleCoords = [(corner + nozzlePos, size)\
                   for corner, size in zip(nozzleCorners[sizes.argsort()], sorted(sizes)) ]

    return nozzleCoords

def calcNozzleCoords_from_pars(pars, t, proj_axis='x'):
    vecX = pars['nozzleVecX']
    vecY = pars['nozzleVecY']
    vecZ = pars['nozzleVecZ']

    norm = np.sqrt(vecX**2+vecY**2+vecZ**2)
    vec0 = np.array([vecX, vecY, vecZ])/norm
    #print vec0

    # See Wikipedia "Rotation formalisms in three dimensions" for the details of the rotation.
    # http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Rotation_matrix_.E2.86.94_Euler_axis.2Fangle

    if 'nozzleAngVelX' in pars:
        angVel = np.array([ pars['nozzleAngVelX'], pars['nozzleAngVelY'], pars['nozzleAngVelZ'] ])
    else:
        angVel = np.zeros(3)
    omega = np.sqrt(sum(angVel*angVel))

    if omega == 0.0:
        zvec = vec0
    else:
        theta = omega*t
        e = angVel/omega    # unit vector of the rotation axis
        e1 = e[0]
        e2 = e[1]
        e3 = e[2]
        cos = np.cos(theta)
        sin = np.sin(theta)
        A = np.array([ [(1-cos)*e1*e1+cos,    (1-cos)*e1*e2-e3*sin, (1-cos)*e1*e3+e2*sin],\
                       [(1-cos)*e2*e1+e3*sin, (1-cos)*e2*e2+cos,    (1-cos)*e2*e3-e1*sin],\
                       [(1-cos)*e3*e1-e2*sin, (1-cos)*e3*e2+e1*sin, (1-cos)*e3*e3+cos] ])
        zvec = np.dot(A, vec0)

    rxvec = np.cross(zvec, np.array([0.,1E-10,1.]))
    rxvec /= np.sqrt(sum(rxvec**2))
    ryvec = np.cross(zvec, rxvec)
    #print rxvec,ryvec,zvec

    r = pars['nozzleRadius']+pars['rFeatherOut']
    l = pars['nozzleHalfL']

    nPnt = 12
    xx = [np.cos(th) for th in np.linspace(0,2*np.pi, nPnt, endpoint=False)]
    yy = [np.sin(th) for th in np.linspace(0,2*np.pi, nPnt, endpoint=False)]
    nozzleCorners = np.zeros([3*nPnt,3])
    for i in range(nPnt):
        nozzleCorners[i,:] =        r*(xx[i]*rxvec + yy[i]*ryvec) + l*zvec
        nozzleCorners[nPnt+i,:] =   r*(xx[i]*rxvec + yy[i]*ryvec)
        nozzleCorners[2*nPnt+i,:] = r*(xx[i]*rxvec + yy[i]*ryvec) - l*zvec

    axisDict = { 'x': 0, 'y': 1, 'z': 2 }
    iaxis = axisDict[proj_axis]
    sizes = nozzleCorners[:,iaxis]/(max(nozzleCorners[:,iaxis])-min(nozzleCorners[:,iaxis]))


    #nozzleEdges = np.zeros([3,2])
    #nozzleCorners = np.zeros([8,3])
    #nozzleEdges[0] = np.array([-nR, nR])
    #nozzleEdges[1] = np.array([-nR, nR])
    #nozzleEdges[2] = np.array([-nL, nL])


    #i = 0
    #for x in nozzleEdges[0]:
    #    for y in nozzleEdges[1]:
    #        for z in nozzleEdges[2]:
    #            nozzleCorners[i] = x*xvec + y*yvec + z*zvec
    #            i += 1

    #for (i, pos) in enumerate(nozzleCorners):
    #    nozzleCorners[i] = np.dot(A, pos)

    #R = np.array([ [vec[2], vec[1]], [-vec[1], vec[2]] ])

    planeDict = { 'x': (1,2), 'y': (2,0), 'z': (0,1) }
    idim = planeDict[proj_axis]

    nozzlePos = [pars['nozzlePosX'] + pars['nozzleLinVelX']*t,\
                 pars['nozzlePosY'] + pars['nozzleLinVelY']*t,\
                 pars['nozzlePosZ'] + pars['nozzleLinVelZ']*t]


    nozzleCoords = [([ corner[idim[0]] + nozzlePos[idim[0]],\
                       corner[idim[1]] + nozzlePos[idim[1]] ], size)\
                   for corner, size in zip(nozzleCorners[sizes.argsort()], sorted(sizes)) ]

    #for corner in nozzleCorners:
    #    print corner
    #for coord in nozzleCoords:
    #    print coord

    return nozzleCoords

def calcDen0_2015(tadd, t1=1.58E13, v=3E9, g=1.33333, r=7.5E20, bf=1.875E20, initM=5, mach=10):
    M = initM+(mach-initM)*np.cos(0.5*np.pi*np.clip((tadd-t1)/t1, -1.0, 0.0))
    den0 = 0.5*1E45/np.pi/v**3/( 0.5*r*r*(1.+1./M**2/(g-1.)) + r*bf*(0.3125+1./M**2/(g-1.))\
                             + bf*bf*(0.06056+0.29736/M**2/(g-1.)) )
    return den0

def calcDen0(data, ptype='io'):
    """
    Calculate the density in the core of the jet according to simulation parameters.
    """
    p = data.ds.parameters
    tadd = data[ptype, 'particle_tadd']
    h = p['sim_helicityjet']
    g = p['sim_gammajet']
    R = p['nozzleradius']
    bf = p['rfeatherout']
    initM = p['sim_initmachjet']
    M = p['sim_machjet']
    t1 = p['sim_duration']/100.0
    x = g/(g-1.0)+h**2/(1.0+h**2)/p['sim_betajet']
    M = initM+(M-initM)*np.cos(0.5*np.pi*np.clip((tadd-t1)/t1, -1.0, 0.0))

    # New simulations in 2016 has 'particle_type' field
    if data.has_key((ptype, 'particle_type')):
        den0 = 0.5*p['sim_powerjet']/np.pi/p['sim_veljet']**3/( R*R*(0.5+x/M**2/g) + R*bf*(0.3125+x/M**2/g)\
                             + bf*bf*(0.06056+0.29736*x/M**2/g) )
    # For simulations in 2015 (did not include magnetic power)
    else:
        den0 = 0.5*p['sim_powerjet']/np.pi/p['sim_veljet']**3/\
               ( 0.5*R*R*(1.+1./M**2/(g-1.)) + R*bf*(0.3125+1./M**2/(g-1.)) + bf*bf*(0.06056+0.29736/M**2/(g-1.)) )

    return den0

def setup_cl(dirs):
    # Set up colors and label names
    colors = {}
    labels = {}
    for dirname in dirs:
        if 'M3_h1' in dirname:
            colors[dirname] = 'pink'
            labels[dirname] = 'low Mach (3)'
        elif 'h1' in dirname:
            colors[dirname] = 'r'   #'#fc8d62'
            labels[dirname] = 'helical'
        elif 'h0' in dirname:
            colors[dirname] = 'b'   #'#8da0cb'
            labels[dirname] = 'poloidal'
        elif 'hinf' in dirname:
            colors[dirname] = 'g'   #'#66c2a5'
            labels[dirname] = 'toroidal'
        elif 'hydro' in dirname:
            colors[dirname] = 'k'
            labels[dirname] = 'hydro'
        elif 'M24_b01' in dirname:
            colors[dirname] = 'purple'
            labels[dirname] = 'low beta (0.01)'
    return colors, labels
