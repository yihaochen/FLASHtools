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
    zvec = np.array([ds.parameters['nozzlevecx'],\
                     ds.parameters['nozzlevecy'],\
                     ds.parameters['nozzlevecz']])
    rxvec = np.cross(zvec, np.array([0.,1E-10,1.]))
    rxvec /= np.sqrt(sum(rxvec**2))
    ryvec = np.cross(zvec, rxvec)

    r = ds.parameters['nozzleradius']+ds.parameters['rfeatherout']
    l = ds.parameters['nozzlehalfl']

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



    planeDict = { 'x': (1,2), 'y': (2,0), 'z': (0,1) }
    idim = planeDict[proj_axis]

    nozzlePos = [ds.parameters['nozzleposx'],\
                 ds.parameters['nozzleposy'],\
                 ds.parameters['nozzleposz']]


    nozzleCoords = [([ corner[idim[0]] + nozzlePos[idim[0]],\
                       corner[idim[1]] + nozzlePos[idim[1]] ], size)\
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
