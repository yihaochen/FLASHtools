import util
import pickle
import h5py
import numpy as np
dir = '/home/ychen/d9/FLASH4/2015_production_runs/0529_L45_M10_b1_h1/'

def rescan(printlist=False):
    files = util.scan_files(dir, '*hdf5_part_[0-9][0-9][0-9][0-9]', printlist=printlist)
    return files
files = rescan(True)

# Velocity threshold
v_thres = 0.01*3E10

# Particles that leave the jet
particles_leave = {}

for f in files[:]:
    h5f = h5py.File(f.fullpath, 'r')
    colname = [item[0].strip() for item in h5f['particle names']]
    tag  = colname.index('tag')
    tadd = colname.index('tadd')
    dens = colname.index('dens')
    tau  = colname.index('tau')
    shok = colname.index('shok')
    velx = colname.index('velx')
    vely = colname.index('vely')
    velz = colname.index('velz')


    # Fields to be written in the output file
    fields = [tadd, tag, dens, tau]

    val = h5f['tracer particles']

    # Go through the list of particles
    for part in val.value:
        vx = part[velx]
        vy = part[vely]
        vz = part[velz]
        vel_magnitude = np.sqrt(vx*vx+vy*vy+vz*vz)
        #if part[shok] > 0.01:
        #    if part[shok] < 1.0:
        #        print par[shok]
        #    next

        # Record new particles that just leave the jet
        if (part[tag], part[tadd]) not in particles_leave and vel_magnitude < v_thres:
            particles_leave[(part[tag], part[tadd])] = [part[field] for field in fields]

    # Collect particle values to a list of dictionaries
    #particles = []
    #for value_particle in val.value:
    #    particle = dict(zip(colname, value_particle))
    #    particles.append(particle)

    # Go through the list of particles
    #for part in particles:
    #    newtag = part['tadd']+part['tag']
    #    if part['shok'] > 0.01:
    #        if part['shok'] < 1.0:
    #            print part['shok']
    #        next
        # Record new particles that just leave the jet
    #    if newtag not in particles_leave and part['velz'] < v_thres:
    #        particles_leave[newtag] = [part[field] for field in fields]
    #print f, '%6i : %6i' % (len(particles), len(particles_leave))

particles_leave = np.array([[key[0] + key[1]] + value for key, value in particles_leave.iteritems()])
pickle.dump(particles_leave, open('%s_particles_leave.pickle' % dir.split('/')[-2], 'wb'))
