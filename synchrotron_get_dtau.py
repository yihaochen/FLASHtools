#! /bin/env python
import util
import pickle
import h5py
import numpy as np
import os
#dir = '/home/ychen/d9/FLASH4/2015_production_runs/0529_L45_M10_b1_h1/'
dir = './'

pickle_path = os.path.join(dir, 'particles_leave_dict.pickle')

if os.path.exists(pickle_path):
    print '%s found, unpickling...' % pickle_path
    read = pickle.load(open(pickle_path, 'r'))
    print 'Loaded %s' % pickle_path

else:

    def rescan(printlist=False):
        files = util.scan_files(dir, '*hdf5_part_[0-9][0-9][0-9][0-9]', printlist=printlist)
        return files
    files = rescan(True)

    # Velocity threshold
    v_thres = 0.05*3E10

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

        sim_time = h5f['real scalars'][0][1]

        # Fields to be written in the output file
        fields = [tadd, tag, dens, tau]

        val = h5f['tracer particles']

        # Go through the list of particles
        for part in val.value:
            # Skip this particle if already recorded
            if (part[tag], part[tadd]) in particles_leave:
                continue
            vx = part[velx]
            vy = part[vely]
            vz = part[velz]
            vel_magnitude = np.sqrt(vx*vx+vy*vy+vz*vz)

            #if part[tag] == 127225.0 and np.isclose(part[tadd], 218857882571403.5):
            #    print sim_time, '%.3e' % (vel_magnitude/v_thres), part[dens], part[tau]

            # Record new particles that just leave the jet
            if vel_magnitude < v_thres:
                # Make sure we have positive density and tau (There might be a better way to do this...)
                if part[dens] > 0.0 and part[tau] > 0.0:
                    particles_leave[(part[tag], part[tadd])] = [part[field] for field in fields] + [sim_time]
        print f.filename, len(particles_leave), val.value.shape[0]


    #particles_leave = np.array([[key[0] + key[1]] + value for key, value in particles_leave.iteritems()])
    #with open('%s_particles_leave_dict.pickle' % dir.split('/')[-2], 'wb') as f:
    with open('particles_leave_dict.pickle', 'wb') as f:
        pickle.dump(particles_leave, f, protocol=pickle.HIGHEST_PROTOCOL)

    read = particles_leave



def rescan(printlist=False):
    files = util.scan_files(dir, '*hdf5_part_[0-9][0-9][0-9]0', printlist=printlist)
    return files
files = rescan(True)

# Write den1 and dtau to the updated particle fields
for f in files[:]:
    h5f  = h5py.File(f.fullpath, 'r')
    h5fr = h5py.File(f.fullpath+'_updated', 'w')

    colname = [item[0].strip() for item in h5f['particle names']]
    itag  = colname.index('tag')
    itadd = colname.index('tadd')
    iden0 = colname.index('den0')
    ivelz = colname.index('velz')
    itau  = colname.index('tau')
    tp = h5f['tracer particles'].value

    dtau = np.zeros(tp.shape[0])
    den1 = tp[:,iden0]

    for i, (tag, tadd)  in enumerate(zip(tp[:,itag], tp[:,itadd])):
        try:
            den1[i] = read[(tag, tadd)][2]
            tau1 = read[(tag, tadd)][3]
            dtau[i] = max(tp[i,itau]-tau1, 1E-100)
        except KeyError:
            print tag, tadd, '%e' % tp[i,ivelz], 'not in pickled data'
        except:
            print den1[i]
            print read
            print tag, tadd
    newcols = np.column_stack((den1,dtau))

    # Go through each dataset in the particle hdf5 file
    # Add new particle fields (den1 and dtau)
    # Copy all other datasets to the new file
    for v in h5f.values():
        if 'particle names' in v.name:
            shape = (v.shape[0]+2, 1)
            newfields = [['den1                    '], ['dtau                    ']]
            data = np.concatenate((v.value, newfields), axis=0)
            h5fr.create_dataset(v.name, shape, v.dtype, data)
        elif 'tracer particles' in v.name:
            shape = (v.shape[0], v.shape[1]+2)
            #newcol = np.expand_dims(v.value[:,0], axis=1)
            #newcols = v.value[:,0:2]
            data = np.concatenate((v.value, newcols), axis=1)
            h5fr.create_dataset(v.name, shape, v.dtype, data)
        else:
            h5fr.create_dataset(v.name, v.shape, v.dtype, v.value)
    h5f.close()
    print 'Saving', h5fr.filename
    h5fr.close()
