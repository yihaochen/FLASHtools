#! /bin/env python
import util
import pickle
import h5py
import numpy as np
import os
#dir = '/home/ychen/d9/FLASH4/2015_production_runs/0529_L45_M10_b1_h1/'
dir = './data/'

force_overwrite = False
pickle_path = os.path.join(dir, 'particles_leave_dict.pickle')

# Read the pickle file if it already exists
if os.path.exists(pickle_path):
    print('%s found, unpickling...' % pickle_path)
    read = pickle.load(open(pickle_path, 'rb'))
    print('Loaded %s' % pickle_path)

else:
################################################################################
# Go through all particle files and record when the velocity of each individual particle drops
# below certain threshold. Use combined particle_tadd and particle_tag to identify particles.
# Save the recorded dictionary to a pickle file.
################################################################################
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
        colname = [item[0].decode().strip() for item in h5f['particle names']]
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

        n_sync_particles = val.value.shape[0]

        if 'type' in colname:
            itype = colname.index('type')
            mask = val.value[:,itype] == 1.0
            n_sync_particles = sum(mask)
        else:
            mask = True

        # Go through the list of particles
        for part in val.value[mask,:]:
            # Skip this particle if already recorded
            if (part[tag], part[tadd]) in particles_leave:
                continue
            vx = part[velx]
            vy = part[vely]
            vz = part[velz]
            vel_magnitude = np.sqrt(vx*vx+vy*vy+vz*vz)

            # Record new particles that just leave the jet
            if vel_magnitude < v_thres:
                # Make sure we have positive density and tau (There might be a better way to do this...)
                if part[dens] > 0.0 and part[tau] > 0.0:
                    particles_leave[(part[tag], part[tadd])] = [part[field] for field in fields] + [sim_time]
        print(f.filename, len(particles_leave), n_sync_particles)
        if len(particles_leave) > 0 and len(particles_leave) == n_sync_particles:
            break


    #particles_leave = np.array([[key[0] + key[1]] + value for key, value in particles_leave.iteritems()])
    #with open('%s_particles_leave_dict.pickle' % dir.split('/')[-2], 'wb') as f:
    with open(pickle_path, 'wb') as f:
        pickle.dump(particles_leave, f, protocol=pickle.HIGHEST_PROTOCOL)

    read = particles_leave



def rescan(printlist=False):
    files = util.scan_files(dir, '*hdf5_part_[0-9][0-9][0-9][0-9]', printlist=printlist)
    return files
files = rescan(True)

################################################################################
# Write den1 and dtau to the updated particle fields as additional HD5 columns.
################################################################################
for f in files[:]:
    h5f  = h5py.File(f.fullpath, 'r')
    if os.path.isfile(f.fullpath+'_updated') and not force_overwrite:
        continue
    h5fr = h5py.File(f.fullpath+'_updated', 'w')

    colname = [item[0].decode().strip() for item in h5f['particle names']]
    itag  = colname.index('tag')
    itadd = colname.index('tadd')
    iden0 = colname.index('den0')
    ivelz = colname.index('velz')
    itau  = colname.index('tau')
    ishok = colname.index('shok')
    if 'type' in colname:
        itype = colname.index('type')
        itype = colname.index('type')
        mask = val.value[:,itype] == 1.0
    else:
        mask = True

    tp = h5f['tracer particles'].value[mask,:]

    dtau = np.zeros(tp.shape[0])
    den1 = tp[:,iden0]

    for i, (tag, tadd)  in enumerate(zip(tp[:,itag], tp[:,itadd])):
        try:
            den1[i] = read[(tag, tadd)][2]
            tau1 = read[(tag, tadd)][3]
            dtau[i] = max(tp[i,itau]-tau1, 1E-100)
        except KeyError:
            print('tag: %6i, tadd: %9.3e, velz: %9.2e, shok: %1i not in pickled data' % (tag, tadd, tp[i,ivelz], tp[i,ishok]))
        except:
            print(den1[i])
            print(read)
            print(tag, tadd)

    newcols = np.column_stack((den1,dtau))

    # Go through each dataset in the particle hdf5 file
    # Add new particle fields (den1 and dtau)
    # Copy all other datasets to the new file
    for v in h5f.values():
        if 'particle names' in v.name:
            shape = (v.shape[0]+2, 1)
            newfields = [[b'den1                    '], [b'dtau                    ']]
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
    print('Saving', h5fr.filename)
    h5fr.close()
