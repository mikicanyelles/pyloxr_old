from MDAnalysis import Merge, Universe
from MDAnalysis.analysis import rms, align
from matplotlib import pyplot as plt
from modules import loader


def ask_mask(u, argsdict):

    rmsddict = {'bb_ref' : None, 'mask' : None, 'subs_ref' : None}

    while True:
        try :
        rmsddict['mask'] = int(input('Do you want to plot the RMSD of the backbone (1), of a/some residue/s (2) or of both (3)? '))
        break
    except TypeError:
        print('Type only the number.')
        continue
    except rmsddict['mask'] not in (1,2,3):
        print('Type only 1, 2 or 3.')
        continue

    return u, argsdict, rmsddict

def ask_resids(u, argsdict):

    while True:
        try :
            index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space or by a colon if a global RMSD is desired). '))
            #mask = u.select_atoms('resid %s' % index)
            for i in range(0, len(index.split(' '))):
                int(index.split()[i])
            break
        except ValueError:
            print('Type only the number.')
            continue

    u_top = Universe(argsdict['parameters'])
    quest = None

    while True:
        print('You have selected this/these residue/s:')
        i = 0; exists = True
        while i < len(index.split()):
            try :
                print(str(u_top.select_atoms('resid %s' % index).residues[i])[1:-1])
                i+=1
            except IndexError:
                print('Residue number %s doesn\'t exist' % index.split()[i])
                exists = False
                break

        if exists == False:
            while True:
                try :
                    index = input('Please, retype the numbers. ')
                    #mask = u.select_atoms('resid %s' % index)
                    for i in range(0, len(index.split())):
                        int(index.split()[i])
                    break
                except ValueError:
                    print('Type only the number.')
                    continue
            continue
        elif exists != False:
            quest = input('Are the selected residues correct ([y]/n)? ')
            if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                break
            if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                while True:
                    try :
                        index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space). '))
                        #mask = u.select_atoms('resid %s' % index)
                        for i in range(0, len(index.split())):
                            int(index.split()[i])
                        break
                    except ValueError:
                        print('Type only the number.')
                        continue
                continue





def mask1(u, argsdict, rmsddict):
    while True:
        try :
            rmsddict['bb_ref'] = int(input('Do you want to set the first frame as the reference (1) or the average structure (2) for the backbone? '))
            break
        except TypeError:
            print('Type only the number.')
            continue
        except rmsddict['bb_ref'] not in (1,2):
            print('Type only 1 or 2.')
            continue

    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)


    time = []
    for i in u.trajectory:
        time.append(int(u.trajectory.time)/1000)
    mask = u.select_atoms('backbone')
    ref = ref_u.select_atoms('backbone')


    if rmsddict['bb_ref'] == 1:
        # creation of the array which has the time, the frame and the RMSD value
        R = rms.RMSD(mask, ref).run()
        rmsd_arr = R.rmsd.T


    elif rmsddict['bb_ref'] == 2:
        # creation of the array which has the time, the frame and the RMSD value
        prealigner = align.AlignTraj(u, reference=ref_u, select='backbone', in_memory=False).run()
        reference_coords = u.trajectory.timeseries(asel=mask).mean(axis=1)
        avg = Merge(mask).load_new(reference_coords[:, None, :], order='afc')

        R = rms.RMSD(mask, avg).run()
        rmsd_arr = R.rmsd.T





def rmsd(u, argsdict):

