'''
Module for creating a proper model for QM/MM calculations by cropping the water box.
It saves the 'pdb' of the croped model, as well as the topology and parameters and initial coordinates.
'''

from MDAnalysis.exceptions import SelectionError

#argsdict=dict({'trajectory': ['3rde.prmtop', None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': 'plots', 'timer': False, 'menu_type' : None, 'u_loaded' : False})

#argsdict=dict({'trajectory': ['3rde.prmtop', None],
#'frame': None, 'latex': False, 'latex_width': None,
#'parallel': False, 'subdir': 'plots', 'timer': False,
#'menu_type' : None, 'u_loaded' : False})

def ask(u, argsdict):

    if argsdict['trajectory'] != None:
        while True:
            try :
                frame = int(input('Of which frame do you want to obtain a model for QM/MM calculations? ')) +1
                print('The selected frame is the %s' % frame)
                from MDAnalysis import Universe
                u_top = Universe(argsdict['trajectory'][0])
                break

            except ValueError:
                print('Type only the number of the frame.')

    elif argsdict['frame'] != None:
        u_top = None; frame = None
        from modules import loader
        u, argsdict = loader.universe_loader_frame(argsdict)


    #while True:
    #    try :
    #        ligand = input('Which is the central residue? (Type the 3-letters code or the number) ')
    #        break
    #    except ValueError:
    #        print('Type only the 3-letters code or the number.')

    while True:
        ligand = input('Which is the central residue? (Type the 3-letters code or the number) ')
        try :
            if u_top != None:
                sel = u_top.select_atoms('resid %s' % ligand)
            elif u_top == None:
                sel = u.select_atoms('resid %s' % ligand)

        except SelectionError:
            try :
                if u_top != None:
                    sel = u_top.select_atoms('resname %s' % ligand)
                elif u_top == None:
                    sel = u.select_atoms('resname %s' % ligand)
            except SelectionError:
                print('The selected residue does not exists or it has not been correctly selected.')
                continue

        locA = str(sel.residues).find('[<') + 2
        locB = str(sel.residues).find('>]>')
        break

    while True:
        print('You have selected this residue: %s' % str(sel.residues)[locA:locB])
        quest = input('Is the residue correct? ([y]/n) ')
        if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
            while True:
                ligand = input('Which is the central residue? (Type the 3-letters code or the number) ')
                try :
                    if u_top != None:
                        sel = u_top.select_atoms('resid %s' % ligand)
                    elif u_top == None:
                        sel = u.select_atoms('resid %s' % ligand)

                except SelectionError:
                    try :
                        if u_top != None:
                            sel = u_top.select_atoms('resname %s' % ligand)
                        elif u_top == None:
                            sel = u.select_atoms('resname %s' % ligand)

                    except SelectionError:
                        print('The selected residue does not exists or it has not been correctly selected.')
                        continue

                locA = str(sel.residues).find('[<') + 2
                locB = str(sel.residues).find('>]>')
                break
            continue
        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            break
        else :
            print('Type yes or not.')
            continue

    while True:
        try :
            radius = float(input('Which radius do you want for the solvent drop? '))
            break
        except ValueError:
            print('Type only the radius (as float).')


    if argsdict['trajectory'] != None and argsdict['u_loaded'] == False:
        from modules import loader
        u, argsdict = loader.universe_loader_traj(argsdict)

    if frame != None and frame > len(u.trajectory):
        print('The selected frame does not exist. Please, select another frame.')
        while True:
            try :
                frame = int(input('Of which frame do you want to obtain a model for QM/MM calculations? ')) +1
                if frame > len(u.trajectory):
                    print('The selected frame does not exist.')
                    continue
                print('The selected frame is the %s' % frame)
                break

            except ValueError:
                print('Type only the number of the frame.')

    return u, argsdict, frame, ligand, radius


def saver(u, argsdict, frame, ligand, radius):
    if 'qmmm_models' not in os.listdir():
        os.mkdir('qmmm_models')
        print('\'qmmm_models\' subdirectory has been created. The created models will be placed here.')
    elif 'qmmm_models' in os.listdir():
        print('Created QM/MM models are going to be placed on \'qmmm_models\' subdirectory.')


    if frame == None:
        selection = u.select_atoms(' not ((resname HOH or resname WAT or resname Na+) and (not around %s resname %s))' % (radius, ligand))
        selection.write('qmmm_models/qmmm_model_%s.pdb' % (str(argsdict['frame'])[:str(argsdict['frame']).find('.pdb')]))

    elif frame != None:
        u.trajectory[frame]
        selection = u.select_atoms(' not ((resname HOH or resname WAT or resname Na+) and (not around %s resname %s))' % (radius, ligand))
        selection.write('qmmm_models/qmmm_model_%s.pdb' % frame)


def topology_modeller_traj(u, argsdict, frame, ligand, radius):
    from parmed import load_file

    topology = load_file(argsdict['trajectory'][0], xyz='qmmm_models/qmmm_model_%s.pdb' % frame)
    cropped = topology[':WAT&:%s<@%s|:Na+&:%s<@%s|:Cl-&:%s<@%s|(!:WAT&!:Na+&!:Cl-)' % (ligand, radius, ligand, radius, ligand, radius)]

    cropped.write_parm('qmmm_models/cropped_%s' % argsdict['trajectory'][0])
    cropped.save('qmmm_models/cropped_%s.inpcrd' % str(argsdict['trajectory'][0])[:str(argsdict['trajectory'][0]).find('.prmtop')])

    print('Cropped topology and coordinates have been saved as \'cropped_%s\' and \'cropped_%s.inpcrd\' in \'qmmm_models\'' % (argsdict['trajectory'][0], str(argsdict['trajectory'][0])[:str(argsdict['trajectory'][0]).find('.prmtop')]))



def topology_modeller_frame(u, argsdict, ligand, radius):
    from parmed import load_file
    print('Since the input structure is a pdb, the topology is not specified.')

    top_file = input('Specify the route to the topology: ')
    topology = load_file(top_file, xyz='qmmm_models/qmmm_model_%s.pdb' % (str(argsdict['frame'])[:str(argsdict['frame']).find(.pdb)]))
    cropped = topology[':WAT&:%s<@%s|:Na+&:%s<@%s|:Cl-&:%s<@%s|(!:WAT&!:Na+&!:Cl-)' % (ligand, radius, ligand, radius, ligand, radius)]

    cropped.write_parm('qmmm_models/cropped_%s' % top_file)
    cropped.save('qmmm_models/cropped_%s.inpcrd' % top_file[:top_file.find('.prmtop')])

    print('Cropped topology and coordinates have been saved as \'cropped_%s\' and \'cropped_%s.inpcrd\' in \'qmmm_models\'' % (top_file, top_file[:top_file.find('.prmtop')]))


def qmmm_modeller(u, argsdict):

    while True:
        u, argsdict, frame, ligand, radius = ask(u, argsdict)
        saver(u, argsdict, frame, ligand, radius)

        while True:
            quest = input('Do you want to generate also the cropped topology and parameters and coordinates files? (y/n) ')
            if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                break
            elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                if argsdict['trajectory'] != None:
                    topology_modeller_traj(u, argsdict, frame, ligand, radius)
                    break
                elif argsdict['frame'] != None:
                    topology_modeller_frame(u, argsdict, ligand, radius)
                    break
                else :
                    print('Answer yes or no.')
                    continue

        if argsdict['trajectory'] != None:
            quest = input('Do you want to save another frame? (y/[n]) ')
            if quest in ('', 'n', 'no', 'N', 'No', 'No', 'nO', '0'):
                break

            elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                continue
