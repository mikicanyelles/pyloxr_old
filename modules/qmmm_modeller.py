'''
Module for creating a proper model for QM/MM calculations by cropping the water box.
It saves the 'pdb' of the cropped model, as well as the topology and parameters and initial coordinates.
Moreover, it can save also the list of active atoms for ChemShell QM/MM calculations.
'''

from MDAnalysis.exceptions import SelectionError
import os
import sys


#argsdict=dict({'trajectory': ['3rde.prmtop', None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': 'plots', 'timer': False, 'menu_type' : None, 'u_loaded' : False})

#argsdict=dict({'trajectory': ['3rde.prmtop', None],
#'frame': None, 'latex': False, 'latex_width': None,
#'parallel': False, 'subdir': 'plots', 'timer': False,
#'menu_type' : None, 'u_loaded' : False})

def ask(u, argsdict):

    if argsdict['trajectory'] != None:
        while True:
            try :
                frame = int(input('Of which frame do you want to obtain a model for QM/MM calculations? '))
                print('The selected frame is the %s' % (frame))
                from MDAnalysis import Universe
                u_top = Universe(argsdict['parameters'])
                break

            except ValueError:
                print('Type only the number of the frame.')

    elif argsdict['frame'] != None:
        u_top = None; frame = None
        from modules import loader
        u, argsdict = loader.universe_loader_frame(argsdict)

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
            radius = float(input('Which radius do you want for the solvent drop (in Å)? '))
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

    stderr_ = open(os.devnull, 'w')
    sys.stderr = stderr_

    if frame == None:
        #selection = u.select_atoms('protein or (byres around %s resname %s)' % (radius, ligand))
        selection = u.select_atoms('byres not ((resname HOH or resname WAT or resname Na+ or resname Cl-) and (not (around %s resname %s)))' % (radius, ligand))
        selection.write('qmmm_models/qmmm_model_%s.pdb' % (str(argsdict['frame'])[:str(argsdict['frame']).find('.pdb')]))

    elif frame != None:
        u.trajectory[frame-1]
#        selection = u.select_atoms('protein or (byres around %s resname %s)' % (radius, ligand))
        selection = u.select_atoms('byres not ((resname HOH or resname WAT or resname Na+ or resname Cl-) and (not (around %s resname %s)))' % (radius, ligand))
        selection.write('qmmm_models/qmmm_model_%s.pdb' % frame)

    stderr_.close()


def topology_modeller_traj(u, argsdict, frame, ligand, radius):
    from parmed import load_file

    stderr_ = open(os.devnull, 'w')
    sys.stderr = stderr_

    u.trajectory[frame-1]
    u.select_atoms('all').write('.temporal.pdb')

    stderr_.close()
    del stderr_


    if len(argsdict['parameters'].split('/')) != 1:
        name_prmtop = 'qmmm_models/cropped_%s_%s_%s.prmtop' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], frame, radius)
        name_inpcrd = 'qmmm_models/cropped_%s_%s_%s.inpcrd' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], frame, radius)
    elif len(argsdict['parameters'].split('\\')) != 1:
        name_prmtop = 'qmmm_models\\cropped_%s_%s_%s.prmtop' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], frame, radius)
        name_inpcrd = 'qmmm_models\\cropped_%s_%s_%s.inpcrd' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], frame, radius)

    if name_prmtop[12:] in os.listdir('qmmm_models') or name_inpcrd[12:] in os.listdir('qmmm_models'):
        while True:
            try :
                quest = input('%s and %s exist. Do you want to overwrite them (1) or to give an alternative name (2) ([1]/2)? ' % (name_prmtop[12:], name_inpcrd[12:]))
                if quest in ('', '1'):
                    try :
                        os.remove(name_prmtop)
                    except FileNotFoundError:
                        pass

                    try :
                        os.remove(name_inpcrd)
                    except FileNotFoundError:
                        pass

                    print('Files have been removed.')
                    break

                elif quest == '2':
                    name = input('Type the alternative name for the \'.prmtop\' and \'.inpcrd\' files (do not include the extension): ')
                    print('Parameters and topology and coordinates will be saved as ' + name)
                    break

                else :
                    print('Type only 1 or 2')
                    continue

            except ValueError:
                print('Type only 1 or 2')
                continue
    else :
        pass

    topology = load_file(argsdict['parameters'], xyz='.temporal.pdb')
    topology.box = None
<<<<<<< HEAD
    topology.strip(':WAT&!(:%s<:%s)' % (ligand, radius))
    topology.strip(':HOH&!(:%s<:%s)' % (ligand, radius))
    topology.strip(':Na+&!(:%s<:%s)' % (ligand, radius))
    topology.strip(':Cl-&!(:%s<:%s)' % (ligand, radius))
=======
    topology.strip(':WAT&!:%s<:%s' % (ligand, radius))
    topology.strip(':Na+&!:%s<:%s' % (ligand, radius))
    topology.strip(':Cl-&!:%s<:%s' % (ligand, radius))
>>>>>>> c47a9c77e8efe2a191871660ed452b3d835c413a

    topology.write_parm(name_prmtop)
    topology.save(name_inpcrd)
    #topology.save(name_inpcrd[:-6] + 'pdb')


    print('Cropped topology and coordinates have been saved as \'%s\' and \'%s\' in \'qmmm_models\'' % (name_prmtop[12:], name_inpcrd[12:]))

    #os.remove('.temporal.pdb')

    return name_prmtop, name_inpcrd


def topology_modeller_frame(u, argsdict, frame, ligand, radius):
    from parmed import load_file
    print('Since the input structure is a pdb, the topology is not specified.')
    top_file = input('Specify the route to the topology: ')


    #    if len(argsdict['trajectory'][0].split('/')) != 1:
    #        name_prmtop = 'qmmm_models/cropped_%s_%s.prmtop' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], radius)
    #        name_inpcrd = 'qmmm_models/cropped_%s_%s.inpcrd' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], radius)
    #        name_prmtop = 'qmmm_models\\cropped_%s_%s.prmtop' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], radius)
    #        name_inpcrd = 'qmmm_models\\cropped_%s_%s.inpcrd' % (str(argsdict['parameters'].split('/')[-1])[:str(argsdict['parameters'].split('/')[-1]).find('.prmtop')], radius)

    if len(top_file.split('/')) != 1:
        name_prmtop = 'qmmm_models/cropped_%s_%s.prmtop' % (str(top_file.split('/')[-1])[:str(top_file.split('/')[-1]).find('.prmtop')], radius)
        name_inpcrd = 'qmmm_models/cropped_%s_%s.inpcrd' % (str(top_file.split('/')[-1])[:str(top_file.split('/')[-1]).find('.prmtop')], radius)
    elif len(top_file.split('\\')) != 1:
        name_prmtop = 'qmmm_models\\cropped_%s_%s.prmtop' % (str(top_file.split('\\')[-1])[:str(top_file.split('\\')[-1]).find('.prmtop')], radius)
        name_inpcrd = 'qmmm_models\\cropped_%s_%s.inpcrd' % (str(top_file.split('\\')[-1])[:str(top_file.split('\\')[-1]).find('.prmtop')], radius)

    elif top_file.find('/') == -1:
        name_prmtop = 'qmmm_models/cropped_%s_%s.prmtop' % (str(top_file)[:top_file.find('.prmtop')], radius)
        name_inpcrd = 'qmmm_models/cropped_%s_%s.inpcrd' % (str(top_file)[:top_file.find('.prmtop')], radius)

    elif top_file.find('\\') == -1:
        name_prmtop = 'qmmm_models\\cropped_%s_%s.prmtop' % (str(top_file)[:top_file.find('.prmtop')], radius)
        name_inpcrd = 'qmmm_models\\cropped_%s_%s.inpcrd' % (str(top_file)[:top_file.find('.prmtop')], radius)


    if name_prmtop[12:] in os.listdir('qmmm_models') or name_inpcrd[12:] in os.listdir('qmmm_models'):
        while True:
            try :
                quest = int(input('%s and %s exist. Do you want to overwrite them (1) or to give an alternative name (2) ([1]/2)? ' % (name_prmtop[12:], name_inpcrd[12:])))
                if quest == 1:
                    try :
                        os.remove(name_prmtop)
                    except FileNotFoundError:
                        pass

                    try :
                        os.remove(name_inpcrd)
                    except FileNotFoundError:
                        pass

                    print('Files have been removed.')
                    break

                elif quest == 2:
                    name = input('Type the alternative name for the \'.prmtop\' and \'.inpcrd\' files (do not include the extension): ')
                    print('Parameters and topology and coordinates will be saved as ' + name)
                    break

                else :
                    print('Type only 1 or 2')
                    continue

            except ValueError:
                print('Type only 1 or 2')
                continue
    else :
        pass

    topology = load_file(top_file, xyz=argsdict['frame'])
    topology.box = None
    topology.strip(':WAT&:%s<@%s' % (ligand, radius))
    topology.strip(':Na+&:%s<@%s' % (ligand, radius))
    topology.strip(':Cl-&:%s<@%s' % (ligand, radius))

    topology.write_parm(name_prmtop)
    topology.save(name_inpcrd)

    print('Cropped topology and coordinates have been saved as \'%s\' and \'%s\' in \'qmmm_models\'' % (name_prmtop[12:], name_inpcrd[12:]))

    return name_prmtop, name_inpcrd


def topology_modeller(u, argsdict, frame, ligand, radius):

    #try:
    #    import parmed as pmd
    #    del pmd
    #except ImportError:
    #    print('ParmEd module is not installed. If you wish to crop also the topology and parameters and coordinates files, please, install it.')
    #    quest = None
    #    name_inpcrd = None
    #    name_prmtop = None
    
    quest = None
    if quest != None:
        while True:
            quest = input('Do you want to generate also the cropped topology and parameters and coordinates files? (y/n) ')
            if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                name_inpcrd = None
                name_prmtop = None
                break

            elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                if argsdict['trajectory'] != argsdict['parameters']:
                    name_prmtop, name_inpcrd = topology_modeller_traj(u, argsdict, frame, ligand, radius)
                    break
                elif argsdict['trajectory'] == argsdict['parameters']:
                    name_prmtop, name_inpcrd = topology_modeller_frame(u, argsdict, frame, ligand, radius)
                    break

            else :
                print('Answer \'yes\' or \'no\'.')

                continue

    return name_prmtop, name_inpcrd


def set_act(name_prmtop, name_inpcrd, radius):
    from MDAnalysis import Universe

    u_set_act = Universe(name_prmtop, name_inpcrd)

    while True:
        try :
            carbon = input("Type the number of the central atom: ")
            sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
            locA = sel_carbon.find('[<') + 2
            locB = sel_carbon.find(' and segid')
            break
        except SelectionError:
            print('The selection does not exists. Please, type an atom that exists.')
            continue



    while True:
        try :
            quest = input('You selected this atom: %s. Is it correct ([y]/n)? ' % sel_carbon[locA:locB])

            if quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                while True:
                    try:
                        carbon = input("Type the number of the central carbon: ")
                        sel_carbon = str(u_set_act.select_atoms("bynum %s" % carbon))
                        locA = sel_carbon.find('[<') + 2
                        locB = sel_carbon.find(' and segid')
                        break
                    except SelectionError:
                        continue
                continue

            elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                break

        except ValueError:
            print("Type just \'yes\' or \'no\'.")

    while True:
        try :
            radius_set_act = float(input("Which is the desired radius (in Å)? "))
            if radius_set_act < radius:
                break
            elif radius_set_act >= radius:
                print('The selected radius for the set_act list is smaller than the radius of solvent. Please, choose a new radius for the active atoms.')
                continue
        except ValueError:
            print("Type just the number, please.")
            continue

    selection = u_set_act.select_atoms(str('byres around %s bynum %s' % (radius_set_act, carbon)))
    txt = open('qmmm_models/set_act_%s_%s' % (name_prmtop[12:-6], radius_set_act), 'w')
    txt.write('set act { ')
    for i in range(0, len(selection)):
        locA = str(selection[i]).find('<Atom ') + 6
        locB = str(selection[i]).find(': ')
        txt.write(str(selection[i])[locA:locB] + " ")
    txt.write("} ")
    txt.close()

    print('%s atoms have been set as active for the QM/MM calculation using ChemShell.' % (len(selection)))



def qmmm_modeller(u, argsdict):

    while True:
        u, argsdict, frame, ligand, radius = ask(u, argsdict)

        saver(u, argsdict, frame, ligand, radius)

        name_prmtop, name_inpcrd = topology_modeller(u, argsdict, frame, ligand, radius)

        if name_inpcrd != None and name_prmtop != None:
            while True:
                quest = input('Do you want to save the tcl list of active atoms (useful for ChemShell QM/MM calculations) (y/n)? ')
                if quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                    break
                elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', 'yeS', '1'):
                    set_act(name_prmtop, name_inpcrd, radius)
                    break
                else :
                    print('Answer just \'yes\' or \'no\'.')
                    continue

        if argsdict['trajectory'] != None:
            quest = input('Do you want to save another frame? (y/[n]) ')
            if quest in ('', 'n', 'no', 'N', 'No', 'No', 'nO', '0'):
                break

            elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                continue

        elif argsdict['frame'] != None:
            break

    return u, argsdict
