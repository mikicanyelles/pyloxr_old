from MDAnalysis import Universe
import MDAnalysis.lib.distances as distanceslib
from numpy import min#, mean, array, max
import sys
import os
from progressbar import ProgressBar, Percentage, Bar, ETA, Timer
import copy
import shutil


def distance_cutoff(argsdict):
    '''
    Function that requires the numbers of the atoms to be measured, the cut-off
    distance and creates a dictionary of type 1 and which includes the list of
    the ref atom, the com atoms and the cut-off.
    The dict will have the following shape:
        dict = {type : 1, features : [[ref_atom], [comp_atoms], [cut-off]]}
    '''

    critdict = {'type' : 1, 'features' : [[], [], []]}

    while True:
        try :
            critdict['features'][0] = [int(input("Input the number of one of the selected atoms for measuring distances. "))]
            break
        except ValueError:
            print('The atom\'s number has not been correctly introduced.\n')
            continue

    while True:
        try :
            subs_atoms_input = input("Input the number of the atom(s) for measuring distances (separate by a comma if you want to consider only the nearest one). ")
            if subs_atoms_input.find(',') != -1:
                critdict['features'][1] = subs_atoms_input.split(',')
                for i in range(0, len(critdict['features'][1])):
                    critdict['features'][1][i] = int(critdict['features'][1][i])
            else :
                critdict['features'][1] = [int(subs_atoms_input)]
            break
        except ValueError:
            print('(Some of) the number(s) of the atom(s) has not been correctly introduced.\n')
            continue

    u_top = Universe(argsdict['parameters'])
    quest = None

    while True:
        print("\nYou have selected those atoms:\n")
        ### Print selected atoms (number, name, type, resname an resid) and save names
        exists = []
        for i in range(0, len(critdict['features'][0:2])):
            name = None
            for j in range(0, len(critdict['features'][0:2][i])):
                a = str(u_top.select_atoms("bynum %s" % critdict['features'][0:2][i][j]))
                exists.append(a.find('AtomGroup []') != -1)
                locA = a.find('[<') +2
                locB = a.find(' and')
                if name == None:
                    name = '\t' + a[locA:locB]
                else :
                    name = name +'; ' + a[locA:locB]
            print(name)
        print("\n")

        if True in exists:
            print('Some of the specified atom doesn\'t exist. Enter the atom numbers again.\n')
            quest = 'no'
            del exists
        else :
            quest = str(input("Are all numbers correct ([y]/n)?"))

        if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
            while True:
                try :
                    critdict['features'][0] = [int(input("Input the number of one of the selected atoms for measuring distances. "))]
                    break
                except ValueError:
                    print('The atom\'s number has not been correctly introduced.\n')
                    continue

            while True:
                try :
                    subs_atoms_input = input("Input the number of the atom(s) for measuring distances (separate by a comma if you want to consider only the nearest one). ")
                    if subs_atoms_input.find(',') != -1:
                        critdict['features'][1] = subs_atoms_input.split(',')
                        for i in range(0, len(critdict['features'][1])):
                            critdict['features'][1][i] = int(critdict['features'][1][i])
                    else :
                        critdict['features'][1] = [int(subs_atoms_input)]
                    break
                except ValueError:
                    print('(Some of) the number(s) of the atom(s) has not been correctly introduced.\n')
                    continue
            continue

        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
        break


    while True:
        try :
            critdict['features'][2] = float(input("Type the cut-off distance (in Å): "))
            break
        except ValueError:
            print("The cut-off distance has not been well specified. Please, retype the cut-off distance.")
            continue

    return critdict


def distance_comparison(argsdict):
    '''
    Function that requires the numbers of the two pairs of atoms, the difference
    (if desired) and creates a dictionary of type 2 which includes the following:
        dict = {type : cut-off, features : [[1st pair], [2nd pair], [difference]]}
    '''

    critdict = {'type' : 2, 'features' : [[], [], []]}

    while True:
        first_inp = input('Type the numbers of the atoms which form to the first bond (which will be the shortest), separated by a comma or a space: ')

        if first_inp.find(',') == -1 and first_inp.find(' ') == -1:
            print('Atoms have not been correctly introduced. Please, introduce them again.')
            continue
        elif len(first_inp.split(',')) > 2 or len(first_inp.split(' ')) > 2:
            print('More than two atoms have been specified. Specify only two.')
            continue
        else :
            if len(first_inp.split(',')) == 2:
                first_inp = first_inp.split(',')
            elif len(first_inp.split(' ')) == 2:
                first_inp = first_inp.split(' ')

            try :
                critdict['features'][0].append(int(first_inp[0]))
                critdict['features'][0].append(int(first_inp[1]))
            except ValueError:
                print('One of the atoms has not been correctly introduced. Type only the corresponding numbers')
                continue

            u_top = Universe(argsdict['parameters'])
            exists = []
            for i in range(0, len(critdict['features'][0])):
                a = str(u_top.select_atoms("bynum %s" % critdict['features'][0][i]))
                exists.append(a.find('AtomGroup []') != -1)
                locA = a.find('[<') +2
                locB = a.find(' and')
                print(a[locA:locB])

            if True in exists:
                print('Some of the specified atom doesn\'t exist. Enter the atom numbers again.\n')
                quest = 'no'
                del exists
            else :
                quest_ = True
                while quest_ == True:
                    quest = str(input("Are all numbers correct ([y]/n)? "))

                    if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0', '', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                        break
                    else :
                        print("Sorry, answer again, please.")
                        quest_ = False
                        continue

                if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                    continue
                elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                    break


    while True:
        second_inp = input('Type the numbers of the atoms which form to the first bond (which will be the longest), separated by a comma or a space: ')

        if second_inp.find(',') == -1 and second_inp.find(' ') == -1:
            print('Atoms have not been correctly introduced. Please, introduce them again.')
            continue
        elif len(second_inp.split(',')) > 2 or len(second_inp.split(' ')) > 2:
            print('More than two atoms have been specified. Specify only two.')
            continue
        else :
            if len(second_inp.split(',')) == 2:
                second_inp = second_inp.split(',')
            elif len(second_inp.split(' ')) == 2:
                second_inp = second_inp.split(' ')

            try :
                critdict['features'][1].append(int(second_inp[0]))
                critdict['features'][1].append(int(second_inp[1]))
            except ValueError:
                print('One of the atoms has not been correctly introduced. Type only the corresponding numbers')
                continue

            u_top = Universe(argsdict['parameters'])
            exists = []
            for i in range(0, len(critdict['features'][1])):
                a = str(u_top.select_atoms("bynum %s" % critdict['features'][1][i]))
                exists.append(a.find('AtomGroup []') != -1)
                locA = a.find('[<') +2
                locB = a.find(' and')
                print(a[locA:locB])

            if True in exists:
                print('Some of the specified atom doesn\'t exist. Enter the atom numbers again.\n')
                quest = 'no'
                del exists
            else :
                quest_ = True
                while quest_ == True:
                    quest = str(input("Are all numbers correct ([y]/n)? "))

                    if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0', '', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                        break
                    else :
                        print("Sorry, answer again, please.")
                        quest_ = False
                        continue

                if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                    continue
                elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                    break

    while True:
        diff_q = input('Do you want to set a minimum distance difference (y/n)? ')

        if diff_q in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
            critdict['features'][2] = None
            break
        elif diff_q in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            diff_q_ = True

            while diff_q_ == True:
                try :
                    critdict['features'][2] = float(input('Which diference do you want to set (in Å)?'))
                    diff_q_ = False
                except ValueError:
                    print('Type only the number.')
            break
        else :
            print("Sorry, answer again, please.")
            continue

    return critdict


def bool_creator(u, argsdict, critdict):#_list):
    '''
    Function for creating a boolean array from a list of critdicts. It takes the type of criteria and 
    returns the final boolean value if all criteria are satisfied or not.
    '''

    if argsdict['u_loaded'] == False:
        from modules import loader
        u, argsdict = loader.universe_loader_traj(argsdict)

    frames_bool = []
    time_ = 0
    print('trajectory length: ' + str(len(u.trajectory)))

    widgets = ['Progress: ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', ETA(), ' | ', Timer(), ' ']

    pbar = ProgressBar(widgets=widgets, maxval=len(u.trajectory))
    pbar.start()

    for ts in u.trajectory:
        bool_ar = []
        for i in range(len(critdict)):

            if critdict[i]['type'] == 1:

                ref_atom = u.select_atoms('bynum ' + str(critdict[i]['features'][0][0]))

                comp_atom_str = 'bynum'
                for j in range(len(critdict[i]['features'][1])):
                    comp_atom_str = comp_atom_str + ' ' + str(critdict[i]['features'][1][j])
                comp_atom = u.select_atoms(comp_atom_str)

                if argsdict['parallel'] == False:
                    dist_ = min(distanceslib.distance_array(ref_atom.positions, comp_atom.positions, backend='serial'))
                elif argsdict['parallel'] == True:
                    dist_ = min(distanceslib.distance_array(ref_atom.positions, comp_atom.positions, backend='OpenMP'))

                if dist_ <= critdict[i]['features'][2]:
                    bool_ar.append(True)
                elif dist_ > critdict[i]['features'][2]:
                    bool_ar.append(False)

            elif critdict[i]['type'] == 2:

                A1 = u.select_atoms('bynum ' + str(critdict[i]['features'][0][0])).positions
                B1 = u.select_atoms('bynum ' + str(critdict[i]['features'][0][1])).positions
                A2 = u.select_atoms('bynum ' + str(critdict[i]['features'][1][0])).positions
                B2 = u.select_atoms('bynum ' + str(critdict[i]['features'][1][1])).positions

                if argsdict['parallel'] == False:
                    dist1 = distanceslib.calc_bonds(A1, B1, backend='serial')
                    dist2 = distanceslib.calc_bonds(A2, B2, backend='serial')
                elif argsdict['parallel'] == True:
                    dist1 = distanceslib.calc_bonds(A1, B1, backend='OpenMP')
                    dist2 = distanceslib.calc_bonds(A2, B2, backend='OpenMP')

                if critdict[i]['features'][2] == None:
                    if dist1 < dist2:
                        bool_ar.append(True)
                    elif dist1 >= dist2:
                        bool_ar.append(False)

                elif critdict[i]['features'][2] != None:
                    if (dist1 < dist2) and (dist1-dist2 < critdict[i]['features'][2]):
                        bool_ar.append(True)
                    elif (dist1 < dist2) and (dist1-dist2 > critdict[i]['features'][2]):
                        bool_ar.append(False)
                    elif dist1 >= dist2:
                        bool_ar.append(False)


        if False in bool_ar:
            frames_bool.append(False)
        elif False not in bool_ar:
            frames_bool.append(True)


        time_+=1
        pbar.update(time_)

    pbar.finish()

    print('frames_bool length: ' + str(len(frames_bool)))

    percent = round(frames_bool.count(True)/len(frames_bool), 2) * 100

    print("The %s %% of the frames (%s frames) of the trajectory satisfy the imposed criterion(a)." % (percent, frames_bool.count(True)))

    return u, argsdict, frames_bool

def txt_saver(frames_bool, critdict):

    txt_name = 'summary_'
    for i in range(len(critdict)):
        if critdict[i]['type'] == 1:
            txt_name = txt_name + str(critdict[i]['features'][0][0]) + '_'
            for j in range(len(critdict[i]['features'][1])):
                txt_name = txt_name + str(critdict[i]['features'][1][j])
                if j != len(critdict[i]['features'][1]) -1:
                    txt_name = txt_name + ','
                #elif j == len(critdict[i]['features'][1]) -1:
                #    txt_name = txt_name + '_'
        if critdict[i]['type'] == 2:
            txt_name = txt_name + str(critdict[i]['features'][0][0]) + ',' + str(critdict[i]['features'][0][1]) + '_' + str(critdict[i]['features'][1][0]) + ',' + str(critdict[i]['features'][1][1])
            
        txt_name = txt_name + '_' + str(critdict[i]['features'][2])

    if txt_name in os.listdir():
        while True:
            quest = input("A previous summary of the frame selection exists. Do you want to overwrite it (1) or to give a new name (2) ([1]/2)? ")
            if quest == '1' or quest == '':
                break
            elif quest == '2':
                while True:
                    txt_name = input('Type the new name for the summary file: ')
                    quest2 = input('Is \'%s\' correct? ([y]/n) ' % txt_name)
                    if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                        break
                    elif quest2 in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                        continue
                break
            else :
                print('Type just \'1\' or \'2\'.')
                continue

    print('Summary file is being saved as \'%s\'.' % txt_name)

    txt_ext_name = txt_name + '_ext'

    if txt_ext_name in os.listdir():
        while True:
            quest = input("A previous extended summary of the frame selection exists. Do you want to overwrite it (1) or to give a new name (2)? ([1]/2) ")
            if quest == '1' or quest == '':
                break
            elif quest == '2':
                while True:
                    txt_ext_name = input('Type the new name for the extended summary file: ')
                    quest2 = input('Is \'%s\' correct? ([y]/n) ' % txt_ext_name)
                    if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                        break
                    elif quest2 in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                        continue
                break
            else :
                print('Type just \'1\' or \'2\'.')
                continue

    print('Extended summary file is being saved as \'%s\'.' % txt_ext_name)

    txt = open(txt_name, 'w')
    txt_ext = open(txt_ext_name, 'w')

    for i in range(len(frames_bool)):
        if frames_bool[i] == True:
            txt.write('Frame %s satisfies the imposed criteria\n' % (i+1))
            txt_ext.write('Frame %s satisfies the imposed criteria\n' % (i+1))
        elif frames_bool[i] == False:
            txt_ext.write('Frame %s does not satisfy the imposed criteria\n' % (i+1))

    T = (frames_bool.count(True)/len(frames_bool)) * 100
    F = (frames_bool.count(False)/len(frames_bool)) * 100

    txt.write('\n\nThe %s %% satisfy the imposed criteria.' % T)
    txt.write('\nThe %s %% do not satisfy the imposed criteria.' % F)

    txt_ext.write('\n\nThe %s %% satisfy the imposed criteria.' % T)
    txt_ext.write('\nThe %s %% do not satisfy the imposed criteria.' % F)

    txt.close()
    txt_ext.close()


def pdb_saver_all(u, frames_bool, critdict):
    '''
    Function for saving the pdbs of all the frames which satisfy the imposed criteria.
    '''

    dir_name = 'frames_'
    for i in range(len(critdict)):
        if critdict[i]['type'] == 1:
            dir_name = dir_name + str(critdict[i]['features'][0][0]) + '_'
            for j in range(len(critdict[i]['features'][1])):
                dir_name = dir_name + str(critdict[i]['features'][1][j])
                if j != len(critdict[i]['features'][1]) -1:
                    dir_name = dir_name + ','
        if critdict[i]['type'] == 2:
            dir_name = dir_name + str(critdict[i]['features'][0][0]) + ',' + str(critdict[i]['features'][0][1]) + '_' + str(critdict[i]['features'][1][0]) + ',' + str(critdict[i]['features'][1][1])
            
        dir_name = dir_name + '_' + str(critdict[i]['features'][2])

    if dir_name not in os.listdir():
        os.mkdir(dir_name)
        #print(dir_name)
    elif dir_name in os.listdir():
        while True:
            quest = input('%s subdirectory aready exists. Do you want to overwrite its content (1) or give an alternative name to the subdirectory (2)? ([1]/2) ' % dir_name)
            if quest in ('1', '2', ''):
                if quest in ('', '1'):
                    shutil.rmtree(dir_name)
                    os.mkdir(dir_name)
                    break
                elif quest == '2':
                    while True:
                        dir_name = input('Type the new name of the subdirectory: ')
                        quest2 = input('Is \'%s\' correct? ([y]/n) ' % dir_name)
                        if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                            break
                        elif quest2 in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                            continue

                    if dir_name not in os.listdir():
                        break
                    elif dir_name in os.listdir():
                        continue

            elif quest not in ('1', '2', ''):
                print('Choose 1 or 2.')
                continue

    time_ = 0
    widgets2 = ['Progress: ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', ETA(), ' | ', Timer(), ' ']

    pbar2 = ProgressBar(widgets=widgets2, maxval=len(frames_bool))
    pbar2.start()

    stderr_ = open(os.devnull, 'w')
    sys.stderr = stderr_

    for i in range(len(frames_bool)):
        time_+=1
        if frames_bool[i] == True:
            u.trajectory[i]
            sel = u.select_atoms('all')
            sel.write('%s/frame_%s.pdb' % (dir_name, (i+1)))

        pbar2.update(time_)

    pbar2.finish()
    stderr_.close()


def pdb_saver_some(u, frames_bool, critdict):
    '''
    Function for saving the pdbs of the selected frames which satisfy the imposed criteria.
    '''

    dir_name = 'frames_'
    for i in range(len(critdict)):
        if critdict[i]['type'] == 1:
            dir_name = dir_name + str(critdict[i]['features'][0][0]) + '_'
            for j in range(len(critdict[i]['features'][1])):
                dir_name = dir_name + str(critdict[i]['features'][1][j])
                if j != len(critdict[i]['features'][1]) -1:
                    dir_name = dir_name + ','
        if critdict[i]['type'] == 2:
            dir_name = dir_name + str(critdict[i]['features'][0][0]) + ',' + str(critdict[i]['features'][0][1]) + '_' + str(critdict[i]['features'][1][0]) + ',' + str(critdict[i]['features'][1][1])
            
        dir_name = dir_name + '_' + str(critdict[i]['features'][2])

    if dir_name not in os.listdir():
        os.mkdir(dir_name)
 
    elif dir_name in os.listdir():
        while True:
            quest = input('%s subdirectory aready exists. Do you want to overwrite its content (1) or give an alternative name to the subdirectory (2) ([1]/2)? ' % dir_name)
            if quest in ('1', '2', ''):
                if quest in ('', '1'):
                    shutil.rmtree(dir_name)
                    os.mkdir(dir_name)
                    break
                elif quest == '2':
                    while True:
                        dir_name = input('Type the new name of the subdirectory: ')
                        quest2 = input('Is \'%s\' correct? ([y]/n) ' % dir_name)
                        if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                            break
                        elif quest2 in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                            continue

                    if dir_name not in os.listdir():
                        break
                    elif dir_name in os.listdir():
                        continue

            elif quest not in ('1', '2', ''):
                print('Choose 1 or 2.')
                continue

    stderr_ = open(os.devnull, 'w')
    sys.stderr = stderr_

    while True:
        try :
            frame_n = str(input('Which frame do you want to save as pdb? (type \'list\' to show a list of frame to select) '))

            if frame_n in ('list', 'LIST', 'List', 'lIST', 'LIst', 'liST', 'LISt', 'lisT'):
                print_ = ''
                for i in range(len(frames_bool)):
                    if frames_bool[i] == True:
                        if print_ == '':
                            print_ = str(int(i+1))
                        else :
                            print_ = print_ + ', ' + str(int(i+1))

                print(print_); del print_
                continue

            elif frame_n == '':
                print('Type the number.')
                continue

            else :
                frame_n = int(frame_n) -1
                if frames_bool[frame_n] == True:
                    u.trajectory[frame_n]
                    sel = u.select_atoms('all')
                    sel.write('%s/frame_%s.pdb' % (dir_name, (frame_n+1)))
                    print('Frame %s saved on %s.' % ((frame_n+1), dir_name))
                    break

                elif frames_bool[frame_n] == False:
                    print('This frame does not satisfy the imposed criteria. Please, select another frame. (Type \'list\' to print the list of frames that satisfy the imposed criteria)')
                    continue

        except TypeError:
            print('Type only the number.')

    while True:
        try :
            frame_n = str(input('Do you want to save another frame as pdb (frame number/[n])? '))

            if frame_n in ('list', 'LIST', 'List', 'lIST', 'LIst', 'liST', 'LISt', 'lisT'):
                print_ = ''
                for i in range(len(frames_bool)):
                    if frames_bool[i] == True:
                        if print_ == '':
                            print_ = str(int(i+1))
                        else :
                            print_ = print_ + ', ' + str(int(i+1))

                print(print_); del print_
                continue

            elif frame_n in ('', 'n', 'no', 'N', 'No', 'NO', 'nO', '0'):
                break

            else :
                frame_n = int(frame_n)-1
                if frames_bool[frame_n] == True:
                    u.trajectory[frame_n]
                    sel = u.select_atoms('all')
                    sel.write('%s/frame_%s.pdb' % (dir_name, (frame_n+1)))
                    print('Frame %s saved on %s.' % ((frame_n+1), dir_name))
                    continue

                elif frames_bool[frame_n] == False:
                    print('This frame does not satisfy the imposed criteria. Please, select another frame. (Type \'list\' to print the list of frames that satisfy the imposed criteria)')
                    continue

        except TypeError:
            print('Type only the number.')

    stderr_.close()


def frame_selector(u, argsdict):

    critdict = []

    while True:
        try :
            quest = int(input("Which type of criteria do you want to impose, cut-off distance (1) or difference of distances between two bonds (2) (1/2)? "))
            if quest == 1:
                critdict.append(distance_cutoff(argsdict))
                break
            elif quest == 2:
                critdict.append(distance_comparison(argsdict))
                break
            else :
                print('Type 1 or 2.')
                continue
        except ValueError:
            print('Type 1 or 2.')
            continue


    while True:
        quest = input("Do you want to add another distance criteria for selecting frames (y/[n])? ")
        if quest in ('', 'n', 'no', 'N', 'No', 'No', 'nO', '0'):
            break
        elif quest in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            a = True
            while a == True:
                try :
                    quest = int(input("Which type of criteria do you want to impose, cut-off distance (1) or difference of distances between two bonds (2) (1/2)? "))
                    if quest == 1:
                        critdict.append(distance_cutoff(argsdict))
                        a = False
                    elif quest == 2:
                        critdict.append(distance_comparison(argsdict))
                        a = False
                    else :
                        print('Type only 1 or 2.')
                except ValueError:
                    print('Type 1 or 2.')
            continue
        else :
            print("Sorry, answer again, please.")
            continue

    # bool_creator has to be converted to the looper
    u, argsdict, frames_bool = bool_creator(u, argsdict, critdict)

    while True:
        quest = input('Do you want to save a summary of the selection results ([y]/n)? ')
        if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            txt_saver(frames_bool, critdict)
            break
        elif quest in ('n', 'no', 'N', 'No', 'NO', 'nO', '0'):
            break
        else :
            print('Please, answer \'yes\' or \'no\'.')
            continue

    while True:
        try :
            quest = int(input('Do you want to save all the frames (1), the selected ones (2) or any of them (3)? '))

            if quest == 1:
                pdb_saver_all(u, frames_bool, critdict)
                break

            elif quest == 2:
                pdb_saver_some(u, frames_bool, critdict)
                break

            elif quest == 3:
                break

            else :
                print('Type 1, 2 or 3')
                continue

        except ValueError:
            print('Type 1, 2 or 3')
            continue

    #del atoms_list; del cut_off; del atoms_list_; del cut_off_
    return u, argsdict
