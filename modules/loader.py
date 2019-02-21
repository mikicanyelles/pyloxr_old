'''
Module for defining all the functions for loading different files or variables (like folders):
    - Trajectory (based on MDAnalysis Universe class)
    - Frame (based on MDAnalysis Universe class)
'''

from MDAnalysis import Universe
import os


def universe_loader_traj(argsdict):
    if argsdict['u_loaded'] == True:
        print("Topology and trajectory(ies) were previously loaded, this will be faster!")
        return True
    else :
        print("Let's load the trajectory!")
        u = Universe(argsdict['trajectory'][0], argsdict['trajectory'][1:])
        argsdict['u_loaded'] = True
        print("Topology and trajectory(ies) loaded!")
        return u, argsdict

def universe_loader_frame(u_loaded, frame):
    global u; global traj#; global top; global dyn
    if u_loaded == True:
        print("Frame was previously loaded, this will be faster!")
        return u_loaded, frame
    else :
        print("Let's load the trajectory!")
        u = Universe(frame)
        u_loaded = True
        print("Frame loaded!")
        return u_loaded, frame

def folders_subfolders(subdir):
    dir_now = os.getcwd()
    while subdir != None:
        subdir_quest = input("Do you want to save the plots in '%s' ([y]/n)? " % subdir)
        if subdir_quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        elif subdir_quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            while True:
                subdir_quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
                if subdir_quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                    subdir = 'plots'
                    break
                elif subdir_quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                    subdir = input("In which folder do you want to save the plots? ")
                    break
                else :
                    print('Sorry, I didn\'t understand you. Please, answer yes or no.')
                    continue
            break
        else:
            print('Sorry, I didn\'t understand you. Please, answer yes or no.')
            continue


    while subdir == None:
        subdir_quest = input("Do you want to save the plots in a subfolder ([y]/n)? " )
        if subdir_quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            while True:
                subdir_quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
                if subdir_quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                    subdir = 'plots'
                    break
                elif subdir_quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                    subdir = input("In which folder do you want to save the plots? ")
                    break
                else :
                    print('Sorry, I didn\'t understand you. Please, answer yes or no.')
                    continue
            break
        elif subdir_quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue
