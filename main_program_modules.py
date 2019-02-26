from modules import loader
from modules import summariser
from modules import distance_plotter_bis
from modules import distance_plotter
from modules import frames_selector
from modules import menus
from modules import rmsd
import modules
import argparse
#import multiprocessing
from sys import exit
import sys
import os
import numpy as np

parser = argparse.ArgumentParser(description="pyLOXr - Python program for LOX reactivity")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    '-t', '--trajectory',
    help='Specify the topology and the coordinates of the trajectory in AMBER format (specify first the topology and then the trajectory/trajectories)',
    nargs='+',
    )
group.add_argument(
    '-f', '--frame',
    help='Specify a structure or a single frame from a trajectory in pdb format. Specifying only one structure limitates the available options.'
    )
parser.add_argument(
    '-c', '--csv',
    help='Save csv files of compluted distances. By default they are not saved.',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-l', '--latex',
    help='Enable LaTeX processing of titles, axis and legends of generated plots',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-lw', '--latex-width',
    help='When LaTeX mode is enabled, specify the width (in cm) of the generated plots. By default, it is set to 10 cm.',
    required=False,
    default=10
    )
parser.add_argument(
    '-m', '--module',
    help='Start the chosen module without showing the menu.',
    required=False,
    type=int
    )
parser.add_argument(
    '-p', '--parallel',
    help='Enable parallel calculation of distances. Requires OpenMP.',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-s', '--subdir',
    help='Specify a subdirectory (its name) for saving the generated plots. \'plots\' folder is set by default. If it is desired to save the plots in the current directory, type \'.\'',
    default='plots'
    )
parser.add_argument(
    '-tm', '--timer',
    help='Enable the timer for the different modules of the program (for testing purposes).',
    required=False,
    action='store_true'
    )


argsdict=vars(parser.parse_args())
if argsdict['trajectory'] != None and argsdict['frame'] == None:
    if len(argsdict['trajectory']) == 1:
        print('Topology or coordinates of the trajectory are not specified. Rerun the script specifying the files using the \'-t\' flag.')
        exit(0)
    if argsdict['trajectory'][0] not in os.listdir()[:]:
        print('The topology file is not in the directory.')
    for i in range(1,len(argsdict['trajectory'])):
        if argsdict['trajectory'][i] not in os.listdir()[:]:
            print('The coordinates file \'%s\' is not in the directory.' % argsdict['trajectory'][i])
    menu_type = 1
elif argsdict['trajectory'] == None and argsdict['frame'] != None:
    menu_type = 2


argsdict['menu_type'] = menu_type
del parser, group, menu_type

if argsdict['trajectory'] != None:
    print('Files %s and %s are being loaded for the analysis.' % (argsdict['trajectory'][0], str(argsdict['trajectory'][1:])))
    argsdict['trajectory'][0] = os.path.abspath(argsdict['trajectory'][0])

    for i in range(1, len(argsdict['trajectory'])):
        argsdict['trajectory'][i] = os.path.abspath(argsdict['trajectory'][i])

if argsdict['frame'] != None:
    print('%s is being loaded as a single structure for the analysis.' % argsdict['frame'])
    argsdict['frame'] = os.path.abspath(argsdict['frame'])

if argsdict['latex'] == True:
    print('LaTeX mode is on. The plots are going to be %s cm of width.' % argsdict['latex_width'])

if argsdict['csv'] == True:
    print('Computed values (RMSD, distances...) will be also saved as csv-type files.')

if argsdict['parallel'] == True:
    print('Parallel calculation for distances is enabled.')

if argsdict['timer'] == True:
    print('Timer is enabled.')

if argsdict['subdir'] != '.':
    if argsdict['subdir'] not in os.listdir():
        os.mkdir(argsdict['subdir'])
    print('Plots will be saved in %s subdirectory.' % argsdict['subdir'])
else :
    print('Plots will be saved in the current directory.')

argsdict['u_loaded'] = False
u = None

if argsdict['module'] == 1:
    print("Summariser module has been selected from the command line")
    u, argsdict = summariser.summarize(u, argsdict)
    argsdict['module'] = None

elif argsdict['module'] == 2:
    print("Distance plots module has been selected from the command line")
    u, argsdict = distance_plotter_bis.general_plotter(u, argsdict)
    argsdict['module'] = None

elif argsdict['module'] == 3:
    u, argsdict = rmsd.rmsd(u, argsdict)
    argsdict['module'] = None

elif argsdict['module'] == 4:
    u, argsdict = frames_selector.frame_selector(u, argsdict)
    argsdict['module'] = None

else :
    pass

if argsdict['module'] == None:

    while True:
        argsdict = menus.menu(argsdict)
        if argsdict['option'] not in ('exit', 'Exit', 'EXIT', 'eXIT', 0, '0') and argsdict['menu_type'] == 1:
            if argsdict['option'] == 1:
                u, argsdict = summariser.summarize(u, argsdict)
                continue
            elif argsdict['option'] == 2:
                u, argsdict = distance_plotter_bis.general_plotter(u, argsdict)
                argsdict['module'] = None
                continue
            elif argsdict['option'] == 3:
                u, argsdict = rmsd.rmsd(u, argsdict)
                argsdict['module'] = None
                continue
            elif argsdict['option'] == 4:
                u, argsdict = frames_selector.frame_selector(u, argsdict)
                argsdict['module'] = None
                continue
        elif argsdict['option'] in ('exit', 'Exit', 'EXIT', 'eXIT', '0', 0):
            exit(0)
            break
