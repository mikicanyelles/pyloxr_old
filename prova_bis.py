import argparse
import multiprocessing
from sys import exit
import os

parser = argparse.ArgumentParser(description="pyLOXr - Python program for LOX reactivity")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument(
    '-t', '--trajectory',
    help='Specify the topology and the coordinates of the trajectory in AMBER format (specify first the topology and then the trajectory/trajectories)',
    nargs='+'
    )
group.add_argument(
    '-f', '--frame',
    help='Specify a structure or a single frame from a trajectory in pdb formats. Specifying only one structure limitates the available options.'
    )
parser.add_argument(
    '-l', '--latex',
    help='Enable LaTeX processing of titles, axis and legends of generated plots',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-tm', '--timer',
    help='Enable the timer for the different modules of the program (for testing purposes)',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-p', '--parallel',
    help='Enable parallel calculation of distances. Requires OpenMP.',
    required=False,
    action='store_true'
    )
parser.add_argument(
    '-pn', '--processors',
    help='Number of processors for parallel computing of distances. Not mandatory when parallelisation is enbabled, by default the number of processors is set to the maximum available.',
    required=False,
    default=multiprocessing.cpu_count()
    )

args = parser.parse_args()
argsdict=vars(args)
#print(argsdict['trajectory'])
#print(argsdict['parallel'])
#print(argsdict['latex'])
#print(argsdict)


if argsdict['trajectory'] != None and argsdict['frame'] == None:
    if len(argsdict['trajectory']) == 1:
        print('Topology or coordinates of the trajectory are not specified. Rerun the script specifying the files using the \'-t\' flap.')
        exit(0)
    if argsdict['trajectory'][0] not in os.listdir()[:]:
        print('The topology file is not in the directory.')
    for i in range(0,len(argsdict['trajectory'])-1):
        if argsdict['trajectory'][i+1] not in os.listdir()[:]:
            print('The coordinates file \'%s\' is not in the directory.' % argsdict['trajectory'][i+1])
    #traj_files = argsdict['trajectory']
    menu_type = 1
    #print(traj_files)
elif argsdict['trajectory'] == None and argsdict['frame'] != None:
    #frame = argsdict['frame']
    menu_type = 2
    print(frame)

if int(argsdict['processors']) > multiprocessing.cpu_count():
    #ppn = multiprocessing.cpu_count()
    #print('The number of used processos for paralellisation for the distances computation is ' + str(ppn))
    print('The maximum number of processors available on this computer is %s but %s are specified. Please, rerun the script specifying %s or less processors.' % (multiprocessing.cpu_count(), argsdict['processors'], multiprocessing.cpu_count()))
    exit(0)

#if argsdict['latex'] == True:
    #plt.rcParams['text.usetex'] = True
    #plt.rcParams['font.family'] = 'serif'
    #plt.rcParams['font.serif'] = 'Computer Modern Roman'
    #plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'
    #latex = True
    #print('LaTeX mode is on!')

argsdict['menu_type'] = menu_type
del parser, group, args, menu_type
#print(argsdict)
