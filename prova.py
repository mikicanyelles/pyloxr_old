import argparse
import multiprocessing
import sys

parser = argparse.ArgumentParser(description="pyLOXr - Python program for LOX reactivity")

parser.add_argument(
    '-t', '--topology',
    help='Specify the topology of the system in AMBER format',
    required=True
    )
parser.add_argument(
    '-c', '--coordinates',
    help='Specify the MD trajectory in AMBER format',
    required=True
    )
parser.add_argument(
    '-f', '--frame',
    help='Specify a structure or a single frame from a trajectory in pdb or mol2 formats. Specifying only one structure limitates the available options.'
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

#**AFEGIR QUE O NECESSITA -t I -c O NECESSITA -f**

args = parser.parse_args()
argsdict=vars(args)
print(argsdict['coordinates'])
print(argsdict['topology'])
print(argsdict['parallel'])
print(argsdict['latex'])
print(argsdict)

if int(argsdict['processors']) > multiprocessing.cpu_count():
    #ppn = multiprocessing.cpu_count()
    #print('The number of used processos for paralellisation for the distances computation is ' + str(ppn))
    print('The maximum number of processors available on this computer is %s but %s are specified. Please, rerun the script specifying %s or less processors.' % (multiprocessing.cpu_count(), argsdict['processors'], multiprocessing.cpu_count()))
    sys.exit(0)

top = argsdict['topology']
dyn = argsdict['coordinates']
if argsdict['latex'] == True:
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'Computer Modern Roman'
    plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'
    print('LaTeX mode is on!')
    latex = True
elif argsdict['latex'] == False:
    latex = False
