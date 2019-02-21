# coding: utf-8

''' 
Description
Using [MDAnalysis](https://www.mdanalysis.org), this program lets the user analize a given dynamics in AMBER. It is able to obtaing plots of distances, plots of RMSd, the number of frames with a certain distance under a given cut-off...
'''

# # Imported packages

#######
# Updates!
#  - [X] Ask for 'width' (in cm) if LaTeX mode before every plotting routine (after asking for how many carbons) and convert the value into inches.
#  - [X] Add LaTeX mode to group1, group2.
#  - [X] Add LaTeX mode to rmsd.
#  - [ ] Add two distances criteria for selecting frames.
#  - [ ] Fix 'serif' family font or hide sterror in LaTeX mode.
#  - [ ] Buscar si 'setted' és correcte o bé és 'set'.
######


import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
#from MDAnalysis import Universe
import MDAnalysis.analysis.distances as distances
import MDAnalysis.lib.distances as distanceslib
#import MDAnalysis.analysis.rms as mdarms
import os
import sys
import time as timer
import pandas as pd

from old_modules import *

u = 0
dir_plots = 0
time_counter = 0


# # Menu

# ### For dynamics files

def menu_dyn():
    print("\n******************************************************************************************")
    print("\nHere you have a list with the options you can choose:")
    print("\t1. Summary of the system and the dynamics.")
    print("\t2. Obtain distance plots.")
    print("\t3. Obtain the plot of RMSD of the the backbone and/or the substrate (using CPPTRAJ).")
    print("\t4. Select frames by H(subs)-protein distance.")
    print("\t5. QM/MM models from frames creation. (Frames have to be saved as pdb by this program (option 4))")
    print("\t6. Create the 'set act' file, where active atoms for ChemShell are specified")
    print("\n******************************************************************************************")


# ### For non dynamics files

def menu_nodyn():
    print("\n******************************************************************************************")
    print("\nYou are in the 'nodynamics' mode. You can just do the following tasks.")
    print("\nIf you want to do other jobs that require the topology and the dynamics, type 'exit' and rerun the script specifing those files (following this sintaxis:\n script.py topology_file_name dynamics_file_name).")
    print("\nHere you have a list with the options you can choose:")
    print("\t1. QM/MM models from frames creation. (Frames have to be saved as pdb by this program (option 4))")
    print("\t2. Create the 'set act' file, where active atoms for ChemShell are specified")
    print("\n******************************************************************************************")


while True:
    try :
        if sys.argv[1] in os.listdir() and sys.argv[2] in os.listdir():
            topology = sys.argv[1]
            dinamica = sys.argv[2]
            print("\nThe topology file is " + str(topology) + ", and the dynamics file is " + str(dinamica) + ".")
            mode = 'dyn'
            break
        elif sys.argv[1] not in os.listdir() or sys.argv[2] not in os.listdir():
            print("Please, check the names of the files and run again the script.")
            sys.exit(1)
            break
    except IndexError: #len(argv[1]) == 0 or len(argv[2]) == 0:
        print("Please, if you want to analyze the trajectory, run the script again introducing the topology and the dynamics files' names following this sintaxis:\n script.py topology_file_name dynamics_file_name")
       #sys.exit(1)
        mode = 'nodyn'
        break
while True:
    try :
        if sys.argv[3] == 'latex':
            plt.rcParams['text.usetex'] = True
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = 'Computer Modern Roman'
            plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'
            print('LaTeX mode is on!')
            latex = True
        else :
            latex = False
        break
    except IndexError:
        latex = False
        break     

        

if mode == 'dyn':
    while True:
        menu_dyn()
        menu_quest = input("\nWhat do you want to do (type the number or \"exit\")? ")
        if menu_quest == '1':
            print("\n")
            old_modules.summarize()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '2':
            print("\n")
            while True:
                c_num = input("How many carbons do you want to analyse in the same plot: 1, 2 or 3? (type \"menu\" if you want to go back) ")
                if c_num == '1':
                    old_modules.group1()
                    break
                elif c_num == '2':
                    old_modules.group2()
                    break
                elif c_num == '3':
                    gold_modules.roup3()
                    break    
                elif c_num in ('menu', 'Menu', 'MENU', 'mENU'):
                    break
                else :
                    print("Please, type 1, 2 or 3.")
                    continue
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '3':
            print("\n")
            old_modules.rmsd_func()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '4':
            print("\n")
            old_modules.frame_sel()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '5':
            print("\n")
            old_modules.model_qmmm()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '6':
            print("\n")
            old_modules.set_act()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest in ('exit', 'Exit', 'EXIT', 'eXIT'):
            print("Bye, bye! See you soon!")
            sys.exit(1)
            break
        elif menu_quest == 'time_off':
            print("\nTime counter turned off.")
            time_counter = 0
            continue
        elif menu_quest == 'time_on':
            print("\nTime counter turned on.")
            time_counter = 1
            continue
        else :
            print('\nSorry, I didn\'t understand you, could you repeat, please? ')
            continue
            
if mode == 'nodyn':
    while True:
        menu_nodyn()
        menu_quest = input("\nWhat do you want to do (type the number or \"exit\")? ")
        if menu_quest == '1':
            print("\n")
            old_modules.model_qmmm()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '2':
            print("\n")
            old_modules.set_act()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest in ('exit', 'Exit', 'EXIT', 'eXIT'):
            print("Bye, bye! See you soon!")
            sys.exit(1)
            break
        elif menu_quest == 'time_off':
            print("\nTime counter turned off.")
            time_counter = 0
            continue
        elif menu_quest == 'time_on':
            print("\nTime counter turned on.")
            time_counter = 1
            continue