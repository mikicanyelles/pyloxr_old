
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
import MDAnalysis.analysis.distances as distances
import MDAnalysis.lib.distances as distanceslib
#import MDAnalysis.analysis.rms as mdarms
import os
import sys
import time as timer
import pandas as pd
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


# # Summary of the system and the dynamics

def summarize():
    global u
    global traj
    ### Time counter starts
    if time_counter == 1:
        time_in = timer.time()

    ### Check if universe is loaded and load it if it's not
    if u != 0: 
        print("Files were previously loaded, this will be faster!")
    else :
        print("Let's load the files!")
        u = mda.Universe(topology, dinamica)
        traj = u.trajectory
        print("Topology and dynamics loaded!")
        
    txt = open("MD_summary.txt", 'w+')
    txt.write("The system has %s atoms, %s residues and %s segments.\n" % (len(u.atoms), len(u.residues), len(u.segments)))
    txt.write("The dynamics starts at %s ns, ends at %s ns, lasts a total of %s ns, has %s frames and a timestep of %s.\n" % (round(traj[0].time/1000,1), traj[-1].time/1000, traj[-1].time/1000-round(traj[0].time/1000,1), len(traj), ((traj[-1].time/1000-round(traj[0].time/1000,1))/(traj[-1].frame-traj[0].frame +1))))
    txt.write("The system has the following shape:\n")
    txt.write("\tEdges length: (%s, %s, %s).\n" % (u.dimensions[0],u.dimensions[1],u.dimensions[2]))
    txt.write("\tAngles: (%s, %s, %s).\n" % (u.dimensions[3],u.dimensions[4],u.dimensions[5]))
    txt.close()
    
    print("Summary saved!")


# # Dynamics analysis

# ### For 1 carbon

def group1():
    global u
    global traj
    global dir_plots
    global time_counter
    global dinamica
    global topology
    dir_now = os.getcwd()
    while dir_plots != 0:
        subplots = input("Do you want to save the plots in '%s' ([y]/n)? " % dir_plots)
        if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            dir_plots = 'plots'
            break
        elif subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue
    while dir_plots == 0:
        subplots = input("Do you want to save the plots in a subfolder ([y]/n)? ")
        if subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        elif subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
            if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                if 'plots' not in list(os.listdir()):
                    dir_plots = 'plots'
                    os.mkdir(dir_plots)
                elif 'plots' in list(os.listdir()):
                    dir_plots = 'plots'
            elif quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                dir_plots = input("In which folder do you want to save the plots? ")
                if dir_plots not in list(os.listdir()):
                    os.mkdir(dir_plots)
            os.chdir(dir_plots)
            topology = '../%s' %(topology)
            dinamica = '../%s' %(dinamica)
            print("Plots will be saved in '%s'" % dir_plots)
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue   
                
    if latex == True:
        width_plots = float(input('Which is the desired width for LaTeX plots?'))*0.39370079

    ### Ask for atom numbers
    while True:
        try :
            o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    
    
    ### Ask if atom numbers are correct
    while True:
        u_top = mda.Universe(topology)
        ### Print selected atoms (number, name, tupe, resname an resid) and save names 
        print("\nYou have selected those atoms:\n")
    
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(u_top.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(list(u_top.select_atoms("bynum %s" % o_prot)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_o_prot = str(a[locA:locB])
    
        a = str(list(u_top.select_atoms("bynum %s" % c9)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c9 = str(a[locA:locB])
        t_h9 = t_c9.replace('C','H')

        quest = str(input("Are all numbers correct ([y]/n)?"))
        if quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            while True:
                try :
                    o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            continue
        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
    
    ### Time counter
    if time_counter == 1:
        time_in = timer.time()
    
    ### Check if universe is loaded and load it if it's not

    if u != 0: 
        print("Files were previously loaded, this will be faster!")
    else :
        print("Let's load the files!")
        u = mda.Universe(topology, dinamica)
        traj = u.trajectory
        print("Topology and dynamics loaded!")

    
    ### Creation of the lists of distances and time
    
    time = []
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o    = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
        time.append(int(traj.time)/1000)
    
    dist_h_9  = []
    for i in range(0, len(traj)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])
    
    ### Avg, min and max distances
    
    avg_9_c = np.mean(dist_c_9)
    min_9_c = np.min(dist_c_9)
    max_9_c = np.max(dist_c_9)
    rng_9_c = np.max(dist_c_9) - np.min(dist_c_9)
    avg_9_h = np.mean(dist_h_9)
    min_9_h = np.min(dist_h_9)
    max_9_h = np.max(dist_h_9)
    rng_9_h = np.max(dist_h_9) - np.min(dist_h_9)
    avg_9_c_ar = np.array([avg_9_c for i in range(0, len(traj))])
    avg_9_h_ar = np.array([avg_9_h for i in range(0, len(traj))])
    print("Lists of distances created")
    
    
    ### Save summary of distances (avg, max and min)
    
    txt = open('summary_of_distances_%s.txt' % t_c9, 'w+')
    
    txt.write('Distances for %s: \n'% t_c9)
    txt.write('\tAverage distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(avg_9_c,3)) + ' Å \n' )
    txt.write('\tMinimum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(min_9_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(max_9_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_c9,t_o_prot) + str(round(rng_9_c,3)) + ' Å \n\n')
    
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n' % t_h9)
    txt.write('\tAverage distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(avg_9_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(min_9_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(max_9_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_h9,t_o_prot) + str(round(rng_9_h,3)) + ' Å \n\n\n')
    
    txt.close()
    
    print("Summary saved")
    
    ### csv files
    csv_array = pd.DataFrame()


    ### Carbon distances plots
    
    ##### Histogram
    plt.hist(dist_c_9, bins=20,range=(0,10), histtype='bar')
    #plt.hist(dist10,bins=40,range=(0,10))
    #plt.hist(dist13,bins=40,range=(0,10))
    if latex == True:
        plt.xlabel('Distance ($\AA$)')
    elif latex == False:
        plt.xlabel('Distance (Å)')    
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s" % (t_c9), y=1.08, loc='center')
    plt.savefig('hist_%s.png' % (t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s.eps' %(t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print('Histogram of carbons saved')
    
    ##### plot w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9)#, label='C_9-OH')
    #ax1.plot(range(0,len(traj)), dist_c_12, color='green')#,  label='C_12-OH')
    #ax1.plot(range(0,len(traj)), dist_c_15, color='coral',  label='C_15-OH')
    #ax1.plot(range(0,len(traj)), avg_9_c_ar, color='purple')
    #ax1.plot(range(0,len(traj)), avg_12_c_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_c_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_9)#,  label='C_15-OH')
    #ax2.plot(time_ar, avg_9_c_ar, color='purple')
    #ax2.plot(time_ar, avg_12_c_ar, color='lime')
    #ax2.plot(time_ar, avg_15_c_ar, color='orangered')
    #ax2.set_xticks(time)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s" % (t_c9), y=1.15, loc='center')
    plt.savefig('plot_%s.png' % (t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s.eps' % (t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9)#, label='C_9-OH')
    #ax1.plot(range(0,len(traj)), dist_c_12, color='green',  label='C_12-OH')
    #ax1.plot(range(0,len(traj)), dist_c_15, color='coral',  label='C_15-OH')
    ax1.plot(range(0,len(traj)), avg_9_c_ar, color='blue')
    #ax1.plot(range(0,len(traj)), avg_12_c_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_c_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    #ax2.plot(time, dist_c_, color='coral',  label='C_15-OH')
    ax2.plot(time, avg_9_c_ar, color='blue')
    #ax2.plot(time, avg_12_c_ar, color='lime')
    #ax2.plot(time, avg_15_c_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s vs. time and frames with average distances" % (t_c9), y=1.15, loc='center')
    plt.savefig('plot_%s_avg.png' % (t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_avg.eps' % (t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print("Plots of carbons saved")
    
    
    ### Hydrogen plots
    
    ##### Histogram
    plt.hist([dist_h_9], bins=20,range=(0,10), histtype='bar')
    #plt.hist(dist10,bins=40,range=(0,10))
    #plt.hist(dist13,bins=40,range=(0,10))
    if latex == True:
        plt.xlabel('Distance ($\AA$)')
    elif latex == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s" % (t_h9), y=1.08, loc='center')
    plt.savefig('hist_%s.png' % (t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s.eps' %(t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print('Histogram of hydrogen saved')
    
    ##### Plot w/o average
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9)#, label='H_9-OH')
    #ax1.plot(range(0,len(traj)), dist_h_12, color='green',  label='H_12-OH')
    #ax1.plot(range(0,len(traj)), dist_h_15, color='coral',  label='H_15-OH')
    #ax1.plot(range(0,len(traj)), avg_9_h_ar, color='purple')
    #ax1.plot(range(0,len(traj)), avg_12_h_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_h_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_9, label='H_15-OH')
    #ax2.plot(time_ar, avg_9_h_ar, color='purple')
    #ax2.plot(time_ar, avg_12_h_ar, color='lime')
    #ax2.plot(time_ar, avg_15_h_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of H-OH distances for %s vs. time and frames" % (t_h9), y=1.15, loc='center')
    plt.savefig('plot_%s.png' % (t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s.eps' % (t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9)#, label='H_9-OH')
    #ax1.plot(range(0,len(traj)), dist_h_12, color='green',  label='H_12-OH')
    #ax1.plot(range(0,len(traj)), dist_h_15, color='coral',  label='H_15-OH')
    ax1.plot(range(0,len(traj)), avg_9_h_ar)
    #ax1.plot(range(0,len(traj)), avg_12_h_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_h_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    #ax2.plot(time, dist_h_15, color='coral',  label='C_15-OH')
    ax2.plot(time, avg_9_h_ar, color='blue')
    #ax2.plot(time, avg_12_h_ar, color='lime')
    #ax2.plot(time, avg_15_h_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of H-OH distances for %s vs. time and frames" % (t_h9), y=1.15, loc='center')
    plt.savefig('plot_%s_avg.png' % (t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_avg.eps' %(t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print("Plots of carbons saved")
    
    
    ### Scattering plot
    
    plt.scatter(dist_h_9,dist_c_9)
    #plt.scatter(dist_h_12,dist_c_12, color='green')
    #plt.scatter(dist_h_15,dist_c_15, color='coral')
    plt.legend(['%s-%s' % (t_c9,t_h9)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.grid(True)
    if latex == True:
        plt.xlabel('O-C distance ($\AA$)')
        plt.ylabel('O-H distance ($\AA$)')
    elif latex == False:
        plt.xlabel('O-C distance (Å)')
        plt.ylabel('O-H distance (Å)')
    plt.savefig('scatter_%s.png' %(t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('scatter_%s.eps' % (t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print("Scatter of carbon vs hydrogen distances saved")
    
    if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'): 
        os.chdir(dir_now)
        
    ### Time counter ends
    if time_counter == 1:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

    topology = topology[3:]
#####################################################################


# ### For 2 carbons

# In[5]:


def group2():
    global u
    global traj
    global dir_plots
    global time_counter
    global dinamica
    global topology
    dir_now = os.getcwd()
    while dir_plots != 0:
        subplots = input("Do you want to save the plots in '%s' ([y]/n)? " % dir_plots)
        if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            dir_plots = 'plots'
            break
        elif subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue
    while dir_plots == 0:
        subplots = input("Do you want to save the plots in a subfolder ([y]/n)? ")
        if subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        elif subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
            if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                if 'plots' not in list(os.listdir()):
                    dir_plots = 'plots'
                    os.mkdir(dir_plots)
                elif 'plots' in list(os.listdir()):
                    dir_plots = 'plots'
            elif quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                dir_plots = input("In which folder do you want to save the plots? ")
                if dir_plots not in list(os.listdir()):
                    os.mkdir(dir_plots)
            os.chdir(dir_plots)
            topology = '../%s' %(topology)
            dinamica = '../%s' %(dinamica)
            print("Plots will be saved in '%s'" % dir_plots)
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue

    if latex == True:
        width_plots = float(input('Which is the desired width for LaTeX plots?'))*0.39370079
    
    ### Ask for atom numbers
    while True:
        try :
            o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c12 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8815
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8816
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    
    
    ### Ask if atom numbers are correct
    while True:
        u_top = mda.Universe(topology)
        ### Print selected atoms (number, name, tupe, resname an resid) and save names 
        print("\nYou have selected those atoms:\n")
    
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(u_top.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(list(u_top.select_atoms("bynum %s" % o_prot)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_o_prot = str(a[locA:locB])
    
        a = str(list(u_top.select_atoms("bynum %s" % c9)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c9 = str(a[locA:locB])
        t_h9 = t_c9.replace('C','H')
        
        a = str(u_top.select_atoms("bynum %s" % c12))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
        
        a = str(list(u_top.select_atoms("bynum %s" % c12)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c12 = str(a[locA:locB])
        t_h12 = str(t_c12.replace('C','H'))
        
        quest = str(input("Are all numbers correct ([y]/n)?"))
        if quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            while True:
                try :
                    o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c12 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8815
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8816
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
    
            
    ### Time counter starts
    if time_counter == 1:
        time_in = timer.time()
    
    
    ### Check if universe is loaded and load it if it's not
    global u
    global traj
    if u != 0: 
        print("Files were previously loaded, this will be faster!")
    else :
        print("Let's load the files!")
        u = mda.Universe(topology, dinamica)
        traj = u.trajectory
        print("Topology and dynamics loaded!")
    
    ### Creation of the lists of distances
    
    ###### For carbon 9
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
    dist_h_9  = []
    for i in range(0, len(traj)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])
    avg_9_c = np.mean(dist_c_9)
    min_9_c = np.min(dist_c_9)
    max_9_c = np.max(dist_c_9)
    rng_9_c = np.max(dist_c_9) - np.min(dist_c_9)
    avg_9_h = np.mean(dist_h_9)
    min_9_h = np.min(dist_h_9)
    max_9_h = np.max(dist_h_9)
    rng_9_h = np.max(dist_h_9) - np.min(dist_h_9)
    avg_9_c_ar = np.array([avg_9_c for i in range(0, len(traj))])#needed for plots, it is just an array of the same value (the average) repeated as many times as frames the dynamics has
    avg_9_h_ar = np.array([avg_9_h for i in range(0, len(traj))])
    
    ##### For carbon 12
    time = []
    dist_c_12  = []
    dist_ha_12 = []
    dist_hb_12 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c12).positions
        pos_ha = u.select_atoms("bynum %s" % h12a).positions
        pos_hb = u.select_atoms("bynum %s" % h12b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
        time.append(traj.time/1000)
    dist_h_12  = []
    for i in range(0, len(traj)):
        if dist_ha_12[i] < dist_hb_12[i]:
            dist_h_12.append(dist_ha_12[i])
        elif dist_ha_12[i] > dist_hb_12[i]:
            dist_h_12.append(dist_hb_12[i])
    avg_12_c = np.mean(dist_c_12)
    min_12_c = np.min(dist_c_12)
    max_12_c = np.max(dist_c_12)
    rng_12_c = np.max(dist_c_12) - np.min(dist_c_12)
    avg_12_h = np.mean(dist_h_12)
    min_12_h = np.min(dist_h_12)
    max_12_h = np.max(dist_h_12)
    rng_12_h = np.max(dist_h_12) - np.min(dist_h_12)
    avg_12_c_ar = np.array([avg_12_c for i in range(0, len(traj))])
    avg_12_h_ar = np.array([avg_12_h for i in range(0, len(traj))])
    
    print("Lists of distances created")
    
    
    ### Summary of distances
    
    txt = open('summary_of_distances_%s_%s.txt' % (t_c9, t_c12), 'w+')
    
    txt.write('Distances for %s: \n'% t_c9)
    txt.write('\tAverage distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(avg_9_c,3)) + ' Å \n' )
    txt.write('\tMinimum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(min_9_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(max_9_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_c9,t_o_prot) + str(round(rng_9_c,3)) + ' Å \n\n')
    
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n' % t_h9)
    txt.write('\tAverage distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(avg_9_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(min_9_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(max_9_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_h9,t_o_prot) + str(round(rng_9_h,3)) + ' Å \n\n\n')
    
    
    txt.write('Distances for %s: \n' % t_c12)
    txt.write('\tAverage distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(avg_12_c,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(min_12_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(max_12_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_c12,t_o_prot) + str(round(rng_12_c,3)) + ' Å \n\n')
    
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n'% t_h12 )
    txt.write('\tAverage distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(avg_12_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(min_12_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(max_12_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_h12,t_o_prot) + str(round(rng_12_h,3)) + ' Å \n\n\n')
    
    txt.close()
    
    print("Summary saved")
    
    
    ### Carbon plots
    
    ##### Histogram
    plt.hist([dist_c_9, dist_c_12], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green'])
    if latex == True:
        plt.xlabel('Distance ($\AA$)')
    elif latex == False:
        plt.xlabel('Distance (Å)')    
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s and %s" % (t_c9,t_c12), y=1.08, loc='center')
    plt.savefig('hist_%s_%s.png' % (t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s_%s.eps' %(t_c9,t_c12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print('Histogram of carbons saved')
    
    ##### Distances w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9, color='indigo')#, label='C_9-OH')
    ax1.plot(range(0,len(traj)), dist_c_12, color='green')#,  label='C_12-OH')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_12, color='green')#,  label='C_15-OH')
    #ax2.plot(time_ar, avg_9_c_ar, color='purple')
    #ax2.plot(time_ar, avg_12_c_ar, color='lime')
    #ax2.plot(time_ar, avg_15_c_ar, color='orangered')
    #ax2.set_xticks(time)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s and %s vs. time and frames" % (t_c9,t_c12), y=1.15, loc='center')
    plt.savefig('plot_%s_%s.png' % (t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s.eps' %(t_c9,t_c12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    ### Distances w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9, color='indigo')#, label='C_9-OH')
    ax1.plot(range(0,len(traj)), dist_c_12, color='green',  label='C_12-OH')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_c_ar, color='purple')
    ax2.plot(time, avg_12_c_ar, color='lime')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s and %s vs. time and frames with average distances" % (t_c9,t_c12), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_avg.png' % (t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_avg.eps' %(t_c9,t_c12,t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print("Plots of carbons saved")
    
    
    ### Hydrogen plots
    ##### Histogram
    plt.hist([dist_h_9, dist_h_12], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green'])
    if latex == True:
        ax1.xlabel('Distance ($\AA$)')
    elif latex == False:
        ax1.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s and %s" % (t_h9,t_h12), y=1.08, loc='center')
    plt.savefig('hist_%s_%s.png' % (t_h9,t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s_%s.eps' %(t_h9,t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print('Histogram of hydrogen saved')
    
    ##### Distance w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9, color='indigo')#, label='H_9-OH')
    ax1.plot(range(0,len(traj)), dist_h_12, color='green')#,  label='H_12-OH')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_12, color='green')#,  label='H_12-OH')
    #ax2.plot(time_ar, avg_9_h_ar, color='purple')
    #ax2.plot(time_ar, avg_12_h_ar, color='lime')
    #ax2.plot(time_ar, avg_15_h_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of H-OH distances for %s and %s vs. time and frames" % (t_h9,t_h12), y=1.15, loc='center')
    plt.savefig('plot_%s_%s.png' % (t_h9,t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s.eps' %(t_h9,t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    ##### Distance w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9, color='indigo')#, label='H_9-OH')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(range(0,len(traj)), dist_h_12, color='green',  label='H_12-OH')
    ax2.plot(time, avg_9_h_ar, color='purple')
    ax2.plot(time, avg_12_h_ar, color='lime')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ps)')
    
    plt.title("Plot of H-OH distances for %s and %s vs. time and frames with averages" % (t_h9,t_h12), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_avg.png' % (t_h9,t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_avg.eps' %(t_h9,t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print("Plots of carbons saved")
    
    ##### Scatter
    
    plt.scatter(dist_h_9,dist_c_9, color='indigo')
    plt.scatter(dist_h_12,dist_c_12, color='green')
    plt.legend(['%s-%s' % (t_c9,t_h9), '%s-%s' % (t_c12,t_h12)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.grid(True)
    if latex == True:
        plt.xlabel('O-C distance ($\AA$)')
        plt.ylabel('O-H distance ($\AA$)')
    elif latex == False:
        plt.xlabel('O-C distance (Å)')
        plt.ylabel('O-H distance (Å)')
    plt.title("Scatter of distances for %s vs. %s and %s vs %s" % (t_c9, t_h9, t_c12, t_h12))
    plt.savefig('scatter_%s_%s.png' %(t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('scatter_%s_%s.eps' %(t_h9,t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print("Scatter of carbon vs hydrogen distances saved")
    
    if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'): 
        os.chdir(dir_now)
    ### Time counter ends
    if time_counter == 1:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

    topology = topology[3:]

#####################################################################

# ### For 3 carbons

# In[6]:


def group3():
    global u
    global traj
    global dir_plots
    global dinamica
    global topology
    dir_now = os.getcwd()
    while dir_plots != 0:
        subplots = input("Do you want to save the plots in '%s' ([y]/n)? " % dir_plots)
        if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            dir_plots = 'plots'
            break
        elif subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue
    while dir_plots == 0:
        subplots = input("Do you want to save the plots in a subfolder ([y]/n)? ")
        if subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        elif subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
            if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                if 'plots' not in list(os.listdir()):
                    dir_plots = 'plots'
                    os.mkdir(dir_plots)
                elif 'plots' in list(os.listdir()):
                    dir_plots = 'plots'
            elif quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                dir_plots = input("In which folder do you want to save the plots? ")
                if dir_plots not in list(os.listdir()):
                    os.mkdir(dir_plots)
            topology = '../%s' %(topology)
            dinamica = '../%s' %(dinamica)
            os.chdir(dir_plots)
            print("Plots will be saved in '%s'" % dir_plots)
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue 

    if latex == True:
        width_plots = float(input('Which is the desired width for LaTeX plots?'))*0.39370079

    ### Ask for atom numbers
    while True:
        try :
            o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c12 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8815
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8816
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    
    while True:
        try :
            c15 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h15a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8822
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h15b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8823
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    
    ### Ask if the numbers are correct
    while True:
        u_top = mda.Universe(topology)
        ### Print selected atoms (number, name, tupe, resname an resid) and save names 
        print("\nYou have selected those atoms:\n")
    
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(u_top.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(list(u_top.select_atoms("bynum %s" % o_prot)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_o_prot = str(a[locA:locB])
    
        a = str(list(u_top.select_atoms("bynum %s" % c9)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c9 = str(a[locA:locB])
        t_h9 = t_c9.replace('C','H')
        
        a = str(u_top.select_atoms("bynum %s" % c12))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
        
        a = str(list(u_top.select_atoms("bynum %s" % c12)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c12 = str(a[locA:locB])
        t_h12 = str(t_c12.replace('C','H'))
        
        a = str(u_top.select_atoms("bynum %s" % c15))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h15a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h15b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
        
        a = str(list(u_top.select_atoms("bynum %s" % c15)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c15 = str(a[locA:locB])
        t_h15 = str(t_c15.replace('C','H'))
        
        quest = str(input("Are all numbers correct ([y]/n)?"))
        if quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            while True:
                try :
                    o_prot  = int(input("Number of the atom which belongs to the protein: ")) #-1 #8784                
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8809
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c12 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8815
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8816
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
    
            while True:
                try :
                    c15 = int(input("\nNumber of one of the carbons: ")) #-1 #8814
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h15a = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8822
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h15b = int(input("Number of the other hydrogen bonded to the previous carbon: ")) #-1 #8823
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            continue
        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
            
    ### Time counter starts
    if time_counter == 1:
        time_in = timer.time()
    
    ### Load universe if it's not loaded
    global u
    global traj
    if u != 0: 
        print("Files were previously loaded, this will be faster!")
    else :
        print("Let's load the files!")
        u = mda.Universe(topology, dinamica)
        traj = u.trajectory
        print("Topology and dynamics loaded!")
        
    ### Creation of the lists of distances
    ##### C9
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
    dist_h_9  = []
    for i in range(0, len(traj)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])
    avg_9_c = np.mean(dist_c_9)
    min_9_c = np.min(dist_c_9)
    max_9_c = np.max(dist_c_9)
    rng_9_c = np.max(dist_c_9) - np.min(dist_c_9)
    avg_9_h = np.mean(dist_h_9)
    min_9_h = np.min(dist_h_9)
    max_9_h = np.max(dist_h_9)
    rng_9_h = np.max(dist_h_9) - np.min(dist_h_9)
    avg_9_c_ar = np.array([avg_9_c for i in range(0, len(traj))])
    avg_9_h_ar = np.array([avg_9_h for i in range(0, len(traj))])
    
    ##### C12
    time = []
    dist_c_12  = []
    dist_ha_12 = []
    dist_hb_12 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c12).positions
        pos_ha = u.select_atoms("bynum %s" % h12a).positions
        pos_hb = u.select_atoms("bynum %s" % h12b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
        time.append(traj.time/1000)
    dist_h_12  = []
    for i in range(0, len(traj)):
        if dist_ha_12[i] < dist_hb_12[i]:
            dist_h_12.append(dist_ha_12[i])
        elif dist_ha_12[i] > dist_hb_12[i]:
            dist_h_12.append(dist_hb_12[i])
    avg_12_c = np.mean(dist_c_12)
    min_12_c = np.min(dist_c_12)
    max_12_c = np.max(dist_c_12)
    rng_12_c = np.max(dist_c_12) - np.min(dist_c_12)
    avg_12_h = np.mean(dist_h_12)
    min_12_h = np.min(dist_h_12)
    max_12_h = np.max(dist_h_12)
    rng_12_h = np.max(dist_h_12) - np.min(dist_h_12)
    avg_12_c_ar = np.array([avg_12_c for i in range(0, len(traj))])
    avg_12_h_ar = np.array([avg_12_h for i in range(0, len(traj))])
    
    ##### C15
    dist_c_15  = []
    dist_ha_15 = []
    dist_hb_15 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c15).positions
        pos_ha = u.select_atoms("bynum %s" % h15a).positions
        pos_hb = u.select_atoms("bynum %s" % h15b).positions
        pos_o    = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_15.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_15.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_15.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
    dist_h_15  = []
    for i in range(0, len(traj)):
        if dist_ha_15[i] < dist_hb_15[i]:
            dist_h_15.append(dist_ha_15[i])
        elif dist_ha_15[i] > dist_hb_15[i]:
            dist_h_15.append(dist_hb_15[i])
    avg_15_c = np.mean(dist_c_15)
    min_15_c = np.min(dist_c_15)
    max_15_c = np.max(dist_c_15)
    rng_15_c = np.max(dist_c_15) - np.min(dist_c_15)
    avg_15_h = np.mean(dist_h_15)
    min_15_h = np.min(dist_h_15)
    max_15_h = np.max(dist_h_15)
    rng_15_h = np.max(dist_h_15) - np.min(dist_h_15)
    avg_15_c_ar = np.array([avg_15_c for i in range(0, len(traj))])
    avg_15_h_ar = np.array([avg_15_h for i in range(0, len(traj))])
    print("Lists of distances created")
    
    
    ### Summary of distances
    
    txt = open('summary_of_distances_%s_%s_%s.txt' % (t_c9, t_c12, t_c15), 'w+')
    txt.write('Distances for %s: \n'% t_c9)
    txt.write('\tAverage distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(avg_9_c,3)) + ' Å \n' )
    txt.write('\tMinimum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(min_9_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_c9,t_o_prot) + str(round(max_9_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_c9,t_o_prot) + str(round(rng_9_c,3)) + ' Å \n\n')
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n' % t_h9)
    txt.write('\tAverage distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(avg_9_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(min_9_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:       ' % (t_h9,t_o_prot) + str(round(max_9_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s:  ' % (t_h9,t_o_prot) + str(round(rng_9_h,3)) + ' Å \n\n\n')
    txt.write('Distances for %s: \n' % t_c12)
    txt.write('\tAverage distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(avg_12_c,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(min_12_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_c12,t_o_prot) + str(round(max_12_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_c12,t_o_prot) + str(round(rng_12_c,3)) + ' Å \n\n')
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n'% t_h12 )
    txt.write('\tAverage distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(avg_12_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(min_12_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_h12,t_o_prot) + str(round(max_12_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_h12,t_o_prot) + str(round(rng_12_h,3)) + ' Å \n\n\n')
    txt.write('Distances for %s: \n' % t_c15 )
    txt.write('\tAverage distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(avg_15_c,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(min_15_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(max_15_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_c15,t_o_prot) + str(round(rng_15_c,3)) + ' Å \n\n')
    txt.write('Distances for the nearest %s, even if it is the A or the B: \n' % t_h15 )
    txt.write('\tAverage distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(avg_15_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(min_15_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(max_15_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_h15,t_o_prot) + str(round(rng_15_h,3)) + ' Å \n\n\n')
    txt.close()
    print("Summary saved")
    
    ### Carbon distances
    ##### Histogram
    
    plt.hist([dist_c_9, dist_c_12, dist_c_15], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green', 'coral'])
    #plt.hist(dist10,bins=40,range=(0,10))
    #plt.hist(dist13,bins=40,range=(0,10))
    if latex == True:
      plt.xlabel('Distance ($\AA$)')
    elif latex == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_c15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s, %s and %s" % (t_c9,t_c12,t_c15), y=1.08, loc='center')
    plt.savefig('hist_%s_%s_%s.png' % (t_c9,t_c12,t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s_%s_%s.eps' % (t_c9,t_c12,t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print('Histogram of carbon saved')
    
    ##### distances w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9, color='indigo')#, label='C_9-OH')
    ax1.plot(range(0,len(traj)), dist_c_12, color='green')#,  label='C_12-OH')
    ax1.plot(range(0,len(traj)), dist_c_15, color='coral',  label='C_15-OH')
    #ax1.plot(range(0,len(traj)), avg_9_c_ar, color='purple')
    #ax1.plot(range(0,len(traj)), avg_12_c_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_c_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_c15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_15, color='coral')#,  label='C_15-OH')
    #ax2.plot(time_ar, avg_9_c_ar, color='purple')
    #ax2.plot(time_ar, avg_12_c_ar, color='lime')
    #ax2.plot(time_ar, avg_15_c_ar, color='orangered')
    #ax2.set_xticks(time)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s, %s and %s vs. time and frames" % (t_c9,t_c12,t_c15), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_%s.png' % (t_c9,t_c12,t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_%s.eps' % (t_c9,t_c12,t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')    
    #plt.show()
    plt.close()
    
    ##### distances w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_c_9, color='indigo')#, label='C_9-OH')
    ax1.plot(range(0,len(traj)), dist_c_12, color='green',  label='C_12-OH')
    ax1.plot(range(0,len(traj)), dist_c_15, color='coral',  label='C_15-OH')
    ax1.plot(range(0,len(traj)), avg_9_c_ar, color='purple')
    ax1.plot(range(0,len(traj)), avg_12_c_ar, color='lime')
    ax1.plot(range(0,len(traj)), avg_15_c_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_c15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_15, color='coral',  label='C_15-OH')
    ax2.plot(time, avg_9_c_ar, color='purple')
    ax2.plot(time, avg_12_c_ar, color='lime')
    ax2.plot(time, avg_15_c_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of C-OH distances for %s, %s and %s vs. time and frames with average distances" % (t_c9,t_c12,t_c15), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_%s_avg.png' % (t_c9,t_c12,t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_%s_avg.eps' % (t_c9,t_c12,t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print("Plots of carbons saved")
    
    
    ### Hydrogen plots
    ##### Histogram
    plt.hist([dist_h_9, dist_h_12,dist_h_15], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green','coral'])
    #plt.hist(dist10,bins=40,range=(0,10))
    #plt.hist(dist13,bins=40,range=(0,10))
    if latex == True:
              plt.xlabel('Distance ($\AA$)')
    elif latex == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s, %s and %s" % (t_h9,t_h12,t_h15), y=1.08, loc='center')
    plt.savefig('hist_%s_%s_%s.png' % (t_h9,t_h12,t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('hist_%s_%s_%s.eps' % (t_h9,t_h12,t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print('Histogram of hydrogen saved')
    
    ##### Distances w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9, color='indigo')#, label='H_9-OH')
    ax1.plot(range(0,len(traj)), dist_h_12, color='green',  label='H_12-OH')
    ax1.plot(range(0,len(traj)), dist_h_15, color='coral',  label='H_15-OH')
    #ax1.plot(range(0,len(traj)), avg_9_h_ar, color='purple')
    #ax1.plot(range(0,len(traj)), avg_12_h_ar, color='lime')
    #ax1.plot(range(0,len(traj)), avg_15_h_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')    
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_15, color='coral',  label='H_15-OH')
    #ax2.plot(time_ar, avg_9_h_ar, color='purple')
    #ax2.plot(time_ar, avg_12_h_ar, color='lime')
    #ax2.plot(time_ar, avg_15_h_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of H-OH distances for %s, %s and %s vs. time and frames" % (t_h9,t_h12,t_h15), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_%s.png' % (t_h9,t_h12,t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_%s.eps' % (t_h9,t_h12,t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    
    ##### Distances w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(traj)), dist_h_9, color='indigo')#, label='H_9-OH')
    ax1.plot(range(0,len(traj)), dist_h_12, color='green',  label='H_12-OH')
    ax1.plot(range(0,len(traj)), dist_h_15, color='coral',  label='H_15-OH')
    ax1.plot(range(0,len(traj)), avg_9_h_ar, color='purple')
    ax1.plot(range(0,len(traj)), avg_12_h_ar, color='lime')
    ax1.plot(range(0,len(traj)), avg_15_h_ar, color='orangered')
    ax1.set_xlabel('Frame')
    if latex == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif latex == False:
        ax1.set_ylabel('Distance (Å)')    
    ax1.set_xticks(range(0,len(traj) +1, round(len(traj)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    
    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_15, color='coral',  label='C_15-OH')
    ax2.plot(time, avg_9_h_ar, color='purple')
    ax2.plot(time, avg_12_h_ar, color='lime')
    ax2.plot(time, avg_15_h_ar, color='orangered')
    #ax2.set_xticks(time_ar)
    ax2.set_xlabel('Time (ns)')
    
    plt.title("Plot of H-OH distances for %s, %s and %s vs. time and frames with averages" % (t_h9,t_h12,t_h15), y=1.15, loc='center')
    plt.savefig('plot_%s_%s_%s_avg.png' % (t_h9,t_h12,t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('plot_%s_%s_%s_avg.eps' % (t_h9,t_h12,t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print("Plots of carbons saved")
    
    
    ### Scatter
    plt.scatter(dist_h_9,dist_c_9, color='indigo')
    plt.scatter(dist_h_12,dist_c_12, color='green')
    plt.scatter(dist_h_15,dist_c_15, color='coral')
    plt.legend(['%s-%s' % (t_c9,t_h9), '%s-%s' % (t_c12,t_h12), '%s-%s' % (t_c15,t_h15)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.grid(True)
    if latex == True:
        plt.xlabel('O-C distance ($\AA$)')
        plt.ylabel('O-H distance ($\AA$)')
    elif latex == False:
        plt.xlabel('O-C distance (Å)')
        plt.ylabel('O-H distance (Å)')
    plt.savefig('scatter_%s_%s_%s.png' %(t_c9,t_c12,t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if latex == True:
        plt.savefig('scatter_%s_%s_%s.eps' %(t_c9,t_c12,t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()
    print("Scatter of carbon vs hydrogen distances saved")
    
    if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'): 
        os.chdir(dir_now)
    ### Time counter ends
    if time_counter == 1:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

    topology = topology[3:]
#####################################################################


# # RMSd

# In[7]:


def rmsd_func():
    global u
    global traj
    global dir_plots
    global time_counter
    global dinamica
    global topology
    dir_now = os.getcwd()
    while dir_plots != 0:
        subplots = input("Do you want to save the plots in '%s' ([y]/n)? " % dir_plots)
        if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            dir_plots = 'plots'
            os.chdir(dir_plots)
            break
        elif subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            dir_plots = 0
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue
    while dir_plots == 0:
        subplots = input("Do you want to save the plots in a subfolder ([y]/n)? ")
        if subplots in ('n', 'no', 'N', 'No', 'No', 'nO'):
            break
        elif subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            quest2 = input("Do you want to save in the \'plots\' folder (recomended) ([y]/n)? ")
            if quest2 in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                if 'plots' not in list(os.listdir()):
                    dir_plots = 'plots'
                    os.mkdir(dir_plots)
                elif 'plots' in list(os.listdir()):
                    dir_plots = 'plots'
            elif quest2 in ('n', 'no', 'N', 'No', 'No', 'nO'):
                dir_plots = input("In which folder do you want to save the plots? ")
                if dir_plots not in list(os.listdir()):
                    os.mkdir(dir_plots)
            os.chdir(dir_plots)
            topology = '../%s' %(topology)
            dinamica = '../%s' %(dinamica)
            print("Plots will be saved in '%s'" % dir_plots)
            break
        else :
            print("Sorry, I didn't understand you. Answer again, please.")
            continue      

    if latex == True:
        width_plots = float(input('Which is the desired width for LaTeX plots?'))*0.39370079

    while True:
        try :
            type_plot = int(input("Do you want to plot the RMSD of the backbone (1), the substrate (2) or both (3)? "))
            if type_plot not in (1, 2, 3):
                print("Type just \'1\', \'2\' or \'3\'.")
                continue
            elif type_plot in (1, 2, 3):
                break
        except ValueError:
            print("Type just \'1\', \'2\' or \'3\'.")
            continue



    ### Plots of backbone
    if type_plot == 1:
        last_resid = int(input("Type the number of the last residue which belongs to the protein: "))
        while True: 
            quest= input("Is '%s' correct ([y]/n)? " % last_resid)
            if quest in ('', 'yes', 'y', 'Y', 'YES', 'Yes', 'yES', 'yEs', 'YeS'):
                break
            elif quest in ('n', 'no', 'N', 'NO', 'No'):
                last_resid = int(input("Type the number of the last residue which belongs to the protein: "))
                continue
            else :
                print("Type 'yes' or 'no'.")
                continue

        while True:
            try :
                rmsd_ref = int(input("Which structure do you want to set as the reference: the first frame of the production (1) or the average structure (2)? "))
                if rmsd_ref not in (1,2):
                    print("Type just '1' or '2'.")
                    continue
                elif rmsd_ref in (1,2):
                    break
            except ValueError:
                print("Type just '1' or '2'.")
                continue

        ### Time counter starts
        if time_counter == 1:
            time_in = timer.time()
        
        if rmsd_ref == 1:
            rmsd_ref_name = "first frame \nof the production"
            rmsd_ref_plot = 'first'
            f = open("rmsd_bb.in", 'w+')
            f.write("trajin %s\n" % dinamica)
            f.write("rmsd rmsd_bb :1-%s@CA,C,N out rmsd_bb.dat\n" % last_resid)
            f.write("run\n")
            f.write("exit\n")
            f.close()
        if rmsd_ref == 2:
            rmsd_ref_name = "average \nstructure"
            rmsd_ref_plot = 'avg'
            in_file = open('rmsd_bb.in', 'w+')
            in_file.write("trajin %s\naverage crdset avg_bb :1-%s@CA,C,N\nrun\nrmsd rmsd_bb :1-%s@CA,C,N out rmsd_bb.dat ref avg_bb\nrun\nexit" % (dinamica, last_resid,last_resid) )
            in_file.close()

        os.system("cpptraj -p %s -i rmsd_bb.in > /dev/null" % topology)
        os.system("sed -i -e 's/^[ ]*//' -e 's/        /, /g' -e 's/       /, /g' -e 's/    /, /g' rmsd_bb.dat")
        os.remove("rmsd_bb.in")
        if rmsd_ref == 1:
            os.rename("rmsd_bb.dat", "rmsd_bb.csv")
            csv_bb = pd.read_csv('rmsd_bb.csv')
        elif rmsd_ref == 2:
            os.rename("rmsd_bb.dat", "rmsd_avg_bb.csv")
            csv_bb = pd.read_csv('rmsd_avg_bb.csv')

        bb_array = np.array(csv_bb)
        
        ## Plot of rmsd
        plt.plot(bb_array[:,0], bb_array[:,1])
        plt.xlabel("Frame")
        if latex == True:
            plt.ylabel("RMSD ($\AA$)")
        elif latex == False:
            plt.ylabel("RMSD (Å)")
        plt.grid(True)
        plt.legend(["Backbone RMSD \ncompared to %s" % rmsd_ref_name], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        plt.title("RMSD of the backbone of the protein")
        plt.savefig('plot_rmsd_bb_%s.png' % rmsd_ref_plot, transparent=False, dpi=300, bbox_inches='tight')
        if latex == True:
            plt.savefig('plot_rmsd_bb_%s.eps' % rmsd_ref_plot, transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
        plt.close()



    ### Plots of substrate
    if type_plot == 2:
        subs_code = input("Type the 3-letters code of the substrate (you can find it in any pdb of the protein): ")
        while True: 
            quest= input("Is '%s' correct ([y]/n)? " % subs_code)
            if quest in ('', 'yes', 'y', 'Y', 'YES', 'Yes', 'yES', 'yEs', 'YeS'):
                break
            elif quest in ('n', 'no', 'N', 'NO', 'No'):
                subs_code = input("Type the 3-letters code of the substrate (you can find it in any pdb of the protein): ")
                continue
            else :
                print("Type 'yes' or 'no'.")
                continue

        while True:
            try :
                rmsd_ref = int(input("Which structure do you want to set as the reference: the first frame of the production (1) or the average structure (2)? "))
                if rmsd_ref not in (1,2):
                    print("Type just '1' or '2'.")
                    continue
                elif rmsd_ref in (1,2):
                    break
            except ValueError:
                print("Type just '1' or '2'.")
                continue

        ### Time counter starts
        if time_counter == 1:
            time_in = timer.time()
        
        if rmsd_ref == 1:
            rmsd_ref_name = "first frame \nof the production"
            rmsd_ref_plot = 'first'
            f = open("rmsd_subs.in", 'w+')
            f.write("trajin %s\n" % dinamica)
            f.write("rmsd rmsd_subs :%s out rmsd_subs.dat\n" % subs_code)
            f.write("run\n")
            f.write("exit\n")
            f.close()
        if rmsd_ref == 2:
            rmsd_ref_name = "average \nstructure"
            rmsd_ref_plot = 'avg'
            in_file = open('rmsd_subs.in', 'w+')
            in_file.write("trajin %s\naverage crdset avg_subs :%s\nrun\nrmsd rmsd_subs :%s out rmsd_subs.dat ref avg_subs\nrun\nexit" % (dinamica, subs_code, subs_code) )
            in_file.close()

        os.system("cpptraj -p %s -i rmsd_subs.in > /dev/null" % topology)
        os.system("sed -i -e 's/^[ ]*//' -e 's/        /, /g' -e 's/       /, /g' -e 's/    /, /g' rmsd_subs.dat")
        os.remove("rmsd_subs.in")
        if rmsd_ref == 1:
            os.rename("rmsd_subs.dat", "rmsd_subs.csv")
            csv_subs = pd.read_csv('rmsd_subs.csv')
        elif rmsd_ref == 2:
            os.rename("rmsd_subs.dat", "rmsd_avg_subs.csv")
            csv_subs = pd.read_csv('rmsd_avg_subs.csv')

        subs_array = np.array(csv_subs)
        
        ## Plot of rmsd
        plt.plot(subs_array[:,0], subs_array[:,1])
        plt.xlabel("Frame")
        if latex == True:
            plt.ylabel("RMSD ($\AA$)")
        elif latex == False:
            plt.ylabel("RMSD (Å)")
        plt.grid(True)
        plt.legend(["Substrate RMSD \ncompared to %s" % rmsd_ref_name], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        plt.title("RMSD of the substrate")
        plt.savefig('plot_rmsd_subs_%s.png' % rmsd_ref_plot, transparent=False, dpi=300, bbox_inches='tight')
        if latex == True:
            plt.savefig('plot_rmsd_subs_%s.eps' % rmsd_ref_plot, transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
        plt.close()




    ### Plots of backbone and the substrate
    if type_plot == 3:
        last_resid = int(input("Type the number of the last residue which belongs to the protein: "))
        while True: 
            quest= input("Is '%s' correct ([y]/n)? " % last_resid)
            if quest in ('', 'yes', 'y', 'Y', 'YES', 'Yes', 'yES', 'yEs', 'YeS'):
                break
            elif quest in ('n', 'no', 'N', 'NO', 'No'):
                last_resid = int(input("Type the number of the last residue which belongs to the protein: "))
                continue
            else :
                print("Type 'yes' or 'no'.")
                continue

        while True:
            try :
                rmsd_ref_bb = int(input("Which structure do you want to set as the reference for the backbone's RMSD: the first frame of the production (1) or the average structure (2)? "))
                if rmsd_ref_bb not in (1,2):
                    print("Type just '1' or '2'.")
                    continue
                elif rmsd_ref_bb in (1,2):
                    break
            except ValueError:
                print("Type just '1' or '2'.")
                continue

        subs_code = input("Type the 3-letters code of the substrate (you can find it in any pdb of the protein): ")
        while True: 
            quest= input("Is '%s' correct ([y]/n)? " % subs_code)
            if quest in ('', 'yes', 'y', 'Y', 'YES', 'Yes', 'yES', 'yEs', 'YeS'):
                break
            elif quest in ('n', 'no', 'N', 'NO', 'No'):
                subs_code = input("Type the 3-letters code of the substrate (you can find it in any pdb of the protein): ")
                continue
            else :
                print("Type 'yes' or 'no'.")
                continue

        while True:
            try :
                rmsd_ref_subs = int(input("Which structure do you want to set as the reference for the substrate's RMSD: the first frame of the production (1) or the average structure (2)? "))
                if rmsd_ref_subs not in (1,2):
                    print("Type just '1' or '2'.")
                    continue
                elif rmsd_ref_subs in (1,2):
                    break
            except ValueError:
                print("Type just '1' or '2'.")
                continue

        ### Time counter starts
        if time_counter == 1:
            time_in = timer.time()
        
        ### .dat for bb
        if rmsd_ref_bb == 1:
            rmsd_ref_name_bb = "first frame \nof the production"
            rmsd_ref_plot_bb = 'first'
            f = open("rmsd_bb.in", 'w+')
            f.write("trajin %s\n" % dinamica)
            f.write("rmsd rmsd_bb :1-%s@CA,C,N out rmsd_bb.dat\n" % last_resid)
            f.write("run\n")
            f.write("exit\n")
            f.close()
        if rmsd_ref_bb == 2:
            rmsd_ref_name_bb = "average \nstructure"
            rmsd_ref_plot_bb = 'avg'
            in_file = open('rmsd_bb.in', 'w+')
            in_file.write("trajin %s\naverage crdset avg_bb :1-%s@CA,C,N\nrun\nrmsd rmsd_bb :1-%s@CA,C,N out rmsd_bb.dat ref avg_bb\nrun\nexit" % (dinamica, last_resid,last_resid) )
            in_file.close()

        os.system("cpptraj -p %s -i rmsd_bb.in > /dev/null" % topology)
        os.system("sed -i -e 's/^[ ]*//' -e 's/        /, /g' -e 's/       /, /g' -e 's/    /, /g' rmsd_bb.dat")
        os.remove("rmsd_bb.in")
        if rmsd_ref_bb == 1:
            os.rename("rmsd_bb.dat", "rmsd_bb.csv")
            csv_bb = pd.read_csv('rmsd_bb.csv')
        elif rmsd_ref_bb == 2:
            os.rename("rmsd_bb.dat", "rmsd_avg_bb.csv")
            csv_bb = pd.read_csv('rmsd_avg_bb.csv')

        bb_array = np.array(csv_bb)
        
        ### .dat for subs
        if rmsd_ref_subs == 1:
            rmsd_ref_name_subs = "first frame \nof the production"
            rmsd_ref_plot_subs = 'first'
            f = open("rmsd_subs.in", 'w+')
            f.write("trajin %s\n" % dinamica)
            f.write("rmsd rmsd_subs :%s out rmsd_subs.dat\n" % subs_code)
            f.write("run\n")
            f.write("exit\n")
            f.close()
        if rmsd_ref_subs == 2:
            rmsd_ref_name_subs = "average \nstructure"
            rmsd_ref_plot_subs = 'avg'
            in_file = open('rmsd_subs.in', 'w+')
            in_file.write("trajin %s\naverage crdset avg_subs :%s\nrun\nrmsd rmsd_subs :%s out rmsd_subs.dat ref avg_subs\nrun\nexit" % (dinamica, subs_code, subs_code) )
            in_file.close()

        os.system("cpptraj -p %s -i rmsd_subs.in > /dev/null" % topology)
        os.system("sed -i -e 's/^[ ]*//' -e 's/        /, /g' -e 's/       /, /g' -e 's/    /, /g' rmsd_subs.dat")
        os.remove("rmsd_subs.in")
        if rmsd_ref_subs == 1:
            os.rename("rmsd_subs.dat", "rmsd_subs.csv")
            csv_subs = pd.read_csv('rmsd_subs.csv')
        elif rmsd_ref_subs == 2:
            os.rename("rmsd_subs.dat", "rmsd_avg_subs.csv")
            csv_subs = pd.read_csv('rmsd_avg_subs.csv')

        subs_array = np.array(csv_subs)
      
        ## Plot of rmsd
        plt.plot(bb_array[:,0], bb_array[:,1], color='darkblue')
        plt.plot(subs_array[:,0], subs_array[:,1], color='red')
        plt.xlabel("Frame")
        if latex == True:
            plt.ylabel("RMSD ($\AA$)")
        elif latex == False:
            plt.ylabel("RMSD (Å)")
        plt.grid(True)
        plt.legend(["Backbone RMSD \ncompared to %s" % rmsd_ref_name_bb, "Substrate RMSD \ncompared to %s" % rmsd_ref_name_subs], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        plt.title("RMSD of the backbone of the protein and of the substrate")
        plt.savefig('plot_rmsd_bb_%s_subs_%s.png' % (rmsd_ref_plot_bb,rmsd_ref_plot_subs), transparent=False, dpi=300, bbox_inches='tight')
        if latex == True:
            plt.savefig('plot_rmsd_bb_%s_subs_%s.eps' % (rmsd_ref_plot_bb,rmsd_ref_plot_subs), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
        plt.close()

    if subplots in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'): 
        os.chdir(dir_now)
        dinamica = dinamica[3:]
        topology = topology[3:]
    
    ### Time counter ends
    if time_counter == 1:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")
#rmsd_func()
#####################################################################


# # Frames selection

# In[8]:


def frame_sel():
    global u
    global traj


    ### Ask for atom numbers
    while True:
        try :
            o_prot  = int(input("Number of the atom which belongs to the protein: "))                
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
            
    while True:
        try :
            c9  = int(input("\nNumber of one of the carbons: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            cut_off  = float(input("Cut-off distance (in Å): "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    
    ### Ask if atom numbers are correct
    while True:
        u_top = mda.Universe(topology)
        ### Print selected atoms (number, name, tupe, resname an resid) and save names 
        print("\nYou have selected those atoms:\n")
    
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(u_top.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")
    
        a = str(list(u_top.select_atoms("bynum %s" % o_prot)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_o_prot = str(a[locA:locB])
    
        a = str(list(u_top.select_atoms("bynum %s" % c9)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c9 = str(a[locA:locB])
        t_h9 = t_c9.replace('C','H')
        print("\nAnd a cut-off distance of %s Å.\n" % cut_off)
        quest = str(input("Are all numbers correct ([y]/n)?"))
        if quest in ('n', 'no', 'N', 'No', 'No', 'nO'):
            while True:
                try :
                    o_prot  = int(input("Number of the atom which belongs to the protein: "))                
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
                    
            while True:
                try :
                    c9  = int(input("\nNumber of one of the carbons: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9a  = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h9b  = int(input("Number of the other hydrogen bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    cut_off  = float(input("Cut-off distance (in Å): "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            continue
        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
    
    
    while True:
        save = str(input("Do you want to save the pdbs (y/[n])?"))
        if save in ('', 'n', 'no', 'N', 'NO', 'No', 'nO'):
            save = 0
            break
        elif save in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
            save = 1
            break
        else :
            print('I didn\'t understand you, could you repeat, please?')
            continue
                  
    ### Time counter starts
    if time_counter == 1:
        time_in = timer.time()

    ### Creates directory for saving pdbs and gives some options if it already exists
    if save == 1:
        a = str(list(u_top.select_atoms("bynum %s" % c9)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c9 = str(a[locA:locB])
        while True:
            try :
                dir_name = "saved_frames_%s_cut_off_%s" % (t_c9, cut_off)
                if dir_name not in os.listdir():
                    os.mkdir(dir_name)
                break
            except FileExistsError :
                print("The folder where the pdb's will be saved already exists.")
                while True:
                    quest = input("Do you want to change the folder's name (1), create a subdirectory (2) or overwrite it (3)?")           
                    if quest == '1':
                        dir_name = str(input("Type the folder's name: "))
                        os.mkdir(dir_name)
                        break
                    elif quest == '2':
                        subdir_name = input("Type the subdirectory name: ")
                        dir_name = "%s/%s" % (dir_name, subdir_name)
                        os.mkdir(dir_name)
                        print("Subdirectory created.")
                        break
                    elif quest == '3':
                        #os.system("rm -rF saved_frames_%s_cut_off_%s/*" % (t_c9, cut_off))
                        #os.remove("%s/*")
                        #os.rmdirs(dir_name)
                        shutil.rmtree(dir_name)
                        os.mkdir(dir_name)
                        print("Directory created again.")
                        break
                    else :
                        print("I didn't understand you. Type 1 or 2, please.")
                        continue       
                break
        print("All selected frames will be saved in %s." % dir_name)
        
        
    ### Check if universe is loaded and load it if it's not
    if u != 0: 
        print("Files were previously loaded, this will be faster!")
    else :
        print("Let's load the files!")
        u = mda.Universe(topology, dinamica)
        traj = u.trajectory
        print("Topology and dynamics loaded!")
    
    ### Creation of the lists of distances and time
    
    #time = []
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o    = u.select_atoms("bynum %s" % o_prot).positions
        dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha)[0])
        dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb)[0])
        dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c)[0])
        #time.append(int(traj.time)/1000)
    
    dist_h_9  = []
    for i in range(0, len(traj)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])
    
    ### Avg, min and max distances
    
    avg_9_c = np.mean(dist_c_9)
    min_9_c = np.min(dist_c_9)
    max_9_c = np.max(dist_c_9)
    rng_9_c = np.max(dist_c_9) - np.min(dist_c_9)
    avg_9_h = np.mean(dist_h_9)
    min_9_h = np.min(dist_h_9)
    max_9_h = np.max(dist_h_9)
    rng_9_h = np.max(dist_h_9) - np.min(dist_h_9)
    avg_9_c_ar = np.array([avg_9_c for i in range(0, len(traj))])
    avg_9_h_ar = np.array([avg_9_h for i in range(0, len(traj))])
    print("Lists of distances created")
    
    ### Loop for changing to True if cut-off is satisfied
    i = 0
    dist_bool = []
    for i in range(0,len(traj)):
        if dist_h_9[i] <= cut_off and dist_h_9[i] < dist_c_9[i]:
            dist_bool.append(True)
        else :
            dist_bool.append(False)

    ### Loop for: Number of frames that satisfy the cut-off (print and save in the txt), save which will or not be saved as pdb,

    if save == 0:
        txt = open('summary_saved_frames_%s.txt' % (t_c9), 'w+')
        txt.write("\nYou have selected those atoms:\n")
        a = str(u.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n" + "\n")
        a = str(u.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        txt.write("\nNumber of saved frames: " + str(dist_bool.count(1)) + "\n")

        txt_ext = open('summary_saved_frames_ext_%s.txt' % (t_c9), 'w+')
        txt_ext.write("\nYou have selected those atoms:\n")
        a = str(u.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n" + "\n")
        a = str(u.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n") 
        txt_ext.write("\nNumber of saved frames: " + str(dist_bool.count(1)) + "\n")
    
        for i in range(0,len(traj)):
            if dist_bool[i] == 1:
                txt_ext.write("Pdb saved for frame " + str(i+1) + "\n")
            elif dist_bool[i] == 0:
                txt_ext.write("Pdb not saved for frame " + str(i+1) + "\n")
                
        txt.close()
        txt_ext.close()
        print("Summaries saved")
    
    if save == 1:
        
        txt = open('summary_saved_frames_%s.txt' % (t_c9), 'w+')
        txt.write("\nYou have selected those atoms:\n")
        a = str(u.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n" + "\n")
        a = str(u.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt.write("\t" + a[locA:locB] + "\n")
        txt.write("\nNumber of saved frames: " + str(dist_bool.count(1)) + "\n")


        while True:
            all = input('Do you want to save all the pdbs (1) or only some of them? (1/[2])? ')
            if all in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                all_quest = 1
                break
            elif all in ('', 'n', 'no', 'N', 'No', 'No', 'nO'):
                all_quest = 0
                break
            else :
                print('Sorry, I did\'t understand you. Could you answer again, please?')
                continue

        txt_ext = open('summary_saved_frames_ext_%s.txt' % (t_c9), 'w+')
        txt_ext.write("\nYou have selected those atoms:\n")
        a = str(u.select_atoms("bynum %s" % o_prot))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n" + "\n")
        a = str(u.select_atoms("bynum %s" % c9))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9a))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n")
        a = str(u.select_atoms("bynum %s" % h9b))
        locA = a.find('[<') +2
        locB = a.find(' and')
        txt_ext.write("\t" + a[locA:locB] + "\n") 
        txt_ext.write("\nNumber of saved frames: " + str(dist_bool.count(1)) + "\n")
        
        if all_quest == 1:
            print("Saving frames in pdb files. This will take a while...")
            f2 = open(os.devnull, 'w')
            sys.stderr = f2
            for i in range(0,len(traj)):
                if dist_bool[i] == 1:
                    u.trajectory[i]
                    sel = u.select_atoms('all')
                    sel.write("%s/frame_%s.pdb" % (dir_name, i+1))
                    txt_ext.write("Pdb saved for frame " + str(i+1) + "\n")
                elif dist_bool[i] == 0:
                    txt_ext.write("Pdb not saved for frame " + str(i+1) + "\n")
            f2.close()
    
            txt.close()
            txt_ext.close()
              
            print("All frames have been saved in " + str(dir_name))

        elif all_quest == 0:
            for i in range(0, len(traj)):
                if dist_bool[i] == 1:
                    txt_ext.write("Pdb saved for frame %s \n" % (i+1))
                elif dist_bool[i] == 0:
                    txt_ext.write("Pdb not saved for frame %s \n" % (i+1))

            while True:
                try :
                    selected_frame = int(input('Which frame do you want to save as a pdb? ')) - 1
                    if dist_bool[selected_frame] == 1:
                        break
                    elif dist_bool[selected_frame] == 0:
                        print('This frame doesn\'t accomplish the required criteria. Please, type another frame (you can check which can be saved in the \'summary_saved_frames_ext_%s.txt file).' % t_c9)
                        continue
                except IndexError:
                    print("This frame doesn't exist.")
                    continue
                except ValueError :
                    print('Type just the number, please.')
                    continue                              

            f2 = open(os.devnull, 'w')
            sys.stderr = f2
            u.trajectory[selected_frame]
            sel = u.select_atoms('all')
            sel.write("%s/frame_%s.pdb" % (dir_name, selected_frame+1))
            txt.write("Pdb saved for frame " + str(selected_frame+1) + "\n")
            f2.close()
            print('Frame %s saved as a pdb' % (selected_frame+1))
            while True:
                more_frames = input('Do you want to save any other frame (y/[n])? ')
                if more_frames in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                    while True:
                        try :
                            selected_frame = int(input('Which frame do you want to save as a pdb? ')) - 1
                            if dist_bool[selected_frame] == 1:
                                break
                            elif dist_bool[selected_frame] == 0:
                                print('This frame doesn\'t accomplish the required criteria. Please, type another frame (you can check which can be saved in the \'summary_saved_frames_ext_%s.txt file).' % t_c9)
                                continue
                        except IndexError:
                            print("This frame doesn't exist.")
                            continue
                        except ValueError :
                            print('Type just the number, please.')
                            continue
                    f2 = open(os.devnull, 'w')
                    sys.stderr = f2
                    u.trajectory[selected_frame]
                    sel = u.select_atoms('all')
                    sel.write("%s/frame_%s.pdb" % (dir_name, selected_frame+1))
                    txt.write("Pdb saved for frame " + str(selected_frame+1) + "\n")
                    f2.close()
                    continue
                elif more_frames in ('','n', 'no', 'N', 'No', 'No', 'nO'):
                    txt.close()
                    break
                else :
                    print('Sorry, I didn\'t understand you. Answer again, please')
                    continue
    ### Time counter ends
    if time_counter == 1:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")


# # QM/MM model as pdb creation

# In[19]:


def model_qmmm():
    exists = []
    del exists[:]
    dir_name_frame = []
    del dir_name_frame[:]
    #lobal dir_name_fin
    for i in range(0, len(os.listdir())):
        if (os.listdir()[i].startswith("saved_frames_") == True) and (os.listdir()[i].find("_cut_off_") != '-1') and int(list(os.listdir()[i])[-1]) % 1 == 0:
            dir_name_frame.append(os.listdir()[i])
            exists.append(1)
    while len(exists) == 0:
        while True:
            dir_name_frame = input('In which directory are the pdbs of the frames placed?')
            if dir_name_frame in os.listdir():
                exists.append(1)
                break
            else :
                print("This directory doesn't exist.")
                continue
            #if dir_name in ('menu', 'Menu', 'MENU'):
    if len(exists) == 1:
        dir_name_fin = dir_name_frame[0]
    if len(exists) > 1:
        print("There are more than one possible folders where pdbs of frames can be placed:")#, which is the desired one?")
        for i in range(0, len(dir_name_frame)):
            print("%s. %s" % (i+1, dir_name_frame[i]))
        while True:
            try :
                quest = int(input("Which is the folder that you want to set as the working directory?"))
                if quest in range(1, len(dir_name_frame)+1):
                    break
                else :
                    print("This number doesn't correspond to any of the options.")
                    continue
            except ValueError:
                print("Type only the number, please.")
                continue
        for i in range(0, len(dir_name_frame)):
            if i +1 == quest:
                dir_name_fin = dir_name_frame[i]
    print("The working directory is setted to \'%s\'." % dir_name_fin)
    #print("\n" + dir_name_fin)

    while True:
        quest = input("Do you want to work with one frame (1) or to work with all the pdb located in the working directory (2)? ")
        if quest in ('1', '2'):
            quest = int(quest)
            break
        else :
            print("Type just 1 or 2, please.")
            continue

    if quest == 1:
        repeat = 0
        while repeat == 0:
            while True:
                try :
                    print('\nThose are the saved frames in the working directory:')
                    for i in range(0, len(os.listdir(dir_name_fin))):
                        if str(os.listdir(dir_name_fin)[i])[:6] == 'frame_' and str(os.listdir(dir_name_fin)[i])[-4:] == '.pdb':
                            print('\t' + os.listdir(dir_name_fin)[i]+ '\n')
                    frame_num = int(input("With which frame do you want to work (type the number)? "))
                    if ("frame_%s.pdb" % (str(frame_num))) in os.listdir(dir_name_fin):
                        frame = 'frame_%s.pdb' % frame_num
                        print("The selected frame is %s" % frame)
                        break
                    else :
                        print("This frame is not saved as a pdb.")
                        continue
                except ValueError:
                    print("Type only the number, please.")
                    continue
            while True :
                try :
                    radius = float(input('Which is the desired radius of solvation in (Å)? '))
                    break
                except ValueError:
                    print("Type just de number.")
                    continue
            u_pdb = mda.Universe('%s/%s' % (dir_name_fin,frame))
            while True :
                try :
                    ligand = str(input('Which is the substrate\'s name (type the three-letters code from the pdb): '))
                    if ligand in str(list(u_pdb.residues)[:]):
                        break
                    elif ligand not in str(list(u_pdb.residues)[:]):
                        print("This residue doesn't exist, retype the code, please.")
                        continue
                except ValueError:
                    print("Type just de three-letters code.")
                    continue
            selection = u_pdb.select_atoms(' not ((resname HOH or resname WAT or resname Na+) and (not around %s resname %s))' % (radius, ligand))
            selection.write('%s/qmmm_model_%s' % (dir_name_fin,frame))
            print("The QM/MM model has been saved.")
            while True:
                save = str(input("Do you want to save any other frame (y/[n])?"))
                if save in ('', 'n', 'no', 'N', 'No', 'No', 'nO'):
                    repeat = 2
                    break
                elif save in ('y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs'):
                    repeat = 0
                    break
                else :
                    print('I didn\'t understand you, could you repeat, please?')
                    continue
    if quest == 2:
        frame = os.listdir(dir_name_fin)[0]
        while True :
            try :
                radius = float(input('Which is the desired radius of solvation (in Å)? '))
                break
            except ValueError:
                print("Type just de number.")
                continue
        u_pdb = mda.Universe('%s/%s' % (dir_name_fin,frame))
        while True :
            try :
                ligand = str(input('Which is the substrate\'s name (type the three-letters code from the pdb): '))
                if ligand in str(list(u_pdb.residues)[:]):
                    break
                elif ligand not in str(list(u_pdb.residues)[:]):
                    print("This residue doesn't exist, retype the code, please.")
                    continue
            except ValueError:
                print("Type just de three-letters code.")
                continue
        print("Be patient, please. This will take a while.")
        while True:
            try :
                for i in range(0, len(os.listdir(dir_name_fin)*2)):
                    if os.listdir(dir_name_fin)[i].startswith('frame_') == True:
                        frame = os.listdir(dir_name_fin)[i]
                        u_pdb = mda.Universe('%s/%s' % (dir_name_fin,frame))
                        selection = u_pdb.select_atoms(' not ((resname HOH or resname WAT or resname Na+) and (not around %s resname %s))' % (radius, ligand))
                        selection.write('%s/qmmm_model_%s' % (dir_name_fin,frame))
                        print(frame)
                break
            except IndexError:
                break
        
        print("The QM/MM models have been saved.")
                
        
     # set act creation

# In[265]:


def set_act():
    dir_now = os.getcwd()
    while True:
        try :
            print("Those are the subdirectories in this folder.")
            subdir_list = list(next(os.walk(os.getcwd()))[1])
            subdir_list.append(0)
            #print("\t1. No, it isn't.")
            if len(next(os.walk(os.getcwd()))[1]) != 0:
                for i in range(0, len(next(os.walk(os.getcwd()))[1])):
                    print("\t%s. " % (i+1) + list(next(os.walk(os.getcwd()))[1])[i])
                print("\t%s. It is in the present folder." % str(i+2))
            elif len(next(os.walk(os.getcwd()))[1]) == 0:
                print("\t%s. It is in the present folder." % str(1))
            quest = int(input("Is the selected pdb in any of those subfolders (type the number)? ")) -1
            if quest in range(0, len(next(os.walk(os.getcwd()))[1]) +1):
                break
            else :
                print("This number doesn't correspond to any of the folders. Please, answer again.")
        except ValueError:
            print("Please, answer just the number of the folder.")
            continue
    dir_name = subdir_list[quest]
    if dir_name != 0:
        os.chdir(dir_name)
    while True:
        try :
            frame = int(input("Which is the selected frame? "))
            if "frame_%s.pdb" % frame in list(os.listdir()):
                frame_pdb = 'frame_%s.pdb' % frame
                break
            elif "qmmm_model_frame_%s.pdb" % frame in list(os.listdir()):
                frame_pdb = 'qmmm_model_frame_%s.pdb' % frame
                break
            else :
                print("This frame doesn't exist. Answer again, please.")
                continue
        except ValueError:
            print("Type just the number, please.")
            continue
    
    while True:
        try :
            radius = float(input("Which is the desired radius (in Å)? "))
            break
        except ValueError:
            print("Type just the number, please.")
            continue
    u_pdb = mda.Universe(frame_pdb)
    while True:
        try :
            carbon = int(input("Type the number of the central carbon: "))
            sel_carbon = str(u_pdb.select_atoms("bynum %s" % carbon))
            locA = sel_carbon.find('[<') + 2
            locB = sel_carbon.find(' and segid')
            print('You selected this atom: %s.' % sel_carbon[locA:locB])
            break
        except ValueError:
            print("Type just the number")
    
    while True:
        quest = input("Is this atom correct ([y]/n)? ")
        if quest in ('', 'y', 'yes', 'YES', 'Yes', 'yES', 'yEs'):
            break
        elif quest in ('no', 'n', 'NO', 'No', 'nO') :
            while True:
                try :
                    carbon = int(input("Type the number of the central carbon."))
                    break
                except ValueError :
                    continue
        break        
    os.chdir(dir_now)
    os.getcwd()
    
    selection = u_pdb.select_atoms(str('around %s bynum %s' % (radius, carbon)))
    txt = open('set_%s_%s' % (radius, frame), 'w+')
    txt.write('set act { ')
    for i in range(0, len(selection)):
        locA = str(selection[i]).find('<Atom ') + 6
        locB = str(selection[i]).find(': ')
        txt.write(str(selection[i])[locA:locB] + " ")
    txt.write("} ")
    txt.close()
        
    print("set_%s_%s file created!" % (radius, frame))    
        


# # Program

# In[43]:


### Loop for checking if arguments are given
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



# dinamica = '5_prod_07.nc'
# topology = '3rde_A_HXA_solvatat_ion.prmtop'
# u = mda.Universe(topology, dinamica)
# traj = u.trajectory
# time_counter = 0

# In[18]:


if mode == 'dyn':
    while True:
        menu_dyn()
        menu_quest = input("\nWhat do you want to do (type the number or \"exit\")? ")
        if menu_quest == '1':
            print("\n")
            summarize()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '2':
            print("\n")
            while True:
                c_num = input("How many carbons do you want to analyse in the same plot: 1, 2 or 3? (type \"menu\" if you want to go back) ")
                if c_num == '1':
                    group1()
                    break
                elif c_num == '2':
                    group2()
                    break
                elif c_num == '3':
                    group3()
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
            rmsd_func()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '4':
            print("\n")
            frame_sel()
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '5':
            print("\n")
            model_qmmm()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '6':
            print("\n")
            set_act()
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
            model_qmmm()
            #remember to say that the pdb have to have the name given by this program
            print("\nLet's go back to the menu!")
            continue
        elif menu_quest == '2':
            print("\n")
            set_act()
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

