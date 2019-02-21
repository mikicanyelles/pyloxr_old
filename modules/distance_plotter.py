'''
Module for plotting the distance between one atom from protein and one, two or three atoms from substrate.
Despite the scripts request three distances (C, Ha and Hb), one can specify the same atom numbers, so the same distance will be computed in the three cases.
'''

#import os
#import sys
from MDAnalysis import Universe
import MDAnalysis.lib.distances as distanceslib
from numpy import max, min, mean, array
from matplotlib import pyplot as plt
from modules import loader


def group1(u, argsdict=dict({'trajectory': [None, None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None})):
#    global u
    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)

    width_plots = None
    if argsdict['latex'] == True:
        width_plots = argsdict['latex_width']*0.39370079
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern Roman'
        plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

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
        u_top = Universe(argsdict['trajectory'][0])
        ### Print selected atoms (number, name, tupe, resname an resid) and save names
        print("\nYou have selected those atoms:\n")

        exists = []
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(u_top.select_atoms("bynum %s" % c9))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        exists.append(a.find('AtomGroup []') != -1)
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

        if True in exists:
            print('Some specified atom doesn\'t exists. Enter the atom numbers again.')
            quest = 'no'
            del exists
        else :
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
        break


    ### Time counter
    if argsdict['timer'] == True:
        import time as timer
        time_in = timer.time()

    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)

    ### Creation of the lists of distances and time
    time = []
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])
        time.append(int(u.trajectory.time)/1000)

    dist_h_9  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])

    ### Avg, min and max distances

    avg_9_c = mean(dist_c_9); min_9_c = min(dist_c_9); max_9_c = max(dist_c_9); rng_9_c = max(dist_c_9) - min(dist_c_9)
    avg_9_h = mean(dist_h_9); min_9_h = min(dist_h_9); max_9_h = max(dist_h_9); rng_9_h = max(dist_h_9) - min(dist_h_9)
    avg_9_c_ar = array([avg_9_c for i in range(0, len(u.trajectory))])
    avg_9_h_ar = array([avg_9_h for i in range(0, len(u.trajectory))])
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

    print("Summary of distances saved")

    ### csv files
    if argsdict['csv'] == True:
        csv = open('%s/distances_%s-%s_%s-%s.csv' % (argsdict['subdir'], t_o_prot, t_c9, t_o_prot, t_h9), 'w')
        csv.write('Frame, Time, %s-%s distance (in Å), %s-%s distance (in Å)\n' %  (t_o_prot, t_c9, t_o_prot, t_h9))
        for i in range(0, len(u.trajectory)):
            csv.write(str(i+1) + ',' + str(time[i]) + ',' + str(dist_c_9[i]) + ',' + str(dist_h_9[i]) + '\n')
        print('csv file saved')

    ### Carbon distances plots

    ##### Histogram
    plt.hist(dist_c_9, bins=20,range=(0,10), histtype='bar')
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s" % (t_c9), y=1.08, loc='center')
    plt.savefig('%s/hist_%s.png' % (argsdict['subdir'],t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s.eps' % (argsdict['subdir'],t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of carbons saved')

    ##### plot w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9)
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_9)
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s" % (t_c9), y=1.15, loc='center')
    plt.savefig('%s/plot_%s.png' % (argsdict['subdir'],t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s.eps' % (argsdict['subdir'],t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9)
    ax1.plot(range(0,len(u.trajectory)), avg_9_c_ar, color='blue')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))#, ['0','1k', '2k', '3k', '4k', '5k', '6k', '7k', '8k', '9k', '10k', '11k', '12k'])
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_c_ar, color='blue')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s vs. time and frames with average distances" % (t_c9), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_avg.png' % (argsdict['subdir'],t_c9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_avg.eps' % (argsdict['subdir'],t_c9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of carbons saved")


    ### Hydrogen plots

    ##### Histogram
    plt.hist([dist_h_9], bins=20,range=(0,10), histtype='bar')
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s" % (t_h9), y=1.08, loc='center')
    plt.savefig('%s/hist_%s.png' % (argsdict['subdir'],t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s.eps' %(argsdict['subdir'],t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of hydrogen saved')

    ##### Plot w/o average
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9)
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_9, label='H_15-OH')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of H-OH distances for %s vs. time and frames" % (t_h9), y=1.15, loc='center')
    plt.savefig('%s/plot_%s.png' % (argsdict['subdir'],t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s.eps' % (argsdict['subdir'],t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9)
    ax1.plot(range(0,len(u.trajectory)), avg_9_h_ar)
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_h_ar, color='blue')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of H-OH distances for %s vs. time and frames" % (t_h9), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_avg.png' % (argsdict['subdir'],t_h9), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_avg.eps' %(argsdict['subdir'],t_h9), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of carbons saved")

    ### Time counter ends
    if argsdict['timer'] == True:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

#####################################################################

def group2(u, argsdict=dict({'trajectory': [None, None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None})):
    if argsdict['latex'] == True:
        width_plots = argsdict['latex_width']*0.39370079
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern Roman'
        plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

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
            c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
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
            c12  = int(input("\nNumber of one of the carbons: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12a  = int(input("Number of one of the hydrogens bonded to the previous carbon: ")) #-1 #8808
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12b  = int(input("Number of the other hydrogen bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break


    ### Ask if atom numbers are correct
    while True:
        u_top = Universe(argsdict['trajectory'][0])
        ### Print selected atoms (number, name, tupe, resname an resid) and save names
        print("\nYou have selected those atoms:\n")

        exists = []
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(u_top.select_atoms("bynum %s" % c9))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        exists.append(a.find('AtomGroup []') != -1)
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
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12b))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(list(u_top.select_atoms("bynum %s" % c12)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c12 = str(a[locA:locB])
        t_h12 = str(t_c12.replace('C','H'))

        if True in exists:
            print('Some specified atom doesn\'t exists. Enter the atom numbers again.')
            quest = 'no'
            del exists
        else :
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
                    c12 = int(input("\nNumber of one of the carbons: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12b = int(input("Number of the other hydrogen bonded to the previous carbon: "))
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
        break


    ### Time counter
    if argsdict['timer'] == True:
        import time as timer
        time_in = timer.time()

    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)

    ### Creation of the lists of distances and time
    ##### For C9 and time
    time = []
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])
        time.append(int(u.trajectory.time)/1000)

    dist_h_9  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])

    ### Avg, min and max distances

    avg_9_c = mean(dist_c_9); min_9_c = min(dist_c_9); max_9_c = max(dist_c_9); rng_9_c = max(dist_c_9) - min(dist_c_9)
    avg_9_h = mean(dist_h_9); min_9_h = min(dist_h_9); max_9_h = max(dist_h_9); rng_9_h = max(dist_h_9) - min(dist_h_9)
    avg_9_c_ar = array([avg_9_c for i in range(0, len(u.trajectory))])
    avg_9_h_ar = array([avg_9_h for i in range(0, len(u.trajectory))])

    ##### For C12
    dist_c_12  = []
    dist_ha_12 = []
    dist_hb_12 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c12).positions
        pos_ha = u.select_atoms("bynum %s" % h12a).positions
        pos_hb = u.select_atoms("bynum %s" % h12b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])

    dist_h_12  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_12[i] < dist_hb_12[i]:
            dist_h_12.append(dist_ha_12[i])
        elif dist_ha_12[i] > dist_hb_12[i]:
            dist_h_12.append(dist_hb_12[i])

    ### Avg, min and max distances

    avg_9_c = mean(dist_c_9); min_9_c = min(dist_c_9); max_9_c = max(dist_c_9); rng_9_c = max(dist_c_9) - min(dist_c_9)
    avg_9_h = mean(dist_h_9); min_9_h = min(dist_h_9); max_9_h = max(dist_h_9); rng_9_h = max(dist_h_9) - min(dist_h_9)
    avg_9_c_ar = array([avg_9_c for i in range(0, len(u.trajectory))])
    avg_9_h_ar = array([avg_9_h for i in range(0, len(u.trajectory))])

    avg_12_c = mean(dist_c_12); min_12_c = min(dist_c_12); max_12_c = max(dist_c_12); rng_12_c = max(dist_c_12) - min(dist_c_12)
    avg_12_h = mean(dist_h_12); min_12_h = min(dist_h_12); max_12_h = max(dist_h_12); rng_12_h = max(dist_h_12) - min(dist_h_12)
    avg_12_c_ar = array([avg_12_c for i in range(0, len(u.trajectory))])
    avg_12_h_ar = array([avg_12_h for i in range(0, len(u.trajectory))])

    print("Lists of distances created")


    ### Save summary of distances (avg, max and min)

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

    print("Summary of distances saved")

    ### csv files
    if argsdict['csv'] == True:
        csv = open('%s/distances_%s-%s_%s-%s_%s-%s_%s-%s.csv' % (argsdict['subdir'], t_o_prot, t_c9, t_o_prot, t_h12, t_o_prot, t_c12, t_o_prot, t_h12), 'w')
        csv.write('Frame, Time, %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å)\n' %  (t_o_prot, t_c9, t_o_prot, t_h9, t_o_prot, t_c12, t_o_prot, t_h12))
        for i in range(0, len(u.trajectory)):
            csv.write(str(i+1) + ',' + str(time[i]) + ',' + str(dist_c_9[i]) + ',' + str(dist_h_9[i]) + ',' + str(dist_c_12[i]) + ',' + str(dist_h_12[i]) + '\n')
        print('csv file saved')

    ### Carbon distances plots

    ##### Histogram
    plt.hist([dist_c_9, dist_c_12], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green'])
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s and %s" % (t_c9,t_c12), y=1.08, loc='center')
    plt.savefig('%s/hist_%s_%s.png' % (argsdict['subdir'], t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s_%s.eps' % (argsdict['subdir'], t_c9,t_c12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of carbons saved')

    ##### plot w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_c_12, color='green')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_12, color='green')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s and %s vs. time and frames" % (t_c9,t_c12), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s.png' % (argsdict['subdir'],t_c9, t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s.eps' % (argsdict['subdir'],t_c9, t_c12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_c_12, color='green')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_c_ar, color='purple')
    ax2.plot(time, avg_12_c_ar, color='lime')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s and %s vs. time and frames with average distances" % (t_c9,t_c12), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_avg.png' % (argsdict['subdir'],t_c9,t_c12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_avg.eps' % (argsdict['subdir'],t_c9,t_c12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of carbons saved")


    ### Hydrogen plots

    ##### Histogram
    plt.hist([dist_h_9, dist_h_12], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green'])
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s and %s" % (t_h9,t_h12), y=1.08, loc='center')
    plt.savefig('%s/hist_%s_%s.png' % (argsdict['subdir'],t_h9, t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s_%s.eps' %(argsdict['subdir'],t_h9, t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of hydrogen saved')

    ##### Plot w/o average
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_h_12, color='green')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_12, color='indigo')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of H-OH distances for %s and %s vs. time and frames" % (t_h9,t_h12), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s.png' % (argsdict['subdir'],t_h9, t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s.eps' % (argsdict['subdir'],t_h9, t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_h_12, color='green')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_h_ar, color='purple')
    ax2.plot(time, avg_12_h_ar, color='lime')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of h-OH distances for %s and %s vs. time and frames with average distances" % (t_h9,t_h12), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_avg.png' % (argsdict['subdir'],t_h9,t_h12), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_avg.eps' % (argsdict['subdir'],t_h9,t_h12), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of hydrogens saved")

    ### Time counter ends
    if argsdict['timer'] == True:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

#####################################################################

def group3(u, argsdict=dict({'trajectory': [None, None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None})):
    if argsdict['latex'] == True:
        width_plots = argsdict['latex_width']*0.39370079
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern Roman'
        plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

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
            c9  = int(input("\nNumber of one of the carbons: ")) #-1 #8807
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
            c12  = int(input("\nNumber of one of the carbons: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12a  = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h12b  = int(input("Number of the other hydrogen bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break

    while True:
        try :
            c15  = int(input("\nNumber of one of the carbons: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h15a  = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break
    while True:
        try :
            h15b  = int(input("Number of the other hydrogen bonded to the previous carbon: "))
        except ValueError:
            print("Enter the number again.")
            continue
        else :
            break

    ### Ask if atom numbers are correct
    while True:
        u_top = Universe(argsdict['trajectory'][0])
        ### Print selected atoms (number, name, tupe, resname an resid) and save names
        print("\nYou have selected those atoms:\n")

        exists = []
        a = str(u_top.select_atoms("bynum %s" % o_prot))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(u_top.select_atoms("bynum %s" % c9))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h9b))
        exists.append(a.find('AtomGroup []') != -1)
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
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h12b))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(list(u_top.select_atoms("bynum %s" % c12)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c12 = str(a[locA:locB])
        t_h12 = str(t_c12.replace('C','H'))

        a = str(u_top.select_atoms("bynum %s" % c15))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h15a))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB])
        a = str(u_top.select_atoms("bynum %s" % h15b))
        exists.append(a.find('AtomGroup []') != -1)
        locA = a.find('[<') +2
        locB = a.find(' and')
        print("\t" + a[locA:locB] + "\n")

        a = str(list(u_top.select_atoms("bynum %s" % c15)))
        locA = a.find(': ') +2
        locB = a.find(' of')
        t_c15 = str(a[locA:locB])
        t_h15 = str(t_c15.replace('C','H'))

        if True in exists:
            print('Some specified atom doesn\'t exists. Enter the atom numbers again.')
            quest = 'no'
            del exists
        else :
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
                    c12 = int(input("\nNumber of one of the carbons: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12a = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h12b = int(input("Number of the other hydrogen bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    c15 = int(input("\nNumber of one of the carbons: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h15a = int(input("Number of one of the hydrogens bonded to the previous carbon: "))
                except ValueError:
                    print("Enter the number again.")
                    continue
                else :
                    break
            while True:
                try :
                    h15b = int(input("Number of the other hydrogen bonded to the previous carbon: "))
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
        break


    ### Time counter
    if argsdict['timer'] == True:
        import time as timer
        time_in = timer.time()

    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)
    ### Creation of the lists of distances and time
    ##### For C9 and time
    time = []
    dist_c_9  = []
    dist_ha_9 = []
    dist_hb_9 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c9).positions
        pos_ha = u.select_atoms("bynum %s" % h9a).positions
        pos_hb = u.select_atoms("bynum %s" % h9b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_9.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_9.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_9.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])
        time.append(int(u.trajectory.time)/1000)

    dist_h_9  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_9[i] < dist_hb_9[i]:
            dist_h_9.append(dist_ha_9[i])
        elif dist_ha_9[i] > dist_hb_9[i]:
            dist_h_9.append(dist_hb_9[i])

    ### Avg, min and max distances

    avg_9_c = mean(dist_c_9); min_9_c = min(dist_c_9); max_9_c = max(dist_c_9); rng_9_c = max(dist_c_9) - min(dist_c_9)
    avg_9_h = mean(dist_h_9); min_9_h = min(dist_h_9); max_9_h = max(dist_h_9); rng_9_h = max(dist_h_9) - min(dist_h_9)
    avg_9_c_ar = array([avg_9_c for i in range(0, len(u.trajectory))])
    avg_9_h_ar = array([avg_9_h for i in range(0, len(u.trajectory))])

    ##### For C12
    dist_c_12  = []
    dist_ha_12 = []
    dist_hb_12 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c12).positions
        pos_ha = u.select_atoms("bynum %s" % h12a).positions
        pos_hb = u.select_atoms("bynum %s" % h12b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_12.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_12.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_12.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])

    dist_h_12  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_12[i] < dist_hb_12[i]:
            dist_h_12.append(dist_ha_12[i])
        elif dist_ha_12[i] > dist_hb_12[i]:
            dist_h_12.append(dist_hb_12[i])

    ##### For C15
    dist_c_15  = []
    dist_ha_15 = []
    dist_hb_15 = []
    for i in u.trajectory:
        pos_c  = u.select_atoms("bynum %s" % c15).positions
        pos_ha = u.select_atoms("bynum %s" % h15a).positions
        pos_hb = u.select_atoms("bynum %s" % h15b).positions
        pos_o  = u.select_atoms("bynum %s" % o_prot).positions
        if argsdict['parallel'] == False:
            dist_ha_15.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='serial')[0])
            dist_hb_15.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='serial')[0])
            dist_c_15.append(distanceslib.calc_bonds(pos_o, pos_c, backend='serial')[0])
        elif argsdict['parallel'] == True:
            dist_ha_15.append(distanceslib.calc_bonds(pos_o, pos_ha, backend='OpenMP')[0])
            dist_hb_15.append(distanceslib.calc_bonds(pos_o, pos_hb, backend='OpenMP')[0])
            dist_c_15.append(distanceslib.calc_bonds(pos_o, pos_c, backend='OpenMP')[0])

    dist_h_15  = []
    for i in range(0, len(u.trajectory)):
        if dist_ha_15[i] < dist_hb_15[i]:
            dist_h_15.append(dist_ha_15[i])
        elif dist_ha_15[i] > dist_hb_15[i]:
            dist_h_15.append(dist_hb_15[i])

    ### Avg, min and max distances

    avg_9_c = mean(dist_c_9); min_9_c = min(dist_c_9); max_9_c = max(dist_c_9); rng_9_c = max(dist_c_9) - min(dist_c_9)
    avg_9_h = mean(dist_h_9); min_9_h = min(dist_h_9); max_9_h = max(dist_h_9); rng_9_h = max(dist_h_9) - min(dist_h_9)
    avg_9_c_ar = array([avg_9_c for i in range(0, len(u.trajectory))])
    avg_9_h_ar = array([avg_9_h for i in range(0, len(u.trajectory))])

    avg_12_c = mean(dist_c_12); min_12_c = min(dist_c_12); max_12_c = max(dist_c_12); rng_12_c = max(dist_c_12) - min(dist_c_12)
    avg_12_h = mean(dist_h_12); min_12_h = min(dist_h_12); max_12_h = max(dist_h_12); rng_12_h = max(dist_h_12) - min(dist_h_12)
    avg_12_c_ar = array([avg_12_c for i in range(0, len(u.trajectory))])
    avg_12_h_ar = array([avg_12_h for i in range(0, len(u.trajectory))])

    avg_15_c = mean(dist_c_15); min_15_c = min(dist_c_15); max_15_c = max(dist_c_15); rng_15_c = max(dist_c_15) - min(dist_c_15)
    avg_15_h = mean(dist_h_15); min_15_h = min(dist_h_15); max_15_h = max(dist_h_15); rng_15_h = max(dist_h_15) - min(dist_h_15)
    avg_15_c_ar = array([avg_15_c for i in range(0, len(u.trajectory))])
    avg_15_h_ar = array([avg_15_h for i in range(0, len(u.trajectory))])

    print("Lists of distances created")


    ### Save summary of distances (avg, max and min)

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

    txt.write('Distances for %s: \n' % t_c15)
    txt.write('\tAverage distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(avg_15_c,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(min_15_c,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_c15,t_o_prot) + str(round(max_15_c,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_c15,t_o_prot) + str(round(rng_15_c,3)) + ' Å \n\n')

    txt.write('Distances for the nearest %s, even if it is the A or the B: \n'% t_h15 )
    txt.write('\tAverage distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(avg_15_h,3)) + ' Å \n')
    txt.write('\tMinimum distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(min_15_h,3)) + ' Å \n')
    txt.write('\tMaximum distance %s - %s:      ' % (t_h15,t_o_prot) + str(round(max_15_h,3)) + ' Å \n')
    txt.write('\tMax-min for distances %s - %s: ' % (t_h15,t_o_prot) + str(round(rng_15_h,3)) + ' Å')

    txt.close()

    print("Summary of distances saved")

    ### csv files
    if argsdict['csv'] == True:
        csv = open('%s/distances_%s-%s_%s-%s_%s-%s_%s-%s.csv' % (argsdict['subdir'], t_o_prot, t_c9, t_o_prot, t_h12, t_o_prot, t_c15, t_o_prot, t_h15), 'w')
        csv.write('Frame, Time, %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å), %s-%s distance (in Å)\n' %  (t_o_prot, t_c9, t_o_prot, t_h9, t_o_prot, t_c12, t_o_prot, t_h12, t_o_prot, t_c15, t_o_prot, t_h15))
        for i in range(0, len(u.trajectory)):
            csv.write(str(i+1) + ',' + str(time[i]) + ',' + str(dist_c_9[i]) + ',' + str(dist_h_9[i]) + ',' + str(dist_c_12[i]) + ',' + str(dist_h_12[i]) + ',' + str(dist_c_15[i]) + ',' + str(dist_h_15[i]) + '\n')
        print('csv file saved')

    ### Carbon distances plots

    ##### Histogram
    plt.hist([dist_c_9, dist_c_12, dist_c_15], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green', 'coral'])
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_c15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of C-OH distances for %s, %s and %s" % (t_c9,t_c12,t_c15), y=1.08, loc='center')
    plt.savefig('%s/hist_%s_%s_%s.png' % (argsdict['subdir'], t_c9,t_c12, t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s_%s_%s.eps' % (argsdict['subdir'], t_c9,t_c12, t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of carbons saved')

    ##### plot w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_c_12, color='green')
    ax1.plot(range(0,len(u.trajectory)), dist_c_15, color='coral')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_c15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_c_15, color='coral')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s, %s and %s vs. time and frames" % (t_c9,t_c12,t_c15), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_%s.png' % (argsdict['subdir'],t_c9, t_c12, t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_%s.eps' % (argsdict['subdir'],t_c9, t_c12, t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_c_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_c_12, color='green')
    ax1.plot(range(0,len(u.trajectory)), dist_c_15, color='coral')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_c9,t_o_prot), '%s-%s' % (t_c12,t_o_prot), '%s-%s' % (t_h12,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_c_ar, color='purple')
    ax2.plot(time, avg_12_c_ar, color='lime')
    ax2.plot(time, avg_15_c_ar, color='orangered')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of C-OH distances for %s, %s and %s vs. time and frames" % (t_c9,t_c12,t_c15), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_%s_avg.png' % (argsdict['subdir'],t_c9,t_c12, t_c15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_%s_avg.eps' % (argsdict['subdir'],t_c9,t_c12, t_c15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of carbons saved")


    ### Hydrogen plots

     ##### Histogram
    plt.hist([dist_c_9, dist_c_12, dist_c_15], bins=20,range=(0,10), histtype='bar', color = ['indigo', 'green', 'coral'])
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,11))
    plt.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of H-OH distances for %s, %s and %s" % (t_h9,t_h12,t_h15), y=1.08, loc='center')
    plt.savefig('%s/hist_%s_%s_%s.png' % (argsdict['subdir'], t_h9,t_h12, t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s_%s_%s.eps' % (argsdict['subdir'], t_h9,t_h12, t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print('Histogram of carbons saved')

    ##### plot w/o avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_h_12, color='green')
    ax1.plot(range(0,len(u.trajectory)), dist_h_15, color='coral')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, dist_h_15, color='coral')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of H-OH distances for %s, %s and %s vs. time and frames" % (t_h9,t_h12,t_h15), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_%s.png' % (argsdict['subdir'],t_h9, t_h12, t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_%s.eps' % (argsdict['subdir'], t_h9, t_h12, t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    ##### plot w/ avg
    fig, ax1 = plt.subplots()
    ax1.plot(range(0,len(u.trajectory)), dist_h_9, color='indigo')
    ax1.plot(range(0,len(u.trajectory)), dist_h_12, color='green')
    ax1.plot(range(0,len(u.trajectory)), dist_h_15, color='coral')
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax1.grid(True)
    ax1.legend(['%s-%s' % (t_h9,t_o_prot), '%s-%s' % (t_h12,t_o_prot), '%s-%s' % (t_h15,t_o_prot)], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

    ax2 = ax1.twiny()
    ax2.plot(time, avg_9_h_ar, color='purple')
    ax2.plot(time, avg_12_h_ar, color='lime')
    ax2.plot(time, avg_15_h_ar, color='orangered')
    ax2.set_xlabel('Time (ns)')

    plt.title("Plot of H-OH distances for %s, %s and %s vs. time and frames" % (t_h9,t_h12,t_h15), y=1.15, loc='center')
    plt.savefig('%s/plot_%s_%s_%s_avg.png' % (argsdict['subdir'],t_h9,t_h12, t_h15), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s_%s_%s_avg.eps' % (argsdict['subdir'],t_h9,t_h12, t_h15), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots of hydrogens saved")

    ### Time counter ends
    if argsdict['timer'] == True:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")

#####################################################################

    return u, argsdict
