from MDAnalysis import Universe
import MDAnalysis.lib.distances as distanceslib
from numpy import max, min, mean, array
#from sys import platform as sys_pf
import sys
if sys.platform == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from progressbar import *

#argsdict=dict({'trajectory': ['3rde.prmtop', None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': 'plots', 'timer': False, 'menu_type' : None, 'u_loaded' : False})
def general_plotter(u, argsdict=dict({'trajectory': [None, None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None, 'u_loaded' : False})):
    #if argsdict['u_loaded'] == False:
    #    from modules import loader
    #    u, argsdict = loader.universe_loader_traj(argsdict)

    if argsdict['latex'] == True:
        width_plots = argsdict['latex_width']*0.39370079
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern Roman'
        plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

    ### Ask for atoms lists
    while True:
        try :
            prot_atoms_input = input("Input the number of the selected atoms for measuring distances (separated by a space). ")
            if prot_atoms_input.find(' ') != -1 and prot_atoms_input.find(',') == -1:
                prot_atoms = prot_atoms_input.split(' ')
                if len(prot_atoms) != len(set(prot_atoms)):
                    print("Some atom numbers are duplicated, this has been fixed.")
                    prot_atoms = list(set(prot_atoms))
                for i in range(0,len(prot_atoms)):
                    prot_atoms[i] = [int(prot_atoms[i])]
                break

            elif prot_atoms_input.find(' ') == -1 and prot_atoms_input.find(',') == -1:
                prot_atoms = [[int(prot_atoms_input)]]
                break

            elif prot_atoms_input.find(',') != -1:
                print('Comparisons have to be specified on the second list of atoms.')
                continue

        except ValueError:
            print('Some of the atom\'s numbers has not been correctly introduced.\n')
            continue

    while True:
        try :
            subs_atoms_input = input("Input the number of the atoms for measuring distances (separated by a space or by a comma if you want that only the nearest gets saved). ")
            if subs_atoms_input.find(' ') != -1 or subs_atoms_input.find(',') != -1:
                if subs_atoms_input.find(',') == -1:
                    comp = False
                    subs_atoms = subs_atoms_input.split(' ')

                    if len(subs_atoms) != len(set(subs_atoms)):
                        print("Some atom numbers are duplicated, this has been fixed.")
                        subs_atoms = list(set(subs_atoms))
                    for i in range(0,len(subs_atoms)):
                        subs_atoms[i] = [int(subs_atoms[i])]

                elif subs_atoms_input.find(',') != -1:
                    comp = True
                    subs_atoms = subs_atoms_input.split(' ')

                    if len(subs_atoms) != len(set(subs_atoms)):
                        print("Some atom numbers are duplicated, this has been fixed.")
                        subs_atoms = list(set(subs_atoms))

                    for i in range(0, len(subs_atoms)):
                        subs_atoms[i] = subs_atoms[i].split(',')
                        if len(subs_atoms[i]) != len(set(subs_atoms[i])):
                            print("Some atom numbers for comparisons are duplicated, this has been fixed.")
                            subs_atoms[i] = list(set(subs_atoms[i]))
                        for j in range(0, len(subs_atoms[i])):
                            subs_atoms[i][j] = int(subs_atoms[i][j])
                break
            elif subs_atoms_input.find(' ') == -1:
                subs_atoms = [[int(subs_atoms_input)]]
                break
        except ValueError:
            print('Some of the atom\'s numbers has not been correctly introduced.\n')
            continue

    u_top = Universe(argsdict['trajectory'][0])
    quest = None

    ### Check if atom numbers are correct
    while True:
        print("\nYou have selected those atoms:\n")
        ### Print selected atoms (number, name, type, resname an resid) and save names
        exists = []
        for i in range(0, len(prot_atoms)):
            a = str(u_top.select_atoms("bynum %s" % prot_atoms[i][0]))
            exists.append(a.find('AtomGroup []') != -1)
            locA = a.find('[<') +2
            locB = a.find(' and')
            print('\t' + a[locA:locB])
        print("\n")

        for i in range(0, len(subs_atoms)):
            name = None
            for j in range(0, len(subs_atoms[i])):
                a = str(u_top.select_atoms("bynum %s" % subs_atoms[i][j]))
                exists.append(a.find('AtomGroup []') != -1)
                locA = a.find('[<') +2
                locB = a.find(' and')
                if name == None:
                    name = '\t' + a[locA:locB]
                else:
                    name = name + '; ' + a[locA:locB]
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
                    prot_atoms_input = input("Input the number of the selected atoms for measuring distances (separated by a space). ")
                    if prot_atoms_input.find(' ') != -1 and prot_atoms_input.find(',') == -1:
                        prot_atoms = prot_atoms_input.split(' ')
                        if len(prot_atoms) != len(set(prot_atoms)):
                            print("Some atom numbers are duplicated, this has been fixed.")
                            prot_atoms = list(set(prot_atoms))
                        for i in range(0,len(prot_atoms)):
                            prot_atoms[i] = [int(prot_atoms[i])]
                        break

                    elif prot_atoms_input.find(' ') == -1 and prot_atoms_input.find(',') == -1:
                        prot_atoms = [[int(prot_atoms_input)]]
                        break

                    elif prot_atoms_input.find(',') != -1:
                        print('Comparisons have to be specified on the second list of atoms.')
                        continue

                except ValueError:
                    print('Some of the atom\'s numbers has not been correctly introduced.\n')
                    continue

            while True:
                try :
                    subs_atoms_input = input("Input the number of the atoms for measuring distances (separated by a space or by a comma if you want that only the nearest gets saved). ")
                    if subs_atoms_input.find(' ') != -1 or subs_atoms_input.find(',') != -1:
                        if subs_atoms_input.find(',') == -1:
                            comp = False
                            subs_atoms = subs_atoms_input.split(' ')

                            if len(subs_atoms) != len(set(subs_atoms)):
                                print("Some atom numbers are duplicated, this has been fixed.")
                                subs_atoms = list(set(subs_atoms))
                            for i in range(0,len(subs_atoms)):
                                subs_atoms[i] = [int(subs_atoms[i])]

                        elif subs_atoms_input.find(',') != -1:
                            comp = True
                            subs_atoms = subs_atoms_input.split(' ')

                            if len(subs_atoms) != len(set(subs_atoms)):
                                print("Some atom numbers are duplicated, this has been fixed.")
                                subs_atoms = list(set(subs_atoms))

                            for i in range(0, len(subs_atoms)):
                                subs_atoms[i] = subs_atoms[i].split(',')
                                if len(subs_atoms[i]) != len(set(subs_atoms[i])):
                                    print("Some atom numbers for comparisons are duplicated, this has been fixed.")
                                    subs_atoms[i] = list(set(subs_atoms[i]))
                                for j in range(0, len(subs_atoms[i])):
                                    subs_atoms[i][j] = int(subs_atoms[i][j])
                        break
                    elif subs_atoms_input.find(' ') == -1:
                        subs_atoms = [[int(subs_atoms_input)]]
                        break
                except ValueError:
                    print('Some of the atom\'s numbers has not been correctly introduced.\n')
                    continue
            continue

        elif quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
            break
        else :
            print("Sorry, answer again, please.")
            continue
        break

    if argsdict['u_loaded'] == False:
        from modules import loader
        u, argsdict = loader.universe_loader_traj(argsdict)

    distances = []; labels = []; l = 0
    time = [[],[]]; time_ = 0

    widgets = ['Progress: ', Percentage(), ' ', Bar(marker='#',left='[',right=']'),
           ' ', ETA(), ' | ', Timer(), ' '] #see docs for other options

    pbar = ProgressBar(widgets=widgets, maxval=len(u.trajectory))
    pbar.start()

    for ts in u.trajectory:
        k = 0
        time[0].append(time_); time_+=1
        time[1].append(int(u.trajectory.time)/1000)
        for i in range(0, len(prot_atoms)):
            pos_prot = u.select_atoms("bynum %s" % prot_atoms[i][0]).positions
            for j in range(0, len(subs_atoms)):
                if len(subs_atoms[j]) == 1:
                    pos_subs = u.select_atoms("bynum %s" % subs_atoms[j][0]).positions
                    if argsdict['parallel'] == False:
                        if l==0:
                            distances.append([])
                            labels.append([])
                            a = str(u_top.select_atoms("bynum %s" % prot_atoms[i][0]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            labels[k] = str(a[locA:locB])
                            a = str(u_top.select_atoms("bynum %s" % subs_atoms[j][0]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            labels[k]=labels[k]+'-'+str(a[locA:locB])
                        else :
                            pass
                        distances[k].append(float(distanceslib.calc_bonds(pos_prot, pos_subs, backend='serial')[0]))

                    elif argsdict['parallel'] == True:
                        if l==0:
                            distances.append([])
                            labels.append([])
                            a = str(u_top.select_atoms("bynum %s" % prot_atoms[i][0]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            labels[k] = str(a[locA:locB])
                            a = str(u_top.select_atoms("bynum %s" % subs_atoms[j][0]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            labels[k]=labels[k]+'-'+str(a[locA:locB])
                        else :
                            pass
                        distances[k].append(float(distanceslib.calc_bonds(pos_prot, pos_subs, backend='OpenMP')[0]))
                    k+=1
                elif len(subs_atoms[j]) > 1:
                    dist = None
                    if l==0:
                        labels.append([])
                        a = str(u_top.select_atoms("bynum %s" % prot_atoms[i][0]))
                        locA = a.find(': ')+2
                        locB = a.find(' of')
                        labels[k] = str(a[locA:locB])
                        m=0
                    for jj in range(len(subs_atoms[j])):
                        pos_subs = u.select_atoms('bynum %s' % subs_atoms[j][jj]).positions
                        if argsdict['parallel'] == False:
                            dist_ = distanceslib.calc_bonds(pos_prot, pos_subs, backend='serial')[0]
                            a = str(u_top.select_atoms("bynum %s" % subs_atoms[j][jj]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            if l==0 and m==0:
                                labels[k]=labels[k]+'-'+str(a[locA:locB])
                                m=1
                            elif l==0 and m!=0:
                                labels[k]=labels[k]+','+str(a[locA:locB])
                        elif argsdict['parallel'] == True:
                            dist_ = distanceslib.calc_bonds(pos_prot, pos_subs, backend='OpenMP')[0]
                            a = str(u_top.select_atoms("bynum %s" % subs_atoms[j][jj]))
                            locA = a.find(': ') +2
                            locB = a.find(' of')
                            if l==0 and m==0:
                                labels[k]=labels[k]+'-'+str(a[locA:locB])
                                m=1
                            elif l==0 and m!=0:
                                labels[k]=labels[k]+','+str(a[locA:locB])
                        if dist != None:
                            if dist_ < dist:
                                dist = dist_
                        elif dist == None:
                            dist = dist_
                    if l==0:
                        distances.append([])
                    else :
                        pass
                    distances[k].append(dist)
                    k+=1
        l=1
        pbar.update(time_)

    pbar.finish()

    label_str = None
    for i in range(0, len(labels)):
        if label_str == None:
            label_str = labels[i]
        else :
            label_str = label_str + ',' + str(labels[i])

    txt = open('summary_of_distances_%s.txt' % label_str, 'w')
    _max = 0
    for i in range(0, len(distances)):
        max_ = max(distances[i])
        min_ = min(distances[i])
        mean_ = mean(distances[i])
        range_ = max(distances[i]) - min(distances[i])
        if max_ > _max:
            _max = max_ #necessary for histograms

        txt.write("Distances for %s: \n" % labels[i])
        txt.write('\tAverage distance:' +'\t'+str(round(mean_,3)) + ' Å \n')
        txt.write('\tMaximum distance:' +'\t'+str(round(max_,3)) + ' Å \n')
        txt.write('\tMinimum distance:' +'\t'+str(round(min_,3)) + ' Å \n')
        txt.write('\tRange of distances:' +'\t'+str(round(range_,3)) + ' Å \n')
        txt.write('\n')
    txt.close(); del txt

    #csv creation
    if argsdict['csv'] == True:
        csv = open('distances_%s.csv' % label_str, 'w')
        csv.write('Frame, Time (ns), ')
        for i in range(0,len(distances)):
            csv.write('Distance ' + str(labels[i])+' (Å), ')
        csv.write('\n')
        for i in range(0,len(time[0])):
            csv.write(str(time[0][i]+1)+', '+ str(time[1][i]))
            for j in range(0,len(distances)):
                csv.write(', '+str(round(distances[j][i],3)))
            csv.write('\n')
        csv.close(); del csv

    #histograms
    bins_ = int((round(_max/10)*10)*2)
    if bins_ < 20:
        bins_ = 20
    elif bins_ % 2 != 0:
        bins_ = bins_ +1
    plt.hist((distances[:]), bins=(bins_), range=(0,bins_/2), histtype='bar')
    if argsdict['latex'] == True:
        plt.xlabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        plt.xlabel('Distance (Å)')
    plt.ylabel('Number of frames')
    plt.xticks(range(0,int(bins_/2+1)))
    plt.legend(labels[:], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
    plt.title("Histogram of %s distances" % (label_str), y=1.08, loc='center')
    plt.savefig('%s/hist_%s.png' % (argsdict['subdir'],label_str), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/hist_%s.eps' %(argsdict['subdir'],label_str), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    #plt.show()
    plt.close()

    #plots
    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()
    ax1.grid(axis='both')
    for i in range(0, len(distances)):
        ax1.plot(time[0], distances[i])
    ax1.legend(labels[:])
    ax1.set_xlabel('Frame')
    if argsdict['latex'] == True:
        ax1.set_ylabel('Distance ($\AA$)')
    elif argsdict['latex'] == False:
        ax1.set_ylabel('Distance (Å)')
    ax1.set_xticks(range(0,len(u.trajectory) +1, round(len(u.trajectory)/10)))
    ax2.set_xlabel('Time (ns)')
    ax2.plot(time[1], distances[-1], c='C%s' % int(len(distances)-1))
    plt.title("Plot of %s distances" % label_str, y=1.15, loc='center')
    plt.savefig('%s/plot_%s.png' % (argsdict['subdir'],label_str), transparent=False, dpi=300, bbox_inches='tight')
    if argsdict['latex'] == True:
        plt.savefig('%s/plot_%s.eps' % (argsdict['subdir'],label_str), width=width_plots, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots have been saved. Let's go back to menu!")
    return u, argsdict
