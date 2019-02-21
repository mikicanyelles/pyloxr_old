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