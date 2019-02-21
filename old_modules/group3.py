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