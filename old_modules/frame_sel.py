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