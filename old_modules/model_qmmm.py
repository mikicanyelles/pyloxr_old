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
        dir_name_fin = dir_name_frame[1]
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