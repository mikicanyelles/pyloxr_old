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