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