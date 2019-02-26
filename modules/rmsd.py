#import MDAnalysis as mda
from MDAnalysis import analysis, Merge, Universe
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import align
from matplotlib import pyplot as plt
from modules import loader



def rmsd(u, argsdict=dict({'trajectory': [None, None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None})):

    width_plots = None
    if argsdict['latex'] == True:
        width_plots = argsdict['latex_width']*0.39370079
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'Computer Modern Roman'
        plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

    ref_u = Universe(argsdict['trajectory'][0], argsdict['trajectory'][1:], ref_frame=0)

    rmsddict = {'bb_ref' : None, 'mask' : None, 'subs_ref' : None}
    while True:
        try :
            rmsddict['mask'] = int(input('Do you want to plot the RMSD of the backbone (1), of the substrate (2) or of both (3)? '))
            break
        except TypeError:
            print('Type only the number.')
            continue
        except rmsddict['mask'] not in (1,2,3):
            print('Type only 1, 2 or 3.')
            continue
        
    if rmsddict['mask'] == 1:


        while True:
            try :
                rmsddict['bb_ref'] = int(input('Do you want to set the first frame as the reference (1) or the average structure (2)? '))
                break
            except TypeError:
                print('Type only the number.')
                continue
            except rmsddict['bb_ref'] not in (1,2):
                print('Type only 1 or 2.')
                continue

        if argsdict['timer'] == True:
            import time as timer
            time_in = timer.time()

        if argsdict['u_loaded'] == False:
            u, argsdict = loader.universe_loader_traj(argsdict)

        time = []
        for i in u.trajectory:
            time.append(int(u.trajectory.time)/1000)
        mask = u.select_atoms('backbone')
        ref = ref_u.select_atoms('backbone')

        if rmsddict['bb_ref'] == 1:
            # creation of the array which has the time, the frame and the RMSD value
            R = rms.RMSD(mask, ref).run()
            array = R.rmsd.T

            ## Plot of rmsd
            fig, ax1 = plt.subplots()
            ax1.plot(array[0, :], array[2, :])
            ax1.set_xlabel('Frame')
            if argsdict['latex'] == True:
                ax1.set_ylabel("RMSD ($\AA$)")
            elif argsdict['latex'] == False:
                ax1.set_ylabel("RMSD (Å)")
            ax1.grid(True)
            ax1.legend(['Backbone RMSD\ncompared to the\nfirst frame of\nthe production'], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
            
        
            ax2 = ax1.twiny()
            ax2.plot(time, array[2, :])
            ax2.set_xlabel('Time (ns)')
        
            plt.title("RMSD of the backbone of the protein compared to the first frame of the production", y=1.15, loc='center')
            plt.savefig('%s/plot_rmsd_backbone_1st_frame.png' % argsdict['subdir'], transparent=False, dpi=300, bbox_inches='tight')
            if argsdict['latex'] == True:
                plt.savefig('%s/plot_rmsd_backbone_1st_frame.eps' % argsdict['subdir'], transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
            plt.close()


        elif rmsddict['bb_ref'] == 2:
            # creation of the array which has the time, the frame and the RMSD value
            prealigner = align.AlignTraj(u,reference=ref_u, select='backbone', in_memory=True).run()
            reference_coords = u.trajectory.timeseries(asel=mask).mean(axis=1)
            avg = Merge(mask).load_new(reference_coords[:, None, :], order='afc')

            R = rms.RMSD(mask, avg).run()
            array = R.rmsd.T

            ## Plot of rmsd
            fig, ax1 = plt.subplots()
            ax1.plot(array[0, :], array[2, :])
            ax1.set_xlabel('Frame')
            if argsdict['latex'] == True:
                ax1.set_ylabel("RMSD ($\AA$)")
            elif argsdict['latex'] == False:
                ax1.set_ylabel("RMSD (Å)")
            ax1.grid(True)
            ax1.legend(['Backbone RMSD\ncompared to the\naverage backbone'], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

            ax2 = ax1.twiny()
            ax2.plot(time, array[2, :])
            ax2.set_xlabel('Time (ns)')
        
            plt.title("RMSD of the backbone of the protein compared to the average backbone", y=1.15, loc='center')
            plt.savefig('%s/plot_rmsd_backbone_avg.png' % argsdict['subdir'], transparent=False, dpi=300, bbox_inches='tight')
            if argsdict['latex'] == True:
                plt.savefig('%s/plot_rmsd_backbone_avg.eps' % argsdict['subdir'], transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
            plt.close()

    if rmsddict['mask'] == 2:
        while True:
            try :
                index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space). '))
                #mask = u.select_atoms('resid %s' % index)
                for i in range(0, len(index.split())):
                    int(index.split()[i])
                break
            except ValueError:
                print('Type only the number.')
                continue
    
        u_top = Universe(argsdict['trajectory'][0])
        quest = None

        while True:
            print('You have selected this/these residue/s:')
            i = 0; exists = True
            while i < len(index.split()):
                try : 
                    print(u_top.select_atoms('resid %s' % index).residues[i])
                    i+=1
                except IndexError:
                    print('Residue number %s doesn\'t exist' % index.split()[i])
                    exists = False
                    break
            
            if exists == False:
                while True:
                    try :
                        index = input('Please, retype the numbers. ')
                        #mask = u.select_atoms('resid %s' % index)
                        for i in range(0, len(index.split())):
                            int(index.split()[i])
                        break
                    except ValueError:
                        print('Type only the number.')
                        continue
                continue     
            elif exists != False:
                quest = input('Are the selected residues correct ([y]/n)? ')
                if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                    break
                if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                    while True:
                        try :
                            index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space). '))
                            #mask = u.select_atoms('resid %s' % index)
                            for i in range(0, len(index.split())):
                                int(index.split()[i])
                            break
                        except ValueError:
                            print('Type only the number.')
                            continue
                    continue


        while True:
            try :
                rmsddict['subs_ref'] = int(input('Do you want to set the first frame as the reference (1) or the average structure (2)? '))
                break
            except TypeError:
                print('Type only the number.')
                continue
            except rmsddict['subs_ref'] not in (1,2):
                print('Type only 1 or 2.')
                continue
        if argsdict['timer'] == True:
            import time as timer
            time_in = timer.time()

        if argsdict['u_loaded'] == False:
            u, argsdict = loader.universe_loader_traj(argsdict)

        time = []
        for i in u.trajectory:
            time.append(int(u.trajectory.time)/1000)
        mask = u.select_atoms('resid %s' % index)
        ref = ref_u.select_atoms('resid %s' % index)

        mask_names = []
        for i in range(0, len(mask.residues)):
            mask_names.append(str(str(mask.residues[i])[9:12]+'-'+str(mask.residues[i])[14:-1]))

        if rmsddict['subs_ref'] == 1:
            # creation of the array which has the time, the frame and the RMSD value
            R = rms.RMSD(mask, ref).run()
            array = R.rmsd.T

            ## Plot of rmsd
            fig, ax1 = plt.subplots()
            ax1.plot(array[0, :], array[2, :])
            ax1.set_xlabel('Frame')
            if argsdict['latex'] == True:
                ax1.set_ylabel("RMSD ($\AA$)")
            elif argsdict['latex'] == False:
                ax1.set_ylabel("RMSD (Å)")
            ax1.grid(True)
            ax1.legend(['RMSD of residue(s) (%s)\ncompared to the\nfirst frame of\nthe production' % (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
            
        
            ax2 = ax1.twiny()
            ax2.plot(time, array[2, :])
            ax2.set_xlabel('Time (ns)')
        
            plt.title("RMSD of the residue(s) (%s) of the protein compared to the first frame of the production" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
            plt.savefig('%s/plot_rmsd_residue_1st_frame.png' % argsdict['subdir'], transparent=False, dpi=300, bbox_inches='tight')
            if argsdict['latex'] == True:
                plt.savefig('%s/plot_rmsd_residue_1st_frame.eps' % argsdict['subdir'], transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
            plt.close()


        elif rmsddict['subs_ref'] == 2:
            # creation of the array which has the time, the frame and the RMSD value
            prealigner = align.AlignTraj(u,reference=ref_u, select='resid %s' % index, in_memory=True).run()
            reference_coords = u.trajectory.timeseries(asel=mask).mean(axis=1)
            avg = Merge(mask).load_new(reference_coords[:, None, :], order='afc')

            R = rms.RMSD(mask, avg).run()
            array = R.rmsd.T

            ## Plot of rmsd
            fig, ax1 = plt.subplots()
            ax1.plot(array[0, :], array[2, :])
            ax1.set_xlabel('Frame')
            if argsdict['latex'] == True:
                ax1.set_ylabel("RMSD ($\AA$)")
            elif argsdict['latex'] == False:
                ax1.set_ylabel("RMSD (Å)")
            ax1.grid(True)
            ax1.legend(['RMSD of the residue(s) (%s)\ncompared to the\nfirst frame of\nthe production'% (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

            ax2 = ax1.twiny()
            ax2.plot(time, array[2, :])
            ax2.set_xlabel('Time (ns)')
        
            plt.title("RMSD of the residue(s) (%s) of the protein compared to the average backbone" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
            plt.savefig('%s/plot_rsidue_backbone_avg.png' % argsdict['subdir'], transparent=False, dpi=300, bbox_inches='tight')
            if argsdict['latex'] == True:
                plt.savefig('%s/plot_residue_backbone_avg.eps' % argsdict['subdir'], transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
            plt.close()

    if rmsddict['mask'] == 3:
        
        while True:
            try :
                index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space). '))
                #mask = u.select_atoms('resid %s' % index)
                for i in range(0, len(index.split())):
                    int(index.split()[i])
                break
            except ValueError:
                print('Type only the number.')
                continue
    
        u_top = Universe(argsdict['trajectory'][0])
        quest = None

        while True:
            print('You have selected this/these residue/s:')
            i = 0; exists = True
            while i < len(index.split()):
                try : 
                    print(u_top.select_atoms('resid %s' % index).residues[i])
                    i+=1
                except IndexError:
                    print('Residue number %s doesn\'t exist' % index.split()[i])
                    exists = False
                    break
            
            if exists == False:
                while True:
                    try :
                        index = (input('Please, retype the numbers. '))
                        #mask = u.select_atoms('resid %s' % index)
                        for i in range(0, len(index.split())):
                            int(index.split()[i])
                        break
                    except ValueError:
                        print('Type only the number.')
                        continue
                continue     
            elif exists != False:
                quest = input('Are the selected residues correct ([y]/n)? ')
                if quest in ('', 'y', 'yes', 'Y', 'YES', 'Yes', 'yES', 'YeS', 'yEs', 'YEs', '1'):
                    break
                if quest in ('n', 'no', 'N', 'No', 'No', 'nO', '0'):
                    while True:
                        try :
                            index = (input('Type the number of the residue corresponding to the substrate or to the desired residue (for more than one residue specify the numbers separated by an space). '))
                            #mask = u.select_atoms('resid %s' % index)
                            for i in range(0, len(index.split())):
                                int(index.split()[i])
                            break
                        except ValueError:
                            print('Type only the number.')
                            continue
                    continue 

        while True:
            try :
                rmsddict['bb_ref'] = int(input('For the backbone, do you want to set the first frame as the reference (1) or the average structure (2)? '))
                break
            except TypeError:
                print('Type only the number.')
                continue
            except rmsddict['bb_ref'] not in (1,2):
                print('Type only 1 or 2.')
                continue

        while True:
            try :
                rmsddict['subs_ref'] = int(input('For the substrate or residues, do you want to set the first frame as the reference (1) or the average structure (2)? '))
                break
            except TypeError:
                print('Type only the number.')
                continue
            except rmsddict['subs_ref'] not in (1,2):
                print('Type only 1 or 2.')
                continue

        if argsdict['timer'] == True:
            import time as timer
            time_in = timer.time()

        if argsdict['u_loaded'] == False:
            u, argsdict = loader.universe_loader_traj(argsdict)

        time = []
        for i in u.trajectory:
            time.append(int(u.trajectory.time)/1000)

        mask_bb = u.select_atoms('backbone')
        mask_subs = u.select_atoms('resid %s' % index)
        ref_bb = ref_u.select_atoms('backbone')
        ref_subs = ref_u.select_atoms('resid %s' % index)

        mask_names = []
        for i in range(0, len(mask_subs.residues)):
            mask_names.append(str(str(mask_subs.residues[i])[9:12]+'-'+str(mask_subs.residues[i])[14:-1]))    

        if rmsddict['bb_ref'] == 1:
            # creation of the array which has the time, the frame and the RMSD value
            R = rms.RMSD(mask_bb, ref_bb).run()
            array_bb = R.rmsd.T

        elif rmsddict['bb_ref'] == 2:
            # creation of the array which has the time, the frame and the RMSD value
            prealigner = align.AlignTraj(u,reference=ref_u, select='backbone', in_memory=True).run()
            reference_coords = u.trajectory.timeseries(asel=mask_bb).mean(axis=1)
            avg = Merge(mask_bb).load_new(reference_coords[:, None, :], order='afc')

            R = rms.RMSD(mask_bb, avg).run()
            array_bb = R.rmsd.T
        
        
        if rmsddict['subs_ref'] == 1:
            # creation of the array which has the time, the frame and the RMSD value
            R = rms.RMSD(mask_subs, ref_subs).run()
            array_subs = R.rmsd.T

        elif rmsddict['subs_ref'] == 2:
            # creation of the array which has the time, the frame and the RMSD value
            prealigner = align.AlignTraj(u,reference=ref_u, select='resid %s' % index, in_memory=True).run()
            reference_coords = u.trajectory.timeseries(asel=mask_subs).mean(axis=1)
            avg = Merge(mask_subs).load_new(reference_coords[:, None, :], order='afc')

            R = rms.RMSD(mask_subs, avg).run()
            array_subs = R.rmsd.T


        ## Plot of rmsd
        fig, ax1 = plt.subplots()
        ax1.plot(array_bb[0, :], array_bb[2, :], color='darkblue')
        ax1.plot(array_subs[0, :], array_subs[2, :], color='red')
        ax1.set_xlabel('Frame')
        if argsdict['latex'] == True:
            ax1.set_ylabel("RMSD ($\AA$)")
        elif argsdict['latex'] == False:
            ax1.set_ylabel("RMSD (Å)")
        ax1.grid(True)
        if rmsddict['bb_ref'] == 1 and rmsddict['subs_ref'] == 1:
            ax1.legend(['Backbone RMSD\ncompared to the\nfirst frame of\nthe production', 'RMSD of the residue(s) (%s)\ncompared to the\nfirst frame of\nthe production'% (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        elif rmsddict['bb_ref'] == 1 and rmsddict['subs_ref'] == 2:
            ax1.legend(['Backbone RMSD\ncompared to the\nfirst frame of\nthe production', 'RMSD of the residue(s)\n(%s) of the\nprotein compared to\nthe average backbone' % (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        elif rmsddict['bb_ref'] == 2 and rmsddict['subs_ref'] == 1:
            ax1.legend(['Backbone RMSD\ncompared to the\naverage backbone', 'RMSD of the residue(s)\n(%s) compared to the\nfirst frame of\nthe production'% (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)
        elif rmsddict['bb_ref'] == 2 and rmsddict['subs_ref'] == 2:
            ax1.legend(['Backbone RMSD\ncompared to the\naverage backbone', 'RMSD of the residue(s)\n(%s)\nof the protein compared to\nthe average backbone' % (str(mask_names[:])[1:-1])], bbox_to_anchor=(1.02,1), loc=2, borderaxespad=0.)

        ax2 = ax1.twiny()
        ax2.plot(time, array_bb[2, :], color='darkblue')
        ax2.set_xlabel('Time (ns)')
        if rmsddict['bb_ref'] == 1 and rmsddict['subs_ref'] == 1:
            plt.title("RMSD of the backbone RMSD of the residue(s) (%s)\nof the protein compared to compared to the first frame of the production" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
        elif rmsddict['bb_ref'] == 1 and rmsddict['subs_ref'] == 2:
            plt.title("RMSD of the backbone compared to the first frame and RMSD of the residue(s) (%s)\nof the protein compared to the average" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
        elif rmsddict['bb_ref'] == 2 and rmsddict['subs_ref'] == 1:
            plt.title("RMSD of the backbone compared to the average and RMSD of the residue(s) (%s)\nof the protein compared to the first frame of the production" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
        elif rmsddict['bb_ref'] == 2 and rmsddict['subs_ref'] == 2:
            plt.title("RMSD of the backbone RMSD of the residue(s) (%s)\n of the protein compared to the average structures" % (str(mask_names[:])[1:-1]), y=1.15, loc='center')
    
        plt.savefig('%s/plot_residue_backbone.png' % argsdict['subdir'], transparent=False, dpi=300, bbox_inches='tight')
        if argsdict['latex'] == True:
            plt.savefig('%s/plot_residue_backbone.eps' % argsdict['subdir'], transparent=False, width=width_plots, dpi=300, bbox_inches='tight')
        plt.close()


    ### Time counter ends
    if argsdict['timer'] == True:
        time_fin = timer.time()
        print("I spent " + str(round((time_fin-time_in)/60,1)) + " min")
        
    return u, argsdict

    
        

