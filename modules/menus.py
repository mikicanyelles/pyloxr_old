'''
Module for loading the menus which reveal the available modules and options.
'''
import sys

# # Menu

# ### For dynamics files

def menu(argsdict=dict({'trajectory': [None, None], 'frame': None, 'csv' : False, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': '.', 'timer': False, 'menu_type' : None})):
    if argsdict['menu_type'] == None:
        while True:
            try :
                argsdict['menu_type'] = int(input('Which menu do you want to display? (Type (1) for the menu for analysing trajectories or (2) for working with a single frame)'))
                break
            except ValueError :
                print('Type only 1 or 2.')
                continue
    def trajectory():
        print("******************************************************************************************")
        print("Here you have a list with the modules you can choose:")
        print("\t1. Summary of the system and the MD trajectory.")
        print("\t2. Obtain distance (between atoms) plots.")
        print("\t3. Obtain the plot of RMSD of the backbone and/or the substrate or selected residues.")
        print("\t4. Select frames by H(subs)-protein distance.")
        #print("\t5. QM/MM models from frames creation. (Frames have to be saved as pdb by this program (option 4))")
        #print("\t6. Create the 'set act' file, where active atoms for ChemShell are specified")
        print("\n\t0. Exit")
        print("\t*. Options")
        print("\t*. Manual")
        print("******************************************************************************************")
    def frame():
        print("\n******************************************************************************************")
        print("\nHere you have a list with the modules you can choose:")
        #print("\nYou are in the 'nodynamics' mode. You can just do the following tasks.")
        print("If you want to do other jobs that require the topology and the trajectory, type 'exit' and rerun the script specifing those files.")
        print("By now, modules regarding one frame are not available")
        #print("\nHere you have a list with the options you can choose:")
        #print("\t1. QM/MM models from frames creation)")
        #print("\t2. Create the 'set act' file, where active atoms for ChemShell are specified")
        print("\n******************************************************************************************")

    def options_menu(argsdict):
        print('\nOptions:')
        print('\tSwitchers: (selecting them will (de)activate the option)')
        if argsdict['csv'] == True:
            print('\t\t-c,  --csv\t\t->\tActivated')
        elif argsdict['csv'] == False:
            print('\t\t-c,  --csv\t\t->\tDeactivated')
        if argsdict['latex'] == True:
            print('\t\t-l,  --latex\t\t->\tActivated')
        elif argsdict['latex'] == False:
            print('\t\t-l,  --latex\t\t->\tDeactivated')
        if argsdict['parallel'] == True:
            print('\t\t-p,  --parallel\t\t->\tActivated')
        elif argsdict['parallel'] == False:
            print('\t\t-p,  --parallel\t\t->\tDeactivated')
        #if argsdict['timer'] == True:
        #    print('\t\t-tm, --timer\t\t->\tActivated')
        #if argsdict['timer'] == False:
        #    print('\t\t-tm, --timer\t\t->\tDeactivated')

        print('\n\tValues: (select them and type the new value which will be stored)')
        if argsdict['latex_width'] == None:
            print('\t\t-lw, --latex-width LATEX_WIDTH (in cm)')
        elif argsdict['latex_width'] != None:
            print('\t\t-lw, --latex-width LATEX_WIDTH (in cm)\t->\tThe currently width for plots is %s cm' % float(argsdict['latex_width']))
        if argsdict['subdir'] == '.':
            print('\t\t-s,  --subdir SUBDIR\t\t\t->\tThe directory for plots is the current directory.')
        elif argsdict['subdir'] != '.':
            print('\t\t-s,  --subdir SUBDIR\t\t\t->\tThe directory for plots is %s.' % argsdict['subdir'])

        print("\n******************************************************************************************\n")

    if argsdict['menu_type'] != None:
        if argsdict['menu_type'] == 1:
            trajectory()
        elif argsdict['menu_type'] == 2:
            frame()

        while True :
            option = input('Type the number of the module, the flag of the option, \'menu\' to reprint de menu, \'options\' to show the options menu or \'exit\'. ')
            if option in ('options', 'option', 'Option', 'Options', 'OPTION', 'OPTIONS', 'oPTIONS', 'oPTION'):
                options_menu(argsdict)

            elif option == 'menu':
                if argsdict['menu_type'] == 1:
                    trajectory()
                elif argsdict['menu_type'] == 2:
                    frameº()

            elif option in ('-c', '--csv', '-l', '--latex', '-p', '--parallel') or (option.find('-lw') == 0 or option.find('--latex-width') == 0 or option.find('-s') == 0 or option.find('--subdir') == 0):
                if option in ('-c', '--csv'):
                    if argsdict['csv'] == True:
                        argsdict['csv'] = False
                        print('csv files saving has been turned off.')
                    elif argsdict['csv'] == False:
                        argsdict['csv'] = True
                        print('csv files saving has been turned on.')
                elif option in ('-l', '--latex'):
                    if argsdict['latex'] == True:
                        argsdict['latex'] = False
                        print('LaTeX processing has been turned off.')
                    elif argsdict['latex'] == False:
                        argsdict['latex'] = True
                        print('LaTeX processing has been turned on.')
                elif option in ('-p', '--parallel'):
                    if argsdict['parallel'] == True:
                        argsdict['parallel'] = False
                        print('Parallel calculation of distances has been turned off.')
                    elif argsdict['parallel'] == False:
                        argsdict['parallel'] = True
                        print('Parallel calculation of distances has been turned on.')
                #elif option in ('-tm', '--timer'):
                #    if argsdict['timer'] == True:
                #        argsdict['timer'] = False
                #        print('Timer has been turned off.')
                #    elif argsdict['timer'] == False:
                #        argsdict['timer'] = True
                #        print('Timer has been turned on.')
                else :
                    pass

                if option.find('-lw') == 0:
                    if len(option) > 4:
                        while True:
                            try :
                                argsdict['latex_width'] = float(option[4:])
                                break
                            except ValueError :
                                while True :
                                    try :
                                        argsdict['latex_width'] = float(input('Input the width (in cm) for the plots generated usin LaTeX processing of text. '))
                                        break
                                    except ValueError:
                                        print('This is not a number.')
                                        continue
                                break
                        print('Plots generated using LaTeX as text processor will have  width of %s cm' % argsdict['latex_width'])
                    elif len(option) == 3 or len(option) == 4:
                        while True :
                            try :
                                argsdict['latex_width'] = float(input('Input the width (in cm) for the plots generated usin LaTeX processing of text. '))
                                break
                            except ValueError:
                                continue
                        print('Plots generated using LaTeX as text processor will have  width of %s cm' % argsdict['latex_width'])

                elif option.find('--latex-width') == 0:
                    if len(option) > 13:
                        while True:
                            try :
                                argsdict['latex_width'] = float(option[14:])
                                break
                            except ValueError :
                                while True :
                                    try :
                                        argsdict['latex_width'] = float(input('Input the width (in cm) for the plots generated usin LaTeX processing of text .'))
                                        break
                                    except ValueError:
                                        print('This is not a number.')
                                        continue
                                break
                        print('Plots generated using LaTeX as text processor will have  width of %s cm' % argsdict['latex_width'])
                    elif len(option) == 12 or len(option) == 13:
                        while True :
                            try :
                                argsdict['latex_width'] = float(input('Input the width (in cm) for the plots generated usin LaTeX processing of text. '))
                                break
                            except ValueError:
                                continue
                        print('Plots generated using LaTeX as text processor will have  width of %s cm' % argsdict['latex_width'])

                elif option.find('-s') == 0:
                    if len(option) > 3:
                        argsdict['subdir'] = option[3:]
                    elif len(option) == 2 or len(option) == 3:
                        argsdict['subdir'] = input('Input the name of the directory where plots will be saved. ')
                    if argsdict['subdir'] != '.':
                        print('Plots will be saved in the \'%s\' subdirectory' % argsdict['subdir'])
                    elif argsdict['subdir'] == '.':
                        print('Plots will be saved in the present directory\t %s' % argsdict['subdir'])

                elif option.find('--subdir') == 0:
                    if len(option) > 3:
                        argsdict['subdir'] = option[9:]
                    elif len(option) == 2 or len(option) == 3:
                        argsdict['subdir'] = input('Input the name of the directory where plots will be saved. ')
                    if argsdict['subdir'] != '.':
                        print('Plots will be saved in the \'%s\' subdirectory' % argsdict['subdir'])
                    elif argsdict['subdir'] == '.':
                        print('Plots will be saved in the present directory')

                options_menu(argsdict)


            elif option in ('1','2','3','4','5','6'):
                if argsdict['menu_type'] == 1:
                    argsdict['option'] = int(option)
                    return argsdict
                #    break
                elif argsdict['menu_type'] == 2:
                    if option in ('1', '2'):
                        argsdict['option'] = int(option)
                        return argsdict
                #        break
                    else :
                        print('This option is not available')
                        continue

            elif option in ('exit', 'Exit', 'EXIT', 'eXIT', '0'):
                argsdict['option'] = option
                return argsdict
                #break

            elif option == '':
                print('Select some module or option or type \'exit\' to leave the program.')
                continue
            else :
                print('This module or option is not available')
                continue
