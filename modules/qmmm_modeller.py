'''
Module for creating a proper model for QM/MM calculations by cropping the water box. 
It saves the 'pdb' of the croped model, as well as the topology and parameters and initial coordinates.
'''
#argsdict=dict({'trajectory': ['3rde.prmtop', None], 'frame': None, 'latex': False, 'latex_width': None, 'parallel': False, 'subdir': 'plots', 'timer': False, 'menu_type' : None, 'u_loaded' : False})

#argsdict=dict({'trajectory': ['3rde.prmtop', None], 
#'frame': None, 'latex': False, 'latex_width': None, 
#'parallel': False, 'subdir': 'plots', 'timer': False, 
#'menu_type' : None, 'u_loaded' : False})

def qmmm_modeller(u, argsdict):

    if argsdict['frame'] != None:
        