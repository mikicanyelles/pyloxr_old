'''
Module for creacting a txt file where all the properties of the MD trajectory are specified
'''

from modules import loader

def summarize(u, argsdict):
    if argsdict['u_loaded'] == False:
        u, argsdict = loader.universe_loader_traj(argsdict)

    txt = open("MD_summary.txt", 'w+')
    txt.write("The system has %s atoms, %s residues and %s segments.\n" % (len(u.atoms), len(u.residues), len(u.segments)))
    txt.write("The dynamics starts at %s ns, ends at %s ns, lasts a total of %s ns, has %s frames and a timestep of %s.\n" % (round(u.trajectory[0].time/1000,1), u.trajectory[-1].time/1000, u.trajectory[-1].time/1000-round(u.trajectory[0].time/1000,1), len(u.trajectory), ((u.trajectory[-1].time/1000-round(u.trajectory[0].time/1000,1))/(u.trajectory[-1].frame-u.trajectory[0].frame +1))))
    txt.write("The system has the following shape:\n")
    txt.write("\tEdges length: (%s, %s, %s).\n" % (u.dimensions[0],u.dimensions[1],u.dimensions[2]))
    txt.write("\tAngles: (%s, %s, %s).\n" % (u.dimensions[3],u.dimensions[4],u.dimensions[5]))
    txt.close()

    print("Summary saved!"); return u, argsdict
