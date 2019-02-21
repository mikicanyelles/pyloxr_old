def summarize():
    global u
    global traj
    ### Time counter starts
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
        
    txt = open("MD_summary.txt", 'w+')
    txt.write("The system has %s atoms, %s residues and %s segments.\n" % (len(u.atoms), len(u.residues), len(u.segments)))
    txt.write("The dynamics starts at %s ns, ends at %s ns, lasts a total of %s ns, has %s frames and a timestep of %s.\n" % (round(traj[0].time/1000,1), traj[-1].time/1000, traj[-1].time/1000-round(traj[0].time/1000,1), len(traj), ((traj[-1].time/1000-round(traj[0].time/1000,1))/(traj[-1].frame-traj[0].frame +1))))
    txt.write("The system has the following shape:\n")
    txt.write("\tEdges length: (%s, %s, %s).\n" % (u.dimensions[0],u.dimensions[1],u.dimensions[2]))
    txt.write("\tAngles: (%s, %s, %s).\n" % (u.dimensions[3],u.dimensions[4],u.dimensions[5]))
    txt.close()
    
    print("Summary saved!")