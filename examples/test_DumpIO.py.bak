from LammpsIO import DumpIO
import numpy as np

# Reads test_in.lammpstrj
dump_in = DumpIO()                 # Creates DumpIO object
dump_in.OpenI('test_in.lammpstrj') # Opens test_in.lammpstrj as input

# Preliminary data
d = dump_in.NumDims()              # Gets number of dimensions (2 or 3)
n_f = dump_in.NumFrames()          # Gets number of frames
n_p = dump_in.NumPartsCurrentf()   # Gets number of particles in the 
                                   #  current frame (frame 0 by default)

# Reading frame = 0
dump_in.LoadNextFrame()            # Loads first frame
dump_in.SortByID()                 # Sorts particles by id
pos_0 = dump_in.GetPos()           # Gets particle positions at f = 0
bb_0 = dump_in.GetBB()             # Gets box bounds at f = 0
t_0 = dump_in.Gett()               # Gets timestep at f = 0
id_0 = dump_in.GetID()             # Gets particle ids at f = 0
data_0 = dump_in.GetDataCol('rand')# Gets data column rand at f = 0

# Skipping frames
for f in range(3):
  dump_in.LoadNextFrame()          # Skips 3 frames to get to f = 3

# Reading frame = 3
dump_in.SortByID()                 # Sorts particles by id
pos_3 = dump_in.GetPos()           # Gets particle positions at f = 3
bb_3 = dump_in.GetBB()             # Gets box bounds at f = 3
t_3 = dump_in.Gett()               # Gets timestep at f = 3
id_3 = dump_in.GetID()             # Gets particle ids at f = 3
data_3 = dump_in.GetDataCol('rand')# Gets data column rand at f = 3

# Closes input
dump_in.CloseI()

# Writes test_out.lammpstrj
dump_out = DumpIO()
dump_out.OpenO('test_out.lammpstrj') # Opens test_out.lammpstrj as
                                     #   output

# Writes frame = 0
dump_out.SetPos(pos_0)              # Sets particle positions at f = 0
dump_out.SetBB(bb_0)                # Sets box bounds at f = 0
dump_out.Sett(t_0)                  # Sets timestep at f = 0
dump_out.SetID(id_0)                # Sets particle ids at f = 0
dump_out.SetDataCol('rand', data_0) # Sets data column rand at f = 0
dump_out.WriteNextFrame()           # writes frame data to test_out

# Writes frame = 1
dump_out.SetPos(pos_3)              # Sets particle positions at f = 3
dump_out.SetBB(bb_3)                # Sets box bounds at f = 3
dump_out.Sett(t_3)                  # Sets timestep at f = 3
dump_out.SetID(id_3)                # Sets particle ids at f = 3
dump_out.SetDataCol('rand', data_3) # Sets data column rand at f = 3
dump_out.WriteNextFrame()           # writes frame data to test_out

# Closes output
dump_out.CloseO()
