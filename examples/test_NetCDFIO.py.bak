from LammpsIO import NetCDFIO
import numpy as np

# Reads test_in.lammpstrj
nc_in = NetCDFIO()                 # Creates NetCDFIO object
nc_in.OpenI('test_in.nc')   # Opens test_in.nc as input

# Preliminary data
d = nc_in.NumDims()                # Gets number of dimensions (2 or 3)
n_f = nc_in.NumFrames()            # Gets number of frames
n_p = nc_in.NumParts()             # Gets number of particles in the 
                                   #  current frame (frame 0 by default)

# Reading frame = 0
f = 0
pos_0 = nc_in.GetPos(f)              # Gets particle positions at f = 0
bb_0 = nc_in.GetBB(f)                # Gets box bounds at f = 0
t_0 = nc_in.Gett(f)                  # Gets timestep at f = 0
data_0 = nc_in.GetDataCol(f, 'rand') # Gets data column rand at f = 0

# Reading frame = 3
f = 3
pos_3 = nc_in.GetPos(f)              # Gets particle positions at f = 3
bb_3 = nc_in.GetBB(f)                # Gets box bounds at f = 3
t_3 = nc_in.Gett(f)                  # Gets timestep at f = 3
data_3 = nc_in.GetDataCol(f, 'rand') # Gets data column rand at f = 3

# Closes input
nc_in.CloseI()

# Writes test_out.lammpstrj
nc_out = NetCDFIO()
nc_out.OpenO('test_out.nc')          # Opens test_out.lammpstrj as
                                     #   output

# Writes frame = 0
f = 0
nc_out.SetPos(f, pos_0)              # Sets particle positions at f = 0
nc_out.SetBB(f, bb_0)                # Sets box bounds at f = 0
nc_out.Sett(f, t_0)                  # Sets timestep at f = 0
nc_out.SetDataCol(f, 'rand', data_0) # Sets data column rand at f = 0

# Writes frame = 1
f = 1
nc_out.SetPos(f, pos_3)              # Sets particle positions at f = 3
nc_out.SetBB(f, bb_3)                # Sets box bounds at f = 3
nc_out.Sett(f, t_3)                  # Sets timestep at f = 3
nc_out.SetDataCol(f, 'rand', data_3) # Sets data column rand at f = 3

# Closes output
nc_out.CloseO()
