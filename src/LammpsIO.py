# =*= coding: utf-8 -*-
"""LammpsIO
Created on 21 September 2018
@author: Robert Ivancic

This module was created to facilitate input and output of LAMMPS
trajectory files in dump and netCDF (AMBER) format.

Notes
-----
  We plan to add a log reading functionality to this module soon.

""" 

import numpy as np
import sys
from netCDF4 import Dataset
from netCDF4 import stringtochar

class NetCDFIO:
  """Input and output of LAMMPS trajectories in NetCDF (AMBER).

  This class inputs and outputs individual frames from an AMBER
  format molecular dynamics trajectory. More on AMBER can be found
  here: http://ambermd.org/netcdf/nctraj.xhtml. More on NetCDF can
  be found here: https://en.wikipedia.org/wiki/NetCDF. There are two
  major advantages to this output format over standard dump files.
  First, these files are in binary format and thus, use about 1/3 of
  the space of standard dump files. Second, reading standard dump files
  must be done sequentially. With NetCDF you may obtain data from any
  frame of interest *f* instantly by simply selecting that frame in the
  data array. The major disadvantage of this format is that it assumes
  a constant number of particles throughout a simulation or experiment.
  Particles are assumed to be ordered by their ids.

  """

  def __init__(self):
    self.__openO = False
    self.__sett = False
    self.__setBB = False
    self.__setPos = False
    self.__setDataCol = False
    self.__closeO = False
    self.__openI = False
    self.__numDims = False
    self.__numFrames = False
    self.__numParts = False
    self.__gett = False
    self.__getBB = False
    self.__getPos = False
    self.__getDataCol = False
    self.__closeI = False
    self.__DataCols = []
    self.__DataCol_vars = []

  def OpenO(self, nc_file_name):
    """Opens NetCDF file for output.
    
    Opens a NetCDF file of nc_file_name for writing (output).

    Parameters
    ----------
    nc_file_name : str
      Name of netCDF file you wish to write to.

    """

    if self.__openI:
      print('ERROR: Cannot set input and output in same instance')
      sys.exit(-1)
    if self.__openO:
      print('ERROR: Cannot open file twice!')
      sys.exit(-1)

    # Conventions
    self.__nc_file = Dataset(nc_file_name, 'w', format='NETCDF4')
    self.__nc_file.Conventions = 'AMBER'
    self.__nc_file.ConventionVersion = '1.0'
    self.__nc_file.program ='LammpsPy.py'
    self.__nc_file.programVersion = '1.0'

    # Dimensions except number of particles
    self.__nc_file.createDimension('frame', 0)
    self.__nc_file.createDimension('spatial', 3)
    self.__nc_file.createDimension('cell_spatial', 3)
    self.__nc_file.createDimension('cell_angular', 3)
    self.__nc_file.createDimension('label', 5)

    # Preliminary variables
    char_sp_var = self.__nc_file.createVariable('char spatial',
      'S1', ('spatial'))
    char_sp_var[:] = ['x','y','z']
    char_csp_var = self.__nc_file.createVariable('char cell_spatial',
      'S1', ('cell_spatial'))
    char_csp_var[:] = ['a','b','c']
    char_sp_var = self.__nc_file.createVariable('char cell_angular',
      'S1',('cell_angular','label'))
    char_sp_var[:] = stringtochar(np.array(['alpha','beta','gamma']))
  
    # notes OpenO has been run
    self.__openO = True

  def Sett(self, f, t):
    """Sets timestep t of frame f
    
    Parameters
    ----------
    f : int
      frame number
    t : int
      timestep at frame f

    """
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    if not self.__sett:
      self.__t_var = self.__nc_file.createVariable('time', 'f4',
        ('frame'))
      self.__t_var.units = 'LJ'
    self.__t_var[f] = t
    self.__sett = True

  def SetBB(self, f, bb_t):
    """Sets box boundaries of frame f. These are formatted as
    [low_x, high_x] for each dimension.
    
    Parameters
    ----------
    f : int
      frame number
    bb_t : numpy array (dimension,2)
      box boundaries at time t

    """

    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    if not self.__setBB:
      self.__cell_o_var = self.__nc_file.createVariable('cell_origin', 
        'f4', ('frame','cell_spatial'))
      self.__cell_o_var.units = 'LJ'
      self.__cell_l_var = self.__nc_file.createVariable('cell_lengths', 
        'f4', ('frame','cell_spatial'))
      self.__cell_l_var.units = 'LJ'  
      self.__cell_a_var = self.__nc_file.createVariable('cell_angles', 
        'f4', ('frame','cell_angular'))
      self.__cell_a_var.units = 'degree'

    # Sets dimension of output
    if not self.__setBB and not self.__setPos:
      if bb_t.shape[0] != 2 and bb_t.shape[0] != 3:
        print('ERROR: bb_t.shape[0] != 2 or 3')
        sys.exit(-1)
      self.__d = bb_t.shape[0]
      self.__nc_file.dimension = self.__d

    # Adds a 3rd dimension for 2D data
    if self.__d != bb_t.shape[0]:
      print('ERROR: dimension != bb_t.shape[0]')
      sys.exit(-1)
    if self.__d == 2:
      bb_t = np.vstack((bb_t, np.array([[-10**(-6),10**(-6)]])))

    self.__cell_o_var[f] = bb_t[:,0]
    self.__cell_l_var[f] = bb_t[:,1]-bb_t[:,0]
    self.__cell_a_var[f] = np.ones(3)*90
    self.__setBB = True

  def SetPos(self, f, pos_t):
    """Sets particle positions (coordinates) at frame f
    
    Parameters
    ----------
    f : int
      frame number
    pos_t : numpy array (number of particles, dimension)
      particle positions (coordinates) at frame f
      
    """ 
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    if not self.__setPos and not self.__setDataCol:
        self.__nc_file.createDimension('atom', pos_t.shape[0])
    if not self.__setPos:
      self.__c_var = self.__nc_file.createVariable('coordinates', 'f4',
        ('frame','atom','spatial'))
      self.__c_var.units = 'LJ'

    # Sets dimension of output
    if not self.__setBB and not self.__setPos:
      if pos_t.shape[1] != 2 and pos_t.shape[1] != 3:
        print('ERROR: pos_t.shape[1] != 2 or 3')
        sys.exit(-1)
      self.__d = pos_t.shape[1]
      self.__nc_file.dimension = self.__d

    # Adds a 3rd dimension for 2D data
    if self.__d != pos_t.shape[1]:
      print('ERROR: dimension != pos_t.shape[1]')
      sys.exit(-1)
    if self.__d == 2:
      pos_t = np.hstack((pos_t, np.array([np.zeros(pos_t.shape[0])]).T))

    self.__c_var[f] = pos_t
    self.__setPos = True

  def SetDataCol(self, f, label, data_col_t):
    """Sets a data column at frame f

    A data column is a way of storing values that are associated with 
    every particle in a simulation or experiment at every timestep.
    Examples include: J2, D2min, p_{hop}, or softness of every particle
    at every timestep. This sets a data column for time t.

    Parameters
    ----------
    f : int
      frame number
    label : str
      data column name
    data_col_t : numpy array (number of particles,)
      data column values for all particles at frame f
      
    """
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    if not self.__setPos and not self.__setDataCol:
        self.__nc_file.createDimension('atom', pos_t.shape[0])
    if label not in self.__DataCols:
      self.__DataCols.append(label)
      self.__DataCol_vars.append(self.__nc_file.createVariable(
        label, 'f4', ('frame','atom')))
      self.__DataCol_vars[-1].units = 'LJ'
    idx_col = (self.__DataCols).index(label)
    (self.__DataCol_vars)[idx_col][f] = data_col_t
    self.__setDataCol = True

  # Closes output file
  def CloseO(self):
    """Closes output NetCDF file
    """
    self.__nc_file.close()
    self.__closeO = True

  # Opens input file
  def OpenI(self, nc_file_name):
    """Opens NetCDF file for input.
    
    Opens a NetCDF file of nc_file_name for reading (input).

    Parameters
    ----------
    nc_file_name : str
      Name of netCDF file you wish to read from.

    """
    if self.__openI:
      print('ERROR: Cannot open file twice!')
      sys.exit(-1)
    if self.__openO:
      print('ERROR: Cannot set input and output in same instance')
      sys.exit(-1)
    self.__nc_file = Dataset(nc_file_name, 'r', format='NETCDF4')

    self.__openI = True

  # returns number of dimensions
  def NumDims(self):
    """Outputs number of dimensions of current simulation being read.

    Returns
    -------
    d : int
      Number of dimensions of simulation you are reading.

    """

    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)
    self.__numDims = True

    # Checks if has attribute dimensions. If it does takes that
    # number to be ground truth. Otherwise assumes 3D.
    if hasattr(self.__nc_file, 'dimension'):
      self.__d = self.__nc_file.dimension
    else:
      self.__d = 3

    return self.__d

  # returns number of frames
  def NumFrames(self):
    """Outputs number of frames of current simulation being read.

    Returns
    -------
    n_f : int
      Number of frames of simulation you are reading.

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)
    self.__numFrames = True
    n_f = self.__nc_file.dimensions['frame'].size
    return n_f

  # returns number of particles
  def NumParts(self):
    """Outputs number of particles of current simulation being read.

    Returns
    -------
    n_p : int
      Number of particles of simulation you are reading.

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)
    self.__numParts = True
    n_p = self.__nc_file.dimensions['atom'].size
    return n_p

  # returns timestep for frame f
  def Gett(self, f):
    """Gets timestep t of frame f
    
    Parameters
    ----------
    f : int
      frame number

    Returns
    -------
    t : int

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__gett = True
    t = self.__nc_file.variables['time'][f]
    return t

  # returns box boundaries for frame f
  def GetBB(self, f):
    """Gets box boundaries (BB) of frame f. These are formatted as
    [low_x, high_x] for each dimension.
    
    Parameters
    ----------
    f : int
      frame number

    Returns
    -------
    bb_t : numpy array (dimension,2)
      box boundaries at time t
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getBB = True
    bb_L = self.__nc_file.variables['cell_origin'][f]
    bb_R = self.__nc_file.variables['cell_lengths'][f] + bb_L

    bb_t = np.vstack((bb_L, bb_R)).T

    # Gets dimensions of array
    d = self.NumDims()

    return bb_t[:d]

  # returns atom coordinates for frame f
  def GetPos(self, f):
    """Gets particle positions (coordinates) at frame f
    
    Parameters
    ----------
    f : int
      frame number

    Returns
    -------
    pos_t : numpy array (number of particles, dimension)
      particle positions (coordinates) at frame f
      
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getPos = True
    pos_t = self.__nc_file.variables['coordinates'][f]

    # Gets dimensions of array
    d = self.NumDims()

    return pos_t[:,:int(d)]

  # returns box boundaries for frame f
  def GetDataCol(self, f, label):
    """Sets a data column at frame f

    A data column is a way of storing values that are associated with 
    every particle in a simulation or experiment at every timestep.
    Examples include: J2, D2min, p_{hop}, or softness of every particle
    at every timestep. This gets a data column for time t.

    Parameters
    ----------
    f : int
      frame number
    label : str
      data column name

    Returns
    -------
    data_col_t : numpy array (number of particles,)
      data column values for all particles at frame f
      
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getDataCol = True
    data_col = self.__nc_file.variables[label][f]
    return data_col

  # Closes output file
  def CloseI(self):
    """Closes input NetCDF file
    """
    self.__nc_file.close()
    self.__closeI = True



class DumpIO:
  """ Input and output LAMMPS trajectories in dump (lammpstrj) format.

  This class inputs and outputs frames of a LAMMPS trajectory. It does
  this *sequentially*, i.e. it reads (writes) in a single frame of a dump
  (lammpstrj) trajectory at a time starting at frame 0. This approach has
  two benefits to reading in an entire trajectory at once both steming 
  from the fact that it is more memory efficient. First (and most 
  importantly), python has memory limits that are easily reached if 
  the number of particles times the number of frames loaded in at one
  time is greater than 25,000,000. Practically, it is very easy to hit
  this limit when running large or long time simulations. Second,
  this approach is generally faster than loading in all frames at once
  as it requires much less time to allocate memory.
  """

  def __init__(self):
    self.__openO = False
    self.__sett = False
    self.__setBB = False
    self.__setPos = False
    self.__setID = False
    self.__setDataCol = False
    self.__id_t = np.array([])
    self.__writeNextFrame = False
    self.__setIDWarning = False
    self.__closeO = False
    self.__openI = False
    self.__loadNextFrame = False
    self.__numDims = False
    self.__numFrames = False
    self.__numParts = False
    self.__gett = False
    self.__getBB = False
    self.__getPos = False
    self.__getDataCol = False
    self.__closeI = False
    self.__DataCols = []
    self.__DataCol_arrs = []

  def OpenO(self, dump_file_name):
    """ Opens dump file for output.

    Opens a dump (lammpstrj) file of name dumpe_file_name for writing
    (output). 
    """

    if self.__openI:
      print('ERROR: Cannot set input and output in same instance')
      sys.exit(-1)
    if self.__openO:
      print('ERROR: Cannot open file twice!')
      sys.exit(-1)

    # Opens file for writing
    self.__dump_file = open(dump_file_name, 'w')

    # notes OpenO has been run
    self.__openO = True

  def Sett(self, t):
    """ Sets timestep t of current frame

    Parameters
    ----------
    t : int
      timestep of current frame
    """
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    self.__t = t
    self.__sett = True

  def SetBB(self, bb_t):
    """ Sets box boundaries (BB) of current frame

    Parameters
    ----------
    bb_t : numpy array (dimension,2)
      box boundaries of current frame
    """

    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)

    # Sets dimension of output
    if not self.__setBB and not self.__setPos:
      if bb_t.shape[0] != 2 and bb_t.shape[0] != 3:
        print('ERROR: bb_t.shape[0] != 2 or 3')
        sys.exit(-1)
      self.__d = bb_t.shape[0]

    # Adds a 3rd dimension for 2D data
    if self.__d != bb_t.shape[0]:
      print('ERROR: dimension != bb_t.shape[0]')
      sys.exit(-1)
    if self.__d == 2:
      bb_t = np.vstack((bb_t, np.array([[-10**(-6),10**(-6)]])))

    self.__bb_t = bb_t
    self.__setBB = True

  def SetID(self, id_t):
    """ Sets particle ids of current frame

    Parameters
    ----------
    id_t : numpy array (number of particles,)
      particle ids of each particle of current frame

    """

    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    self.__id_t = id_t
    self.__setID = True

  def SetPos(self, pos_t):
    """ Sets particle positions of current frame

    Parameters
    ----------
    pos_t : numpy array (number of particles, dimensions)
      particle positions (coordinates) of current frame in order of 
      id_t

    """
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)

    # Sets dimension of output
    if not self.__setBB and not self.__setPos:
      if pos_t.shape[1] != 2 and pos_t.shape[1] != 3:
        print('ERROR: pos_t.shape[1] != 2 or 3')
        sys.exit(-1)
      self.__d = pos_t.shape[1]

    # Adds a 3rd dimension for 2D data
    if self.__d != pos_t.shape[1]:
      print('ERROR: dimension != pos_t.shape[1]')
      sys.exit(-1)

    self.__pos_t = pos_t
    self.__setPos = True

  def SetDataCol(self, label, data_col_t):
    """ Sets particle positions of current frame
    
    Parameters
    ----------
    label : str
      name of the data column
    data_col_t : numpy array (number of particles,)
      particle positions (coordinates) of current frame in order of 
      id_t

    """
    if not self.__openO:
      print('ERROR: Must use OpenO before set functions')
      sys.exit(-1)
    if label not in self.__DataCols:
      self.__DataCols.append(label)

    # Finds index of data column
    idx_col = (self.__DataCols).index(label)
    while(idx_col >= len(self.__DataCol_arrs)):
      self.__DataCol_arrs.append([])
    self.__DataCol_arrs[idx_col] = data_col_t
    self.__setDataCol = True

  def WriteNextFrame(self):
    """ Writes the next frame to the dump file.
    
    Writes all timestep, box boundary, particle id, particle position,
    and data column information given in frame to file. Must use OpenO
    and Set functions before use.

    Warnings
    --------
    Will run without having set particle ids. Assumes particles in order
    from 1,...,(number of particles).
    """
    error_message = 'ERROR: Incorrect input formatting. Ensure:\n\t'
    if not self.__openO:
      print('ERROR: Must use OpenO before Write functions')
      sys.exit(-1)
    if not self.__sett:
      print(error_message + 'Sett is used')
      sys.exit(-1)
    if not self.__setBB:
      print(error_message + 'SetBB is used')
      sys.exit(-1)
    if not self.__setPos:
      print(error_message + 'SetPos is used')
      sys.exit(-1)
    if len(self.__DataCols) != len(self.__DataCol_arrs):
      print(error_message+'All data columns are stored in each frame')
      sys.exit(-1)
    if self.__id_t.shape[0] != self.__pos_t.shape[0] and self.__setID:
      print(error_message+'id_t.shape[0] == pos_t.shape[0]')
      sys.exit(-1)
    elif not self.__setID:
      # Automatically assumes particles in order if ids not set
      if not self.__setIDWarning:
        print('WARNING: SetID function not used. Assuming particle ID')
        print('\tis in order particles are in pos_t')
        self.__setIDWarning = True
      self.__id_t = np.arange(self.__pos_t.shape[0])+1
    for col in range(len(self.__DataCols)):
      if self.__DataCol_arrs[col].shape[0] != self.__pos_t.shape[0]:
        print(error_message+'For data column '+\
              str(self.__DataCols[col])+', '+\
              'col.shape[0] == pos_t.shape[0]')
        sys.exit(-1)

    # Finds number of particles and dimension
    n_p = self.__pos_t.shape[0]
    bb_t = self.__bb_t
    if self.__d == 2:
      dim_str = 'x y '
    elif self.__d==3:
      dim_str = 'x y z '
    else:
      print('ERROR: Must have 2 or 3 dimensions')
      sys.exit(-1)

    # Writes frame header
    self.__dump_file.write('ITEM: TIMESTEP\n'+str(self.__t)+'\n')
    self.__dump_file.write('ITEM: NUMBER OF ATOMS\n'+str(n_p)+'\n')
    self.__dump_file.write('ITEM: BOX BOUNDS\n')
    bb_str_t = ' \n'.join([' '.join(lin) for lin \
          in bb_t.astype(str)])+'\n'
    self.__dump_file.write(bb_str_t)
    dataH_str = 'ITEM: ATOMS '+'id '+dim_str+\
          ' '.join(self.__DataCols)+' \n'
    self.__dump_file.write(dataH_str)

    # Writes frame data
    ids = self.__id_t.astype(str)
    pos_t = self.__pos_t.astype(str)
    if len(self.__DataCol_arrs) > 0:
      dataCols_t = np.array(self.__DataCol_arrs).astype(str)
      data_t = np.vstack((ids, pos_t.T, dataCols_t)).T
    else:
      data_t = np.vstack((ids, pos_t.T)).T
    data_str_t = ' \n'.join([' '.join(lin) for lin in data_t])+'\n'
    self.__dump_file.write(data_str_t)
    

  def CloseO(self):
    """ Closes output file for dump
    """
    self.__dump_file.close()
    self.__closeO = True

  def OpenI(self, dump_file_name):
    """ Opens lammpstrj file for input

    Opens a dump (lammpstrj) file of name dump_file_name for reading
    (input)
    """
    if self.__openI:
      print('ERROR: Cannot open file twice!')
      sys.exit(-1)
    if self.__openO:
      print('ERROR: Cannot set input and output in same instance')
      sys.exit(-1)

    # Opens file for writing
    self.__dump_file = open(dump_file_name, 'r')

    self.__openI = True

  def LoadNextFrame(self):
    """ Loads next frame in the dump (lammpstrj) file.

    Reads all timestep, box boundary, particle id, particle position, 
    and data column information given in frame to file. Must use OpenI
    before use. Get functions may be used after this function.
    """
    # Reads the header of dump file
    line = ''
    while line[:11] != 'ITEM: ATOMS':
      line = self.__dump_file.readline()
      if line[:16] == 'ITEM: BOX BOUNDS':
        bb_t = []
        line = self.__dump_file.readline()
        while line[:5] != 'ITEM:':
          last_pos = self.__dump_file.tell()
          bb_t.append(np.array(line.split()).astype(float))
          line = self.__dump_file.readline()
        self.__dump_file.seek(last_pos)
        self.__bb_t = np.array(bb_t)
        line = ''
      if line[:14] == 'ITEM: TIMESTEP':
        line = self.__dump_file.readline()
        self.__t = int(line)
      if line[:21] == 'ITEM: NUMBER OF ATOMS':
        line = self.__dump_file.readline()
        n_p = int(line)

    # Obtains data columns of data
    DataCols = np.array(line[11:].split())
    col_id = np.where(DataCols=='id')[0][0]
    col_x = np.where(DataCols=='x')[0][0]
    col_y = np.where(DataCols=='y')[0][0]
    if 'z' in DataCols:
      self.__d = 3
      col_z = np.where(DataCols=='z')[0][0]
      col_coordinates = np.array([col_x, col_y, col_z, col_id])
    else:
      self.__d = 2
      col_coordinates = np.array([col_x, col_y, col_id])
    remove_cols = np.array(['id','x','y','z'])
    DataCols_no_coordinates = np.setdiff1d(DataCols,remove_cols)
    self.__DataCols = DataCols_no_coordinates
    DataCols_idx = np.array([np.where(DataCols==col)[0][0] \
          for col in self.__DataCols])
    n_DataCols = len(self.__DataCols)

    # Obtains data from frame
    self.__pos_t = np.zeros((n_p, self.__d))
    self.__id_t = np.zeros(n_p)
    self.__DataCol_arrs = [np.zeros(n_p) for ii in range(n_DataCols)]
    for pp in range(n_p):
      line = self.__dump_file.readline()
      data_pp = line.split()
      
      # Stores data
      self.__pos_t[pp,0] = float(data_pp[col_x])
      self.__pos_t[pp,1] = float(data_pp[col_y])
      if self.__d != 2:
        self.__pos_t[pp,2] = float(data_pp[col_z])
      self.__id_t[pp] = int(data_pp[col_id])
      for idx, col in enumerate(DataCols_idx):
        self.__DataCol_arrs[idx][pp] = float(data_pp[col])

    self.__loadNextFrame = True

  def NumDims(self):
    """Outputs number of dimensions of current simulation being read.

    Returns
    -------
    d : int
      Number of dimensions of simulation you are reading.
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)
    if not self.__loadNextFrame:
      last_pos = self.__dump_file.tell()
      self.LoadNextFrame()
      self.__dump_file.seek(last_pos)

    self.__numDims = True

    return self.__d

  def NumFrames(self):
    """Outputs number of frames of current simulation being read.

    Returns
    -------
    n_f : int
      Number of frames of the simulation you are reading.
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)

    # stores current position in file and restarts read
    # loops through all lines in file to find ITEM: ATOMS
    n_f = 0
    last_pos = self.__dump_file.tell()
    self.__dump_file.seek(0)
    for line in self.__dump_file:
      if line[:11] == 'ITEM: ATOMS':
        n_f += 1
    self.__dump_file.seek(last_pos)
    
    self.__numFrames = True

    return n_f

  def NumPartsCurrentf(self):
    """Outputs number of particles in current frame.

    Returns
    -------
    n_p : int
      Number of particles in the current frame.

    Warnings
    --------
    If used without loading any frames, gives number of particles of
    first frame.

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before num functions')
      sys.exit(-1)
    if not self.__loadNextFrame:
      last_pos = self.__dump_file.tell()
      self.LoadNextFrame()
      self.__dump_file.seek(last_pos)

    self.__numParts = True

    return self.__pos_t.shape[0]

  def Gett(self):
    """ Gets timestep t of current frame.

    Returns
    -------
    t : int
      timestep of the current frame.
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__gett = True
    t = self.__t
    return t

  def GetBB(self):
    """ Gets box boundary (BB) of current frame. These are formatted as
    [low_x, high_x] for each dimension.

    Returns
    -------
    bb_t : numpy array (dimension,2)
      box boundaries of current frame
    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getBB = True
    bb_t = self.__bb_t
    return bb_t[:self.__d]

  def GetID(self):
    """ Gets particle ids of current frame

    Returns
    ----------
    id_t : numpy array (number of particles,)
      particle ids of each particle of current frame

    """

    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getBB = True
    id_t = self.__id_t.astype(int)
    return id_t

  def GetPos(self):
    """ Gets particle positions (coordinates) of current frame

    Returns
    ----------
    pos_t : numpy array (number of particles, dimension)
      particle positions (coordinates) of each particle of current 
      frame

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getPos = True
    pos_t = self.__pos_t
    return pos_t

  # returns box boundaries for current frame
  def GetDataCol(self, label):
    """ Gets data column of current frame

    A data column is a way of storing values that are associated with 
    every particle in a simulation or experiment at every timestep.
    Examples include: J2, D2min, p_{hop}, or softness of every 
    particle at every timestep. This gets a data column for the current
    frame in the dump file.

    Parameters
    ----------
    label : str
      name of data column label

    Returns
    ----------
    id_t : numpy array (number of particles,)
      particle ids of each particle of current frame

    """
    if not self.__openI:
      print('ERROR: Must use OpenI before get functions')
      sys.exit(-1)
    self.__getDataCol = True
    try:
      data_col = np.where(self.__DataCols==label)[0][0]
      data = self.__DataCol_arrs[data_col]
    except:
      print('ERROR: label for data column not in dump file')
      sys.exit(-1)
    return data

  # Closes output file
  def CloseI(self):
    """ Closes input (reading) dump file.
    """
    self.__dump_file.close()
    self.__closeO = True

  def SortByID(self):
    """ Sorts particles of all arrays by their id

    Sorts pos_t, id_t, and data_col_t arrays so that 
    id_t = range(n_p)
    """

    # ranking of id values
    ranks = self.__id_t.argsort()

    # Performs sorting
    self.__id_t = self.__id_t[ranks]
    self.__pos_t = self.__pos_t[ranks]
    for col in range(len(self.__DataCol_arrs)):
      self.__DataCol_arrs[col] = self.__DataCol_arrs[col][ranks]

