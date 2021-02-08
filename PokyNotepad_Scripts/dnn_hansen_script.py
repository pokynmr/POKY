
logo="""
  Script to reconstruct sparsely sampled 2D NMR spectra using
  Deep neural networks.

  Please cite:
  D. Flemming Hansen
  "Using Deep Neural Networks to Reconstruct Non-uniformly Sampled NMR Spectra"
  J. Biomol. NMR (2019)

  DNN NUS reconstruction code originally written by D. Flemming Hansen
  March 2019

  POKY implementation with GUI and NMRGlue processing
  by Woonghee Lee
  January 2021

  You will be asked three times (or more).
  First, where is (or will be) your DeepNUS directory?
  Second, where will be your processing directory?
  Third, what will be the name of sparse time-domain data?

  You will have "_td" and "_fd" reconstructed files.
"""

print(logo)
import os

NP=256      # Number of complex points in the indirect dimension
Verbose=2   # Verbose level
LoadData=True # Load data after completion

# POKY libraries
import __main__
s = __main__.main_session

# Path you have DeepNUS files
# Download and unpac from here: http://www.ucl.ac.uk/hansen-lab/DeepNUS/DeepNUS.tar.gz
# Make sure you have run MakeSpare.py before executing this script.

sampling_methods = ('256_pg_32_002', '256_rand_32_000', 'Cancel')
samp = s.show_selectionexdialog('Sampling', 'Select sampling method', sampling_methods)
if samp in [-1, 2]:
  raise SystemExit

deep_nus_path = s.open_directorydialog('Set DeepNUS home directory by Hansen D.F.', '')
if deep_nus_path == '':
  raise SystemExit

out_path = s.open_directorydialog('Set your output directory.', deep_nus_path)
if out_path == '':
  raise SystemExit

# Sparse TD pipe file name (prepared by using modified DFH's MakeSparse.py)
# You can either generate a new one or use the existing file.
sparse_specpath = s.save_filedialog('Select a sparse time-domain data to save.',
            'NMRPipe (*.ft1 *.ft2 *.ft3 *.dat);; Any (*)', out_path)
#sparse_specpath = s.open_filedialog('Select a sparse time-domain data to save.',
#            'NMRPipe (*.ft1 *.ft2 *.ft3 *.dat);; Any (*)', out_path)

# Reconstructed TD spectrum name (will be created)
td_outname = os.path.splitext(os.path.basename(sparse_specpath))[0] + '_td'
# Reconstructed FD spectrum name (will be created)
fd_outname = os.path.splitext(os.path.basename(sparse_specpath))[0] + '_fd'

# Sampling Schedule
SampleFile= os.path.join(deep_nus_path, 'SamplingSchedules',
                                         sampling_methods[samp] + '.sched')
# Optimised Parameter tensors
ParamFile= os.path.join(deep_nus_path, 'ParameterTensors',
                                          sampling_methods[samp] + '.h5')

# The network graph
ModelFile= os.path.join(deep_nus_path,
            'NetGraphs',
            'model_lstm_3_np256_ss32.json')

# download if DeepNUS is not in place.
if not os.path.exists(SampleFile) or not os.path.exists(ParamFile) \
    or not os.path.exists(ModelFile):
  print('DeepNUS package not found. Download from web.')
  if not s.show_message_yes_no('DeepNUS package not found.',
        'Do you want to download from D. Flemming Hansen\'s homepage? ' +
        'It will take a while and your POKY will not response. ' +
        'You can visit his homepage and download manually instead.'):
    raise SystemExit
  from wlutil import download_from_web, gunzip, move_path
  download_from_web('http://www.ucl.ac.uk/hansen-lab/DeepNUS/DeepNUS.tar.gz',
                    os.path.join(deep_nus_path, 'DeepNUS.tar.gz'))
  gunzip(os.path.join(deep_nus_path, 'DeepNUS.tar.gz'))
  move_path(os.path.join(deep_nus_path, 'DeepNUS'), deep_nus_path)

# make sparse td file if not exists.
# load modules here because we will use anyway
import sys
import numpy as np
import nmrglue as ng

if not os.path.exists(sparse_specpath):
  #fully sampled spectrum
  infile = os.path.join(deep_nus_path, 'Uniform', 'test.ft1')
  dic, data = ng.pipe.read(infile)

  #
  # let's read the sampling schedule
  ss=[]
  for l in open(SampleFile,'r').readlines():
    ss.append(int(l.strip() ))

  #
  # make the small file (only sampled points)
  presize=dic['FDF1TDSIZE']
  dic['FDF1TDSIZE']=len(ss)
  dic['FDF1APOD']=len(ss)
  dic['FDSLICECOUNT']=len(ss)
  dic['FDSPECNUM']=len(ss)

  ndata=np.zeros( (2*len(ss), data.shape[1]) , dtype=np.float32 )
  for n1 in range( len(ss) ):
    for n2 in range(ndata.shape[1]):
      ndata[n1*2  ,n2]=data[ss[n1]*2,  n2]
      ndata[n1*2+1,n2]=data[ss[n1]*2+1,n2]

  ng.pipe.write(sparse_specpath,dic,ndata,overwrite=True)

##### No need to change below here ######
# If you need to change processing parameters, go to the bottom.
#

# POKY libraries
import __main__
s = __main__.main_session

infile = sparse_specpath
inpath = os.path.dirname(infile)
td_outname = os.path.join(inpath, td_outname+'.ft1')
fd_outname = os.path.join(inpath, fd_outname+'.ft2')

if Verbose > 1:
  print('Input path: ' + infile)
  print('Reconstructed TD path: ' + td_outname)
  print('Reconstructed FD path: ' + fd_outname)

# Tensorflow
if Verbose>1:
  print(' -- Loading tensorflow libraries ... ')

import tensorflow as tf
from tensorflow import keras

if Verbose>1:
  print(' DONE ')

#
#Reduce output from tensorflow
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

SamplingSchedule=[]
for l in open(SampleFile,'r').readlines():
  SamplingSchedule.append( int(l.split()[0]))
SamplingSchedule=np.array(SamplingSchedule)
SampledPoints=len(SamplingSchedule)

model_json=open(ModelFile,'r').read()

if Verbose>1:
  print(' -- Loading parameter tensors  ... ')
#with tf.device('/device:GPU:1'):  # << Include this to use specific device
model=keras.models.model_from_json(model_json,
        custom_objects={'cos': tf.math.cos} )
model.load_weights(ParamFile)

if Verbose>1:
  print(' DONE ')

#
# Get spectrum
dic, data = ng.pipe.read(infile)

#
# Allocate space for input for neural network
inp   = np.zeros( (data.shape[1],2,SampledPoints))
x2_train= np.zeros( (data.shape[1],SampledPoints))

#
# Allocate space for output data
ndata = np.zeros( (NP*2, data.shape[1]), dtype=np.float32 )

if Verbose>1:
  print(' -- Start making input       ... ')

for n2 in range(data.shape[1]):
  #
  x2_train[n2,:]=np.copy(SamplingSchedule/NP)
  #
  norm=np.sqrt( data[0,n2]*data[0,n2] + data[1,n2]*data[1,n2] )
  #
  for n in range( int(data.shape[0]/2) ):
    inp[n2,0,n] = data[n*2,  n2]/norm
    inp[n2,1,n] = data[n*2+1,n2]/norm

if Verbose>1:
  print(' DONE ')
if Verbose>1:
  print(' -- Start predicting spectrum  ... ')

# Predict spectrum
#with tf.device('/device:GPU:1'): # << Include this to use specific device
out=model.predict([ inp,x2_train ] )

if Verbose>1:
  print(' DONE ')
if Verbose>1:
  print(' -- Start saving spectrum    ... ')

#
# Untangle output and put back
for n2 in range(data.shape[1]):
  norm=np.sqrt( data[0,n2]*data[0,n2] + data[1,n2]*data[1,n2] )

  for n in range(len(out[n2,0,:])):
    ndata[ n*2  , n2 ] = out[n2,0,n]*norm
    ndata[ n*2+1, n2 ] = out[n2,1,n]*norm
  #
  # Put back data data were already sampled
  for n in range( len(SamplingSchedule) ):
    ndata[ 2*SamplingSchedule[n], n2 ]   = data[2*n,n2]
    ndata[ 2*SamplingSchedule[n]+1, n2 ] = data[2*n+1,n2]

dic['FDF1TDSIZE']  =NP
dic['FDF1APOD']    =NP
dic['FDSLICECOUNT']=NP
dic['FDSPECNUM']   =NP

ng.pipe.write(td_outname, dic, ndata, overwrite=True)

if Verbose>1:
  print('Reconstruction finished.')

########################################################################
# Processing from reconstructed data
# Using NMRGlue
#
dic, data = ng.pipe.read(td_outname)
dic, data = ng.pipe_proc.tp(dic, data)
dic, data = ng.pipe_proc.sp(dic, data, off=0.42, end=0.98, pow=2, c=0.5)
dic, data = ng.pipe_proc.zf(dic, data, auto=True)
dic, data = ng.pipe_proc.zf(dic, data, zf=1)
dic, data = ng.pipe_proc.ft(dic, data, alt=True)
dic, data = ng.pipe_proc.ps(dic, data, p0=0.0, p1=0.0)
dic, data = ng.pipe_proc.di(dic, data)
dic, data = ng.pipe_proc.tp(dic, data)
# polynomial baseline correction has not been implemented in NMRGlue yet
#dic, data = ng.pipe_proc.poly(dic, data)
# median baseline correction is used instead
dic,data = ng.pipe_proc.med(dic, data)
ng.pipe.write(fd_outname,dic,data,overwrite=True)

# open files. this will generate ucsf files in the same directory as well.
if LoadData:
  s.open_spectrum(infile)
  s.open_spectrum(td_outname)
  s.open_spectrum(fd_outname)

if Verbose>1:
  print('Fourier Transform finished.')

