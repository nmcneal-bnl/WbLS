#!/usr/bin/env python2
'''
These classes handle data from the PTI spectrometer.
TEXT data are assumed (not the .gx* nonsense).
'''
import copy
import os
from enum import Enum
import time
import pandas
import numpy

class PTIData(object):
    '''PTI spectrometer data class.'''
    run_types = Enum('RunType', 'Unknown Emission Excitation Synchronous')
    file_types = Enum('FileType', 'Unknown Session Trace Group')
    print_initialize= False
    def __init__(self, fname):
        if self.print_initialize:
            print("Initializing PTI_Data at {0}".format(time.asctime(time.localtime())))
        
        ## Member variables for reference ##
        self.file_path = str() # The file name to be read in
        self.file_type = str() # Possibilities: Session, Trace, Group

        self.acq_start = None
        self.num_samples = -1
        self.step_size = -1
        self.PMT_mode = str()
        self.ex_range = list([-2,-1])
        self.em_range = list([-2,-1])

        self.wavelengths = None
        self.raw_data = None
        self.cor_data = None
        self.diode   = None
        
        self.baseline_incpt = None
        self.baseline_slope = None
        self.baseline_incpt_se = None
        self.baseline_slope_se = None

        self.ex_monochromator_offset = -1
        self.em_monochromator_offset = -1


        ## Reading in the file ##
        # Take in the given parameter
        self.file_path = fname
        
        # Checking for file's existence and opening if possible
        if not os.path.exists(fname):
            print("ERROR!! File does not exist.")
            self.read_success = False            
            return
        

        with open(self.file_path, 'r') as target_file:
            firstline = target_file.readline()
            if '<Session>' in firstline:
                self.file_type = self.file_types.Session
            elif '<Trace>' in firstline:
                self.file_type = self.file_types.Trace
            elif '<Group>' in firstline:
                self.file_type = self.file_types.Group
            else:
                print("ERROR!! Unknown file format.")
                self.file_type = self.file_types.Unknown
                self.read_success = False
                return

        self.read_success = self.ReadHeaderInfo()
        self.WL = [0]*self.num_samples
        if self.file_type == self.file_types.Session:
            self.Spec = [0]*self.num_samples
            self.SpecRaw = [0]*self.num_samples
            self.USpecRaw = [0]*self.num_samples
            self.ExCorr = [0]*self.num_samples #Note ExCorr here is the photodiode signal.
            self.FileSpecCorrected = [0]*self.num_samples
            self.UFileSpecCorrected = [0]*self.num_samples
        elif self.file_type.value > 1:
            self.Trace = [0]*self.num_samples
            self.UTrace = [0]*self.num_samples
        
        self.ReadSpecData()
        self.SpecCorrected = None
        self.USpecCorrected = None
        return

    def RegisterCorrSpec(self, CorrSpec, UCorrSpec):
        '''
        Define the SpecCorrected and USpecCorrected members.
        '''
        self.SpecCorrected = CorrSpec
        self.USpecCorrected = UCorrSpec
        return

    def ReadHeaderInfo(self):
        '''
        Read the header (first 7 lines) and extract useful info.
        This is called in initialization.
        
        Useful info that is set:
        - Start date and time (as a time struct).
        - Number of data points acquired.
        - Excitation Wavelength range.
        - Emission Wavelength Range.
        - Run type (from run_types enum).
        '''
        #Read the header info to determine the run type
        with open(self.file_path, 'r') as thefile:
            if self.file_type == self.file_types.Session:
                success = self._ReadHdrSession(thefile)
            elif self.file_type == self.file_types.Trace:
                success = self._ReadHdrTrace(thefile)
            elif self.file_type == self.file_types.Group:
                success = self._ReadHdrGroup(thefile)
        return success
        
    def _ReadHdrSession(self, thefile):
        for i, line in enumerate(thefile):

            # Reading in the second line to get the start time
            if i==1:
                wrds = line.split()
                self.acq_start = time.strptime(wrds[-2] + ' ' + wrds[-1], 
                                              '%Y-%m-%d %H:%M:%S')
            
            # Reading in the number of samples from line 6
            elif i==5:
                words = line.split()
                self.num_samples = int(words[0])

            elif i==6:
                success = self._ReadWLRangeLine(line)
                break
        return success

    def _ReadHdrTrace(self, thefile):
        for i, line in enumerate(thefile):
            if i==1:
                #No acquisition time in trace files, just use file creation time as an estimate.
                self.acq_start = time.localtime(os.path.getctime(self.file_path))
                self.num_samples = int(line)
            elif i==2:
                success = self._ReadWLRangeLine(line)
                break
        return success

    def _ReadHdrGroup(self, thefile):
        #No acquisition time in trace files, just use file creation time as an estimate.
        success = True
        self.acq_start = time.localtime(os.path.getctime(self.file_path))
        self.PMT_mode = 'Correction'
        if 'excorr' in self.file_path:
            self.RunType = self.run_types.Excitation
        elif 'emcorr' in self.file_path:
            self.RunType = self.run_types.Emission
        self.num_samples = -100
        for i, line in enumerate(thefile):
            if i==3:
                self.num_samples = int(line)
            elif i==6:
                wrds = line.split()
                if self.RunType == self.run_types.Excitation:
                    self.ex_range = [float(wrds[0])]
                elif self.RunType == self.run_types.Emission:
                    self.em_range = [float(wrds[0])]
                else:
                    print("ERROR!! Bad correction file (the names need excorr or emcorr).")
                    success = False
                    break
            elif i==6+self.num_samples-1:
                wrds = line.split()
                if self.RunType == self.run_types.Excitation:
                    self.ex_range += [float(wrds[0])]
                elif self.RunType == self.run_types.Emission:
                    self.em_range += [float(wrds[0])]
                break
        return success
    
    def _ReadWLRangeLine(self, line):
        success = True
       
        # Get PMT mode 
        if line[0]=='D':
            self.PMT_mode = 'Digital'
        elif line[0] == 'A':
            self.PMT_mode = 'Analogue'
        else:
            self.PMT_mode = 'Unknown'
            success = False
        
        # Get the ranges for exitation and emission
        words = line.split()
        self.ex_range = [float(val) for val in words[1].split(':')[0].split('-')]
        self.em_range = [float(val) for val in words[1].split(':')[1].split('-')]
        
        
        if len(self.ex_range)>1 and len(self.em_range)>1:
            self.RunType = self.run_types.Synchronous
        elif len(self.ex_range)>1:
            self.RunType = self.run_types.Excitation
        elif len(self.em_range)>1:
            self.RunType = self.run_types.Emission
        else:
            self.RunType = self.run_types.Unknown
            success = False
        return success

    def ReadSpecData(self):
        '''
        Read the data from the file.
        
        If the file is a trace, this will read:
        - WL (list of wavelengths)
        - Trace (the trace - the caller is expected to know what it is)
        If the file is a session, this will also read:
        - Spec (the spectrum)
        - SpecRaw (uncorrected for excitation/emission)
        - ExCorr (the excitation correction data from the photodiode)
        '''
        if self.file_type == self.file_types.Session:
            self._ReadSessionData()
        elif self.file_type == self.file_types.Trace:
            self._ReadTraceData()
        elif self.file_type == self.file_types.Group:
            self._ReadGroupData()
        return
        
    def _ReadSessionData(self):
        NoCorr = False
        
        self.wavelengths = numpy.genfromtxt(self.file_path,
                                      skip_header=8,
                                      max_rows  = self.num_samples,
                                      usecols = [0])

        self.raw_data = numpy.genfromtxt(self.file_path,
                                      skip_header=8,
                                      max_rows  = self.num_samples,
                                      usecols = [1])

        self.cor_data = numpy.genfromtxt(self.file_path,
                                      skip_header=8,
                                      max_rows  = self.num_samples,
                                      usecols = [3])
        self.step_size = self.wavelengths[1] - self.wavelengths[0]

        excorr_header_line_num = -1
        f = open(self.file_path, 'r')
        for tup in enumerate(f):   
            if 'excorr' in tup[1].lower():
                excorr_header_line_num = tup[0] + 1
                break
        
        self.diode =  numpy.genfromtxt(self.file_path,
                                      skip_header=excorr_header_line_num+1,
                                      max_rows  = self.num_samples,
                                      usecols = [1])
        

        f.close()
    
    def _ReadTraceData(self):
        self.wavelengths = numpy.genfromtxt(self.file_path,
                                      skip_header=4,
                                      max_rows  = self.num_samples,
                                      usecols = [0])        
            
        read_data = numpy.genfromtxt(self.file_path,
                                      skip_header=4,
                                      max_rows  = self.num_samples,
                                      usecols = [1]) 
        with open(self.file_path,'r') as thefile:
            for i, line in enumerate(thefile):
                if i == 2:
                    if 'COR' in line:
                        self.cor_data = read_data
                    else:
                        self.raw_data = read_data

        return

    def _ReadGroupData(self):
        with open(self.file_path, 'r') as thefile:
            for i, line in enumerate(thefile):
                if i > 5 and i < (6 + self.num_samples):
                    wrds = line.split()
                    self.WL[i-6] = float(wrds[0])
                    self.Trace[i-6] = float(wrds[1])
                    self.UTrace[i-6] = abs(float(wrds[1]))**(0.5)
        return

    def get_date(self, space="     "):
        date = "%d-%d-%d%s%d:%02d" %(self.acq_start.tm_year,
                        self.acq_start.tm_mon,
                        self.acq_start.tm_mday,
                        space,
                        self.acq_start.tm_hour,
                        self.acq_start.tm_min)
        return date

    def plot(self, fig_size = (10,10)):
        from matplotlib import pyplot as plt

        fig = plt.figure(figsize = fig_size)
        ax1 = plt.subplot2grid((2,2), (0,0))
        ax2 = plt.subplot2grid((2,2), (0,1))
        ax3 = plt.subplot2grid((2,2), (1,0), colspan = 2)

        ax1.plot(self.wavelengths, self.raw_data)
        ax2.plot(self.wavelengths, self.diode)
        ax3.plot(self.wavelengths, self.cor_data)

        

        date = self.get_date()

        fig.suptitle(self.file_path + '\n' + date)
        ax1.set_title("Raw Data", fontsize = 15)
        ax2.set_title("Diode Signal", fontsize = 15)
        ax3.set_title("Fully Corrected Data", fontsize = 15)

        for ax in [ax1, ax2, ax3]:
            ax.grid()
            ax.set_xlim([self.wavelengths[0], self.wavelengths[-1]])
        return fig
       

    def __add__(self, other):
        new = copy.deepcopy(self)
        new.wavelengths= (new.wavelengths + other.wavelengths) / 2.0
        new.raw_data = (new.raw_data + other.raw_data) / 2.0
        new.cor_data = (new.cor_data + other.cor_data) / 2.0
        new.diode = (new.diode + other.diode) / 2.0
        new.file_path = "Merged"
        return new
            
        
