#!/usr/bin/python

'''
Created on Oct 7, 2014
@author: ralonso
04/09/17: 16:15:00 Se cambio a process_analyze_paralellel> para probar carpetas disintas.

'''

import matplotlib
matplotlib.use('Agg')
import os
import glob
from matplotlib import dates
import datetime
import time
import traceback
from optparse import OptionParser
import sys 
import copy

#import plot_results_ric as ric_plt
#import process_analyze
#import process_analyze_b as process_analyze#original
import analyze as process_analyze
#import process_analyze_recuperacion as process_analyze
#import process_plot
import shutil

#
# Interferometry:
# colorsys.hsv_to_rgb(hue,sat,val)
# colorsys.hls_to_rgb
#

class FreqSettings():
    """
    I use this class to describe the options of a specific operating frequency
    """
    def __init__(self, id, label, subdirectory, channel_dirs=[],  value=0, procdata_folder="", end_flag=''):
        self.id = id
        self.label = label
        self.subdirectory = subdirectory
        self.channel_dirs = channel_dirs 
        self.value = value
        self.procdata_folder = procdata_folder
        self.end_flag = end_flag

class AutoMode():
    
    def __init__(self, procdata_folder="", graphics_folder="", rawdata_from_doy_folder="", frequency="", N_max=30, n_max_bucles=10, station_name="", incoherent_int = 0,send_graphs=0, auto_delete_rawdata=0,mode_campaign=1,profiles_number=600,stationtx_codes='0,2'):
        
        self.procdata_folder = procdata_folder #
        self.graphics_folder = graphics_folder
        self.rawdata_from_doy_folder = rawdata_from_doy_folder #
        self.frequency = frequency
        self.N_max = N_max
        self.n_max_bucles = n_max_bucles
        self.station_name=station_name
	self.incoherent_int      = incoherent_int
        self.send_graphs         = send_graphs
        self.convert_images      = 1
        self.auto_delete_rawdata = auto_delete_rawdata
	self.mode_campaign       = mode_campaign
	self.profiles_number     = profiles_number
	self.stationtx_codes     = stationtx_codes
        

	#filter
	search_station= self.stationtx_codes.split(',')
	for i in range(len(search_station)):
		print search_station[i]
		try:
			int(search_station[i])
		except:
			print "Corregir el numero de codificacion  %d Tx que es equivalente a %s"%(i,search_station[i])
			sys.exit(1)

        print "Auto Mode v 1.0.2.1_beta"
        print "Number to analyze"+str(self.N_max)
        
        f0 = FreqSettings(id="f0", label="2.72MHz", subdirectory="0", channel_dirs=['0','2'],
                          value=2.72e6, procdata_folder="sp01_f0", end_flag='procdata_end_f0_flag')
        f1 = FreqSettings(id="f1", label="3.64MHz", subdirectory="1", channel_dirs=['1','3'],
                          value=3.64e6, procdata_folder="sp01_f1", end_flag='procdata_end_f1_flag')
        
        if self.frequency == "":
            print "Choose a frequency"
            sys.exit(10)
            
        if self.frequency == "f0":
            print "Choosing only f0"
            self.frequencies=[f0]
            
        if self.frequency == "f1":
            print "Choosing only f1"
            self.frequencies=[f1]
        
        
    def _run(self):
        
        nth_bucles=0
        send_web=0
        while nth_bucles < self.n_max_bucles:
            
            nth_bucles = nth_bucles + 1            
            print "**==> Bucle: "+str(nth_bucles)+" of "+str(self.n_max_bucles)
        
            print "===== Analyzing DOY folders ====="      
            freq_i = self.frequencies[0]      

            if self.procdata_folder == "":
                print "PROCDATA FOLDER DOES NOT EXIST"
                sys.exit(10)
            
            if self.rawdata_from_doy_folder == "":
                
                print "There's not a RAWDATA FOLDER, ONLY PROCESS "
                procdata_doy_folder = self.procdata_folder
                graphics_doy_folder = self.graphics_folder
            
            elif self.rawdata_from_doy_folder != "":
                print "==> Option RAWDATA FROM DOY FOLDER <=="
                print self.rawdata_from_doy_folder 
                #Removing the las slash if exists!
                if self.rawdata_from_doy_folder[-1] == '/':
                    #print "Removing last slash /"
                    self.rawdata_from_doy_folder = self.rawdata_from_doy_folder[:-1]
            
                rawdata_parent_folder, first_doy_folder = os.path.split(self.rawdata_from_doy_folder)
                dirns = []
                dirns = dirns + glob.glob("%s/d*"%(rawdata_parent_folder))
                dirns.sort()
                
                #Get the sublist of folders starting on my initial folder
		try:	
	                dirns = dirns[dirns.index(self.rawdata_from_doy_folder):]
        	except:
                        dirns = dirns[dirns.index(dirns[0]):]
    
                #Remove folder that contain the procdata end flag
                temp = copy.deepcopy(dirns)
                for dir in temp:
                    #print "Checking dir: %s"%(dir
			
                    #if os.path.isfile("%s/procdata_end_flag"%(dir)):
                    if os.path.isfile("%s/%s"%(dir, freq_i.end_flag)):
                        #print "removing folder: %s"%(dir)
                        print "The folder %s have been already processed ..."%(dir)
                        dirns.remove(dir)

                        for channel_i in freq_i.channel_dirs:
                            rawdir_ch_i = "%s/%s"%(dir, channel_i)			   
                            
                            #If folder exists
                            if os.path.isdir(rawdir_ch_i):
                                print " ------  > Dir exists: %s"%(rawdir_ch_i)
                                if self.auto_delete_rawdata == 1:
                                    
                                    #Wait 10 seconds if the deleting flag exists
                                    while os.path.isfile("%s/deleting_flag"%(dir)):
                                        print "Other process is deleting files, I'm going to wait for 10 seconds ..."
                                        time.sleep(10)
                                    
                                    
                                    os.system("touch %s/deleting_flag"%(dir))
                                    print "I'm about to delete folder in 10 more seconds: %s"%(rawdir_ch_i)
                                    time.sleep(10)
                                    print " ------  > Dir to delete: %s"%(rawdir_ch_i)
                                    
                                    list_files_to_delete = open("%s/timestamps.log"%(rawdir_ch_i))
                                    
                                    var_i = 0
                                    for file_i in list_files_to_delete:
					print "que esta pasando 0"
					print file_i
				        #rint file_i,file_i.split()[0]
					print "que esta pasando 1"
                                        var_i = var_i + 1 
					try: 
                                         if os.path.isfile(file_i.split()[0]):
                                             print "borrando"
                                             os.remove(file_i.split()[0])
					except:
                                         print "except"
					 continue
                                        if var_i > 1000:
                                            print "Waiting 1 second :D"
                                            time.sleep(1)
                                            var_i=0
                                    del var_i
                                    
                                    shutil.rmtree(rawdir_ch_i)                                    
                                    os.system("rm %s/deleting_flag"%(dir))
                                    print "Delete was done..."    

                            else:
                                print " ------  > Dir DOESN'T exists: %s"%(rawdir_ch_i)
                        
                del temp
                
                #Check if the
                if len(dirns) == 0:
                    print "There's no folders to analyze"
                    sys.exit(0)
                
                #Get the folder to analyze and if rawdata_end_flag exists!
                dirns = [dirns[0]]
                rawdata_end_flag = os.path.isfile("%s/rawdata_end_flag"%(dirns[0]))
                rawdata_parent_folder, new_first_doy_folder = os.path.split(dirns[0])
                procdata_doy_folder = "%s/%s"%(self.procdata_folder, new_first_doy_folder[:8])
                graphics_doy_folder = "%s/%s"%(self.graphics_folder, new_first_doy_folder[:8])
                
                print "raw data folder to analyze: %s "%(dirns[0])
                print "rawdata_end_flag: %s"%(rawdata_end_flag)
                print "procdata_doy_folder: %s"%(procdata_doy_folder)
                
                """
                At this point I have the folder to work (rawdata DOY folder)
                And if rawdata_end_flag_exist
                And procdata_doy_folder where procdata going to be save
                
                3 vars:
                    rawdata_end_flag
                    dirns[0]
                    procdata_doy_folder
                """
            
                folder = dirns[0]
                print "=> Analyzing folder: %s"%(folder)
                    
                print "    => Frequency %s"%(freq_i.label)
                print "procdata_doy_folder: %s"%(procdata_doy_folder)
		if self.mode_campaign == 1:			
			print "MODE CAMPAING ON"
			self.incoherent_int=0
			self.profiles_number=600

		else:
			print "Mode CAMPAING OFF"
			self.profiles_number=100	

                print "Profiles Number value : %s"%(self.profiles_number)
                print "Incoh Int. value : %s"%(self.incoherent_int)

		
                #ric_plt.analyze_dirs_ric(dirn=folder, proc_folder=procdata_doy_folder, freq=freq_i, max_N_analyze=self.N_max)
                process_analyze.analyze_dirs_ric(dirn=folder, proc_folder=procdata_doy_folder, freq=freq_i, max_N_analyze=self.N_max,inc_int=self.incoherent_int,profiles_number=self.profiles_number,stationtx_codes=self.stationtx_codes)


            #print "    => Generating graphics ..."
            #ric_plt.plot_ric(dirn=folder, freq=freq_i.id, label_freq=freq_i.label)
            #ric_plt.plot_ric(procdata_doy_folder=procdata_doy_folder, freq=freq_i, graphics_doy_folder = graphics_doy_folder)
            #process_plot.plot_ric(procdata_doy_folder=procdata_doy_folder, freq=freq_i, graphics_doy_folder = graphics_doy_folder)
            #print "Graphics generated successfully!!! (I hope)"

	    self.convert_images=0            
            if self.convert_images == 1:
                print "Converting Images ... "
                #Creating directory 'figures'
                os.system("mkdir -p %s/%s/figures"%(graphics_doy_folder, freq_i.procdata_folder))
                
                #print "List of PNG files"
                png_files = glob.glob("%s/%s/*.png"%(graphics_doy_folder, freq_i.procdata_folder))
                png_files.sort()
                for i in png_files:
                    print i
                    tmp_file = os.path.basename(i)
                    new_jpeg_file = tmp_file[:-4]+".jpeg"
                    #new_jpeg_file = tmp_file[:-4]+"_N.jpeg"
                    #print "convert %s %s/%s/cspec/figures/%s"%(i, folder, freq_i.subdirectory, new_jpeg_file)
                    os.system("convert %s %s/%s/figures/%s"%(i, graphics_doy_folder, freq_i.procdata_folder, new_jpeg_file))

            #print "......... Continue NEXT ..........."
            #continue

            if self.send_graphs == 1:
                """ Sending images """
                # Make remote folder
                send_web= send_web+1
                print "Web counter",send_web
                if send_web==10:
		  send_web=0
		  file_1 = os.path.basename(png_files[0])
		  thisyear = file_1[14+8:18+8]
		  thismonth = file_1[19+8:21+8]
		  thisday = file_1[22+8:24+8]
		  
		  remote_folder="/home/wmaster/web2/web_signalchain/data/JRO/%s/%s/%s/%s/figures/"%(self.station_name, thisyear, thismonth, thisday)
		  print "Remote_folder: %s"%(remote_folder)
		  command_1="mkdir -p %s"%(remote_folder)
		  os.system("ssh wmaster@181.177.232.125 %s"%(command_1))

		  #os.system("ls %s/%s/cspec/figures/*.jpeg"%(folder, freq_i.subdirectory))
		  temp_command = "scp %s/%s/figures/*.jpeg wmaster@181.177.232.125:%s"%(graphics_doy_folder, freq_i.procdata_folder, remote_folder)
		  print "#==> Temp_Command:"
		  print temp_command 
		  os.system(temp_command)
		  #os.system("ls -f %s"%(folder))
                
            else:
                print "..... No Sending Images ..... "

            print "\n"
                 
            #print "\n============================"
            #print "Taking 10 seconds of break"
            #time.sleep(8)
            #print "Left only 2 seconds Don't stop"
            #time.sleep(2)
            
            print "\n============================"
            print "Starting again"
            

if __name__ == '__main__':
    """
    Testing Auto Mode
    
    python hfrx2_auto_mode_ric_v2.py --freq=f0 -N 90 -x 10 
    --rawdata-from-doy="/full/path/rawdata/dYYYYDOY/"  
    --procdata="/full/path/procdata/" 
    --graphics="/full/path/graphics/"
    
    """
    
    parser = OptionParser(usage="%prog: [options]")

    parser.add_option("", "--procdata",
                      dest="procdata",
                      type="string",
                      action="store",
                      default="",
                      help="procdata specific directory.")
    
    parser.add_option("", "--graphics",
                      dest="graphics",
                      type="string",
                      action="store",
                      default="",
                      help="graphics specific directory.")    

    parser.add_option("", "--rawdata-from-doy",
                      dest="rawdata_from_doy",
                      type="string",
                      action="store",
                      default="",
                      help="rawdata specific directory to start, after that follows the next days.")

    parser.add_option("-f", "--freq",
                      dest="frequency",
                      type="string",
                      action="store",
                      default="",
                      help="frequency")

    parser.add_option("-N", "",
                      dest="N_max",
                      type="int",
                      action="store",
                      default=30,
                      help="Max number of hdf5 files to create each bucle")

    parser.add_option("-I", "--incoherent-integration",
                      dest="incoherent_integration",
                      type="int",
                      action="store",
                      default=0,
                      help="Integrate N blocks of 10 sec each one, and save data in spectra mode.")

    parser.add_option("-x", "--max-bucles",
                      dest="max_bucles",
                      type="int",
                      action="store",
                      default=10,
                      help="Max number of bucles before program dies")

    parser.add_option("", "--station-name",
                      dest="station_name",
                      type="string",
                      action="store",
                      default="HF",
                      help="Name to identify this station (Ex: HF, HFT)")
    
    parser.add_option("", "--send-graphs",
                      dest="send_graphs",
                      type="int",
                      action="store",
                      default=1,
                      help="Option to send graphs to the web page (1: YES, 0: NO)")    

    parser.add_option("", "--auto-delete-rawdata",
                      dest="auto_delete_rawdata",
                      type="int",
                      action="store",
                      default=0,
                      help="Option to auto delete rawdata (1: YES, 0: NO)")    
 
    parser.add_option("", "--mode-campaign",
                      dest="mode_campaign",
                      type="int",
                      action="store",
                      default=1,
                      help="Option to Select Mode Campaign  (1: 600 profiles and 0 incoherent integration, 0: 100 profiles and incoherent integration depend of user)")

##### ULTIMA ADICION 20/02/2018
    parser.add_option("", "--stationtx-codes",
                      dest="stationtx_codes",
                      type="str",
                      action="store",
                      default="0,2",
                      help="Option to process Data in Rx from any Tx Station, the input is the number of code used in Tx for example: {0}, {0,1} , {0,2} o {0,1,2}")





    (op, args) = parser.parse_args()
    
    if op.station_name == "":
        print "ERROR: Choose a valid name for the station"
        sys.exit(1)
  
    currentdoy=(datetime.datetime.now()-datetime.timedelta(minutes = 10)).strftime("%j")
    #Let to adquisition process be 10 minutes forward
    currentyear=datetime.datetime.now().strftime("%Y")
    daystring = 'd'+ currentyear + currentdoy
    #print 'THE CODE IS HARDCODED!!!!'
    #daystring = 'd2018008_1722'
    #------------new--------------------------------  
    if op.mode_campaign == 1:
        print "INSIDE MODE CAMPAIGN",op.procdata
        try:
                op.procdata= op.procdata+"CAMPAIGN/"
        except  OSError:
                if not os.path.isdir(path):
                        raise

        print "NEW", op.procdata
        time.sleep(4)
    #----------------------------------------------- 
    print daystring
    my_auto = AutoMode(procdata_folder=op.procdata, graphics_folder=op.graphics,
                       rawdata_from_doy_folder=(op.rawdata_from_doy+daystring), frequency=op.frequency, 
                       N_max=op.N_max, n_max_bucles=op.max_bucles, 
                       station_name=op.station_name, send_graphs=op.send_graphs, auto_delete_rawdata=op.auto_delete_rawdata,incoherent_int=op.incoherent_integration,mode_campaign = op.mode_campaign,stationtx_codes=op.stationtx_codes)
    my_auto._run()
    print "... This program have finished"
