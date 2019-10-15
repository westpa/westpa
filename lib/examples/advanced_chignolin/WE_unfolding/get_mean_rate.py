import sys
import h5py
import numpy


####### READ IN TAU FROM COMMAND LINE  #############
if len(sys.argv) == 3 :
   tau = float(sys.argv[1])
   window_time = float(sys.argv[2])
   window_frames = int(window_time/tau)
else :
   print "\n\nPLEASE, PROVIDE THE WESTPA RESAMPLING TIME, TAU, AND THE WIDTH FOR THE WINDOW-AVERAGING  (BOTH IN SECONDS) AS COMMAND LINE ARGUMENTS, E.G.:\n\n\t'python get_mean_rate.py  20e-12  1e-9'\n\n\n"
   sys.exit()
###################################################################################



####### READ IN THE FLUX AT EVERY ITERATION #############
fluxanl              = h5py.File('fluxanl.h5')
flux                 = numpy.array(fluxanl['target_flux']['target_0']['flux'])
###########################################################################################



###### GET THE AVERAGE, WRITE OUT THE RATE AS A FUNCTION OF MOLECULAR TIME ####################
FileOut =  open("mean_rate.dat", "w")
print>>FileOut, "#molecular time [ns]\trate [1/s]\n"

for i in range(len(flux)) :
   if i <= window_frames+1 :
      flux_mean = sum(flux[0 : (i+1)])  /  (i+1)
   else :
      flux_mean  = sum(flux[i-window_frames+1  : i+1])  / window_frames

   print>>FileOut, "{0:.3f}".format((i+1)/(1e-9/tau)), "\t\t", "{0:.8e}".format(flux_mean/tau)
###########################################################################################


