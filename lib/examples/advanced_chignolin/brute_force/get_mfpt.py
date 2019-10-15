import sys
import re
from bootstrap import get_CI, get_CR




#####################################
########## LOAD INPUT DATA ##########
#####################################

if len(sys.argv) != 5 :
   print "\n\nPLEASE, PROVIDE THE RMSD FILE, ITS TIME RESOLUTION (IN SECONDS), THE LOWER AND THE UPPER TARGET STATE RMSD VALUES AS COMMAND LINE ARGUMENTS, E.G.:\n\n\t'python get_mfpt.py  rmsd.dat  20e-12  0.5  4.0'\n\n\n"
   sys.exit()
else :
   FileIn = sys.argv[1]
   dt = float(sys.argv[2])
   MinVal = float(sys.argv[3])
   MaxVal = float(sys.argv[4])

#####################################






#####################################
####### DEFINE & LOAD STATES ########
#####################################

states = []
data = []
for line in open(FileIn).readlines() :
  Words = line.split()
  if (len(Words) > 0) and (re.search("\d", Words[0])) :
     data.append(float(Words[1]))

for i in data:
    if i <= MinVal:
       states.append(1)		# 'state 1' = folded
    elif i >= MaxVal:
       states.append(3)		# 'state 3' = unfolded
    else:
       states.append(2)		# 'state 2' = transition state

#####################################






#####################################
#### CALCULATE THE PASSAGE TIMES ####
#####################################

state_F = [1] 		# folded   state definition
state_U = [3]		# unfolded state definition
PT_FU = [] 		# passage time F to U
PT_UF = [] 		# passage time U to F
count = 0  		# first passage time counter

P_state = "PREVIOUS" 	# previous state

# loop through all state frames, determine the current state, calculate the passage time
for frame in states :	
   if frame in state_F:
      state = "F"
   elif frame in state_U:
      state = "U"
   else:
      state = P_state

   if (state == "F") or (state == "U"):
      count += 1
   if P_state == "F" and state == "U":
      PT_FU.append(count)
      count = 0
   elif P_state == "U" and state == "F":
      PT_UF.append(count)
      count = 0
   elif P_state == "PREVIOUS" and (state == "F" or state == "U"):
      count = 0

   P_state = state

#####################################






#####################################
######## CALCULATE THE MFPTs ########
#####################################

# check the existence of passage times
if len(PT_FU) == 0 :
   print "WARNING: No F->U events observed!"
   sys.exit()
if len(PT_UF) == 0 :
   print "WARNING: No U->F events observed!"
   sys.exit() 


# mfpts in seconds
mfpt_FU = float(sum(PT_FU)) / len(PT_FU) * dt
mfpt_UF = float(sum(PT_UF)) / len(PT_UF) * dt


# Bayesian bootstrapping for uncertainty estimation
[PT_FU_CRmin,PT_FU_CRmax] = get_CR(PT_FU, 10000) 
[PT_UF_CRmin,PT_UF_CRmax] = get_CR(PT_UF, 10000)

#####################################






#####################################
######### PRINTING RESULTS ##########
#####################################

print "MFTS (in seconds):" 
print "F-->U: mean = ", mfpt_FU, "\t", "confidence interval =  [", PT_FU_CRmin*dt, ",", PT_FU_CRmax*dt, "]"
print "U-->F: mean = ", mfpt_UF, "\t", "confidence interval =  [", PT_UF_CRmin*dt, ",", PT_UF_CRmax*dt, "]"
 
#####################################



