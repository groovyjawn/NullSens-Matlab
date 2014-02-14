
# This script is used to generate parameter sweep experiments for NullSens.
# The output file is formated for running the experiments on servers running Condor:
# http://research.cs.wisc.edu/htcondor/ 
# Copyright (C) 2014  Steven D. Essinger

# This file is part of NullSens-Matlab.

# NullSens-Matlab is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


outFile = open('experiments_to_run.txt', 'w') # Output file containing all experiemnts to run via Condor
sites_to_run = [20,40,60,80,100,200,500,1000] 

# False Positive Experiments
for Type in range(0, 3, 1): # Type of community covariation (0 Negative, 1 Positive, 2 Mixed)
    for n in sites_to_run: # Number of sites to include in the community
        for p in range(20, 21, 20): # Number of species to include in the community
            for q in range(1,5,1): # Number of environmental factors (gradients) to include
                for N in range(1,13,1): # Noise parameter selection - integer values 1 to 12
					M = 1; CD = 0 # False positive experiments so not needed
					name = str(n) + '_' + str(p) + '_' + str(q) + '_' + str(N) + '_' + str(M)+ '_' + str(CD) + '_' + str(Type) # Name of file to save each experiments
					outFile.write('executable = /usr/local/bin/matlab' + '\n') # Use Matlab to call the program
					outFile.write('arguments = -nodisplay -singleCompThread -nojvm -r NullSens('+ str(n) + ',' + str(p) + ',' + str(q) + ',' + str(N) + ',' + str(M) + ',' + str(CD) + ',' + str(Type) + ',\'RESULTS_' + str(name) + '\')' + '\n') # Input arguments for NullSens.m
					outFile.write('error = Null_Sens_error_' + str(name) + '.txt' + '\n')
					outFile.write('Requirements = Memory > 1024 \n') # Requirements for job submission - see Condor documentation for more information http://research.cs.wisc.edu/htcondor/ 
					outFile.write('queue' + '\n\n') # Add the job to Condor queue.

# True Positive Experiments
covary_pairs_to_run = [1,3,10] # Number of covaring pairs to include in the community
for Type in range(0, 3, 1):
    for n in  list_n:
        for p in range(60, 61, 20):
            for q in range(1, 5, 1):
                for N in range(1,13,1):
                    for M in range(1,6,1): # Covariation magnitude parameter selection - integer values 1 to 5
                        for CD in covary_pairs_to_run: 
                        	name = str(n) + '_' + str(p) + '_' + str(q) + '_' + str(N) + '_' + str(M)+ '_' + str(CD) + '_' + str(Type)
                        	outFile.write('executable = /usr/local/bin/matlab' + '\n')
                        	outFile.write('arguments = -nodisplay -singleCompThread -nojvm -r NullSens('+ str(n) + ',' + str(p) + ',' + str(q) + ',' + str(N) + ',' + str(M) + ',' + str(CD) + ',' + str(Type) + ',\'RESULTS_' + str(name) + '\')' + '\n')
                        	outFile.write('error = Null_Sens_error_' + str(name) + '.txt' + '\n')
                        	outFile.write('Requirements = Memory > 1024 \n')
                        	outFile.write('queue' + '\n\n')

outFile.close()
