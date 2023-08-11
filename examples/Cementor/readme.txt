1. phase1_generateHostSample.py is an example code for generating a host sample with a target porosity and void ratio. 
The host sample is saved after the execution of phase1_generateHostSample.py.

2. phase2_addFinesInTargetMicrostructure.py is for generating fines with a specifical distribution pattern in the given host sample. 
The host sample generated from phase1_generateHostSample.py is imported. 
Therefore, phase2_addFinesInTargetMicrostructure.py should be executed after phase1_generateHostSample.py.

You can also use your own sample as the host sample.

3. Execute phase1_generateHostSample.py and phase2_addFinesInTargetMicrostructure.py by:
yade phase1_generateHostSample.py
yade phase2_addFinesInTargetMicrostructure.py
