1. phase1_generateHostSample.py is an example code for generating a host sample with a target porosity and void ratio. 
The host sample is saved after the execution of phase1_generateHostSample.py.

2. phase2_Cementor.py is for generating fines with a specifical distribution pattern in the given host sample. 
The host sample generated from phase1_generateHostSample.py is imported. 
Therefore, phase2_Cementor.py should be executed after phase1_generateHostSample.py.

You can also use your own sample as the host sample.

3. Execute phase1_generateHostSample.py and phase2_Cementor.py by:
yade phase1_generateHostSample.py
yade phase2_Cementor.py

Note that using "yade phase2_Cementor.py", you work on one set of parameters at a time.
To try multiple sets of parameters, you may use" yade-batch batch_params.txt phase2_Cementor.py"

For more information on batch mode, you may refer to YADE examples:
https://yade-dev.gitlab.io/trunk/tutorial-examples.html?highlight=batch#oedometric-test