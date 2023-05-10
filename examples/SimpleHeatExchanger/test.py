sh = SimpleHeatExchanger()

sh.bodyIds = [-2,-1,0,1,111]
sh.clumpIds = [-1,-1,111,111,111]
sh.mass = [1,1,1,1,2]
sh.cap = [1,2,3,4,5]
sh.cond = [1,2,3,4,5]
sh.T = [100,10,10,10,10]
sh.dummyIntA = [0.01]
sh.dummyIntId1 = [-2]
sh.dummyIntId2 = [-1]
sh.iterPeriod = 1

O.engines += [sh]

O.dt = 0.0001
O.run(2, wait=True)

sh.dict()


"""
# touch clump member
sh.needsInit = True # init energy based on temp
sh.dummyIntId1 = [-2]
sh.dummyIntId2 = [1]

"""
