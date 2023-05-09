sh = SimpleHeatExchanger()

sh.bodyIds = [-2,-1,0,1]
sh.clumpIds = [-1,-1,111,111]
sh.mass = [1,1,1,1]
sh.cap = [1,2,3,4]
sh.T = [10,10,10,10]
sh.iterPeriod = 1

O.engines += [sh]

O.run(2, wait=True)

sh.dict()


"""
# change

sh.bodyIds = [-2,-1,0,1,2]
sh.clumpIds = [-1,-1,111,111,11]
sh.mass = [1,1,1,1,1]
sh.cap = [1,2,3,4,5]
sh.T = [10,10,10,10,10]

"""
