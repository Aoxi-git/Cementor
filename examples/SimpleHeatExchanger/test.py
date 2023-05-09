sh = SimpleHeatExchanger()

sh.bodyIds = [-2,-1,0,1,2,6,7,8,9,10,11,12]
sh.clumpIds = [-1,-1,10,111,111,10,10,30,30,30,30,111]
sh.iterPeriod = 1

O.engines += [sh]

O.run(2, wait=True)

sh.dict()
