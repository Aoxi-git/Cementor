# -*- encoding=utf-8 -*-
O.bodies.append(geom.facetBox((0, 0, 0), (10, 5, 2)))
O.bodies.append(geom.facetBox((0, 0, 10), (10, 5, 2), wallMask=0b110011, wire=False))
O.bodies.append(geom.facetBox((0, 10, 5), (10, 5, 2), wallMask=0b101010, wire=False))
