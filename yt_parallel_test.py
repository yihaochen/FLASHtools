#!/usr/bin/env python
import yt
yt.enable_parallelism()

ds = yt.load("c")
v, c = ds.find_max("density")
print v, c
p = yt.ProjectionPlot(ds, "x", "density")
p.save()
