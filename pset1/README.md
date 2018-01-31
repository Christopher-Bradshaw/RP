# Question 3, Pset 1

Run with `python3 q3.py` (requires numpy, scipy and matplotlib). This will print some results and pop up some graphs.

Integration uses scipy integrate (with some slightly clever (maybe?) limit finding). Peak finding uses the `find_planck_peaks` which relies on the fact that the functions we care about increase, peak, decrease. We can start low, increase until we start to decrease, and then find the peak.
