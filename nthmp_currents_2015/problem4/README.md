
# Problem 4

See http://coastal.usc.edu/currents_workshop/problems/prob4.html
for the description.

Download `all_data.zip` from that page and unzip into `../all_data`
for some of the scripts in this directory to work.

* Run `save_bathy.m` in Matlab to create text files `x_onshore.txt`,
  `y_onshore.txt`, `z_onshore.txt` that have the onshore topography
  in a form readable in Python.

* Execute

    $ python maketopo.py

  to create `pwlinear2.topotype1` and `seaside_onshore.tt2` to be read by
  Geoclaw.

* Execute

    $ python make_fgmax.py

  to create `fgmax_grid.txt` used to specify region where maximum 
  amplitude, speed, etc. should be monitored.

* Run GeoClaw and plot the results:

    $ make .plots


