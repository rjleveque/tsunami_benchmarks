
Run using the changes to GeoClaw for proposed version 5.4.0, and
changes to riemann solver proposed by Wenwen Li, on riemann branch
geoclaw_sLsR.

Also reduced resolution and added Manning coefficient on 
flat part of domain to better match experiments.

===========
amrclaw
===========
/Users/rjl/git/clawpack/amrclaw

--- last commit ---
826b95c Merge pull request #154 from mjberger/bcseq

--- branch and status ---
## master...origin/master
 M examples/run_tests.py


===========
clawutil
===========
/Users/rjl/git/clawpack/clawutil

--- last commit ---
0465a1f Merge pull request #88 from ortegagingrich/library

--- branch and status ---
## master...origin/master


===========
riemann
===========
/Users/rjl/git/clawpack/riemann

--- last commit ---
1ef063f modication of geoclaw sLsR solver to use one-sided speeds rather than Roe averages

--- branch and status ---
## geoclaw_sLsR


===========
geoclaw
===========
/Users/rjl/git/clawpack/geoclaw

--- last commit ---
8e9bca2 Merge pull request #194 from mandli/fix-template-make

--- branch and status ---
## master...origin/master
 M examples/multi-layer/plane_wave/setplot.py
 M examples/run_tests.py
 M examples/tsunami/chile2010/setplot.py
 M src/python/geoclaw/dtopotools.py
 M src/python/geoclaw/fgmax_tools.py
 M src/python/geoclaw/topotools.py
