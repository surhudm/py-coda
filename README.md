# py-coda
Python routines for MCMC diagnostics using the excellent MCMC diagnostic tools
from R.

Installation:
=============

1) Install R (https://www.r-project.org/) and rpy2

On Mac with Homebrew
```bash
brew tap homebrew/science
brew install r
pip install rpy2
```

or Ubuntu
```bash
sudo apt-get install r-base r-base-dev python-rpy2
```

2) Install the coda package

```bash
R CMD INSTALL mypkg -l $HOME/Rpackages/
cat R_LIBS_USER="$HOME/Rpackages" > ~/.Renviron
```

3) Install py-coda

```bash
git clone https://github.com/surhudm/py-coda.git
cd py-coda
python setup.py install
```

Run MCMC diagnostics:
=====================

```python
from py-coda import 
mcmcobj = read_pycoda("chainfile.trimmed.out", "chainnew.ind", thin=1)
mcmcobj.geweke()
mcmcobj.get_stats()
mcmcobj.plot_traces()
mcmcobj.plot_autocorr()
```
