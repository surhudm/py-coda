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
mkdir $HOME/Rpackages
echo R_LIBS_USER="$HOME/Rpackages" > ~/.Renviron
Rscript -e 'install.packages("coda", repos="http://cran.us.r-project.org")'
```

3) Install py-coda

```bash
git clone https://github.com/surhudm/py-coda.git
cd py-coda
python setup.py install
```

Run MCMC diagnostics:
=====================

Prepare MCMC file say chainfile.out where every line is the parameter set at
every step of the chain. Prepare another file where each line lists the
parameter name, write this parameter name in latex form, as it will get set as a
label on the plots. Then you could run the diagnostics as follows:

```python
from py_coda import read_pycoda
mcmcobj = read_pycoda("chainfile.out", "chainnew.ind", thin=1)
mcmcobj.geweke()
mcmcobj.get_stats()
mcmcobj.plot_traces()
mcmcobj.plot_autocorr()
mcmcobj.heidelberger_welch()
```

In your python code you could even use the python bindings for coda directly.
For example, if you have a numpy array with shape (Nchain, Nparam) which stores the
Monte Carlo Markov Chain (MCMC) where Nchain is the number of iterations of the
MCMC and Nparam is the number of parameters, and an array of labels for your
parameters, you could simply do
```python
from py_coda import mcmc
mcmcobj = mcmc(labels, data, thin=1)
mcmcobj.geweke()
mcmcobj.get_stats()
mcmcobj.plot_traces()
mcmcobj.plot_autocorr()
mcmcobj.heidelberger_welch()
```

Please file bug reports or issues if you find any.
