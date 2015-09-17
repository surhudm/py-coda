Installation:
=============

1) Install R (https://www.r-project.org/) and rpy2

On Mac with Homebrew

.. sourcecode:: bash           
 
    $ brew tap homebrew/science
    $ brew install r
    $ pip install rpy2

or Ubuntu

.. sourcecode:: bash           

    $ sudo apt-get install r-base r-base-dev python-rpy2

2) Install the coda package

.. sourcecode:: bash           

    $ mkdir $HOME/Rpackages
    $ echo R_LIBS_USER="$HOME/Rpackages" > ~/.Renviron
    $ Rscript -e 'install.packages("coda", repos="http://cran.us.r-project.org")'

3) Install py-coda

.. sourcecode:: bash           

    $ git clone https://github.com/surhudm/py-coda.git
    $ cd py-coda
    $ python setup.py install
