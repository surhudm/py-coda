import numpy as np
import sys as sys
import pandas

class IOError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class mcmc():
    def __init__(self, labels, data, thin):
        '''
        Initializes mcmc object

        :Parameters:

        -   labels : An array of strings containing the latex labels for the
            parameters
        -   data : An array containing the output of the Monte Carlo Markov
            chain. The array should have shape (Niterations, Nparam)
        -   thin : Specify the thinning parameter, this will use only one
            sample every thin iterations.

        :Returns:

        -   MCMC object

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> help(a)
        '''

        from rpy2 import robjects
        from rpy2.robjects.packages import importr
        self.coda = importr("coda")

        self.labels = labels

        self.nparam = len(self.labels)
        data_columns = data[0].size

        if data_columns < self.nparam:
            raise IOError("The number of columns in chainfile are less than equal to number of labels in indexfile")

        if data_columns > self.nparam:
            data = data[:, :self.nparam-data_columns]

        data = np.transpose(data)
        data = data[:, ::thin]

        self.chain_elements = data[0].size

        self.codamcmc = self.coda.mcmc(robjects.r['matrix'](robjects.FloatVector(data.flatten()), nrow=self.chain_elements))

    def plot_traces(self, n_at_time=5, backend="TkAgg", savefig=0):
        '''
        Generate trace plots

        :Parameters:

        -   n_at_time : Optional argument to specify the number of parameters to
            be plotted on one page, the subplots are arranged in n_at_time rows,
            with default value equal to 5
        -   backend : Optional argument to specify the matplotlib backend to be
            used. By default uses the TkAgg backend.
        -   savefig : Optional string argument which is appended to the figure
            names, if non-empty, the plots are not displayed. By default figures
            are not saved, but just displayed.

        :Returns:

        -   0

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> a.plot_traces(n_at_time=5, backend=\"TkAgg\", savefig=\"MCMC\")

        '''
        import matplotlib
        # matplotlib.use(backend)
        import pylab as pl

        if not savefig:
            pl.ion()

        sys.stdout.write("Computing traces...\n")

        fignumber = 0
        for i in range(self.nparam):
            ax = pl.subplot(n_at_time, 1, i%n_at_time+1)
            array = np.array(self.codamcmc[i*self.chain_elements:(i+1)*self.chain_elements])
            ax.plot(range(self.chain_elements),
                    array,
                    color='grey')
            ax.set_ylabel(r"%s"%self.labels[i])

            if i%n_at_time==n_at_time-1 or i==self.nparam-1:
                ax.set_xlabel("Iteration number")
                print("Iteration number : ", i)

                if savefig:
                    pl.savefig(savefig+"Trace_%d.png" % fignumber)
                    fignumber += 1
                else:
                    input("press any key to continue")
                pl.clf()
            else:
                ax.set_xticklabels([])

        if not savefig:
            pl.ioff()

        return 0

    def get_stats(self, fp=sys.stdout):
        '''
        Display statistics of the chain, such as the 2.5, 16, 50, 84 and 97.5
        percentiles

        :Parameters:

        -   fp : Optional python file object, default is sys.stdout

        :Returns:

        -   0

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> a.get_stats()

        '''
        sys.stdout.write("Computing stats...\n")

        i025 = int(0.025 * self.chain_elements)
        i16 = int(0.16 * self.chain_elements)
        i50 = int(0.50 * self.chain_elements)
        i84 = int(0.84 * self.chain_elements)
        i975 = int(0.975 * self.chain_elements)

        # Let us compute the mean and percentiles
        list_to_print = ["Param", "0.025%", "16%", "50%", "84%", "97.5%"]
        fp.write('\t'.join([x for x in list_to_print])+"\n")
        for i in range(self.nparam):
            array = self.codamcmc[i*self.chain_elements:(i+1)*self.chain_elements]
            array = np.sort(array)
            list_to_print = [self.labels[i], "%.5e" % array[i025], "%.5e" %
                    array[i16], "%.5e" % array[i50], "%.5e" % array[i84], "%.5e" %
                    array[i975] ]
            fp.write('\t'.join([x for x in list_to_print])+"\n")

    def plot_autocorr(self, n_at_time=6, backend="TkAgg", savefig=0):
        '''
        Generate plots with autocorrelations

        :Parameters:

        -   n_at_time : Optional argument to specify the number of parameters to
            be plotted on one page, the subplots are arranged in n_at_time rows,
            with default value equal to 5
        -   backend : Optional argument to specify the matplotlib backend to be
            used. By default uses the TkAgg backend.
        -   savefig : Optional string argument which is appended to the figure
            names, if non-empty, the plots are not displayed. By default figures
            are not saved, but just displayed.

        :Returns:

        -   0

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> a.plot_autocorr(n_at_time=5, backend=\"TkAgg\", savefig=\"MCMC\")

        '''

        import matplotlib
        # matplotlib.use(backend)
        import pylab as pl

        if not savefig:
            pl.ion()

        sys.stdout.write("Computing autocorrelations, please be patient...\n")
        fignumber = 0
        for i in range(self.nparam):
            ax = pl.subplot(n_at_time, 1, i%n_at_time+1)

            array = self.codamcmc[i*self.chain_elements:(i+1)*self.chain_elements]
            f = np.fft.fft(array-np.mean(array), n=2*self.chain_elements-1)
            acf = np.real(np.fft.ifft(f * np.conjugate(f)))
            acf = acf[:self.chain_elements]
            acf /= acf[0]
            ax.plot(np.arange(self.chain_elements), acf, color='grey')
            ax.set_ylabel(r"%s"%self.labels[i])

            if i%n_at_time==n_at_time-1 or i==self.nparam-1:
                ax.set_xlabel("Iteration number")

                if savefig:
                    pl.savefig(savefig+"autocorr_%d.png" % fignumber)
                    fignumber += 1
                else:
                    input("press any key to continue")
                pl.clf()
            else:
                ax.set_xticklabels([])

        if not savefig:
            pl.ioff()

        return 0

    def heidelberger_welch(self, pvalue=0.05, eps=0.1, fp=sys.stdout):
        '''
        Performs the Heidelberger Welch test for stationarity and half width.

        :Parameters:

        -   pvalue : Optional pvalue to be used to decide the success or failure
            of the stationarity test (default=0.05)
        -   eps : Optional argument which sets a tolerance level to decide the
            success of the halfwidth test (default=0.1)

        :Returns:

        -   Prints diagnostic text output, can optionally return the
            quantitative results of the test in an array as well

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> a.heidelberger_welch()
        '''

        hw_output = self.coda.heidel_diag(self.codamcmc, eps=eps, pvalue=pvalue)
        sys.stdout.write("Computing Heidelberger Welch stationarity and the half width test, please be patient...\n")
        df = pandas.DataFrame(np.array(hw_output), columns=hw_output.colnames, index=self.labels)
        df.to_csv(fp, sep=" ", na_rep="NaN", index_label="Variable")

        return df

    def geweke(self, nbins=20, n_at_time=6, backend="TkAgg", savefig=0, fp=sys.stdout):
        '''
        Plots the z-scores of the Geweke test

        :Parameters:

        -   nbins : Optional argument which sets the number of different
            starting points in the chain to use to perform the Geweke test
            (default=20)
        -   n_at_time : Optional argument to specify the number of parameters to
            be plotted on one page, the subplots are arranged in n_at_time rows,
            with default value equal to 5
        -   backend : Optional argument to specify the matplotlib backend to be
            used. By default uses the TkAgg backend.
        -   savefig : Optional string argument which is appended to the figure
            names, if non-empty, the plots are not displayed. By default figures
            are not saved, but just displayed.
        -   fp: Optional file handle, else print to stdout

        :Returns:

        -   Generates plots with Geweke z-scores

        :Examples:

        >>> from py_coda import mcmc
        >>> a = mcmc(label_array, mcmc_array, 1)
        >>> a.geweke()
        '''

        import matplotlib
        # matplotlib.use(backend)
        import pylab as pl

        if not savefig:
            pl.ion()

        sys.stdout.write("Computing Geweke test, please be patient...\n")

        bins = np.linspace(1., self.chain_elements/2-1, nbins)
        geweke_output = np.zeros(nbins*self.nparam).reshape(nbins, self.nparam)
        for i in range(nbins):
            tmp_geweke = self.coda.geweke_diag(self.coda.window_mcmc(self.codamcmc, start=bins[i]))
            geweke_output[i, :] = tmp_geweke[0]

        fignumber = 0
        for i in range(self.nparam):
            ax = pl.subplot(n_at_time, 1, i%n_at_time+1)
            ax.fill_between(bins, -1.96, 1.96, color="grey", alpha=0.3)
            ax.scatter(bins, geweke_output[:, i], color='k', label=r"%s" % self.labels[i])
            ax.set_ylabel(r"Z-score")
            ax.legend()

            if i%n_at_time==n_at_time-1 or i==self.nparam-1:
                ax.set_xlabel("Starting iteration number")

                if savefig:
                    pl.savefig(savefig+"geweke_%d.png" % fignumber)
                    fignumber += 1
                else:
                    input("press any key to continue")
                pl.clf()
            else:
                ax.set_xticklabels([])

        if not savefig:
            pl.ioff()

        df = pandas.DataFrame(geweke_output.T, columns=bins.astype(int), index=self.labels)
        df.to_csv(fp, sep=" ", na_rep="NaN", index_label="Variable")

        return 0


def read_pycoda(chainfile, indexfile, thin=1, weight_col=None, **kwargs):
    '''
    Reads in flat text file MCMC output

    :Parameters:

    -   chainfile : Name of the chainfile to be read. Should contain the
        parameters of the model at every step on every line.
    -   indexfile : Name of the file that contains latex labels for every
        parameter
    -   thin : Optional argument to specify the thinning interval (default=1)
    -   weight_col : Optional argument to specify the thinning interval (default=1)
    -   kwargs : Optional keyword arguments that can be passed to the
        pandas.read_csv command for reading the chainfile (e.g., read only
        specific columns using usecols)

    :Returns:

    -   An py_coda mcmc object

    :Examples:

    >>> from py_coda import read_pycoda
    >>> a = read_pycoda(\"chainfile.out\", \"chainfile.labels\", thin=1,
    usecols=(0, 1, 2))
    >>> help(a)
    '''

    labels = [line.rstrip("\n") for line in open(indexfile)]
    #print (labels)

    #data = np.loadtxt(chainfile, **kwargs)
    df = pandas.read_csv(chainfile, header=None, names=labels, **kwargs)

    if weight_col is not None:
        df = df.iloc[np.repeat(np.arange(len(df)), df[weight_col].values.astype(int))]
        df.reset_index(drop=True, inplace=True)
        df.drop(weight_col, inplace=True, axis=1)
        labels.remove(weight_col)

    data = df.values

    return mcmc(labels, data, thin)


if __name__ == "__main__":
    mcmcobj = read_pycoda("jkdir/chainfile.trimmed.out", "jkdir/chainnew.ind", thin=1)
    mcmcobj.geweke()
    mcmcobj.get_stats()
    mcmcobj.plot_traces()
    mcmcobj.plot_autocorr()
    mcmcobj.heidelberger_welch()
