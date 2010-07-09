"""The metascientific statistical package."""

import scipy.stats as sps

##########################
### Contingency Tables ###
##########################

class CTScaleError(Exception):
    """Contingency Table Scale Error.
    """

    def __init__(self, scale_type, scale_str):
        """
        Args:
           * `scale_type` (str): Variable that is being scaled, such as x, y, or R.
           * `scale_str`  (str): Flag for scale, such as 'linear', 'log', or 'nines'."""

        self.scale_type = scale_type
        self.scale_str  = scale_str

    def __str__(self):
        return "Contingency table scale flag wrong! {0} scale = {1}".format(self.scale_type, self.scale_str)


#####################
### General Stats ###
#####################

def LeastSquaresSummary(lsout):
    """
    SciPy contains a very nice optimization package for fitting various curves.
    However unlike pure statistical packages (at present), the `scipy.optimize.leastsq` 
    function does not have an option to return a summary table of commonly derived results.

    LeastSquaresSummary() mimics the `summary()` function found in R though it acts on the 
    output of a scipy optimization/curve-fit.  Note that the `scipy.optimize.leastsq`
    must have the `full_output` keyword set non-zero for this function to work::

       #Generating summary info from scipy.optimize.leastsq()
       lsout = scipy.optimize.leastsq(..., full_output=1)
       lssum = metasci.stats.LeastSquaresSummary(lsout)

    This function was inspired by `Ondrej Certik on the SciPy-dev list 
    <http://mail.scipy.org/pipermail/scipy-dev/2009-March/011527.html>`_.

    Args:
       * `lsout` (tuple): Least squares optimization full output.  Tuple has the 
         following members: (x, cov_x, infodict, mesg, ier).
         See `the documentation <http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html>`_ 
         for more information.

    Returns:
       * `summary` (dict): Dictionary of summary values.  Dictionary contains the following 
         keys and values.  

         Copied from output:
            * `params`    (1D numpy array):     Fit parameter solution values :math:`p`. Copy of lsout[0].
            * `covar`     (2D numpy array):     Fit parameter covariance matrix :math:`\Sigma`. Copy of lsout[1].
            * `nfev`      (int):                Number of function calls in least-squares calculation. Copy of lsout[2]['nfev'].
            * `residuals` (1D numpy array):     Residuals :math:`\\varepsilon`. Copy of lsout[2]['fvec'] times -1.
            * `fjac`      (2D numpy array):     Jacobian approximation at output :math:`J`. Copy of lsout[2]['fjac'].
            * `ipvt`      (1D int numpy array): Permutation matrix definition. Copy of lsout[2]['ipvt'].
            * `qtf`       (1D numpy array):     Vector defining q. Copy of lsout[2]['qtf'].
            * `message`   (str):                If optimization failed, string describing failure. Copy of lsout[3]
            * `errcode`   (int):                Error code. If 1, 2, 3, or 4, then no error. Copy of lsout[4]

         Calculated:
            * `chi2` (float): Chi-Squared statistic.
              Call :math:`\\varepsilon_n` the nth point of the residual, 

              .. math:: \chi^2 = \sum_n^N \\varepsilon_n^2

            * `dof` (int): Degrees of freedom::

                 summary['dof'] = len(summary['residuals']) - len(summary['params'])

            * `stderr` (1D numpy array): Standard error for each fit parameter. 
              With :math:`\Sigma` as the covariance matrix and i as the index of the fit parameters,

              .. math:: \mbox{SE}_i = s_{e_i} = \Sigma_{ii} \cdot \sqrt{\\frac{\chi^2}{\mbox{dof}}}

            * `tscore` (1D numpy array): The t-score or t-statistic for each fit parameter.
              With :math:`p_i` as the ith fit parameter,

              .. math:: t_i = \\frac{p_i}{\mbox{SE}_i}

            * `pvalue` (1D numpy array): The P-value each fit parameter.
              Here :math:`t` is the Student's t-distribution for dof Degrees of Freedom,

              .. math:: P_i = t(t_i, \mbox{dof})

            * `cor` (2D numpy array): Correlation matrix :math:`\\rho_{ij}` for the fit parameters. 
              Calculated from the covariance matrix :math:`\Sigma`, 

              .. math:: \\rho_{ij} = \\frac{\Sigma_{ij}}{\sqrt{\Sigma_{ii}\cdot\Sigma_{jj}}}
    """

    summary = {}	
    
    summary["params"]    = lsout[0].copy()
    summary["covar"]     = lsout[1].copy()
    summary["nfev"]      = lsout[2]['nfev']
    summary["residuals"] = - lsout[2]['fvec'].copy()
    summary["fjac"]      = lsout[2]['fjac'].copy()
    summary["ipvt"]      = lsout[2]['ipvt'].copy()
    summary["qtf" ]      = lsout[2]['qtf'].copy()
    summary["message"]   = lsout[3]
    summary["errcode"]   = lsout[4]
    
    summary["chi2"]   = sum(summary["residuals"]*summary["residuals"])
    summary["dof"]    = len(summary['residuals']) - len(summary['params'])
    summary["stderr"] = np.array([np.sqrt(summary["covar"][i,i])*np.sqrt(summary["chi2"]/summary["dof"]) for i in range(len(summary["params"]))])
    summary["tscore"] = np.array([summary["params"][i]/summary["stderr"][i] for i in range(len(summary["params"]))])
    summary["pvalue"] = sps.t.pdf(summary["tscore"], summary["dof"])
    summary["cor"]    = np.array([[summary["covar"][i][j] / np.sqrt(summary["covar"][i][i] * summary["covar"][j][j]) for j in range(len(summary["covar"][i]))] for i in range(len(summary["covar"]))])

    return summary

def SummaryStr(summary, paramlabels=None, sinc="rsvc", sbuf=20, sform=".6G"):
    """
    Takes a least-squares summary dictionary and returns its string representation.

    Args:
       * `summary` (dict): output from LeastSquaresSummary().

    Keyword Args:
       * `paramlabels` (list): String labels for params, length len(summary["params"]). Default ["p0", "p1", ...].
       * `sinc`        (str):  String include flags for which summary data to print. For example, 'sv' would only 
         print the summary table and the covariance matrix.  Default prints data.  Flags have the following meanings:
            * `r`: Residual summary table.
            * `s`: Regression parameter summary.
            * `v`: Parameter covariance matrix.
            * `c`: Parameter correlation matrix.
       * `sbuf`        (int):  
       * `sform`       (str):  Data format string specification.

    Returns:
       * `s` (str): LaTeX-valid string.
    """

    salign = "^" + str(sbuf)

    I = len(summary["params"])

    if paramlabels == None:
        paramlabels = []
        for i in range(I):
            paramlabels.append("p" + str(i))

    s = ""

    if "r" in sinc:
        s = s + "Regression Residuals:\n"
        s = s + "{0:{salign}}{1:{salign}}{2:{salign}}\n".format("Min", "Median", "Max", salign=salign)
        s = s + "{0:{salign}{sform}}{1:{salign}{sform}}{2:{salign}{sform}}\n".format(summary["residuals"].min(), 
            np.median(summary["residuals"]), summary["residuals"].max(), salign=salign, sform=sform)
        s = s + "\n\n"

    if "s" in sinc:
        s = s + "Regression Parameter Summary:\n"
        s = s + "{0:{salign}}{1:{salign}}{2:{salign}}{3:{salign}}{4:{salign}}\n".format("", 
            "Value", "Std. Err.", "t", "P", salign=salign)
        for i in range(I):
            s = s + "{0:{salign}}".format(paramlabels[i], salign=salign)
            s = s + "{0:{salign}{sform}}".format(summary["params"][i], salign=salign, sform=sform)
            s = s + "{0:{salign}{sform}}".format(summary["stderr"][i], salign=salign, sform=sform)
            s = s + "{0:{salign}{sform}}".format(summary["tscore"][i], salign=salign, sform=sform)
            s = s + "{0:{salign}{sform}}\n".format(summary["pvalue"][i], salign=salign, sform=sform)
        s = s + "\n\n"

    if "v" in sinc:
        s = s + "Parameter Covariance Matrix:\n"
        for i in range(I):
            s = s + "{0:{1}}|".format("", sbuf-1)
            for j in range(I-1):
                s = s + "{0:{salign}{sform}}".format(summary["covar"][i][j], salign=salign, sform=sform)
            s = s + "{0:{salign}{sform}}|\n".format(summary["covar"][i][I-1], salign=salign, sform=sform) 
        s = s + "\n\n"

    if "c" in sinc:
        s = s + "Parameter Correlation Matrix:\n"
        for i in range(I):
            s = s + "{0:{1}}|".format("", sbuf-1)
            for j in range(I-1):
                s = s + "{0:{salign}{sform}}".format(summary["cor"][i][j], salign=salign, sform=sform)
            s = s + "{0:{salign}{sform}}|\n".format(summary["cor"][i][I-1], salign=salign, sform=sform)
        s = s + "\n\n"

    return s 

def SummaryLaTeX(summary, paramlabels=None, sinc="rsvc", sform=".6G"):
    """
    Takes a least-squares summary dictionary and returns its LaTeX string representation.

    Args:
       * `summary` (dict): output from LeastSquaresSummary().

    Keyword Args:
       * `paramlabels` (list): String labels for params, length len(summary["params"]). Default ["p0", "p1", ...].
       * `sinc`        (str):  String include flags for which summary data to print. For example, 'sv' would only 
         print the summary table and the covariance matrix.  Default prints data.  Flags have the following meanings:
            * `r`: Residual summary table.
            * `s`: Regression parameter summary.
            * `v`: Parameter covariance matrix.
            * `c`: Parameter correlation matrix.
       * `sform`       (str):  Data format string specification.

    Returns:
       * `s` (str): LaTeX-valid string.  Compiles to the following for sample y = p0 + p1 x:

         .. math::
            \\begin{tabular}{|c|c|c|}
            \hline
            \\textbf{Min}&\\textbf{Median}&\\textbf{Max}\\\\
            \hline
            -2.41757&0.00465894&2.51009\\\\
            \hline
            \end{tabular}

            \\begin{tabular}{|l|c|c|c|c|}
            \hline
            & \\textbf{Value} & \\textbf{Std. Err.} & \\textbf{t} & \\textbf{P} \\\\
            \hline
            p0 & 10.0298 & 0.0512374 & 195.753 & 0\\\\
            \hline
            p1 & 2.28153 & 0.0887236 & 25.715 & 1.97358E-111\\\\
            \hline
            \end{tabular}

            \Sigma = \left| \\begin{array}{cc}
            0.00399002 & -0.00598204\\\\
            -0.00598204 & 0.0119641\\\\
            \end{array} \\right|

            \\rho = \left| \\begin{array}{cc}
            1 & -0.865809\\\\
            -0.865809 & 1\\\\
            \end{array} \\right|
    """

    #LaTeX formatting shortcuts
    endrow  = "\\\\\n"
    hline   = "\\hline\n"
    I = len(summary["params"])

    if paramlabels == None:
        paramlabels = []
        for i in range(I):
            paramlabels.append("p" + str(i))

    s = ""

    if "r" in sinc:
        s = s + "\\begin{center}\n"
        s = s + "\\begin{table}[htbp]\n"
        s = s + "\\caption{Regression Residuals:}\n"
        s = s + "\\begin{center}\n"
        s = s + "\\begin{tabular}{|c|c|c|}\n"
        s = s + hline
        s = s + "\\textbf{Min}&\\textbf{Median}&\\textbf{Max}" + endrow
        s = s + hline
        s = s + "{0:{sform}}&{1:{sform}}&{2:{sform}}".format(summary["residuals"].min(), 
            np.median(summary["residuals"]), summary["residuals"].max(),sform=sform) + endrow
        s = s + hline
        s = s + "\\end{tabular}\n"
        s = s + "\\end{center}\n"
        s = s + "\\end{table}\n"
        s = s + "\\end{center}\n\n"

    if "s" in sinc:
        s = s + "\\begin{center}\n"
        s = s + "\\begin{table}[htbp]\n"
        s = s + "\\caption{Regression Parameter Summary:}\n"
        s = s + "\\begin{center}\n"
        s = s + "\\begin{tabular}{|l|c|c|c|c|}\n"
        s = s + hline
        s = s + "& \\textbf{Value} & \\textbf{Std. Err.} & \\textbf{t} & \\textbf{P} " + endrow
        s = s + hline
        for i in range(I):
            s = s + "{0}".format(paramlabels[i])
            s = s + " & {0:{sform}}".format(summary["params"][i], sform=sform)
            s = s + " & {0:{sform}}".format(summary["stderr"][i], sform=sform)
            s = s + " & {0:{sform}}".format(summary["tscore"][i], sform=sform)
            s = s + " & {0:{sform}}".format(summary["pvalue"][i], sform=sform) + endrow + hline
        s = s + "\\end{tabular}\n"
        s = s + "\\end{center}\n"
        s = s + "\\end{table}\n"
        s = s + "\\end{center}\n\n"

    if "v" in sinc:
        s = s + "\\[ \\Sigma = \\left| \\begin{array}{"
        for i in range(I):
            s = s + "c"
        s = s + "}\n"
        for i in range(I):
            for j in range(I-1):
                s = s + "{0:{sform}} & ".format(summary["covar"][i][j], sform=sform)
            s = s + "{0:{sform}}".format(summary["covar"][i][I-1], sform=sform) + endrow
        s = s + "\\end{array} \\right| \\]\n\n\n"

    if "c" in sinc:
        s = s + "\\[ \\rho = \\left| \\begin{array}{"
        for i in range(I):
            s = s + "c"
        s = s + "}\n"
        for i in range(I):
            for j in range(I-1):
                s = s + "{0:{sform}} & ".format(summary["cor"][i][j], sform=sform)
            s = s + "{0:{sform}}".format(summary["cor"][i][I-1], sform=sform) + endrow
        s = s + "\\end{array} \\right| \\]\n"

    return s 


from ContingencyTable2D import *
from ContingencyTable3D import *
