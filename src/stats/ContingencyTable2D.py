"""Module containing a 2D contingency table class."""

import numpy as np
import tables as tb

from .. import LinearUniformBins, LogUniformBins, NinesUniformBins
from . import CTScaleError

class ContingencyTable2D(object):
    """Two Dimensional Contingency Table Object
    """

    def __init__(self, x, R, I=2, J=2, xscale="linear", Rscale="linear", 
        xlabel="x", Rlabel="R", title=None, xrange=[], Rrange=[], datarange=[], 
        xdscrt=False, Rdscrt=False, sbuf=40, sform=".6G", extracols={}):
        """Args:
           * `x` (list): Independent variable dataset, length N.
           * `R` (list): Response variable dataset, length N.

        Keyword Arguments:
           * `I`         (int):  Number of response bins. 
           * `J`         (int):  Number of independent variable bins. 
           * `xscale`    (str):  Scale of independent variable binning.
           * `Rscale`    (str):  Scale of response variable binning.
           * `xlabel`    (str):  String label for independent variable.
           * `Rlabel`    (str):  String label for response variable.
           * `title`     (str):  String label for the contingency table itself.
           * `xrange`    (list): Range on which the the independent variable is defined, length 2.
           * `Rrange`    (list): Range on which the the response variable is defined, length 2.
           * `datarange` (list): Concatenation of xrange + Rrange, length 4.
           * `xdscrt`    (bool): Flag for x being a discrete (or continuous) variable.
           * `Rdscrt`    (bool): Flag for R being a discrete (or continuous) variable.
           * `sbuf`      (str):  Column width for string representation of tables.
           * `sform`     (str):  Data format for string representation of tables.
           * `extracols` (dict): Dictionary of extra columns to include in the HDF5 model. See CT2D_TableModel().
        """

        self.x = np.array(x)	#Independent variable (parameter)
        self.R = np.array(R)	#Response (dependnent) variable

        self.N = len(self.x)	#Size of Contingency Table Data
        self.J = int(J)         #Number of Parameter Bins
        self.I = int(I)         #Number of Response Bins

        self.xdiscrete = xdscrt #Parameter a discrete variable?
        self.Rdiscrete = Rdscrt #Response a discrete variable?

        #Set datarange, if not given
        if not (len(datarange) == 4):
            datarange = []

            if len(xrange) == 2:
                datarange.append(xrange[0])
                datarange.append(xrange[1])
            else:
                datarange.append(min(self.x))
                datarange.append(max(self.x))

            if len(Rrange) == 2:
                datarange.append(Rrange[0])
                datarange.append(Rrange[1])
            else:
                datarange.append(min(self.R))
                datarange.append(max(self.R))

        #Set x-scale
        if self.xdiscrete:
            xbinnum = self.J - 1
        else:
            xbinnum = self.J

        if xscale == "linear":
            self.xbounds = LinearUniformBins(datarange[0], datarange[1], xbinnum)
        elif xscale == "log":
            self.xbounds = LogUniformBins(datarange[0], datarange[1], xbinnum)
        elif xscale == "nines":
            self.xbounds = NinesUniformBins(datarange[0], datarange[1], xbinnum)
        else:
            raise CTScaleError(xlabel, xscale)

        #Set R-scale
        if self.Rdiscrete:
            Rbinnum = self.I - 1
        else:
            Rbinnum = self.I

        if Rscale == "linear":
            self.Rbounds = LinearUniformBins(datarange[2], datarange[3], Rbinnum)
        elif Rscale == "log":
            self.Rbounds = LogUniformBins(datarange[2], datarange[3], Rbinnum)
        elif Rscale == "nines":
            self.Rbounds = NinesUniformBins(datarange[2], datarange[3], Rbinnum)
        else:
            raise CTScaleError(Rlabel, Rscale)

        self.N_idot = np.zeros(self.I, dtype=np.int)
        self.N_dotj = np.zeros(self.J, dtype=np.int)

        self.N_ij = np.zeros( (self.I, self.J), dtype=np.int)

        #Now fill the data sets!
        for n in range(self.N):
            i = 0
            j = 0

            in_i = False
            in_j = False

            while (not in_i):
                if (not self.Rdiscrete) and (self.Rbounds[i] <= self.R[n]) and (self.R[n] <= self.Rbounds[i+1]):
                    in_i = True
                elif (self.Rdiscrete) and (self.Rbounds[i] == self.R[n]):
                    in_i = True
                else:
                    i = i + 1

            while (not in_j):
                if (not self.xdiscrete) and (self.xbounds[j] <= self.x[n]) and (self.x[n] <= self.xbounds[j+1]):
                    in_j = True
                elif (self.xdiscrete) and (self.xbounds[j] == self.x[n]):
                    in_j = True
                else:
                    j = j + 1

            self.N_idot[i] = self.N_idot[i] + 1
            self.N_dotj[j] = self.N_dotj[j] + 1
            self.N_ij[i][j] = self.N_ij[i][j] + 1

        #Calculates the Chi-Squared value for the contingency table!
        self.E_ij = np.zeros( (self.I, self.J), dtype=np.float)
        self.chi2 = 0.0

        for i in range(self.I):	
            for j in range(self.J):
                self.E_ij[i][j] = float(self.N_idot[i]) * float(self.N_dotj[j]) / float(self.N)
                #Expected value QA
                if self.E_ij[i][j] == 0.0:
                    self.E_ij[i][j] = 10.0**-300
                #Chi-Squared calculation
                self.chi2 = self.chi2 + ( ((float(self.N_ij[i][j]) - self.E_ij[i][j])**2) / self.E_ij[i][j] )

        #Calculate Cramer's V
        mindim = min(self.I -1, self.J - 1)
        self.V = np.sqrt(self.chi2 / (self.N * mindim) )

        #Calculate the Contingency Coefficient C
        self.C = np.sqrt( self.chi2 / (self.chi2 + self.N) )

        #Calculate the probabilities
        self.p_ij   = (self.N_ij.astype(float)) / float(self.N)
        for i in range(self.I):
            for j in range(self.J):
                #QA on p_ij zero entries
                if self.p_ij[i][j] == 0.0:
                    self.p_ij[i][j] = 10.0**-300

        self.p_idot = (self.N_idot.astype(float)) / float(self.N)
        for i in range(self.I):
            if self.p_idot[i] == 0.0:
                self.p_idot[i] = 10.0**-300
        
        self.p_dotj = (self.N_dotj.astype(float)) / float(self.N)
        for j in range(self.J):
            if self.p_dotj[j] == 0.0:
                self.p_dotj[j] = 10.0**-300


        #Calulate Entropy H
        self.H_x = 0.0
        for j in range(self.J):
            self.H_x = self.H_x - (self.p_dotj[j] * np.log(self.p_dotj[j]))

        self.H_R = 0.0
        for i in range(self.I):
            self.H_R = self.H_R - (self.p_idot[i] * np.log(self.p_idot[i]))

        self.H_Rx = 0.0
        self.H_R_x = 0.0
        self.H_x_R = 0.0
        self.I_Rx = 0.0 	#Mutual Information
        for i in range(self.I):
            for j in range(self.J):
                self.H_Rx  = self.H_Rx  - (self.p_ij[i][j] * np.log(self.p_ij[i][j]))
                self.H_R_x = self.H_R_x - (self.p_ij[i][j] * np.log(self.p_ij[i][j] / self.p_dotj[j]))
                self.H_x_R = self.H_x_R - (self.p_ij[i][j] * np.log(self.p_ij[i][j] / self.p_idot[i]))

                if (10.0**-300 == self.p_ij[i][j] == self.p_idot[i] == self.p_dotj[j]):
                    self.I_Rx = self.I_Rx + 0.0
                else:
                    self.I_Rx = self.I_Rx + (self.p_ij[i][j] * np.log(self.p_ij[i][j] / (self.p_idot[i] * self.p_dotj[j])) )

        #Calculate the R-statistic
        self.R_stat = np.sqrt(1.0 - np.exp(-2.0*self.I_Rx))

        #Normalized Metrics
        self.U_R_x = self.I_Rx / self.H_R
        self.U_x_R = self.I_Rx / self.H_x
        self.U_Rx  = 2.0 * self.I_Rx / (self.H_R + self.H_x)

        self.M1 = self.I_Rx / np.sqrt(self.H_R * self.H_x)
        self.M2 = self.I_Rx / min(self.H_R, self.H_x)
        self.M3 = self.I_Rx / max(self.H_R, self.H_x)
        self.M4 = self.I_Rx / self.H_Rx
        self.M5 = 1.0 - (self.I_Rx / self.H_Rx)
        self.S = 1.0

        #Finally, coppy over some options to the CT itself
        self.sbuf   = int(sbuf)			#String buffer for contingency table printing
        self.salign = "^" + str(self.sbuf)	#String buffer for contingency table printing
        self.sform  = sform			#String formatting for contingency table data
        self.xlabel = xlabel			#Parameter label
        self.Rlabel = Rlabel			#Response label
        self.title  = title			#Table title / Label
        if self.title == None:
            self.title = "{xlabel} to {Rlabel}".format(**self.__dict__)

        return 

    def __call__(self, x, R, **kwargs):
        """Updates the contingency table, works by calling __init__()."""

        self.__dict__.update(kwargs)
        self.__init__(x, R, **self.__dict__)
        return

    def __str__(self):
        """Returns string representation of the contingency table."""

        s = "Contingency Table for {title}:\n".format(**self.__dict__)
        #Top Row
        s = s + format("", self.salign)
        for j in range(self.J):
            if self.xdiscrete:
                s = s + format("{xlabel} = {bounds[0]:{sform}}".format(bounds=self.xbounds[j:j+1], **self.__dict__), self.salign)
            else:
                s = s + format("{bounds[0]:{sform}} < {xlabel} < {bounds[1]:{sform}}".format(bounds=self.xbounds[j:j+2], **self.__dict__), self.salign)
        s = s + format("", self.salign)
        s = s + "\n"
        #Middle, data rows
        for i in range(self.I):
            if self.Rdiscrete:
                s = s + format("{Rlabel} = {bounds[0]:{sform}}".format(bounds=self.Rbounds[i:i+1], **self.__dict__), self.salign)
            else:
                s = s + format("{bounds[0]:{sform}} < {Rlabel} < {bounds[1]:{sform}}".format(bounds=self.Rbounds[i:i+2], **self.__dict__), self.salign)
            for j in range(self.J):
                s = s + "{0:{1}{2}}".format(self.N_ij[i][j], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.N_idot[i], self.salign, self.sform)
            s = s + "\n"
        #Bottom row
        s = s + format("", self.salign)
        for j in range(self.J):
            s = s + "{0:{1}{2}}".format(self.N_dotj[j], self.salign, self.sform)
        s = s + "{0:{1}{2}}".format(self.N, self.salign, self.sform)
        s = s + "\n"
        s = s + "\n"

        s = s + "Chi-Squared               = " + str(self.chi2) + "\n"
        s = s + "Cramer's V                = " + str(self.V) + "\n"
        s = s + "Contingency Coefficient C = " + str(self.C) + "\n"
        s = s + "R-value R                 = " + str(self.R_stat) + "\n"
        s = s + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Probability Table for {title}:\n".format(**self.__dict__)
        #Top Row
        s = s + format("", self.salign)
        for j in range(self.J):
            if self.xdiscrete:
                s = s + format("{xlabel} = {bounds[0]:{sform}}".format(bounds=self.xbounds[j:j+1], **self.__dict__), self.salign)
            else:
                s = s + format("{bounds[0]:{sform}} < {xlabel} < {bounds[1]:{sform}}".format(bounds=self.xbounds[j:j+2], **self.__dict__), self.salign)
        s = s + format("", self.salign)
        s = s + "\n"
        #Middle, data rows
        for i in range(self.I):
            if self.Rdiscrete:
                s = s + format("{Rlabel} = {bounds[0]:{sform}}".format(bounds=self.Rbounds[i:i+1], **self.__dict__), self.salign)
            else:
                s = s + format("{bounds[0]:{sform}} < {Rlabel} < {bounds[1]:{sform}}".format(bounds=self.Rbounds[i:i+2], **self.__dict__), self.salign)
            for j in range(self.J):
                s = s + "{0:{1}{2}}".format(self.p_ij[i][j], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.p_idot[i], self.salign, self.sform)
            s = s + "\n"
        #Bottom row
        s = s + format("", self.salign)
        for j in range(self.J):
            s = s + "{0:{1}{2}}".format(self.p_dotj[j], self.salign, self.sform)
        s = s + "{0:{1}{2}}".format(1.0, self.salign, self.sform)
        s = s + "\n"
        s = s + "\n"

        s = s + "Entropies:\n"
        s = s + "H(R)   = " + str(self.H_R)   + "\n"
        s = s + "H(x)   = " + str(self.H_x)   + "\n"
        s = s + "H(R,x) = " + str(self.H_Rx)  + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Conditional Entropies:\n"
        s = s + "H(R|x) = " + str(self.H_R_x)  + "\n"
        s = s + "H(x|R) = " + str(self.H_x_R)  + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Mutual Information:\n"
        s = s + "I(R,x)   = " + str(self.I_Rx)   + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Nomralized Measures:\n"
        s = s + "U(R|x) = " + str(self.U_R_x) + "\n"
        s = s + "U(x|R) = " + str(self.U_x_R) + "\n"
        s = s + "U(R,x) = " + str(self.U_Rx)  + "\n"
        s = s + "M1     = " + str(self.M1)    + "\n"
        s = s + "M2     = " + str(self.M2)    + "\n"
        s = s + "M3     = " + str(self.M3)    + "\n"
        s = s + "M4     = " + str(self.M4)    + "\n"
        s = s + "M5     = " + str(self.M5)    + "\n"
        s = s + "S      = " + str(self.S)    + "\n"

        
        return s


    def LaTeX(self, sinc="cprm"):
        """Return a LaTeX string represnetation of the contingency table.

        Keyword Args:
           * `sinc` (str):  String include flags that dictates which summary data to print. For example,
             'cm' would only print the contingency table and the normalized metrics. Default prints all data available.
             Flags have the following meanings:
                * `c`: The contingency table itself.
                * `p`: The probability table itself.
                * `r`: Basic contingency table results.
                * `m`: Normalized Measures.
                
        Returns:
           * LaTeX valid string that represnets the contingency table and its results.
        """

        #begin tabular command shortcut!
        #easy way to set the columns
        begin_tabular = "\\begin{tabular}{|c|"
        for j in range(self.J):
            begin_tabular = begin_tabular + "|c"
        begin_tabular = begin_tabular + "||c|}\n"

        #LaTeX formatting shortcuts
        self.endrow  = "\\\\\n"
        self.hline   = "\\hline\n"
        self.columns = self.J + 2

        s = ""

        if 'c' in sinc:
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table for {title}.}}\n".format(**self.__dict__) 
            s = s + "\\begin{center}\n"
            s = s + begin_tabular
            s = s + self.hline

            #Top Row
            for j in range(self.J):
                if self.xdiscrete:
                    s = s + "&{xlabel} $= {bounds[0]:{sform}}$".format(bounds=self.xbounds[j:j+1], **self.__dict__)
                else:
                    s = s + "&${bounds[0]:{sform}} <$ {xlabel} $< {bounds[1]:{sform}}$".format(bounds=self.xbounds[j:j+2], **self.__dict__)
            s = s + "&"
            s = s + self.endrow
            s = s + self.hline
            #Middle, data rows
            for i in range(self.I):
                if self.Rdiscrete:
                    s = s + "{Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                else:
                    s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.N_ij[i][j], self.sform)
                s = s + "&${0:{1}}$".format(self.N_idot[i], self.sform)
                s = s + self.endrow
                s = s + self.hline
            #Bottom row
            for j in range(self.J):
                s = s + "&${0:{1}}$".format(self.N_dotj[j], self.sform)
            s = s + "&${0:{1}}$".format(self.N, self.sform)
            s = s + self.endrow
            s = s + self.hline

            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"

        if 'p' in sinc:
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Probability Table for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + begin_tabular
            s = s + self.hline

            #Top Row
            for j in range(self.J):
                if self.xdiscrete:
                    s = s + "&{xlabel} $= {bounds[0]:{sform}}$".format(bounds=self.xbounds[j:j+1], **self.__dict__)
                else:
                    s = s + "&${bounds[0]:{sform}} <$ {xlabel} $< {bounds[1]:{sform}}$".format(bounds=self.xbounds[j:j+2], **self.__dict__)
            s = s + "&"
            s = s + self.endrow
            s = s + self.hline
            #Middle, data rows
            for i in range(self.I):
                if self.Rdiscrete:
                    s = s + "{Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                else:
                    s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.p_ij[i][j], self.sform)
                s = s + "&${0:{1}}$".format(self.p_idot[i], self.sform)
                s = s + self.endrow
                s = s + self.hline
            #Bottom row
            for j in range(self.J):
                s = s + "&${0:{1}}$".format(self.p_dotj[j], self.sform)
            s = s + "&${0:{1}}$".format(1.0, self.sform)
            s = s + self.endrow
            s = s + self.hline

            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"

        if 'r' in sinc:
            #Contingency Table Results
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Results for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "$\\chi^2$&${chi2:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "Cramer\'s V&${V:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "Contingency Coefficient C&${C:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "R-value Statistic&${R_stat:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R)&${H_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x)&${H_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,x)&${H_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R$|$x)&${H_R_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x$|$R)&${H_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(R,x)&${I_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"

        if 'm' in sinc:
            #Normalized Metrics
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Normalized Measures for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "U(R$|$x)&${U_R_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(x$|$R)&${U_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,x)&${U_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M1: $\\frac{{I(R,x)}}{{\\sqrt{{H(R) H(x)}}}}$&${M1:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M2: $\\frac{{I(R,x)}}{{\\min(H(R),H(x))}}$&${M2:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M3: $\\frac{{I(R,x)}}{{\\max(H(R),H(x))}}$&${M3:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M4: $\\frac{{I(R,x)}}{{H(R,x))}}$&${M4:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M5: $1 - \\frac{{I(R,x)}}{{H(R,x))}}$&${M5:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "S: $\\frac{{U(R|x)}}{{U(R|x))}}$&${S:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"

        return s

    def FillTableRow(self, row, extracols={}):
        """Fills an HDF5 table row with this contingency table's data and results."""

        if self.H5table == None:
            #create a copy of an HDF5 table model for this CT
            self.H5table = CT2D_TableModel(self.I, self.J, extracols, self.xdiscrete, self.Rdiscrete)

        for key in self.H5table:
            row[key] = getattr(self, key)
        row.append()
        return

    #Document Strings for class attributes
    x = None
    """Independent variable (parameter) dataset."""
    R = None
    """Response (dependent) variable dataset."""

    N = None
    """Number of data points in contingency table."""
    I = None
    """Number of R Response bins."""
    J = None
    """Number of x independent variable bins."""


    Rdiscrete = None
    """Flag (bool) for if the response is a discrete (or continuous) variable."""
    xdiscrete = None
    """Flag (bool) for if the indpendent parameter is a discrete (or continuous) variable."""


    datarange = None
    """List of ranges that the table data and bins are defined on.
       datarange = [x-lower, x-upper,  R-lower, R-upper]
    """

    xscale = None
    """Scale of the independent variable x. The following are valid flags:
       * 'linear': Linearly uniform bins.
       * 'log':    Log-uniform bins. 
       * 'nines':  One-minus-log-uniform bins.
    """
    Rscale = None
    """Scale of the response variable R. The following are valid flags:
       * 'linear': Linearly uniform bins.
       * 'log':    Log-uniform bins. 
       * 'nines':  One-minus-log-uniform bins.
    """

    xbounds = None
    """List of independent variable x bounds, length J+1.  Denoted :math:`b^x`."""
    Rbounds = None
    """List of response variable R bounds, length I+1.  Denoted :math:`b^x`."""
    

    N_idot = None
    """List of sums over the jth index. Array of length I.

    .. math::
       N_{i\cdot} = \sum_{j} N_{ij}
    """

    N_dotj = None
    """List of sums over the ith index. Array of length J.

    .. math::
       N_{\cdot j} = \sum_{i} N_{ij}
    """

    N_ij = None
    """Raw Contingency Table (Matrix) of size IxJ.

    .. math::
       N_{ij} = \sum_n^N \Pi_{b_i^R, b_{i+1}^R}(R_n) \cdot \Pi_{b_j^x, b_{j+1}^x}(x_n)

    Here, :math:`\Pi_{a,b}(x)` is the `Boxcar Function <http://mathworld.wolfram.com/BoxcarFunction.html>`_
    and :math:`b^x` are the bounds in x (self.xbounds), etc.
    """

    E_ij = None
    """Expected Frequency Matrix of size IxJ.

    .. math::
       E_{ij} = \\frac{N_{i\cdot} \cdot N_{\cdot j}}{N}	
    """

    chi2 = None
    """Chi-Squared Statistic.

    .. math::
       \chi^2 =  \sum_{ij} \\frac{\left(N_{ij} - E_{ij}\\right)^2}{E_{ij}}
    """

    V = None
    """Cramer's V Statistic. Calculated with Chi-Squared.

    .. math::
       V = \sqrt{\\frac{\chi^2}{N \cdot (\min[I, J] - 1)}}
    """

    C = None
    """Contingency Coefficient Statistic.

    .. math::
       C = \sqrt{\\frac{\chi^2}{\chi^2 + N}}
    """ 

    #Probabilities
    p_idot = None
    """List of ith bin probabilities.

    .. math::
       p_{i\cdot} = \\frac{N_{i\cdot}}{N}	
    """

    p_dotj = None
    """List of jth bin probabilities.

    .. math::
       p_{\cdot j} = \\frac{N_{\cdot j}}{N}
    """

    p_ij = None
    """Raw Probability Table (Matrix) of size IxJ.

    .. math::
       p_{ij} = \\frac{N_{ij}}{N}	
    """

    #Entropies
    H_x = None
    """Entropy H(x) of independent variable x.

    .. math::
       H(x) = - \sum_j p_{\cdot j} \ln(p_{\cdot j})	
    """
    H_R = None
    """Entropy H(R) of response variable R.

    .. math::
       H(R) = - \sum_i p_{i\cdot} \ln(p_{i\cdot})	
    """

    H_Rx = None
    """Joint Entropy H(R,x) of response variable R and independent variable x.

    .. math::
       H(R,x) = - \sum_{ij} p_{ij} \ln(p_{ij})	
    """

    #Conditional Entropies
    H_R_x  = None
    """Conditional Entropy H(R|x).

    .. math::
       H(R|x) = - \sum_{i,j} p_{ij} \ln\left(\\frac{p_{ij}}{p_{\cdot j}}\\right)
    """

    H_x_R  = None
    """Conditional Entropy H(x|R).

    .. math::
       H(x|R) = - \sum_{i,j} p_{ij} \ln\left(\\frac{p_{ij}}{p_{i\cdot}}\\right)
    """

    #Mutual Information
    I_Rx = None
    """Mutual Information I(R,x).

    .. math::
       I(R,x) = - \sum_{i,j} p_{ij} \ln\left(\\frac{p_{ij}}{p_{i\cdot} \cdot p_{\cdot j}}\\right)
    """

    #Calculate the R-statistic
    R_stat = None
    """R-Value Statistic.

    .. math::
       r = \sqrt{1 - e^{-2 I(R,x)}}
    """

    #Normalized Metrics
    U_R_x = None
    """Uncertainty coefficient (coefficient of constraint) U(R|x).  Calculated from,

    .. math::

       U(R|x) = \\frac{I(R,x)}{H(R)}	   
    """

    U_x_R = None
    """Uncertainty coefficient (coefficient of constraint) U(x|R).  Calculated from,

    .. math::

       U(x|R) = \\frac{I(R,x)}{H(x)}	   
    """
    U_Rx  = None
    """Symetric Uncertainty U(R,x).  Calculated from,

    .. math::

       U(R,x) = 2\\frac{I(R,x)}{H(R) + H(x)}	   
    """
    
    M1 = None
    """Normalized Measure 1, from Strehl, Alexander; Joydeep Ghosh (2002). *"Cluster ensembles -- a knowledge reuse 
    framework for combining multiple partitions"*. Journal of Machine Learning Research 3: 583-617. 
    doi:10.1162/153244303321897735.  Equation (2):

    .. math::

       M_1 = \\frac{I(R,x)}{\sqrt{H(R) H(x)}}
    """
    M2 = None
    """Normalized Measure 2, from Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data mining, 
    in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.34),

    .. math::

       M_2 = \\frac{I(R,x)}{\min[H(R), H(x)]}
    """
    M3 = None
    """Normalized Measure 3, from Kvalseth, T.O. *"Entropy and correlation: some comments"*, IEEE Transactions On Systems, Man, 
    and Cybernetics, SMC-17, 517-519, 1987.  

    .. math::

       M_3 = \\frac{I(R,x)}{\max[H(R), H(x)]}
    """
    M4 = None
    """Normalized Measure 4, from Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data mining, 
    in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.35),

    .. math::

       M_4 = \\frac{I(R,x)}{H(R,x)}
    """
    M5 = None
    """Normalized Measure 5, from Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data mining, 
    in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.35),

    .. math::

       M_5 = 1 - \\frac{I(R,x)}{H(R,x)}
    """

    S = None
    """Normalized Measure S, from the minds of E. A. Schnider and A. M. Scopatz.  Measures whether the response is determined
    by only one of the indpendent variables or if the dependence is based on many indepenent parameters.  The interpretation
    of this normalized metric is as follows:
       * S = 1: Response soley determined by one variable.  
       * S < 0: Respose by combination of variables, though one independent parameter shows a weaker dependence.
       * S = 0: Response is perfectly determined by all variables in tandem.

    Note that this metric is trivial in the case of one indepenent parameter.

    .. math::

       S = \\frac{U(R|x)}{U(R|x)} = 1
    """

    #Finally, copy over some options to the CT itself
    sbuf = None
    """String buffer for contingency table printing; integer number spaces per column, default = 40."""
    salign = None
    """String alignment and buffer for contingency table printing; default = '^sbuf' [centered]."""
    sform = None
    """String formatting for contingency table data, default = '.6G'."""
    xlabel = None
    """Independent variable label, default = 'x'."""
    Rlabel = None
    """Response variable label, default = 'R'."""
    title = None
    """Contingency table label, default = '{xlabel} to {Rlabel}'."""

    #create a copy of an HDF5 table model for this CT
    H5table = None
    """HDF5 table model for this contingency table and its results."""





#2D Contingency Table Data Model
def CT2D_TableModel(I, J, extracols={}, xdscrt=False, Rdscrt=False): 
    """Two-Dimensional Contingency Table Model for HDF5

    Args:
        * `I` (int): Number of response parameter bins.
        * `J` (int): Number of independent parameter bins.

    Keyword Args:
        * `extracols` (dict): Dictionary of extra columns that 
          should be included in the table.  The keys are the 
          name of the extra columns while the values specify the
          type of column in the usual PyTables way.  For example, 
          You may wish to tag contingeny tables as time series data::

            extracols['time'] = tables.Time64Col(pos=-1)

          The standard colums range from positions 0 to 27.  
          To prepend data, start from pos=-1 and work down.
          To append columns, start from pos=28 and work up!

          Note that to have ContingencyTable2D.FillTableRow(row) automatically
          fill in any extra column data you'll need to set the approriate data 
          as an attribute of the Contingency Table itself first::

            ct = ContingencyTable2D(...)
            setattr(ct, 'time', 1.0)
            ct.FillTableRow(mytablerow)
           
        * `xdscrt` (bool): Flag for x being a discrete (or continuous) variable.
        * `Rdscrt` (bool): Flag for x being a discrete (or continuous) variable.
    """
    CT2D_TM = {
        'xlabel': tb.StringCol(50, pos=0),
        'Rlabel': tb.StringCol(50, pos=1),

        'H_R':    tb.Float64Col(pos=2),
        'H_x':    tb.Float64Col(pos=3),

        'H_Rx':   tb.Float64Col(pos=4),
        'H_R_x':  tb.Float64Col(pos=5),
        'H_x_R':  tb.Float64Col(pos=6),

        'I_Rx':   tb.Float64Col(pos=7),

        'U_Rx':   tb.Float64Col(pos=8),
        'U_R_x':  tb.Float64Col(pos=9),
        'U_x_R':  tb.Float64Col(pos=10),

        'M1':     tb.Float64Col(pos=11),
        'M2':     tb.Float64Col(pos=12),
        'M3':     tb.Float64Col(pos=13),
        'M4':     tb.Float64Col(pos=14),
        'M5':     tb.Float64Col(pos=15),
        'S':      tb.Float64Col(pos=16),

        'chi2':   tb.Float64Col(pos=17),
        'C':      tb.Float64Col(pos=18),
        'V':      tb.Float64Col(pos=19),
        'R_stat': tb.Float64Col(pos=20),

        'N':      tb.Int32Col(pos=21),

        'N_dotj': tb.Int32Col(pos=22, shape = J),
        'N_idot': tb.Int32Col(pos=23, shape = I),
        'N_ij':   tb.Int32Col(pos=24, shape = (I, J)),
        'E_ij':   tb.Float64Col(pos=25, shape = (I, J)),
        }

    if Rdscrt:
        CT2D_TM["Rbounds"] = tb.Float64Col(pos=26, shape = I)
    else:
        CT2D_TM["Rbounds"] = tb.Float64Col(pos=26, shape = I+1)

    if xdscrt:
        CT2D_TM["xbounds"] = tb.Float64Col(pos=27, shape = J)
    else:
        CT2D_TM["xbounds"] = tb.Float64Col(pos=27, shape = J+1)

    CT2D_TM.update(extracols)

    return CT2D_TM
