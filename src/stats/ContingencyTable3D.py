"""Module containing a 3D contingency table class."""

import numpy as np 
import tables as tb

from .  import CTScaleError
from .. import LinearUniformBins, LogUniformBins, NinesUniformBins

class ContingencyTable3D(object):
    "Three Dimensional Contingency Table Object"

    def __init__(self, x, y, R, I=2, J=2, K=2, xscale="linear", yscale="linear", Rscale="linear", 
        xlabel="x", ylabel="y", title=None, Rlabel="R", xrange=[], yrange=[], Rrange=[], datarange=[], 
        xdscrt=False, ydscrt=False, Rdscrt=False, sbuf=40, sform=".6G", extracols={}):
        """Args:
           * `x` (list): First independent variable dataset, length N.
           * `y` (list): Second independent variable dataset, length N.
           * `R` (list): Response variable dataset, length N.

        Keyword Args:
           * `I`         (int):  Number of response bins. 
           * `J`         (int):  Number of first independent variable bins. 
           * `K`         (int):  Number of second independent variable bins. 
           * `xscale`    (str):  Scale of first independent variable binning.
           * `yscale`    (str):  Scale of second independent variable binning.
           * `Rscale`    (str):  Scale of response variable binning.
           * `xlabel`    (str):  Label for first independent variable.
           * `ylabel`    (str):  Label for second independent variable.
           * `Rlabel`    (str):  Label for response variable.
           * `title`     (str):  Label for contingency table itself.
           * `xrange`    (list): Range on which the the first independent variable is defined, length 2.
           * `yrange`    (list): Range on which the the second independent variable is defined, length 2.
           * `Rrange`    (list): Range on which the the response variable is defined, length 2.
           * `datarange` (list): Concatenation of xrange + yrange + Rrange, length 6.
           * `xdscrt`    (bool): Flag for x being a discrete (or continuous) variable.
           * `ydscrt`    (bool): Flag for y being a discrete (or continuous) variable.
           * `Rdscrt`    (bool): Flag for R being a discrete (or continuous) variable.
           * `sbuf`      (str):  Column width for string representation of tables.
           * `sform`     (str):  Data format for string representation of tables.
           * `extracols` (dict): Dictionary of extra columns to include in the HDF5 model. See CT3D_TableModel().
        """

        self.x = np.array(x)    #First Independent variable (parameter)
        self.y = np.array(y)    #Second Independent variable (parameter)
        self.R = np.array(R)    #Response (dependent) variable

        self.N = len(self.x)    #Size of Contingency Table Data
        self.I = int(I)         #Number of Response Bins
        self.J = int(J)         #Number of Parameter Bins
        self.K = int(K)         #Number of Parameter Bins

        self.xdiscrete = xdscrt #First Parameter a discrete variable?
        self.ydiscrete = ydscrt #Seconf Parameter a discrete variable?
        self.Rdiscrete = Rdscrt #Response a discrete variable?


        #Set datarange, if not given
        if not (len(datarange) == 6):
            datarange = []

            if len(xrange) == 2:
                datarange.append(xrange[0])
                datarange.append(xrange[1])
            else:
                datarange.append(min(self.x))
                datarange.append(max(self.x))

            if len(yrange) == 2:
                datarange.append(yrange[0])
                datarange.append(yrange[1])
            else:
                datarange.append(min(self.y))
                datarange.append(max(self.y))

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
        elif xscale =="log":
            self.xbounds = LogUniformBins(datarange[0], datarange[1], xbinnum)
        elif xscale == "nines":
            self.xbounds = NinesUniformBins(datarange[0], datarange[1], xbinnum)
        else:
            raise CTScaleError(xlabel, xscale)

        #Set y-scale
        if self.ydiscrete:
            ybinnum = self.K - 1
        else:
            ybinnum = self.K

        if yscale == "linear":
            self.ybounds = LinearUniformBins(datarange[2], datarange[3], ybinnum)
        elif yscale == "log":
            self.ybounds = LogUniformBins(datarange[2], datarange[3], ybinnum)
        elif yscale == "nines":
            self.ybounds = NinesUniformBins(datarange[2], datarange[3], ybinnum)
        else:
            raise CTScaleError(ylabel, yscale)

        #Set R-scale
        if self.Rdiscrete:
            Rbinnum = self.I - 1
        else:
            Rbinnum = self.I

        if Rscale == "linear":
            self.Rbounds = LinearUniformBins(datarange[4], datarange[5], Rbinnum)
        elif Rscale == "log":
            self.Rbounds = LogUniformBins(datarange[4], datarange[5], Rbinnum)
        elif Rscale == "nines":
            self.Rbounds = NinesUniformBins(datarange[4], datarange[5], Rbinnum)
        else:
            raise CTScaleError(Rlabel, Rscale)

        self.N_idotdot = np.zeros(self.I, dtype=np.int)
        self.N_dotjdot = np.zeros(self.J, dtype=np.int)
        self.N_dotdotk = np.zeros(self.K, dtype=np.int)

        self.N_ijdot = np.zeros((self.I, self.J), dtype=np.int)
        self.N_idotk = np.zeros((self.I, self.K), dtype=np.int)
        self.N_dotjk = np.zeros((self.J, self.K), dtype=np.int)

        self.N_ijk = np.zeros((self.I, self.J, self.K), dtype=np.int)

        #Now fill the data sets!
        for n in range(self.N):
            i = 0
            j = 0
            k = 0

            in_i = False
            in_j = False
            in_k = False

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

            while (not in_k):
                if (not self.ydiscrete) and (self.ybounds[k] <= self.y[n] <= self.ybounds[k+1]):
                    in_k = True
                elif (self.ydiscrete) and (self.ybounds[k] == self.y[n]):
                    in_k = True
                else:
                    k = k + 1

            self.N_idotdot[i] = self.N_idotdot[i] + 1
            self.N_dotjdot[j] = self.N_dotjdot[j] + 1
            self.N_dotdotk[k] = self.N_dotdotk[k] + 1

            self.N_ijdot[i][j] = self.N_ijdot[i][j] + 1
            self.N_idotk[i][k] = self.N_idotk[i][k] + 1
            self.N_dotjk[j][k] = self.N_dotjk[j][k] + 1

            self.N_ijk[i][j][k] = self.N_ijk[i][j][k] + 1

        #Calculates the Chi-Squared & G-Statistic value for the contingency table!
        self.E_ijk = np.zeros((self.I, self.J, self.K), dtype=np.float)
        for i in range(self.I):	
            for j in range(self.J):
                for k in range(self.K):
                    self.E_ijk[i][j][k] = float(self.N_idotdot[i]) * float(self.N_dotjdot[j]) * float(self.N_dotdotk[k])/ (float(self.N)**2)

        self.chi2 = ( ((self.N_ijk.astype(float) - self.E_ijk)**2) / self.E_ijk ).sum()
        self.G = (2.0 * self.N_ijk.astype(float) * np.log(self.N_ijk.astype(float) / self.E_ijk)).sum()

        #Calculate Cramer's V
        mindim = min(self.I - 1, self.J - 1, self.K - 1)
        self.V = np.sqrt(self.G / (self.N * mindim) )

        #Calculate the Contingency Coefficient C
        self.C = np.sqrt( self.G / (self.G + self.N) )

        #Calculate Probabilities
        self.p_ijk = (self.N_ijk.astype(float)) / float(self.N)
        for i in range(self.I):
            for j in range(self.J):
                for k in range(self.K):
                    #QA on p_ij zero entries
                    if self.p_ijk[i][j][k] == 0.0:
                        self.p_ijk[i][j][k] = 10.0**-300

        self.p_idotdot = (self.N_idotdot.astype(float)) / float(self.N)
        self.p_dotjdot = (self.N_dotjdot.astype(float)) / float(self.N)
        self.p_dotdotk = (self.N_dotdotk.astype(float)) / float(self.N)

        self.p_ijdot = (self.N_ijdot.astype(float)) / float(self.N)
        self.p_idotk = (self.N_idotk.astype(float)) / float(self.N)
        self.p_dotjk = (self.N_dotjk.astype(float)) / float(self.N)

        for i in range(self.I):
            if self.p_idotdot[i] == 0.0:
                self.p_idotdot[i] = 10.0**-300
            for j in range(self.J):
                if self.p_ijdot[i][j] == 0.0:
                    self.p_ijdot[i][j] = 10.0**-300
        
        for j in range(self.J):
            if self.p_dotjdot[j] == 0.0:
                self.p_dotjdot[j] = 10.0**-300
            for k in range(self.K):
                if self.p_dotjk[j][k] == 0.0:
                    self.p_dotjk[j][k] = 10.0**-300

        for k in range(self.K):
            if self.p_dotdotk[k] == 0.0:
                self.p_dotdotk[k] = 10.0**-300
            for i in range(self.I):
                if self.p_idotk[i][k] == 0.0:
                    self.p_idotk[i][k] = 10.0**-300

        #Calculate Entropy H
        self.H_x = -1.0 * (self.p_dotjdot * np.log(self.p_dotjdot)).sum()
        self.H_y = -1.0 * (self.p_dotdotk * np.log(self.p_dotdotk)).sum()
        self.H_R = -1.0 * (self.p_idotdot * np.log(self.p_idotdot)).sum()

        self.H_Rx = -1.0 * (self.p_ijdot * np.log(self.p_ijdot)).sum()
        self.H_Ry = -1.0 * (self.p_idotk * np.log(self.p_idotk)).sum()
        self.H_xy = -1.0 * (self.p_dotjk * np.log(self.p_dotjk)).sum()

        self.H_Rxy = -1.0 * (self.p_ijk * np.log(self.p_ijk)).sum()

        #Conditional Entropies and Mutual Information
        self.H_R_x = 0.0 
        self.H_x_R = 0.0 
        self.I_Rx = 0.0
        for i in range(self.I):
            for j in range(self.J):
                self.H_R_x = self.H_R_x - (self.p_ijdot[i][j] * np.log(self.p_ijdot[i][j] / self.p_dotjdot[j]))
                self.H_x_R = self.H_x_R - (self.p_ijdot[i][j] * np.log(self.p_ijdot[i][j] / self.p_idotdot[i]))
                self.I_Rx = self.I_Rx + (self.p_ijdot[i][j] * np.log(self.p_ijdot[i][j] / (self.p_idotdot[i] * self.p_dotjdot[j])))

        self.H_R_y = 0.0 
        self.H_y_R = 0.0 
        self.I_Ry = 0.0
        for i in range(self.I):
            for k in range(self.K):
                self.H_R_y = self.H_R_y - (self.p_idotk[i][k] * np.log(self.p_idotk[i][k] / self.p_dotdotk[k]))
                self.H_y_R = self.H_y_R - (self.p_idotk[i][k] * np.log(self.p_idotk[i][k] / self.p_idotdot[i]))
                self.I_Ry = self.I_Ry + (self.p_idotk[i][k] * np.log(self.p_idotk[i][k] / (self.p_idotdot[i] * self.p_dotdotk[k])))

        self.H_x_y = 0.0
        self.H_y_x = 0.0 
        self.I_xy = 0.0
        for j in range(self.J):
            for k in range(self.K):
                self.H_x_y = self.H_x_y - (self.p_dotjk[j][k] * np.log(self.p_dotjk[j][k] / self.p_dotdotk[k]))
                self.H_y_x = self.H_y_x - (self.p_dotjk[j][k] * np.log(self.p_dotjk[j][k] / self.p_dotjdot[j]))
                self.I_xy = self.I_xy + (self.p_dotjk[j][k] * np.log(self.p_dotjk[j][k] / (self.p_dotjdot[j] * self.p_dotdotk[k])))

        self.H_R_xy = 0.0 
        self.H_x_Ry = 0.0
        self.H_y_Rx = 0.0

        self.H_Rx_y = 0.0 
        self.H_Ry_x = 0.0
        self.H_xy_R = 0.0

        self.I_Rxy  = 0.0
        self.I_Rx_y = 0.0
        self.I_Ry_x = 0.0
        self.I_xy_R = 0.0
        for i in range(self.I):
            for j in range(self.J):
                for k in range(self.K):
                    self.H_R_xy = self.H_R_xy - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_dotjk[j][k]))
                    self.H_x_Ry = self.H_x_Ry - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_idotk[i][k]))
                    self.H_y_Rx = self.H_y_Rx - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_ijdot[i][j]))

                    self.H_Rx_y = self.H_Rx_y - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_dotdotk[k]))
                    self.H_Ry_x = self.H_Ry_x - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_dotjdot[j]))
                    self.H_xy_R = self.H_xy_R - (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / self.p_idotdot[i]))

                    self.I_Rxy  = self.I_Rxy  + (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / (self.p_idotdot[i] * self.p_dotjdot[j] * self.p_dotdotk[k])))
                    self.I_Rx_y = self.I_Rx_y + (self.p_ijk[i][j][k] * np.log(self.p_dotdotk[k] * self.p_ijk[i][j][k] / (self.p_idotk[i][k] * self.p_dotjk[j][k])))
                    self.I_Ry_x = self.I_Ry_x + (self.p_ijk[i][j][k] * np.log(self.p_dotjdot[j] * self.p_ijk[i][j][k] / (self.p_ijdot[i][j] * self.p_dotjk[j][k])))
                    self.I_xy_R = self.I_xy_R + (self.p_ijk[i][j][k] * np.log(self.p_idotdot[i] * self.p_ijk[i][j][k] / (self.p_ijdot[i][j] * self.p_idotk[i][k])))

        #Calculate the R-statistic
        self.R_stat = np.sqrt(1.0 - np.exp(-2.0*self.I_Rxy))

        #Normalized Metrics
        self.U_R_x  = self.I_Rx / self.H_R
        self.U_x_R  = self.I_Rx / self.H_x

        self.U_R_y  = self.I_Ry / self.H_R
        self.U_y_R  = self.I_Ry / self.H_y

        self.U_x_y  = self.I_xy / self.H_x
        self.U_y_x  = self.I_xy / self.H_y

        self.U_Rx_y = self.I_Rxy / self.H_Rx
        self.U_Ry_x = self.I_Rxy / self.H_Ry
        self.U_xy_R = self.I_Rxy / self.H_xy

        self.U_R_xy = self.I_Rxy / self.H_R
        self.U_y_Rx = self.I_Rxy / self.H_y
        self.U_x_Ry = self.I_Rxy / self.H_x

        self.U_Rx   = 2.0 * self.I_Rx  / (self.H_R + self.H_x)
        self.U_Ry   = 2.0 * self.I_Ry  / (self.H_R + self.H_y)
        self.U_xy   = 2.0 * self.I_xy  / (self.H_x + self.H_y)
        self.U_Rxy  = 3.0 * self.I_Rxy / (self.H_R + self.H_x + self.H_y)

        self.M1 = self.I_Rxy / ((self.H_R * self.H_x * self.H_y)**(1.0/3.0))
        self.M2 = self.I_Rxy / min(self.H_R, self.H_x, self.H_y)
        self.M3 = self.I_Rxy / max(self.H_R, self.H_x, self.H_y)
        self.M4 = self.I_Rxy / self.H_Rxy
        self.M5 = 1.0 - (self.I_Rxy / self.H_Rxy)
        self.S  = (self.U_x_R + self.U_y_R) / (2.0 * self.U_xy_R)

        # Try Jun's method
        # Cut along y
        U_x_R_y = []
        for k in range(self.K):
            H_x_k = -1.0 * (self.p_dotjk[k] * np.log(self.p_dotjk[k])).sum()
            I_Rx_k = 0.0
            for i in range(self.I):
                for j in range(self.J):
                    I_Rx_k = I_Rx_k + (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / (self.p_idotk[i][k] * self.p_dotjk[j][k])))
            U_x_R_k  = I_Rx_k / H_x_k
            U_x_R_y.append(U_x_R_k)

        self.mu_U_x_R = np.mean(U_x_R_y)
        self.sigma_U_x_R = np.std(U_x_R_y)
        self.se_U_x_R = self.sigma_U_x_R / np.sqrt(self.K)
        self.cv_U_x_R = self.sigma_U_x_R / self.mu_U_x_R

        # Cut along x
        U_y_R_x = []
        for j in range(self.J):
            H_y_j = -1.0 * (self.p_dotjk[j] * np.log(self.p_dotjk[j])).sum()
            I_Ry_j = 0.0
            for i in range(self.I):
                for k in range(self.K):
                    I_Ry_j = I_Ry_j + (self.p_ijk[i][j][k] * np.log(self.p_ijk[i][j][k] / (self.p_ijdot[i][j] * self.p_dotjk[j][k])))
            U_y_R_j  = I_Ry_j / H_y_j
            U_y_R_x.append(U_y_R_j)

        self.mu_U_y_R = np.mean(U_y_R_x)
        self.sigma_U_y_R = np.std(U_y_R_x)
        self.se_U_y_R = self.sigma_U_y_R / np.sqrt(self.J)
        self.cv_U_y_R = self.sigma_U_y_R / self.mu_U_y_R

        # Make Jun's Method Symmetric
        self.cv_U_x_y_R = 0.5 * (self.cv_U_x_R + self.cv_U_y_R)

        # Try Entropy metric as sensitivity of sensitivity
        self.eta_R_x_y = np.zeros((self.J, self.K))
        for j in range(self.J):
            for k in range(self.K):
                H_R = -1.0 * (self.p_ijk[:,j,k] * np.log(self.p_ijk[:,j,k])).sum()
                eta_R = H_R / np.log(self.I)
                self.eta_R_x_y[i][j] = eta_R

        self.mu_eta_R_x_y = self.eta_R_x_y.mean()
        self.sigma_eta_R_x_y = self.eta_R_x_y.std()
        self.se_eta_R_x_y = self.sigma_eta_R_x_y / np.sqrt(self.J * self.K)
        self.cv_eta_R_x_y = self.sigma_eta_R_x_y / self.mu_eta_R_x_y

        #Finally, copy over some options to the CT itself
        self.sbuf   = int(sbuf) 		#String buffer for contingency table printing
        self.salign = "^" + str(self.sbuf)	#String buffer for contingency table printing
        self.sform  = sform 			#String formatting for contingency table data
        self.xlabel = xlabel			#Parameter label
        self.ylabel = ylabel			#Parameter label
        self.Rlabel = Rlabel			#Response label
        self.title  = title			#table label
        if self.title == None:
            self.title = "{xlabel}, {ylabel} to {Rlabel}".format(**self.__dict__)

        return

    def __call__(self, x, y, R, **kwargs):
        "Updates the contingency table, works by calling __init__()."
        self.__dict__.update(kwargs)
        self.__init__(x, y, R, **self.__dict__)
        return

    def __str__(self):
        """Returns string representation of the contingency table."""

        s = "Contingency Table for {title}:\n".format(**self.__dict__)
        for k in range(self.K):
            if self.ydiscrete:
                s = s + format("Slice for {ylabel} = {bounds[0]:{sform}}:".format(bounds=self.ybounds[k:k+1], **self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
            else:
                s = s + format("Slice for {bounds[0]:{sform}} < {ylabel} < {bounds[1]:{sform}}:".format(bounds=self.ybounds[k:k+2], **self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
            s = s + "\n"
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
                    s = s + "{0:{1}{2}}".format(self.N_ijk[i][j][k], self.salign, self.sform)
                s = s + "{0:{1}{2}}".format(self.N_idotk[i][k], self.salign, self.sform)
                s = s + "\n"
            #Bottom row
            s = s + format("", self.salign)
            for j in range(self.J):
                s = s + "{0:{1}{2}}".format(self.N_dotjk[j][k], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.N_dotdotk[k], self.salign, self.sform)
            s = s + "\n"
            s = s + "\n"

        s = s + format("Summary table for sum over {ylabel}:".format(**self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
        s = s + "\n"
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
                s = s + "{0:{1}{2}}".format(self.N_ijdot[i][j], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.N_idotdot[i], self.salign, self.sform)
            s = s + "\n"
        #Bottom row
        s = s + format("", self.salign)
        for j in range(self.J):
            s = s + "{0:{1}{2}}".format(self.N_dotjdot[j], self.salign, self.sform)
        s = s + "{0:{1}{2}}".format(self.N, self.salign, self.sform)
        s = s + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Chi-Squared (poorly defined)               = " + str(self.chi2) + "\n"
        s = s + "G-test                                     = " + str(self.G) + "\n"	
        s = s + "Cramer's V                                 = " + str(self.V) + "\n"
        s = s + "Contingency Coefficient C (poorly defined) = " + str(self.C) + "\n"
        s = s + "R-value R                                  = " + str(self.R_stat) + "\n"
        s = s + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Probability Table for {title}:\n".format(**self.__dict__)
        for k in range(self.K):
            if self.ydiscrete:
                s = s + format("Slice for {ylabel} = {bounds[0]:{sform}}:".format(bounds=self.ybounds[k:k+1], **self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
            else:
                s = s + format("Slice for {bounds[0]:{sform}} < {ylabel} < {bounds[1]:{sform}}:".format(bounds=self.ybounds[k:k+2], **self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
            s = s + "\n"
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
                    s = s + "{0:{1}{2}}".format(self.p_ijk[i][j][k], self.salign, self.sform)
                s = s + "{0:{1}{2}}".format(self.p_idotk[i][k], self.salign, self.sform)
                s = s + "\n"
            #Bottom row
            s = s + format("", self.salign)
            for j in range(self.J):
                s = s + "{0:{1}{2}}".format(self.p_dotjk[j][k], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.p_dotdotk[k], self.salign, self.sform)
            s = s + "\n"
            s = s + "\n"

        s = s + format("Probability summary table for sum over {ylabel}:".format(**self.__dict__), "^" + str(self.sbuf*(self.J+2)) )
        s = s + "\n"
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
                s = s + "{0:{1}{2}}".format(self.p_ijdot[i][j], self.salign, self.sform)
            s = s + "{0:{1}{2}}".format(self.p_idotdot[i], self.salign, self.sform)
            s = s + "\n"
        #Bottom row
        s = s + format("", self.salign)
        for j in range(self.J):
            s = s + "{0:{1}{2}}".format(self.p_dotjdot[j], self.salign, self.sform)
        s = s + "{0:{1}{2}}".format(1.0, self.salign, self.sform)
        s = s + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Entropies:\n"
        s = s + "H(R)     = " + str(self.H_R)   + "\n"
        s = s + "H(x)     = " + str(self.H_x)   + "\n"
        s = s + "H(y)     = " + str(self.H_y)   + "\n"
        s = s + "H(R,x)   = " + str(self.H_Rx)  + "\n"
        s = s + "H(R,y)   = " + str(self.H_Ry)  + "\n"
        s = s + "H(x,y)   = " + str(self.H_xy)  + "\n"
        s = s + "H(R,x,y) = " + str(self.H_Rxy) + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Conditional Entropies:\n"
        s = s + "H(R|x)   = " + str(self.H_R_x)  + "\n"
        s = s + "H(R|y)   = " + str(self.H_R_y)  + "\n"
        s = s + "H(x|R)   = " + str(self.H_x_R)  + "\n"
        s = s + "H(y|R)   = " + str(self.H_y_R)  + "\n"
        s = s + "H(x|y)   = " + str(self.H_x_y)  + "\n"
        s = s + "H(y|x)   = " + str(self.H_y_x)  + "\n"
        s = s + "H(R|x,y) = " + str(self.H_R_xy) + "\n"
        s = s + "H(x|R,y) = " + str(self.H_x_Ry) + "\n"
        s = s + "H(y|R,x) = " + str(self.H_y_Rx) + "\n"
        s = s + "H(R,x|y) = " + str(self.H_Rx_y) + "\n"
        s = s + "H(R,y|x) = " + str(self.H_Ry_x) + "\n"
        s = s + "H(x,y|R) = " + str(self.H_xy_R) + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Mutual Information:\n"
        s = s + "I(R,x)   = " + str(self.I_Rx)   + "\n"
        s = s + "I(R,y)   = " + str(self.I_Ry)   + "\n"
        s = s + "I(x,y)   = " + str(self.I_xy)   + "\n"
        s = s + "I(R,x,y) = " + str(self.I_Rxy)  + "\n"
        s = s + "\n"
        s = s + "\n"

        s = s + "Conditional Mutual Information:\n"
        s = s + "I(R,x|y) = " + str(self.I_Rx_y)  + "\n"
        s = s + "I(R,y|x) = " + str(self.I_Ry_x)  + "\n"
        s = s + "I(x,y|R) = " + str(self.I_xy_R)  + "\n"

        s = s + "Nomralized Measures:\n"
        s = s + "U(R|x)   = " + str(self.U_R_x)  + "\n"
        s = s + "U(x|R)   = " + str(self.U_x_R)  + "\n"
        s = s + "U(R|y)   = " + str(self.U_R_y)  + "\n"
        s = s + "U(y|R)   = " + str(self.U_y_R)  + "\n"
        s = s + "U(x|y)   = " + str(self.U_x_y)  + "\n"
        s = s + "U(y|x)   = " + str(self.U_y_x)  + "\n"
        s = s + "U(R,x|y) = " + str(self.U_Rx_y) + "\n"
        s = s + "U(R,y|x) = " + str(self.U_Ry_x) + "\n"
        s = s + "U(x,y|R) = " + str(self.U_xy_R) + "\n"
        s = s + "U(R|x,y) = " + str(self.U_R_xy) + "\n"
        s = s + "U(x|R,y) = " + str(self.U_x_Ry) + "\n"
        s = s + "U(y|R,x) = " + str(self.U_y_Rx) + "\n"
        s = s + "U(R,x)   = " + str(self.U_Rx)   + "\n"
        s = s + "U(R,y)   = " + str(self.U_Ry)   + "\n"
        s = s + "U(x,y)   = " + str(self.U_xy)   + "\n"
        s = s + "U(R,x,y) = " + str(self.U_Rxy)  + "\n"
        s = s + "M1       = " + str(self.M1)     + "\n"
        s = s + "M2       = " + str(self.M2)     + "\n"
        s = s + "M3       = " + str(self.M3)     + "\n"
        s = s + "M4       = " + str(self.M4)     + "\n"
        s = s + "M5       = " + str(self.M5)     + "\n"
        s = s + "S        = " + str(self.S)      + "\n"

        s = s + "mu_U_x_R =    " + str(self.mu_U_x_R)    + "\n"
        s = s + "sigma_U_x_R = " + str(self.sigma_U_x_R) + "\n"
        s = s + "se_U_x_R =    " + str(self.se_U_x_R)    + "\n"
        s = s + "cv_U_x_R =    " + str(self.cv_U_x_R)    + "\n"

        s = s + "mu_U_y_R =    " + str(self.mu_U_y_R)    + "\n"
        s = s + "sigma_U_y_R = " + str(self.sigma_U_y_R) + "\n"
        s = s + "se_U_y_R =    " + str(self.se_U_y_R)    + "\n"
        s = s + "cv_U_y_R =    " + str(self.cv_U_y_R)    + "\n"

        return s

    def LaTeX(self, SingleTable=False, sinc="cprheim"):
        """
        LaTeX representation of the contingency table.

        Keyword Args:		
           * `SingleTable` (True or False): Flag for the number of output tablulars for the contingency table.
           * `sinc`        (str):  String include flags that dictates which summary data to print. For example, 
             'si' would only print the probability table and the mutual information. Default prints all data available.  
             Flags have the following meanings:
                * `c`: The contingency table itself.
                * `p`: The probability table itself.
                * `r`: Basic contingency table results.
                * `h`: Entropies.
                * `e`: Conditional Entropies.
                * `i`: Mutual Information & Conditional Mutual Information.
                * `m`: Normalized Measures.

        Returns:
           * LaTeX valid string that represnets the contingency table and its results.
           * If `SingleTable` is True, then LaTeX output of the contingency table will be in a single, large table 
             rather than being split into separate tables for each slice.
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
            #Contingency Table(s)
            if SingleTable:
                s = s + "\\begin{center}\n"
                s = s + "\\begin{table}[htbp]\n"
                s = s + "\\caption{{Contingency Tables for {title}.}}\n".format(**self.__dict__)
                s = s + "\\begin{center}\n"
                s = s + begin_tabular
                s = s + self.hline
            else:
                s = s + "Contingency Tables for {title}.\n".format(**self.__dict__)

            for k in range(self.K):
                if SingleTable:
                    if self.ydiscrete:
                        s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Slice for {ylabel} $=$ {bounds[0]:{sform}}}}{endrow}".format(bounds=self.ybounds[k:k+1], **self.__dict__) 
                    else:
                        s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Slice for {bounds[0]:{sform}} $<$ {ylabel} $<$ {bounds[1]:{sform}}}}{endrow}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
                    s = s + self.hline
                else:
                    s = s + "\\begin{center}\n"
                    s = s + "\\begin{table}[htbp]\n"
                    if self.ydiscrete:
                        s = s + "\\caption{{Slice for {ylabel} $=$ {bounds[0]:{sform}}}}".format(bounds=self.ybounds[k:k+1], **self.__dict__) 
                    else:
                        s = s + "\\caption{{Slice for {bounds[0]:{sform}} $<$ {ylabel} $<$ {bounds[1]:{sform}}}}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
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
                        s = s + "${Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                    else:
                        s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                    for j in range(self.J):
                        s = s + "&${0:{1}}$".format(self.N_ijk[i][j][k], self.sform)
                    s = s + "&${0:{1}}$".format(self.N_idotk[i][k], self.sform)
                    s = s + self.endrow
                    s = s + self.hline
                #Bottom row
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.N_dotjk[j][k], self.sform)
                s = s + "&${0:{1}}$".format(self.N_dotdotk[k], self.sform)
                s = s + self.endrow
                s = s + self.hline

                if SingleTable:
                    s = s + self.hline
                else:
                    s = s + "\\end{tabular}\n"
                    s = s + "\\end{center}\n"
                    s = s + "\\end{table}\n"
                    s = s + "\\end{center}\n\n"

            if SingleTable:
                s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Summary for summation over {ylabel}}}{endrow}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
                s = s + self.hline
            else:
                s = s + "\\begin{center}\n"
                s = s + "\\begin{table}[htbp]\n"
                s = s + "\\caption{{Summary for summation over {ylabel}}}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
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
                    s = s + "${Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                else:
                    s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.N_ijdot[i][j], self.sform)
                s = s + "&${0:{1}}$".format(self.N_idotdot[i], self.sform)
                s = s + self.endrow
                s = s + self.hline
            #Bottom row
            for j in range(self.J):
                s = s + "&${0:{1}}$".format(self.N_dotjdot[j], self.sform)
            s = s + "&${0:{1}}$".format(self.N, self.sform)
            s = s + self.endrow
            s = s + self.hline

            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"
            s = s + "\n"
            s = s + "\n"


        if 'p' in sinc:
            #Probability Table(s)
            if SingleTable:
                s = s + "\\begin{center}\n"
                s = s + "\\begin{table}[htbp]\n"
                s = s + "\\caption{{Probability Tables for {title}.}}\n".format(**self.__dict__)
                s = s + "\\begin{center}\n"
                s = s + begin_tabular
                s = s + self.hline
            else:
                s = s + "Probability Tables for {title}.\n".format(**self.__dict__)

            for k in range(self.K):
                if SingleTable:
                    if self.ydiscrete:
                        s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Slice for {ylabel} $=$ {bounds[0]:{sform}}}}{endrow}".format(bounds=self.ybounds[k:k+1], **self.__dict__) 
                    else:
                        s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Slice for {bounds[0]:{sform}} $<$ {ylabel} $<$ {bounds[1]:{sform}}}}{endrow}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
                    s = s + self.hline
                else:
                    s = s + "\\begin{center}\n"
                    s = s + "\\begin{table}[htbp]\n"
                    if self.ydiscrete:
                        s = s + "\\caption{{Slice for {ylabel} $=$ {bounds[0]:{sform}}}}".format(bounds=self.ybounds[k:k+1], **self.__dict__) 
                    else:
                        s = s + "\\caption{{Slice for {bounds[0]:{sform}} $<$ {ylabel} $<$ {bounds[1]:{sform}}}}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
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
                        s = s + "${Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                    else:
                        s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                    for j in range(self.J):
                        s = s + "&${0:{1}}$".format(self.p_ijk[i][j][k], self.sform)
                    s = s + "&${0:{1}}$".format(self.p_idotk[i][k], self.sform)
                    s = s + self.endrow
                    s = s + self.hline
                #Bottom row
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.p_dotjk[j][k], self.sform)
                s = s + "&${0:{1}}$".format(self.p_dotdotk[k], self.sform)
                s = s + self.endrow
                s = s + self.hline

                if SingleTable:
                    s = s + self.hline
                else:
                    s = s + "\\end{tabular}\n"
                    s = s + "\\end{center}\n"
                    s = s + "\\end{table}\n"
                    s = s + "\\end{center}\n\n"

            if SingleTable:
                s = s + "\\multicolumn{{{columns}}}{{|c|}}{{Probability summary for summation over {ylabel}}}{endrow}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
                s = s + self.hline
            else:
                s = s + "\\begin{center}\n"
                s = s + "\\begin{table}[htbp]\n"
                s = s + "\\caption{{Probability summary for summation over {ylabel}}}".format(bounds=self.ybounds[k:k+2], **self.__dict__) 
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
                    s = s + "${Rlabel} $= {bounds[0]:{sform}}$".format(bounds=self.Rbounds[i:i+1], **self.__dict__)
                else:
                    s = s + "${bounds[0]:{sform}} <$ {Rlabel} $< {bounds[1]:{sform}}$".format(bounds=self.Rbounds[i:i+2], **self.__dict__)
                for j in range(self.J):
                    s = s + "&${0:{1}}$".format(self.p_ijdot[i][j], self.sform)
                s = s + "&${0:{1}}$".format(self.p_idotdot[i], self.sform)
                s = s + self.endrow
                s = s + self.hline
            #Bottom row
            for j in range(self.J):
                s = s + "&${0:{1}}$".format(self.p_dotjdot[j], self.sform)
            s = s + "&${0:{1}}$".format(1.0, self.sform)
            s = s + self.endrow
            s = s + self.hline

            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"
            s = s + "\n"
            s = s + "\n"

        if 'r' in sinc:
            #Contingency Table Results
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Results for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "$\\chi^2$ (poorly defined)&${chi2:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "G-test&${G:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "Cramer\'s V&${V:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "Contingency Coefficient C (poorly defined)&${C:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "R-value Statistic&${R_stat:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"
            s = s + "\n"
            s = s + "\n"

        if 'h' in sinc:
            #Entropies
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Entropies for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "H(R)&${H_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x)&${H_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(y)&${H_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,x)&${H_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,y)&${H_Ry:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x,y)&${H_xy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,x,y)&${H_Rxy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"
            s = s + "\n"
            s = s + "\n"

        if 'e' in sinc:
            #Conditional Entropies
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Conditional Entropies for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "H(R$|$x)&${H_R_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R$|$y)&${H_R_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x$|$R)&${H_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(y$|$R)&${H_y_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x$|$y)&${H_x_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(y$|$x)&${H_y_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R$|$x,y)&${H_R_xy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x$|$R,y)&${H_x_Ry:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(y$|$R,x)&${H_y_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,x$|$y)&${H_Rx_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R,y$|$x)&${H_Ry_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(x,y$|$R)&${H_xy_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "H(R$|$x$|$x)&${H_R_x_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "\\end{tabular}\n"
            s = s + "\\end{center}\n"
            s = s + "\\end{table}\n"
            s = s + "\\end{center}\n\n"
            s = s + "\n"
            s = s + "\n"


        if 'i' in sinc:
            #Mutual and Conditional Mutual Information
            s = s + "\\begin{center}\n"
            s = s + "\\begin{table}[htbp]\n"
            s = s + "\\caption{{Contingency Table Mutual Information for {title}.}}\n".format(**self.__dict__)
            s = s + "\\begin{center}\n"
            s = s + "\\begin{tabular}{|l|c|}\n"
            s = s + self.hline
            s = s + "I(R,x)&${I_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(R,y)&${I_Ry:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(x,y)&${I_xy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(R,x,y)&${I_Rxy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(R,x$|$y)&${I_Rx_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(R,y$|$x)&${I_Ry_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "I(x,y$|$R)&${I_xy_R:{sform}}${endrow}".format(**self.__dict__) 
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
            s = s + "U(R$|$y)&${U_R_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(y$|$R)&${U_y_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(x$|$y)&${U_x_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(y$|$x)&${U_y_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,x$|$y)&${U_Rx_y:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,y$|$x)&${U_Ry_x:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(x,y$|$R)&${U_xy_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R$|$x,y)&${U_R_xy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(x$|$R,y)&${U_x_Ry:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(y$|$R,x)&${U_y_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,x)&${U_Rx:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,y)&${U_Ry:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(x,y)&${U_xy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "U(R,x,y)&${U_Rxy:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M1: $\\frac{{I(R,x,y)}}{{\\sqrt[3]{{H(R) H(x) H(y)}}}}$&${M1:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M2: $\\frac{{I(R,x,y)}}{{\\min[H(R), H(x), H(y)]}}$&${M2:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M3: $\\frac{{I(R,x,y)}}{{\\max[H(R), H(x), H(y)]}}$&${M3:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M4: $\\frac{{I(R,x,y)}}{{H(R,x,y))}}$&${M4:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "M5: $1 - \\frac{{I(R,x,y)}}{{H(R,x,y))}}$&${M5:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "S: $\\frac{{U(x|R) + U(y|R)}}{{2 U(x,y|R))}}$&${S:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline

            s = s + "$\mu(U(x|R|y))$:&${mu_U_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\sigma(U(x|R|y))$:&${sigma_U_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\mbox{{Std. Error}}(U(x|R|y))$:&${se_U_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\mbox{{Coef. of Var.}}(U(x|R|y))$:&${cv_U_x_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline

            s = s + "$\mu(U(y|R|x))$:&${mu_U_y_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\sigma(U(y|R|x))$:&${sigma_U_y_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\mbox{{Std. Error}}(U(y|R|x))$:&${se_U_y_R:{sform}}${endrow}".format(**self.__dict__) 
            s = s + self.hline
            s = s + "$\mbox{{Coef. of Var.}}(U(y|R|x))$:&${cv_U_y_R:{sform}}${endrow}".format(**self.__dict__) 
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
            self.H5table = CT3D_TableModel(self.I, self.J, self.K, extracols, 
                self.xdiscrete, self.ydiscrete, self.Rdiscrete)

        for key in self.H5table:
            row[key] = getattr(self, key)
        row.append()
        return

    #Document Strings for class attributes
    x = None
    """First Independent variable (parameter) dataset."""
    y = None
    """Second Independent variable (parameter) dataset."""
    R = None
    """Response (dependent) variable dataset."""

    N = None
    """Number of data points in contingency table."""
    I = None
    """Number of R response bins."""
    J = None
    """Number of x independent variable bins."""
    K = None
    """Number of y independent variable bins."""


    Rdiscrete = None
    """Flag (bool) for if the response is a discrete (or continuous) variable."""
    xdiscrete = None
    """Flag (bool) for if the first indpendent parameter is a discrete (or continuous) variable."""
    ydiscrete = None
    """Flag (bool) for if the second indpendent parameter is a discrete (or continuous) variable."""


    datarange = None
    """List of ranges that the table data and bins are defined on.
       datarange = [x-lower, x-upper, y-lower, y-upper, R-lower, R-upper]
    """

    xscale = None
    """Scale of the first independent variable x. The following are valid flags:
       * 'linear': Linearly uniform bins.
       * 'log':    Log-uniform bins. 
       * 'nines':  One-minus-log-uniform bins.
    """
    yscale = None
    """Scale of the second independent variable y. The following are valid flags:
       * 'linear': Linearly uniform bins.
       * 'log':    Log-uniform bins. 
       * 'nines':  One-minus-log-uniform bins.
    """
    Rscale = None
    """Scale of the Response variable R. The following are valid flags:
       * 'linear': Linearly uniform bins.
       * 'log':    Log-uniform bins. 
       * 'nines':  One-minus-log-uniform bins.
    """

    xbounds = None
    """List of first independent variable x bounds, length J+1.  Denoted :math:`b^x`."""
    ybounds = None
    """List of second independent variable y bounds, length K+1.  Denoted :math:`b^y`."""
    Rbounds = None
    """List of response variable R bounds, length I+1.  Denoted :math:`b^R`."""
    

    N_idotdot = None
    """List of sums over j- and k-indices.

    .. math::
       N_{i\cdot\cdot} = \sum_{jk} N_{ijk} 
    """
    N_dotjdot = None
    """List of sums over i- and k-indices.

    .. math::
       N_{\cdot j \cdot} = \sum_{ik} N_{ijk} 
    """
    N_dotdotk = None
    """List of sums over i- and j-indices.

    .. math::
       N_{\cdot\cdot k} = \sum_{ij} N_{ijk} 
    """
    N_ijdot = None
    """Matrix of sums over kth index. Matrix is of size IxJ.

    .. math::
       N_{ij\cdot} = \sum_{k} N_{ijk} 
    """
    N_idotk = None
    """Matrix of sums over jth index. Matrix is of size IxK.

    .. math::
       N_{i\cdot k} = \sum_{j} N_{ijk} 
    """
    N_dotjk = None
    """Matrix of sums over ith index. Matrix is of size JxK.

    .. math::
       N_{\cdot jk} = \sum_{i} N_{ijk} 
    """

    N_ijk = None
    """Raw Contingency Table (Matrix) of size IxJxK.

    .. math::
       N_{ijk} = \sum_n^N \Pi_{b_i^R, b_{i+1}^R}(R_n) \cdot \Pi_{b_j^x, b_{j+1}^x}(x_n) \cdot \Pi_{b_k^y, b_{k+}^y}(y_n)

    Here, :math:`\Pi_{a,b}(x)` is the `Boxcar Function <http://mathworld.wolfram.com/BoxcarFunction.html>`_ 
    and :math:`b^x` are the bounds in x (self.xbounds), etc.
    """


    E_ijk = None
    """Expected Frequency Matrix of size IxJxK.  
    Calculated via Equation (4.7) in *"The Analysis of Contingency Tables"* by B.S. Everitt.

    .. math::
       E_{ijk} = \\frac{N_{i\cdot\cdot} \cdot N_{\cdot j\cdot} \cdot N_{\cdot\cdot k}}{N^2}
    """
    chi2 = None
    """Chi-Squared Statistic.  Calculated via Equation (4.8) in *"The Analysis of Contingency Tables"* by B.S. Everitt.

    .. math::	
       \chi^2 =  \sum_{ijk} \\frac{\left(N_{ijk} - E_{ijk}\\right)^2}{E_{ijk}}
    """
    G = None
    """Log-likelihood ratio G-test statistic. Better than Chi-Squared.

    .. math::
       G = 2 \sum_{ijk} N_{ijk} \cdot \ln\left(\\frac{N_{ijk}}{E_{ijk}}\\right)
    """
    V = None
    """Cramer's V Statistic. Calculated with G and not Chi-Squared.

    .. math::
       V = \sqrt{\\frac{G}{N \cdot (\min[I, J, K] - 1)}}	
    """
    C = None
    """Contingency Coefficient Statistic.  Calculated with *G* and not Chi-Squared. BEWARE of use in 3D context!

    .. math::
       C = \sqrt{\\frac{G}{G + N}}
    """

    #Probabilities
    p_idotdot = None
    """List of ith bin probabilities.
    
    .. math::
       p_{i\cdot\cdot} = \\frac{N_{i\cdot\cdot}}{N}
    """
    p_dotjdot = None
    """List of jth bin probabilities.
    
    .. math::
       p_{\cdot j \cdot} = \\frac{N_{\cdot j \cdot}}{N}
    """
    p_dotdotk = None
    """List of kth bin probabilities.
    
    .. math::
       p_{\cdot\cdot k} = \\frac{N_{\cdot\cdot k}}{N}
    """

    p_ijdot = None
    """Matrix of ith- and jth-bin probabilities. Matrix is of size IxJ.

    .. math::
       p_{ij\cdot} = \\frac{N_{ij\cdot}}{N}
    """
    p_idotk = None
    """Matrix of ith- and kth-bin probabilities. Matrix is of size IxK.

    .. math::
       p_{i\cdot k} = \\frac{N_{i\cdot k}}{N}
    """
    p_dotjk = None
    """Matrix of jth- and kth-bin probabilities. Matrix is of size JxK.

    .. math::
       p_{\cdot jk} = \\frac{N_{\cdot jk}}{N}
    """

    p_ijk = None
    """Raw Probability Table (Matrix) of size IxJxK.

    .. math::
       p_{ijk} = \\frac{N_{ijk}}{N}
    """

    #Entropies
    H_x = None
    """Entropy H(x) of first independent variable x.

    .. math::
       H(x) = - \sum_j p_{\cdot j \cdot} \ln(p_{\cdot j \cdot})	
    """
    H_y = None
    """Entropy H(y) of second independent variable y.

    .. math::
       H(y) = - \sum_j p_{\cdot\cdot k} \ln(p_{\cdot\cdot k})	
    """
    H_R = None
    """Entropy H(R) of response variable R.

    .. math::
       H(R) = - \sum_j p_{i\cdot\cdot} \ln(p_{i\cdot\cdot})	
    """

    H_Rx = None
    """Joint Entropy H(R,x) of response variable R and first independent variable x.

    .. math::
       H(R,x) = - \sum_{i,j} p_{ij\cdot} \ln(p_{ij\cdot})	
    """
    H_Ry = None
    """Joint Entropy H(R,y) of response variable R and second independent variable y.

    .. math::
       H(R,y) = - \sum_{i,j} p_{i\cdot k} \ln(p_{i\cdot k})	
    """
    H_xy = None
    """Joint Entropy H(x,y) of independent variables x and y.

    .. math::
       H(x,y) = - \sum_{j,k} p_{\cdot jk} \ln(p_{\cdot jk})	
    """

    H_Rxy = None
    """Joint Entropy H(R,x,y) of response R and independent variables x and y.

    .. math::
       H(R,x,y) = - \sum_{i,j,k} p_{ijk} \ln(p_{ijk})	
    """

    #Conditional Entropies
    H_R_x  = None
    """Conditional Entropy H(R|x).

    .. math::
       H(R|x) = - \sum_{i,j} p_{ij\cdot} \ln\left(\\frac{p_{ij\cdot}}{p_{\cdot j \cdot}}\\right)	
    """
    H_R_y  = None
    """Conditional Entropy H(R|y).

    .. math::
       H(R|y) = - \sum_{i,k} p_{i\cdot k} \ln\left(\\frac{p_{i\cdot k}}{p_{\cdot\cdot k}}\\right)	
    """

    H_x_R  = None
    """Conditional Entropy H(x|R).

    .. math::
       H(x|R) = - \sum_{i,j} p_{ij\cdot} \ln\left(\\frac{p_{ij\cdot}}{p_{i\cdot\cdot}}\\right)	
    """
    H_x_y  = None
    """Conditional Entropy H(x|y).

    .. math::
       H(x|y) = - \sum_{j,k} p_{\cdot jk} \ln\left(\\frac{p_{\cdot jk}}{p_{\cdot\cdot k}}\\right)	
    """

    H_y_R  = None
    """Conditional Entropy H(y|R).

    .. math::
       H(y|R) = - \sum_{i,k} p_{i\cdot k} \ln\left(\\frac{p_{i\cdot k}}{p_{i\cdot\cdot}}\\right)	
    """
    H_y_x  = None
    """Conditional Entropy H(y|x).

    .. math::
       H(y|x) = - \sum_{j,k} p_{\cdot jk} \ln\left(\\frac{p_{\cdot jk}}{p_{\cdot j \cdot}}\\right)	
    """

    H_R_xy = None
    """Conditional Entropy H(R|x,y).

    .. math::
       H(R|x,y) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{\cdot jk}}\\right)	
    """
    H_x_Ry = None
    """Conditional Entropy H(x|R,y).

    .. math::
       H(x|R,y) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{i\cdot k}}\\right)	
    """
    H_y_Rx = None
    """Conditional Entropy H(y|R,x).

    .. math::
       H(y|R,x) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{ij\cdot}}\\right)	
    """

    H_Rx_y = None
    """Conditional Entropy H(R,x|y).

    .. math::
       H(R,x|y) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{\cdot\cdot k}}\\right)	
    """
    H_Ry_x = None
    """Conditional Entropy H(R,y|x).

    .. math::
       H(R,y|x) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{\cdot j \cdot}}\\right)	
    """
    H_xy_R = None
    """Conditional Entropy H(x,y|R).

    .. math::
       H(x,y|R) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{i\cdot\cdot}}\\right)	
    """

    #Mutual Information
    I_Rx = None
    """Mutual Information I(R,x).

    .. math::
       I(R,x) = - \sum_{i,j} p_{ij\cdot} \ln\left(\\frac{p_{ij\cdot}}{p_{i\cdot\cdot} \cdot p_{\cdot j \cdot}}\\right)	
    """
    I_Ry = None
    """Mutual Information I(R,y).

    .. math::
       I(R,y) = - \sum_{i,j} p_{i\cdot k} \ln\left(\\frac{p_{i\cdot k}}{p_{i\cdot\cdot} \cdot p_{\cdot\cdot k}}\\right)	
    """
    I_xy = None
    """Mutual Information I(x,y).

    .. math::
       I(x,y) = - \sum_{j,k} p_{\cdot jk} \ln\left(\\frac{p_{\cdot jk}}{p_{\cdot j \cdot} \cdot p_{\cdot\cdot k}}\\right)	
    """

    I_Rxy = None
    """Mutual Information I(R,x,y).

    .. math::
       I(R,x,y) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{ijk}}{p_{i\cdot\cdot} \cdot p_{\cdot j \cdot} \cdot p_{\cdot\cdot k}}\\right)
    """

    #Conditional Mutual Information
    I_Rx_y = None
    """Conditional Mutual Information I(R,x|y).

    .. math::
       I(R,x|y) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{\cdot\cdot k} \cdot p_{ijk}}{p_{i\cdot k} \cdot p_{\cdot jk}}\\right)
    """
    I_Ry_x = None
    """Conditional Mutual Information I(R,y|x).

    .. math::
       I(R,y|x) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{\cdot j \cdot} \cdot p_{ijk}}{p_{ij\cdot} \cdot p_{\cdot jk}}\\right)
    """
    I_xy_R = None
    """Conditional Mutual Information I(x,y|R)

    .. math::
       I(x,y|R) = - \sum_{i,j,k} p_{ijk} \ln\left(\\frac{p_{i\cdot\cdot} \cdot p_{ijk}}{p_{ij\cdot} \cdot p_{i\cdot k}}\\right)
    """

    #Calculate the R-statistic
    R_stat = None
    """R-Value Statistic.

    .. math::
       r = \sqrt{1 - e^{-2 I(R,x,y)}}
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
    U_R_y = None
    """Uncertainty coefficient (coefficient of constraint) U(R|y).  Calculated from,

    .. math::

       U(R|y) = \\frac{I(R,y)}{H(R)}	   
    """
    U_y_R = None
    """ Uncertainty coefficient (coefficient of constraint) U(y|R).  Calculated from,

    .. math::

       U(y|R) = \\frac{I(R,y)}{H(y)}	   
    """
    U_x_y = None
    """Uncertainty coefficient (coefficient of constraint) U(x|y).  Calculated from,

    .. math::

       U(x|y) = \\frac{I(x,y)}{H(x)}	   
    """
    U_y_x = None
    """Uncertainty coefficient (coefficient of constraint) U(y|x).  Calculated from,

    .. math::

       U(y|x) = \\frac{I(x,y)}{H(y)}	   
    """
    U_Rx_y = None
    """Uncertainty coefficient (coefficient of constraint) U(R,x|y).  Calculated from,

    .. math::

       U(R,x|y) = \\frac{I(R,x,y)}{H(R,x)}	   
    """
    U_Ry_x = None
    """Uncertainty coefficient (coefficient of constraint) U(R,y|x).  Calculated from,

    .. math::

       U(R,y|x) = \\frac{I(R,x,y)}{H(R,y)}	   
    """
    U_xy_R = None
    """Uncertainty coefficient (coefficient of constraint) U(x,y|R).  Calculated from,

    .. math::

       U(x,y|R) = \\frac{I(R,x,y)}{H(x,y)}	   
    """

    U_R_xy = None
    """Uncertainty coefficient (coefficient of constraint) U(R|x,y).  Calculated from,

    .. math::

       U(R|x,y) = \\frac{I(R,x,y)}{H(R)}
    """
    U_x_Ry = None
    """Uncertainty coefficient (coefficient of constraint) U(x|R,y).  Calculated from,

    .. math::

       U(x|R,y) = \\frac{I(R,x,y)}{H(x)}
    """

    U_y_Rx = None
    """Uncertainty coefficient (coefficient of constraint) U(y|R,x).  Calculated from,

    .. math::

       U(y|R,x) = \\frac{I(R,x,y)}{H(y)}
    """

    U_Rx  = None
    """Symetric Uncertainty U(R,x).  Calculated from,

    .. math::

       U(R,x) = 2\\frac{I(R,x)}{H(R) + H(x)}	   
    """

    U_Ry  = None
    """Symetric Uncertainty U(R,y).  Calculated from,

    .. math::

       U(R,y) = 2\\frac{I(R,y)}{H(R) + H(y)}	   
    """

    U_xy  = None
    """Symetric Uncertainty U(x,y).  Calculated from,

    .. math::

       U(x,y) = 2\\frac{I(x,y)}{H(R) + H(y)}	   
    """
    U_Rxy  = None
    """Symetric Uncertainty U(R,x,y).  Calculated from,

    .. math::

       U(R,x,y) = 3\\frac{I(R,x,y)}{H(R) + H(x) + H(y)}	   
    """
    
    M1 = None
    """Normalized Measure 1, in analogy to Strehl, Alexander; Joydeep Ghosh (2002). *"Cluster ensembles -- a knowledge reuse 
    framework for combining multiple partitions"*. Journal of Machine Learning Research 3: 583-617. 
    doi:10.1162/153244303321897735.  Equation (2):

    .. math::

       M_1 = \\frac{I(R,x,y)}{\sqrt[3]{H(R) H(x) H(y)}}
    """

    M2 = None
    """Normalized Measure 2, in analogy to Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data 
    mining, in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.34),

    .. math::

       M_2 = \\frac{I(R,x,y)}{\min[H(R), H(x), H(y)]}
    """

    M3 = None
    """Normalized Measure 3, in analogy to Kvalseth, T.O. *"Entropy and correlation: some comments"*, IEEE Transactions On Systems, Man, 
    and Cybernetics, SMC-17, 517-519, 1987.  

    .. math::

       M_3 = \\frac{I(R,x,y)}{\max[H(R), H(x), H(y)]}
    """

    M4 = None
    """Normalized Measure 4, in analogy to Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data 
    mining, in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.35),

    .. math::

       M_4 = \\frac{I(R,x,y)}{H(R,x,y)}
    """

    M5 = None
    """Normalized Measure 5, in analogy to Yao, Y. Y. (2003) *"Information-theoretic measures for knowledge discovery and data 
    mining, in Entropy Measures, Maximum Entropy Principle and Emerging Applications"*, Karmeshu (ed.), Springer, pp. 115-136.
    Equation (6.35),

    .. math::

       M_5 = 1 - \\frac{I(R,x,y)}{H(R,x)}
    """

    S = None
    """Normalized Measure S, from the minds of E. A. Schnider and A. M. Scopatz.  Measures whether the response is determined 
    by only one of the indpendent variables or if the dependence is based on both indepenent parameters.  The interpretation
    of this normalized metric is as follows:
       * S = 1: Response soley determined by one variable.  Second varibale contribution irrelevant.
       * S < 0: Respose by combination of both variables, though one independent parameter shows a weaker dependence.
       * S = 0: Response is perfectly determined by both variables in tandem.

    .. math::

       S = \\frac{U(x|R) + U(y|R)}{2 U(x,y|R)}
    """
    mu_U_x_R = None
    """Dr. Jun Li's method:  Take the mean of U(x|R) of each kth slice of the table."""
    sigma_U_x_R = None
    """Dr. Jun Li's method:  Take the standard deviation of U(x|R) of each kth slice of the table."""
    se_U_x_R = None
    """Dr. Jun Li's method:  Take the standard error of U(x|R) of each kth slice of the table."""
    cv_U_x_R = None
    """Dr. Jun Li's method:  Take the coefficient of variation of U(x|R) of each kth slice of the table."""

    mu_U_y_R = None
    """Dr. Jun Li's method:  Take the mean of U(y|R) of each jth slice of the table."""
    sigma_U_y_R = None
    """Dr. Jun Li's method:  Take the standard deviation of U(y|R) of each jth slice of the table."""
    se_U_y_R = None
    """Dr. Jun Li's method:  Take the stamdard error of U(y|R) of each jth slice of the table."""
    cv_U_y_R = None
    """Dr. Jun Li's method:  Take the coefficient of variation of U(y|R) of each jth slice of the table."""

    cv_U_x_y_R = None
    """Dr. Jun Li's method:  Average of cv_U_x_R and cv_U_y_R."""

    eta_R_x_y = None
    """The entopy of R for given (x, y), normalized by ln(I)."""
    mu_eta_R_x_y = None
    """The mean of all values in eta(R|x|y)"""
    sigma_eta_R_x_y = None
    """The standard deviation of all values in eta(R|x|y)"""
    se_eta_R_x_y = None
    """The standard error of all values in eta(R|x|y)"""
    cv_eta_R_x_y = None
    """The coefficent of variation of all values in eta(R|x|y)"""

    #Finally, copy over some options to the CT itself
    sbuf = None
    """String buffer for contingency table printing; integer number spaces per column, default = 40."""
    salign = None
    """String alignment and buffer for contingency table printing; default = '^sbuf' [centered]."""
    sform = None
    """String formatting for contingency table data, default = '.6G'."""
    xlabel = None
    """First independent variable label, default = 'x'."""
    ylabel = None
    """Second independent variable label, default = 'y'."""
    Rlabel = None
    """Response variable label, default = 'R'."""
    title = None
    """Contingency table label, default = '{xlabel}, {ylabel} to {Rlabel}'."""

    #create a copy of an HDF5 table model for this CT
    H5table = None
    """HDF5 table model for this contingency table and its results."""

#3D Contingency Table Data Model
def CT3D_TableModel(I, J, K, extracols={}, xdscrt=False, ydscrt=False, Rdscrt=False):
    """Three-Dimensional Contingency Table Model for HDF5.

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

          The standard colums range from positions 0 to 67.  
          To prepend data, start from pos=-1 and work down.
          To append columns, start from pos=28 and work up!

          Note that to have ContingencyTable3D.FillTableRow(row) automatically
          fill in any extra column data you'll need to set the approriate data 
          as an attribute of the Contingency Table itself first::

            ct = ContingencyTable3D(...)
            setattr(ct, 'time', 1.0)
            ct.FillTableRow(mytablerow)

        * `xdscrt` (bool): Flag for x being a discrete (or continuous) variable.
        * `ydscrt` (bool): Flag for y being a discrete (or continuous) variable.
        * `Rdscrt` (bool): Flag for R being a discrete (or continuous) variable.
    """

    CT3D_TM = {
        'xlabel':  tb.StringCol(50, pos=0),
        'ylabel':  tb.StringCol(50, pos=1),
        'Rlabel':  tb.StringCol(50, pos=2),

        'H_x':     tb.Float64Col(pos=3),
        'H_y':     tb.Float64Col(pos=4),
        'H_R':     tb.Float64Col(pos=5),

        'H_Rx':    tb.Float64Col(pos=6),
        'H_Ry':    tb.Float64Col(pos=7),
        'H_xy':    tb.Float64Col(pos=8),

        'H_Rxy':   tb.Float64Col(pos=9),

        'H_R_x':   tb.Float64Col(pos=10),
        'H_x_R':   tb.Float64Col(pos=11),
        'I_Rx':    tb.Float64Col(pos=12),

        'H_R_y':   tb.Float64Col(pos=13),
        'H_y_R':   tb.Float64Col(pos=14),
        'I_Ry':    tb.Float64Col(pos=15),

        'H_x_y':   tb.Float64Col(pos=16),
        'H_y_x':   tb.Float64Col(pos=17),
        'I_xy':    tb.Float64Col(pos=18),

        'H_R_xy':   tb.Float64Col(pos=19),
        'H_x_Ry':   tb.Float64Col(pos=20),
        'H_y_Rx':   tb.Float64Col(pos=21),

        'H_Rx_y':   tb.Float64Col(pos=22),
        'H_Ry_x':   tb.Float64Col(pos=23),
        'H_xy_R':   tb.Float64Col(pos=24),

        'I_Rxy':    tb.Float64Col(pos=25),
        'I_Rx_y':   tb.Float64Col(pos=26),
        'I_Ry_x':   tb.Float64Col(pos=27),
        'I_xy_R':   tb.Float64Col(pos=28),

        'U_R_x':   tb.Float64Col(pos=29),
        'U_x_R':   tb.Float64Col(pos=30),

        'U_R_y':   tb.Float64Col(pos=31),
        'U_y_R':   tb.Float64Col(pos=32),

        'U_x_y':   tb.Float64Col(pos=33),
        'U_y_x':   tb.Float64Col(pos=34),

        'U_Rx_y':  tb.Float64Col(pos=35),
        'U_Ry_x':  tb.Float64Col(pos=36),
        'U_xy_R':  tb.Float64Col(pos=37),

        'U_R_xy':  tb.Float64Col(pos=38),
        'U_y_Rx':  tb.Float64Col(pos=39),
        'U_x_Ry':  tb.Float64Col(pos=40),

        'U_Rx':    tb.Float64Col(pos=41),
        'U_Ry':    tb.Float64Col(pos=42),
        'U_xy':    tb.Float64Col(pos=43),
        'U_Rxy':   tb.Float64Col(pos=44),

        'M1':      tb.Float64Col(pos=45),
        'M2':      tb.Float64Col(pos=46),
        'M3':      tb.Float64Col(pos=47),
        'M4':      tb.Float64Col(pos=48),
        'M5':      tb.Float64Col(pos=49),
        'S':       tb.Float64Col(pos=50),

        'mu_U_x_R':    tb.Float64Col(pos=51),
        'sigma_U_x_R': tb.Float64Col(pos=52),
        'se_U_x_R':    tb.Float64Col(pos=53),
        'cv_U_x_R':    tb.Float64Col(pos=54),

        'mu_U_y_R':    tb.Float64Col(pos=55),
        'sigma_U_y_R': tb.Float64Col(pos=56),
        'se_U_y_R':    tb.Float64Col(pos=57),
        'cv_U_y_R':    tb.Float64Col(pos=58),

        'cv_U_x_y_R': tb.Float64Col(pos=59),

        'eta_R_x_y':       tb.Float64Col(pos=60, shape=(J, K)),
        'mu_eta_R_x_y':    tb.Float64Col(pos=61),
        'sigma_eta_R_x_y': tb.Float64Col(pos=62),
        'se_eta_R_x_y':    tb.Float64Col(pos=63),
        'cv_eta_R_x_y':    tb.Float64Col(pos=64),

        'chi2':    tb.Float64Col(pos=65),
        'G':       tb.Float64Col(pos=66),
        'C':       tb.Float64Col(pos=67),
        'V':       tb.Float64Col(pos=68),
        'R_stat':  tb.Float64Col(pos=69),

        'N':       tb.Int32Col(pos=70),

        'N_idotdot': tb.Int32Col(pos=71, shape=I),
        'N_dotjdot': tb.Int32Col(pos=72, shape=J),
        'N_dotdotk': tb.Int32Col(pos=73, shape=K),

        'N_ijdot':   tb.Int32Col(pos=74, shape=(I, J)),
        'N_idotk':   tb.Int32Col(pos=75, shape=(I, K)),
        'N_dotjk':   tb.Int32Col(pos=76, shape=(J, K)),

        'N_ijk':   tb.Int32Col(pos=77, shape=(I, J, K)),
        'E_ijk':   tb.Float64Col(pos=78, shape=(I, J, K)),
        }

    if Rdscrt:
        CT3D_TM["Rbounds"] = tb.Float64Col(pos=79, shape=I)
    else:
        CT3D_TM["Rbounds"] = tb.Float64Col(pos=79, shape=I+1)

    if xdscrt:
        CT3D_TM["xbounds"] = tb.Float64Col(pos=80, shape=J)
    else:
        CT3D_TM["xbounds"] = tb.Float64Col(pos=80, shape=J+1)

    if ydscrt:
        CT3D_TM["ybounds"] = tb.Float64Col(pos=81, shape=K)
    else:
        CT3D_TM["ybounds"] = tb.Float64Col(pos=81, shape=K+1)

    CT3D_TM.update(extracols)

    return CT3D_TM

