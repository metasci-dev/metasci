from . import isoname

###################################
### ORIGEN Input Deck Functions ###
###################################

def write_tape4(isovec, name="TAPE4.INP"):
    """Writes a TAPE4.INP-type ORIGEN inputfile for an arbitrary isotopic vector.

    Args:
        * `isovec` (dict): Isotopic vector of {iso: mass} form with mass in units of grams.

    Keyword Args:
        * `name` (str): Name of tape 4 file.
    """
    tape4=open(name, 'w')

    for ikey in isovec.keys():
        iso = isoname.mixed_2_zzaaam(ikey)

        #Actinide or Fission Product or Other Identifier
        if  89 < (iso/10000):
            typeflag = "2"
        else:
            typeflag = "1"

        tape4.write("{0} {1} {2:.10E}   0 0   0 0   0 0\n".format(typeflag, iso, isovec[ikey]))

    tape4.write("0 0 0 0\n")
    tape4.close()
    return


    def make_input_tape5(self, t, p = ""):
        """Writes a TAPE5 files to the current directory + p.
        Only use this function from within the libs/ORIGEN/ directory."""

        # Grab important data from library file.
        libfile = tb.openFile('../{0}.h5'.format(reactor), 'r')
        lfr = libfile.root
        n_t = FineTime.index(t)
        IRFtime = FineTime[n_t] - FineTime[n_t-1]
        IRFflux = lfr.Fine.flux[n_t]
        libfile.close()

        # Grab the library numbers used in the original TAPE9 file...
        NLB = []
        tape9_orig = open("../../{0}.tape9".format(reactor), 'r')
        for line in tape9_orig:
            ls = line.split()
            if ls == []:
                continue
            elif (3 < int(ls[0])) and not (ls[0] in NLB):
                NLB.append(ls[0])
        tape9_orig.close()

        # Make template file fill-value dictionary
        tape5_kw = {
            'CUT_OFF': self.cut_off,
            'IRFtime': '{0:.10E}'.format(IRFtime),
            'IRFflux': '{0:.10E}'.format(IRFflux),
            'NLB1':    NLB[0],
            'NLB2':    NLB[1],
            'NLB3':    NLB[2],
            }

        # Fill the template
        with open("../../../templates/{0}.tape5.origen.template".format(reactor), 'r') as f:
            template_file = f.read()

        with open("{0}{1}_T{2}.tape5".format(p, reactor, t), 'w') as f:
            f.write(template_file.format(**tape5_kw))

        return

    def run(self):
        """Runs the ORIGEN Burnup Calculations."""
        os.chdir('libs/ORIGEN/')

        # Grab General Data from the HDF5 File
        libfile = tb.openFile("../{0}.h5".format(reactor), 'r')
        CoreLoadIsos = list(libfile.root.CoreLoad_zzaaam)
        libfile.close()

        if 0 < verbosity:
            print(message("Preping the ORIGEN Directories..."))
        t1 = time.time()
        for t in FineTime[1:]:
            self.make_input_tape5(t)
            self.make_input_tape9(t)
        for iso in CoreLoadIsos: 
            os.mkdir("{0}".format(iso))
        t2 = time.time()
        if 0 < verbosity:
            print(message("...Done!  That only took {0:time} min.\n", "{0:.3G}".format((t2-t1)/60.0) ))

        if 0 < verbosity:
            print(message("  ~~~~~  Starting ORIGEN Runs  ~~~~~  "))
        orit1 = time.time()

        # Initialize the data structures
        self.BU  = {}
        self.k   = {}
        self.Pro = {}
        self.Des = {}
        self.Tij = {}

        for iso in CoreLoadIsos:
            isoLL = isoname.zzaaam_2_LLAAAM(iso)
            if 0 < verbosity:
                print(message("  ~~~~~  Now on {0:iso}  ~~~~~  \n", "Isotope {0}".format(isoLL)))
            isot1 = time.time()

            # Initilize iso data, for t = 0
            self.BU[iso]  = [0.0]
            self.k[iso]   = [0.0]
            self.Pro[iso] = [0.0]
            self.Des[iso] = [0.0]
            self.Tij[iso] = [{iso: 1000.0}]

            for t in FineTime[1:]:
                if 0 < verbosity:
                   print(message("Starting ORIGEN run for {0:iso} at {1:time}...", isoLL, "Time {0}".format(t)))
                t1 = time.time()

                os.chdir("{0}".format(iso))

                # Make/Get Input Decks
                self.make_input_tape4(Tij[iso][-1])
                shutil.copy("../{0}_T{1}.tape5".format(reactor, t), "TAPE5.INP")
                shutil.copy("../{0}_T{1}.tape9".format(reactor, t), "TAPE9.INP")

                # Run ORIGEN
                subprocess.call(self.run_str, shell=True)

                # Parse Output
                self.parse_tape6()
                self.BU[iso].append(  BU[iso][-1] + self.tape6_BU )
                self.k[iso].append(   self.tape6_k )
                self.Pro[iso].append( self.tape6_Pro )
                self.Des[iso].append( self.tape6_Des )
                self.Tij[iso].append( self.tape6_outvec )

                # Clean up the directory
                for f in os.listdir('.'):
                    if f[-4:] in ['.INP', '.OUT']:
                        metasci.SafeRemove(f)
                os.chdir('../') # Back to ORIGEN Directory

                t2 = time.time()
                if 0 < verbosity:
                    print(message("ORIGEN run completed in {0:time} min!", "{0:.3G} min".format((t2-t1)/60.0) ))
    
            isot2 = time.time()
            if 0 < verbosity:
                print(message("  ~~~~~  Isotope {0:iso} took {1:time} min!  ~~~~~  \n", isoLL, "{0:.3G} min".format((isot2-isot1)/60.0) ))


        # Kludge to put Tij in the right units and form
        allORIGENisoList = []
        for iso in CoreLoadIsos:
            for t in Tij[iso]:
                for j in t.keys():
                    if (j not in allORIGENisoList):
                        allORIGENisoList.append(j)
        for iso in CoreLoadIsos:
            for n_t in range(len(Tij[iso])):
                for j in allORIGENisoList:
                    if j in Tij[iso][n_t].keys():
                        Tij[iso][n_t][j] = Tij[iso][n_t][j] / (10.0**3)
                    else:
                        Tij[iso][n_t][j] = 0.0
    
        orit2 = time.time()
        if 0 < verbosity:
            print(message("  ~~~~~  ORIGEN took {0:time} to run!  ~~~~~  ", "{0:.3G} min".format((orit2-orit1)/60.0) ))

        os.chdir('../../') #Back to 'reactor' root
        return 

    def parse_tape6(self, p = ""):
        """Parses an ORIGEN TAPE6.OUT file that is in the current directory + path p."""
        InTable5 = False

        # (Re-)Initialize data structures
        if hasattr(self, "tape6_BU"):
            del self.tape6_BU, self.tape6_k, self.tape6_Pro, self.tape6_Des, self.tape6_outvec
        self.tape6_BU     = -1.0
        self.tape6_k      = -1.0
        self.tape6_Pro    = -1.0
        self.tape6_Des    = -1.0
        self.tape6_outvec = {}

        tape6 = open("%sTAPE6.OUT"%p, 'r')
        for line in tape6:
            if "BURNUP,MWD" in line:
                ls = line.split()
                self.tape6_BU = float(ls[-1])
                continue
            elif "K INFINITY" in line:
                ls = line.split()
                self.tape6_k = float(ls[-1])
                continue
            elif "NEUT PRODN" in line:
                ls = line.split()
                self.tape6_Pro = float(ls[-1])
                continue
            elif "NEUT DESTN" in line:
                ls = line.split()
                self.tape6_Des = float(ls[-1])
                continue
            elif "5 SUMMARY TABLE:  CONCENTRATIONS, GRAMS" in line:
                InTable5 = True
                continue
            elif InTable5 and ("OUTPUT UNIT =  6" in line):
                InTable5 = False
                continue
            elif InTable5:
                ls = line.split()
                try:
                    iso = isoname.LLAAAM_2_zzaaam(ls[0])
                except:
                    try:
                        iso = isoname.LLAAAM_2_zzaaam(ls[0] + ls[1])
                    except:
                        continue
                self.tape6_outvec[iso] = float(ls[-1])
            else:
                continue
        tape6.close()
        return 




def write_tape9(name_new="TAPE9.INP", name_org="original.tape9", SNG={}, SN2N={}, 
    SN3N={}, SNF={}, SNA={}, SNP={}, SNGX={}, SN2NX={}, iso_list=[], sigfig=3):
    """Writes a new ORIGEN TAPE9.INP file based on an original TAPE9 template 
    file and given cross section values.

    Keyword Args:
    	* `name_new` (str): Path to the new tape9 file.
    	* `name_org` (str): Path to the tape9 template file.
        * `SNG` (dict): (n, gamma) cross-section isotopic vector, of form {iso: xs_(n,g) [barns]}.
        * `SN2N` (dict): (n, 2n) cross-section isotopic vector, of form {iso: xs_(n,2n) [barns]}.
        * `SN3N` (dict): (n, 3n) cross-section isotopic vector, of form {iso: xs_(n,3n) [barns]}.
        * `SNF` (dict): (n, fission) cross-section isotopic vector, of form {iso: xs_(n,f) [barns]}.
        * `SNA` (dict): (n, alpha) cross-section isotopic vector, of form {iso: xs_(n,alpha) [barns]}.
        * `SNP` (dict): (n, proton) cross-section isotopic vector, of form {iso: xs_(n,p) [barns]}.
        * `SNGX` (dict): (n, gamma*) excited state cross-section isotopic vector, of form {iso: xs_(n,g*) [barns]}.
        * `SN2NX` (dict): (n, 2n*) excited state cross-section isotopic vector, of form {iso: xs_(n,2n*) [barns]}.
        * `iso_list` (list): List of isotopes to write over.
        * `sigfig` (int): Ensures that all overwritten data is given to this many digits beyond the 
          decimal point via "{xs:.{sigfig}E}".format().
    """

    #perform some basic set-up
    zero_space = "0.0       "

    SNG   = isoname.isovec_keys_2_zzaaam(SNG)
    SN2N  = isoname.isovec_keys_2_zzaaam(SN2N)
    SN3N  = isoname.isovec_keys_2_zzaaam(SN3N)
    SNF   = isoname.isovec_keys_2_zzaaam(SNF)
    SNA   = isoname.isovec_keys_2_zzaaam(SNA)
    SNP   = isoname.isovec_keys_2_zzaaam(SNP)
    SNGX  = isoname.isovec_keys_2_zzaaam(SNGX)
    SN2NX = isoname.isovec_keys_2_zzaaam(SN2NX)

    if iso_list == []:
        iso_set = set()
        iso_set.update(SNG.keys())
        iso_set.update(SN2N.keys())
        iso_set.update(SN3N.keys())
        iso_set.update(SNF.keys())
        iso_set.update(SNA.keys())
        iso_set.update(SNP.keys())
        iso_set.update(SNGX.keys())
        iso_set.update(SN2NX.keys())
        iso_list = list(iso_set)
    else:
        iso_list = isoname.mixed_2_zzaaam_List(iso_list)

    #Open the first tape9 for reading and the new tape9 for writing.
    tape9_o = open(name_org, 'r')   #Original File
    tape9_n = open(name_new, 'w')   #New File

    #Write the new file...
    for line in tape9_o:
        ls = line.split()
        if ls == []:
            #Rewrites blank lines
            tape9_n.write(line)
            continue
        elif int(ls[0]) <= 3:
            #Rewrites Decay data and '-1' spacer lines
            tape9_n.write(line)
            continue

        #Rewrites Title Lines
        try:
            int(ls[1])
            iso = isoname.mixed_2_zzaaam(ls[1])
        except ValueError:
            tape9_n.write(line)
            continue

        #If we don't care about the isotope, just rewrite the line...
        if iso not in iso_list:
            tape9_n.write(line)
            continue

        #Fixes messed up exponential data in the original file...
        orig_data = []
        SkipNext = False
        for n in range(len(ls)):
            if SkipNext:
                SkipNext = False
            elif 'E' in ls[n]:
                try:
                    orig_data.append(float(ls[n]))
                except ValueError:
                    orig_data.append(float(ls[n] + ls[n+1]))
                    SkipNext = True
            else:
                orig_data.append(float(ls[n]))
        ls = orig_data	#This is what ls was suppossed to look like! Stupid ORIGEN...			 

        newline = line[:13]

        #(n, gamma) XS
        if iso in SNG.keys(): 
            if SNG[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SNG[iso], sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[2], sigfig)

        #(n, 2n) XS
        if iso in SN2N.keys():
            if SN2N[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SN2N[iso], sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[3], sigfig)

        #Check if in actinide set, because positional arguments mean different things.
        if (iso/10000) in isoname.act:
            #(n, 3n) XS
            if iso in SN3N.keys():
                if SN3N[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SN3N[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[4], sigfig)

            #(n, fission) XS
            if iso in SNF.keys():
                if SNF[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNF[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[5], sigfig)

        #If not an actinide, then...
        else:
            #(n, alpha) XS
            if iso in SNA.keys():
                if SNA[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNA[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[4], sigfig)

            #(n, proton) XS
            if iso in SNP.keys():
                if SNP[iso] == 0.0:
                    newline = newline + zero_space
                else:
                    newline = newline + "{0:.{1}E} ".format(SNP[iso], sigfig)
            else:
                newline = newline + "{0:.{1}E} ".format(ls[5], sigfig)

        #(n, g*) XS
        if iso in SNGX.keys(): 
            if SNGX[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SNGX[iso], sigfig)
        elif ls[2] == 0.0:
            newline = newline + zero_space
        elif iso in SNG.keys():
            sngx = SNG[iso] * ls[6] / ls[2]
            if sngx == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(sngx, sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[6], sigfig)

        #(n, 2n*) XS
        if iso in SN2NX.keys(): 
            if SN2NX[iso] == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(SN2NX[iso], sigfig)
        elif ls[3] == 0.0:
            newline = newline + zero_space
        elif iso in SN2N.keys():
            sn2nx = SN2N[iso] * ls[7] / ls[3]
            if sn2nx == 0.0:
                newline = newline + zero_space
            else:
                newline = newline + "{0:.{1}E} ".format(sn2nx, sigfig)
        else:
            newline = newline + "{0:.{1}E} ".format(ls[7], sigfig)

        #End the line
        newline = newline + line[-8:]
        tape9_n.write(newline)

    #Close out the files and return
    tape9_o.close()
    tape9_n.close()
    return
