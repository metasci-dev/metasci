from . import isoname

###################################
### ORIGEN Input Deck Functions ###
###################################

def writeTAPE4(isovec, name="TAPE4.INP"):
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

def writeTAPE9(name_new="TAPE9.INP", name_org="original.tape9", SNG={}, SN2N={}, 
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
