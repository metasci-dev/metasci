"""A Pure Python implementation of an isotopic naming package.

Converts between naming conventions for nuclides:
    * zzaaam (int) is for numerals only (923350).
    * LLAAAM (str) is for letters as well (U-235).
    * MCNP (int) is for numerals without the meta-stable flag (92235), as used in MCNP.
"""

from .. import ReverseDic

LLzz = {
    "AC": 89,  "AL": 13, "AM": 95,  "SB": 51,  "AR": 18,  "AS": 33,  "AT": 85, 
    "BA": 56,  "BK": 97, "BE": 4,   "BI": 83,  "BH": 107, "B":  5,   "BR": 35, 
    "CD": 48,  "CS": 55, "CA": 20,  "CF": 98,  "C":  6,   "CE": 58,  "CL": 17, 
    "CR": 24,  "CO": 27, "CU": 29,  "CM": 96,  "DS": 110, "DB": 105, "DY": 66, 
    "ES": 99,  "ER": 68, "EU": 63,  "FM": 100, "F":  9,   "FR": 87,  "GD": 64, 
    "GA": 31,  "GE": 32, "AU": 79,  "HF": 72,  "HS": 108, "HE": 2,   "HO": 67,
    "H":  1,   "IN": 49, "I":  53,  "IR": 77,  "FE": 26,  "KR": 36,  "LA": 57, 
    "LR": 103, "PB": 82, "LI": 3,   "LU": 71,  "MG": 12,  "MN": 25,  "MT": 109,
    "MD": 101, "HG": 80, "MO": 42,  "ND": 60,  "NE": 10,  "NP": 93,  "NI": 28, 
    "NB": 41,  "N":  7,  "NO": 102, "OS": 76,  "O":  8,   "PD": 46,  "P":  15, 
    "PT": 78,  "PU": 94, "PO": 84,  "K":  19,  "PR": 59,  "PM": 61,  "PA": 91,
    "RN": 86,  "RE": 75, "RH": 45,  "RG": 111, "RB": 37,  "RU": 44,  "RF": 104,
    "SM": 62,  "SC": 21, "SG": 106, "SE": 34,  "SI": 14,  "AG": 47,  "NA": 11, 
    "SR": 38,  "S":  16, "TA": 73,  "TC": 43,  "TE": 52,  "TB": 65,  "TL": 81, 
    "TH": 90,  "TM": 69, "SN": 50,  "TI": 22,  "W":  74,  "U":  92,  "V":  23, 
    "XE": 54,  "YB": 70, "Y":  39,  "ZN": 30,  "ZR": 40,  "RA": 88,
    }

zzLL = ReverseDic(LLzz)

LAN = ["CE", "DY", "ER", "EU", "GD", "HO", "LA", "LU", "ND", "PM", "PR", "SM", 
    "TB", "TM", "YB"]
ACT = ["AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
    "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG"]
TRU = ["NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", 
    "DB", "SG", "BH", "HS", "MT", "DS", "RG"]
MA = ["NP", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", 
    "SG", "BH", "HS", "MT", "DS", "RG"]
FP = []
for key in LLzz.keys():
    if not (key in ACT):
        FP.append(key)

lan = [LLzz[iso] for iso in LAN]
act = [LLzz[iso] for iso in ACT]
tru = [LLzz[iso] for iso in TRU]
ma  = [LLzz[iso] for iso in MA]
fp  = [LLzz[iso] for iso in FP]

class NotANuclide(Exception):
    """Error for value that connat be converted to a nuclide."""
    def __init__(self, nucwas='', nucnow=''):
        """Keyword Args:
            * `nucwas` (str or int): Nuclide that tried to be converted.
            * `nownow` (str or int): Resultant, failed nuclide.
        """
        self.nucwas = nucwas
        self.nucnow = nucnow
        return
    def __str__(self):
        s = "Not a Nuclide!" 
        if self.nucwas:
            s = s + " {0}".format(self.nucwas)
        if self.nucnow:
            s = s + " --> {0}".format(self.nucnow)
        return s

class IndeterminateNuclideForm(Exception):
    """Error for when a value does not have a distinct nuclide form."""
    def __init__(self, nuc):
        """Args:
            * `nuc` (str or int): Nuclide whose form tried to be identified.
        """
        self.nuc = nuc
        return
    def __str__(self):
        s = "Nuclide Form Could Not Be Determined: {0}.\n".format(self.nuc)
        s = s + "Please Ensure Nuclide is of zzaaam, LLAAAM, or MCNP form."
        return s

def CurrentForm(nuc):
    """Returns current form of a nuclide.

    Args:
        * `nuc` (str or int): Nuclide whose form to determine.

    Returns:
        * `form` (str): Current form flag.  May be of the following values:
          'zzaaam', 'LLAAAM', or 'MCNP'.
    """
    nuc = str(nuc)
    nuc = nuc.upper()
    nuc = nuc.strip('-')

    if nuc[0] in list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'):
        return "LLAAAM"
    else:
        if len(nuc) == 7:
            return "zzaaam"
        elif nuc[-1] in list(['23456789']):
            return "MCNP"
        if int(nuc[:-4]) <= int(nuc[-4:-1]) <= int(nuc[:-4]) * 5:
            return "zzaaam"
        elif int(nuc[:-3]) <= int(nuc[-3:]) <= int(nuc[:-3]) * 5:
            return "MCNP"
        else:
            raise IndeterminateNuclideForm(nuc)
         

############################
### LLAAAM_2_* Functions ###
############################

def LLAAAM_2_zzaaam(nuc):
    """Converts nuclide from LLAAAM form to zzaaam form.
    
    Args:
        * `nuc` (str): LLAAAM nuclide.

    Returns:
        * `newnuc` (int): zzaaam nuclide.
    """
    newnuc = ""
    nucstr = nuc.upper()
    nucstr = nucstr.strip('-')

    anum = nucstr.strip('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')
    anum = int(anum)

    if nucstr[-1] == "M":
        newnuc = (10*anum) + 1
    elif nucstr[-1] in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]:
        newnuc = (10*anum)
    else:
        raise NotANuclide(nucstr, newnuc)

    LL = nucstr[:-1].strip('0123456789-')

    if LL in LLzz.keys():
        newnuc = (LLzz[LL]*10000) + newnuc
    else:
        newnuc = "zz{0}".format(newnuc)
        raise NotANuclide(nucstr, newnuc)

    return newnuc

def LLAAAM_2_MCNP(nuc):
    """Converts nuclide from LLAAAM form to MCNP form.

    Args:
        * `nuc` (str): LLAAAM nuclide.

    Returns:
        * `newnuc` (int): MCNP nuclide.
    """
    return zzaaam_2_MCNP( LLAAAM_2_zzaaam(nuc) )

############################
### zzaaam_2_* Functions ###
############################

def zzaaam_2_LLAAAM(nuc):
    """Converts nuclide from azzzm form to LLAAAM form.

    Args:
        * `nuc` (int): zzaaam nuclide.

    Returns:
        * `newnuc` (str): LLAAAM nuclide.
    """
    newnuc = ""
    
    #grab number values
    try:
        znum = nuc/10000
        anum = (nuc%10000)/10
        mnum = nuc%10
    except:
        raise NotANuclide(nuc, newnuc)

    #grab metastable flag
    if 0 == mnum:
        mflag = ""
    elif 0 < mnum :
        mflag = "M"
    else:
        raise NotANuclide(nuc, newnuc)

    #Grab LL and make newnuc
    if znum in zzLL.keys():
        newnuc = "{0}{1}{2}".format(zzLL[znum], anum, mflag)
    else:
        newnuc = "{0}{1}{2}".format('LL', anum, mflag)
        raise NotANuclide(nuc, newnuc)
    return newnuc

def zzaaam_2_MCNP(nuc):
    """Converts nuclide from zzaaam form to MCNP form.

    Args:
        * `nuc` (int): zzaaam nuclide.

    Returns:
        * `newnuc` (int): MCNP nuclide.
    """

    if (nuc%10) == 0:
        newnuc = nuc / 10
    else:
        znum = nuc / 10000
        anum = (nuc/10) - (znum*1000) + 300
        anum = anum + ((nuc%10)*100)
        newnuc = (znum*1000) + anum

    return newnuc    


##########################
### MCNP_2_* Functions ###
##########################

def MCNP_2_zzaaam(nuc):
    """Converts nuclide from MCNP form to zzaaam form.

    Args:
        * `nuc` (int): MCNP nuclide.

    Returns:
        * `newnuc` (int): zzaaam nuclide.
    """

    if (nuc%1000)-400 < 0:
        newnuc = nuc * 10
    else:
        #Please make more general so that more that the first metastable state is returned...
        newnuc = (nuc - 400)*10 + 1

    return newnuc

def MCNP_2_LLAAAM(nuc):
    """Converts nuclide from MCNP form to LLAAAM form.

    Args:
        * `nuc` (int): MCNP nuclide.

    Returns:
        * `newnuc` (str): LLAAAM nuclide.
    """
    return zzaaam_2_LLAAAM( MCNP_2_zzaaam(nuc) )

############################
### mixed_2_*_ Functions ###
############################

def mixed_2_zzaaam(nuc):
    """Converts nuclide of unknown/mixed form to zzaaam form.

    Args:
        * `nuc` (str or int): Nuclide of unknown form.

    Returns:
        * `newnuc` (int): zzaaam nuclide.
    """
    currentform = CurrentForm(nuc)
    if currentform == "zzaaam":
        return nuc
    elif currentform == "LLAAAM":
        return LLAAAM_2_zzaaam(nuc)
    elif currentform == "MCNP":
        return MCNP_2_zzaaam(nuc)
    else:
        raise IndeterminateNuclideForm(nuc)

def mixed_2_LLAAAM(nuc):
    """Converts nuclide from mixed form to LLAAAM form.

    Args:
        * `nuc` (str or int): Nuclide of unknown form.

    Returns:
        * `newnuc` (str): LLAAAM nuclide.
    """
    currentform = CurrentForm(nuc)
    if currentform == "zzaaam":
        return zzaaam_2_LLAAAM(nuc)
    elif currentform == "LLAAAM":
        return nuc
    elif currentform == "MCNP":
        return MCNP_2_LLAAAM(nuc)
    else:
        raise IndeterminateNuclideForm(nuc)

def mixed_2_MCNP(nuc):
    """Converts nuclide from mixed form to MCNP form.

    Args:
        * `nuc` (str or int): Nuclide of unknown form.

    Returns:
        * `newnuc` (int): MCNP nuclide.
    """
    currentform = CurrentForm(nuc)
    if currentform == "zzaaam":
        return zzaaam_2_MCNP(nuc)
    elif currentform == "LLAAAM":
        return LLAAAM_2_MCNP(nuc)
    elif currentform == "MCNP":
        return nuc
    else:
        raise IndeterminateNuclideForm(nuc)

##################################
### (a*)_2_(b*)_List Functions ###
##################################

def LLAAAM_2_zzaaam_List(nuclist):
    """Converts a list of LLAAAM form to a list of zzaaam form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [LLAAAM_2_zzaaam(nuc) for nuc in nuclist]

def LLAAAM_2_MCNP_List(nuclist):
    """Converts a list of LLAAAM form to a list of MCNP form.

    Args:
        * `nuclist` (str list): List of LLAAAM nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [LLAAAM_2_MCNP(nuc) for nuc in nuclist]

def zzaaam_2_LLAAAM_List(nuclist):
    """Converts a list of zzaaam form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [zzaaam_2_LLAAAM(nuc) for nuc in nuclist]

def zzaaam_2_MCNP_List(nuclist):
    """Converts a list of zzaaam form to a list of MCNP form.

    Args:
        * `nuclist` (int list): List of zzaaam nuclides.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return [zzaaam_2_MCNP(nuc) for nuc in nuclist]

def MCNP_2_LLAAAM_List(nuclist):
    """Converts a list of MCNP form to a list of LLAAAM form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return [MCNP_2_LLAAAM(nuc) for nuc in nuclist]

def MCNP_2_zzaaam_List(nuclist):
    """Converts a list of MCNP form to a list of zzaaam form.

    Args:
        * `nuclist` (int list): List of MCNP nuclides.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return [MCNP_2_zzaaam(nuc) for nuc in nuclist]


################################
### mixed_2_*_List Functions ###
################################
def RearRemoveDuplicates(l):
    """Removes duplicate entries from list l, starting from the back.

    Args:
        * `l` (list): List with possible duplicates.

    Returns:
        * `l` (list): List with no duplicates.
    """
    for n in range(len(l)-1, -1, -1):
        if 1 < l.count(l[n]):
            l.pop(n)
    return l
    
def mixed_2_zzaaam_List(nuclist):
    """Converts a list of mixed form to a list of zzaaam form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of zzaaam nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_zzaaam(nuc) for nuc in nuclist] )

def mixed_2_LLAAAM_List(nuclist):
    """Converts a list of mixed form to a list of LLAAAM form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (str list): List of LLAAAM nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_LLAAAM(nuc) for nuc in nuclist] )

def mixed_2_MCNP_List(nuclist):
    """Converts a list of mixed form to a list of MCNP form.

    Args:
        * `nuclist` (str or int list): List of nuclides of mixed form.

    Returns:
        * `newnuclist` (int list): List of MCNP nuclides.
    """
    return RearRemoveDuplicates( [mixed_2_MCNP(nuc) for nuc in nuclist] )


#################################
### isovec_keys_2_* Functions ###
#################################
def isovec_keys_2_zzaaam(isovec):
    """Converts all keys of an isotopic vector dictionary to zzaaam form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of zzaaam (int) form.
    """
    newvec = {}
   
    for iso in isovec.keys():
        newvec[mixed_2_zzaaam(iso)] = isovec[iso]
    
    return newvec

def isovec_keys_2_LLAAAM(isovec):
    """Converts all keys of an isotopic vector dictionary to LLAAAM form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of LLAAAM (str) form.
    """
    newvec = {}
   
    for iso in isovec.keys():
        newvec[mixed_2_LLAAAM(iso)] = isovec[iso]
    
    return newvec

def isovec_keys_2_MCNP(isovec):
    """Converts all keys of an isotopic vector dictionary to MCNP form.

    Args:
        * `isovec` (dict): isotopic vector with keys of unknown/mixed form.

    Returns:
        * `newvec` (dict): isotopic vector with keys of MCNP (int) form.
    """
    newvec = {}
   
    for iso in isovec.keys():
        newvec[mixed_2_MCNP(iso)] = isovec[iso]
    
    return newvec
