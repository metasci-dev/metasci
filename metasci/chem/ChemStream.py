###############################
### Chemical Form Functions ###
###############################

def ChemicalForm(s):
    """Takes string s of appropriate specification and makes a Chemical Form dictionary out of it."""
    if s == None or s == '':
        return None
    elif type(s) == dict:
        return s
    s = s.upper().replace(' ', '').replace('[', '').replace(']', '')
    chemform = {}

    l = AlphanumericSplit(s)
    for n in range(0, len(l), 2):
        l[n] = l[n].replace('$ACT', str(ACT).replace(' ', '').replace('[', '').replace(']', '').replace('\'', '') )
        l[n] = l[n].replace('$TRU', str(TRU).replace(' ', '').replace('[', '').replace(']', '').replace('\'', '') )
        l[n] = l[n].replace('$MA',  str(MA).replace(' ', '').replace('[', '').replace(']', '').replace('\'', '') )
        l[n] = l[n].replace('$FP',  str(FP).replace(' ', '').replace('[', '').replace(']', '').replace('\'', '') )
        chemform[ l[n] ] = float( l[n+1] )
    return chemform

def AtomsPerMolecule(chemform):
    """Returns the number of atoms per molecule (as a float) of a chemical form.
    The chemical form may be in either string or dictionary format."""
    if type(chemform) == str:
        chemform = ChemicalForm(chemform)

    apm = 0.0
    for key in chemform.keys():
        apm = apm + chemform[key]
    return apm

