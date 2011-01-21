"""This module provides a way to grab the raw data from various web sources, since providing a mirror of the original 
data is not acceptable in some cases.  However, the output of these functions is a new expression of the underlying
data and is therefore not subject to the original copyright.  These new datasets are provided here."""

import isoname
import urllib2

from urllib import urlopen
from urllib import urlencode


# Note that since ground state and meta-stable isotopes are of the same atomic weight, 
# the meta-stables have been discluded from the following data sets.

def grab_kaeri_atomic_weights(file_out='atomic_weight.txt'):
    """Makes the atomic weight library.
    Library rows have the the following form:

    iso	AW	AW_sig	Abund

    where:
        iso	= Isotope in LLZZZM format
        AW	= Atomic Weight [amu]
        AW_sig	= Atomic Weight Uncertainty [amu]
        Abund	= Natural fractional atomic abundance [unitless]

    Not to be used under normal circumstances.
    More like an embedded script, in case the librrary file is lost and unrecoverable.

    FIXME: This could use a rewrite such that it doesn't have to grab them all at once.
    """

    isolist = []
    
    for key in isoname.LLaadic.keys():
        NucFetched = False

        while not NucFetched:
            try:
                print key 
                kaeri = urllib2.urlopen( 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=%s'%(key) )
                NucFetched = True
            except:
                print "Failed to grab, retrying",

        for line in kaeri:
            if 0 < line.count("/cgi-bin/nuclide?nuc="):
                nuc = line.partition("/cgi-bin/nuclide?nuc=")[2].partition("\"")[0].upper()
                if not (nuc == key):
                    isolist.append(nuc)
        kaeri.close()

    print "\n~~~~~~~~\n"

    isotab = []

    for key in isolist:
        AW = 0.0
        AW_sig = 0.0
        Abund = 0.0			

        NucFetched = False

        while not NucFetched:
            try:
                print key
                kaeri = urllib2.urlopen( 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=%s'%(key) )
                NucFetched = True
            except:
                print "Failed to grab, retrying",

        for line in kaeri:
            if 0 < line.count("Atomic Mass:"):
                ls = line.split()
                AW = ls[2]
                AW_sig = ls[4]
            elif 0 < line.count("Atomic Percent Abundance:"):
                ls = line.split()
                abund_try = ls[-1]
                while Abund == 0.0:
                    try:
                        Abund = float(abund_try) / 100.0
                    except:
                        abund_try = abund_try[:-1]
        kaeri.close()

        if AW == 0.0:
            continue

        isotab.append([isoname.LLZZZM_2_aazzzm(key), AW, AW_sig, '%G'%Abund])


    isotab = sorted(isotab)

    libfile = open(file_out, 'w')
    for row in isotab:
        new_row = '{0:<6}  {1:<11}  {2:<9}  {3}\n'.format(isoname.aazzzm_2_LLZZZM(row[0]), row[1], row[2], row[3])
        libfile.write(new_row)
    libfile.close()

    return



def grab_kaeri_neutron_xs(nuclist, dir_out='xs_html/'):
    """Grapbs the neutron cross-section summary webpages from the KAERI website.

    Args:
        * nuclist (list): list of nuclides to grab in LLAAAM form.  It is 
          a good idea to pipe the first column from the atomic weight library in 
          as this value.

    Keyword Args:
        * dir_out (str): Path to output directory. 
    """
    param_dict = {'nuc': 'XX', 'n': 2}

    for nuc in nuclist:
        nuc_fetched = False
        aram_dict['nuc'] = nuc
        params = urlencode(param_dict)

        while not nuc_fetched:
            try:
                print "Grabbing " + nuc 
                kaeri = urlopen("http://atom.kaeri.re.kr/cgi-bin/nuclide?%s"%params)
                nuc_fetched = True
            except:
                print "Failed to grab; retrying.",

        with open(dir_out + nuc + '.html', 'w') as f:
            f.write(kaeri.read())


def grab_scattering_lengths(file_out='scattering_lengths.html'):
    """Grapbs the scattering cross-section lengths for neutrons from the NIST website.

    Keyword Args:
        * file_out (str): Path to output file. 
    """
    nist = urlopen("http://www.ncnr.nist.gov/resources/n-lengths/list.html")

    with open(file_out, 'w') as f:
        f.write(nist.read())


# The following is what I used to grab the decay library, 
# but good lord does it need a rewrite! I am too ashamed to 
# expose it in the module, yet I feel that it needs to be 
# under version control.
"""\
import os
import fpformat
import isoname

def delimit(mystr, mylc):
	tempout = []

	newlist = mystr.split(mylc[0])
	for newel in newlist:
		if len(mylc) == 1:
			return newlist
		else:
			tempout.extend( delimit(newel, mylc[1:]) )
	for el in tempout:
		if el == "":
			tempout.remove("")
	return tempout
	
def time2sec(time, unit):
	unit = unit.lower()
	if time == "inf":
		return time + "\t"
	elif unit in ["fs", "femtosec", "femtosecond", "femtoseconds"]:
		return fpformat.sci(float(time)*(10.0**(-15)), 6)
	elif unit in ["ps", "picosec", "picosecond", "picoseconds"]:
		return fpformat.sci(float(time)*(10.0**(-12)), 6)
	elif unit in ["ns", "nanosec", "nanosecond", "nanoseconds"]:
		return fpformat.sci(float(time)*(10.0**(-9)), 6)
	elif unit in ["us", "microsec", "microsecond", "microseconds"]:
		return fpformat.sci(float(time)*(10.0**(-6)), 6)
	elif unit in ["ms", "millisec", "millisecond", "milliseconds"]:
		return fpformat.sci(float(time)*(10.0**(-3)), 6)
	elif unit in ["s", "sec", "second", "seconds"]:
		return fpformat.sci(float(time), 6)
	elif unit in ["m", "min", "minute", "minutes"]:
		return fpformat.sci(float(time)*60.0, 6)
	elif unit in ["h", "hour", "hours"]:
		return fpformat.sci(float(time)*3600.0, 6)
	elif unit in ["d", "day", "days"]:
		return fpformat.sci(float(time)*3600.0*24.0, 6)
	elif unit in ["y", "year", "years"]:
		return fpformat.sci(float(time)*3600.0*24.0*365.0, 6)
	elif unit ==  "mev":
		return fpformat.sci(float(time) * (10.0**(-15)) * float(7.6e-08)/6.03, 6)
	else:
		print "Time given is not a number!"
		print "Time =", time, "Unit = ", unit
		input()
		return "nan"

def write2lib(fromiso, time, unit, toiso, branchratio):
	decay = open("Decay.LIB", 'a')
	timeNsec = time2sec(time, unit)
	decay.write(fromiso + "\t" + timeNsec + "\t" + toiso + "\t" + fpformat.sci(float(branchratio) * 0.01, 6) + "\n")
	decay.close()
	return

try:
	os.remove("Decay.LIB")
except:
	pass

listfile = open("isolist.txt", 'r')

modeofdecaycount = 0
tempbranchratio = "100"

for isoline in listfile:
#for isoline in ["H1", "MN52", "AM242"]:
#for isoline in ["U234"]:
	iso2parse = isoline.split()[0]
	print "Now grabing data for " + iso2parse

	put, get = os.popen2("wget -O iso2parse.txt http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=" + iso2parse)
	for gline in get:
		print gline

	i2p = open("iso2parse.txt", 'r')
	spincount = 0
	metacount = 0
	temptime = tempunit = temptoiso = ""
	tempbranchratio = "100"
	for i2pline in i2p:
		if i2pline.count("Spin") == 1:
			spincount = spincount + 1
			modeofdecaycount = 0
			temptime = tempunit = ""
			temptoiso = ""
			tempbranchratio = "100"
		if i2pline.count("Meta state") == 1:
			metacount = metacount + 1

		if metacount == 0:
			if i2pline.count("Stable Isotope") == 1:
				write2lib(iso2parse, "inf", "s", "None", "100") 
			elif i2pline.count("Half life") == 1:
				i2pls = i2pline.split()
				temptime = i2pls[-2].strip("<(~=?)>").split("(")[0]
				tempunit = i2pls[-1].strip("<(~=?)>")
				print 
			elif i2pline.count("Mode of decay") == 1:
				modeofdecaycount = modeofdecaycount + 1
				tempmodedont = False
				if 0 < i2pline.count("+") or 0 < i2pline.count("-XN")  or (tempbranchratio == "100" and 1 < modeofdecaycount):
					tempmodedont = True
				tempmode = ""
				tempmode = delimit(i2pline, ["<", ">", "/a", "ul", " ", "\n"])[-1]
				if tempmode in ["Alpha", "ECF", "Ne", "Mg"] or tempmodedont:
					temptoiso = iso2parse
				else:
					temptoiso = isoname.aazzzm_2_LLZZZM( isoname.LLZZZM_2_aazzzm(tempmode))
			elif i2pline.count("Branch ratio") == 1:
				tempbranchratio = i2pline.split()[-2].strip("<(~=?)>")
			elif i2pline.count("Decay energy") == 1:
				if tempunit == "Unknown" or tempmode in ["Alpha", "Ne", "Mg"] or tempunit[-2:].upper() == "EV" or tempbranchratio[-2:] == "E-" or tempmodedont:
					pass
				else:
					write2lib(iso2parse, temptime, tempunit, temptoiso, tempbranchratio) 
				temptoiso = ""
#				tempbranchratio = "100"

		elif metacount == 1:
			if i2pline.count("Half life") == 1:
				i2pls = i2pline.split()
				temptime = i2pls[-2].strip("<(~=?)>").split("(")[0]
				tempunit = i2pls[-1].strip("<(~=?)>")
			elif i2pline.count("Mode of decay") == 1:
				modeofdecaycount = modeofdecaycount + 1
				tempmodedont = False
				if 0 < i2pline.count("+") or  0 < i2pline.count("-XN") or (tempbranchratio == "100" and 1 < modeofdecaycount):
					tempmodedont = True
				tempmode = ""
				tempmode = delimit(i2pline, ["<", ">", "/a", "ul", " ", "\n"])[-1]
				if tempmode in ["IT", "SF", "Alpha", "ECF", "Ne", "Mg"] or tempmodedont:
					temptoiso = iso2parse
				else:
					temptoiso = isoname.aazzzm_2_LLZZZM( isoname.LLZZZM_2_aazzzm(tempmode))
			elif i2pline.count("Branch ratio") == 1:
				tempbranchratio = i2pline.split()[-2].strip("<(~=?)>")
			elif i2pline.count("Decay energy") == 1:
				if tempunit == "Unknown" or tempmode in ["Alpha", "Ne", "Mg"] or tempunit[-2:].upper() == "EV" or tempbranchratio[-2:] == "E-" or tempmodedont:
					pass
				else:
					write2lib(iso2parse + "M", temptime, tempunit, temptoiso, tempbranchratio) 
				temptoiso = ""
#				tempbranchratio = "100"

		else:
			pass


listfile.close()
"""
