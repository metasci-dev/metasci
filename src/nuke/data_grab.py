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
	    param_dict['nuc'] = nuc
    	params = urlencode(param_dict)

	    while not nuc_fetched:
		    try:
			    print("Grabbing " + nuc)
			    kaeri = urlopen("http://atom.kaeri.re.kr/cgi-bin/nuclide?%s"%params)
			    nuc_fetched = True
		    except:
			    print("Failed to grab; retrying.", end="  ")


	    with open(dir_out + nuc + '.html', 'w') as f:
		    f.write(kaeri.read())
