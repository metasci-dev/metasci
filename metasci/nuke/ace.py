"""This module helps grab information out of continuous-energy ACE files (used in both MCNP and Serpent)."""
import re
import linecache

import numpy as np

import isoname


def parse_zaid(zaid):
    """Parses a ZAID into its constituent parts.  This function may also be used to check the
    validity of a ZAID string.

    Args: 
        * zaid (str): A 10 character isotope specificaton string (HZ).  For example, "92235.06c" 
          represents continuous-energy cross-sections for U-235 at 600 K.

    Returns:
        * iso (int): Isotope in zzaaam form.
        * temp (int): Temperature for this ZAID.
        * type (str): The flag at the end of the ZAID that determines the cross-section type.
    """

    zaid_pattern = "(\d{1,6}).(\d{2})([A-Za-z])"

    # Match the pattern
    m = re.match(zaid_pattern, zaid)

    # Check that a match was found
    if m is None:
        raise ValueError("The ZAID '{0}' is not valid.".format(zaid))

    # perform conversions
    iso = isoname.MCNP_2_zzaaam(int(m.group(1)))
    temp = int(m.group(2)) * 100
    type = m.group(3)

    return iso, temp, type


def zaids(ace_file):
    """Grabs the ZAIDs available from an ACE file.

    Args: 
        * ace_file (str): The path the an ACE file that zaid lives in.

    Returns:
        * zaids (dict): A dictioray whose keys are the available ZAIDs and
          whose values are the line numbers on which they start.
    """
    zaids = {}

    current_line_num = 1

    while 0 < current_line_num: 
        # Try to get the ZAID line
        zaid_line = linecache.getline(ace_file, current_line_num)

        if zaid_line == '':
            # break out of loop if line number not available
            current_line_num = -1
            break

        # Get and Verify ZAID
        zaid = zaid_line.split()[0] 
        parse_zaid(zaid)

        # Add this ZAID to dictionary
        zaids[zaid] = current_line_num

        # Increment current line number
        number_of_entries_line = linecache.getline(ace_file, current_line_num - 1 + 7)
        number_of_table_entries = int(number_of_entries_line.split()[0])
        lines_for_table = (number_of_table_entries + 3) / 4 
        lines_for_zaid = lines_for_table + 13

        current_line_num  = current_line_num + lines_for_zaid - 1

    return zaids


def mt(zaid, ace_file):
    """Grabs the MT numbers for a specific isotope out of an ACE file.

    Args: 
        * zaid (str): A 10 character isotope specificaton string (HZ).  For example, "92235.06c" 
          represents continuous-energy cross-sections for U-235 at 600 K.
        * ace_file (str): The path the an ACE file that zaid lives in.

    Returns:
        * mt (set): A set of all valid MT numbers for this zaid.
    """
    # Grab the line number that the ZAIDS of this ACE file start on
    zs = zaids(ace_file)

    # Verify that the ZAID is in the file
    if zaid not in zs:
        raise ValueError("The ZAID {0} is not in the file {1}".format(zaid, ace_file))

    # Get the number of MTs for this ZAID
    number_of_mts_line = linecache.getline(ace_file, zs[zaid] - 1 + 7)
    number_of_mts = int(number_of_mts_line.split()[3])

    # Get the index of the first MT number in the table
    mts_index_line = linecache.getline(ace_file, zs[zaid] - 1 + 9)
    mts_index = int(mts_index_line.split()[2])

    # Get the MT lines
    mts_first_line_num = zs[zaid] - 1 + 13 + ((mts_index - 1) / 4)
    mts_lines = [linecache.getline(ace_file, mts_first_line_num + n) for n in range(number_of_mts/4 + 2)]

    # Convert the lines to a list
    mts = []
    for line in mts_lines:
        mts.extend(line.split())

    offset = (mts_index - 1)%4
    mts = mts[offset:offset + number_of_mts]
    mts = np.array(mts, dtype=int)

    # Covert the list to a set and verify it is the right length
    mt = set(mts)
    assert len(mt) == number_of_mts

    return mt

if __name__ == "__main__":
    ace = "/usr/share/serpent/xsdata/ENDF7/95242ENDF7.ace"
    zaid = "95242.09c"
    
    z = zaids(ace)
    mts = mt(zaid, ace)

    print z 
    print mts
