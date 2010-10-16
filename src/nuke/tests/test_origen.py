import os

from nose.tools import assert_equal

from metasci.nuke import origen as msno


def test_write_tape4():
    isovec = {"U235": 0.95, 80160: 0.05}
    msno.write_tape4(isovec, "test.tape4")

    observed_file = open("test.tape4", 'r')
    observed = observed_file.read()

    expected = ("1 80160 5.0000000000E-02   0 0   0 0   0 0\n"
                "2 922350 9.5000000000E-01   0 0   0 0   0 0\n"
                "0 0 0 0\n")

    assert_equal(observed, expected)

    observed_file.close()
    os.remove("test.tape4")


def test_out_table_string1():
    obs = msno.out_table_string(None, None)
    exp = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
    assert_equal(obs, exp)

def test_out_table_string2():
    obs = msno.out_table_string((False, False, True), None)
    exp = "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1"
    assert_equal(obs, exp)

def test_out_table_string3():
    obs = msno.out_table_string((False, False, True), range(1, 25))
    exp = "7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7"
    assert_equal(obs, exp)

def test_out_table_string4():
    obs = msno.out_table_string((False, False, True), [10, 5])
    exp = "8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
    assert_equal(obs, exp)

def test_out_table_string5():
    obs = msno.out_table_string((True, False, True), [10, 5])
    exp = "8 8 8 8 3 8 8 8 8 3 8 8 8 8 8 8 8 8 8 8 8 8 8 8"
    assert_equal(obs, exp)

def test_write_tape5_irradiation():
    msno.write_tape5_irradiation("IRP", 100, 0.550, [204, 205, 206], name="test.tape5",
                                 out_table_nes=(False, False, True), out_table_laf=(True,  False,  True),  out_table_num=[5, 10])

    observed_file = open("test.tape5", 'r')
    observed = observed_file.read()

    expected = ("  -1\n"
                "  -1\n"  
                "  -1\n"
                "  CUT     5 1.000E-10 -1\n"
                "  RDA     FIND CROSS SECTION LIBRARY IDENTIFIER NUMBERS IN YOUR LIBRARY FILE\n"
                "  LIB     0 1 2 3 204 205 206 9 3 0 3 0\n"
                "  OPTL    8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  OPTA    8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  OPTF    8 8 8 8 7 8 8 8 8 7 8 8 8 8 8 8 8 8 8 8 8 8 8 8\n"
                "  INP     1 -1  0  -1  4  4\n"
                "  HED     1     IN FUEL\n"
                "  RDA     ALL IRRADIATION (IRF and IRP) CARDS MUST TAKE PLACE IN BETWEEN BURNUP (BUP) CARDS\n"
                "  BUP\n"
                "  IRP     1.0000000000E+02  5.5000000000E-01   1   2   4  2\n"
                "  BUP\n"
                "  OUT     2  1 1 0\n"
                "  END\n")

    assert_equal(observed, expected)

    observed_file.close()
    os.remove("test.tape5")
