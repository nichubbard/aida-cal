#!/usr/bin/env python3

import datetime
import numpy as np
import matplotlib.pyplot as plt
import progress.bar as p
from progress.spinner import Spinner
import scipy.optimize as opt
import lmfit
from optparse import OptionParser

"""
Options propsed:
    --save-spn : Save Spn calculation for later
    --use-spn  : Used saved Spn calculations

    --use-gain : Used saved intrinsic gains (for absolute gain)
    --absolute=A=E : Gain matched "A" should equal "E" keV
                   :  Default is 1=0.7

    --minimum-energy=N : Minimum energy for Spn prediction (defaults 0 keV)
                       : Can be a tuple for DSSDs

    --dssds=X : X Number of DSSDs
    --wide    : Wide DSSDs
"""

print("AIDA Intrinstic Calibrator")
print(" based on Reese et al. 2015")
print("")

parser = OptionParser()
parser.add_option("--save-spn", dest="save", action="store_true",
        help="Save Spn calculation for later")
parser.add_option("--use-spn", dest="use", action="store_true",
        help="Use saved Spn calculations")
parser.add_option("--use-gains", dest="gains", action="store_true",
        help="Reuse saved instrinsic gains file")
parser.add_option("--absolute", dest="absolute", metavar="A=E",
        help="Scale so gain matched \"A\" should equal \"E\" keV.\
            Default is 1=0.7\
            Can be a tuple for per-DSSD (e.g. (A,B,C)=E")
parser.add_option("--minimum-energy", dest="minimum", metavar="E",
        help="Minimum energy for Spn prediction (default: 0)\
        Can be a tuple for per-DSSD (e.g. (0,0,10000))")
parser.add_option("--dssds", dest="dssds",
        help="Number of DSSDs")
parser.add_option("--wide", dest="wide", action="store_true",
        help="Wide DSSDs")

(options, args) = parser.parse_args()

scale_factor = np.array([0.7])

n_dssds = int(options.dssds)
n_x = 386 if options.wide else 128
n_y = 128

def scale_write_gains(final):
    if (len(scale_factor) > 1):
        for i in range(n_dssds):
            final[i] *= scale_factor[i]
            print("Average gain is {0}".format(np.mean(final[i])))
    else:
        final *= scale_factor
        print("Average gain is {0}".format(np.mean(final)))
    f = open("aida_gains.txt", "w")
    f.write("#AIDA Gains generated on {0}\n".format(datetime.datetime.now()))
    f.write("#!DSSD\n")
    f.write("# DSSD X/Y STRIP GAIN\n")
    it = np.nditer(final, flags=['multi_index'])
    for e in it:
        # print("{0} {1}".format(it.multi_index, e))
        dssd, strip = it.multi_index
        xy = "X" if strip < n_x else "Y"
        strip = strip if strip < n_x else strip - n_x
        f.write("{0} {1} {2} {3}\n".format(dssd, xy, strip, e))



if options.absolute:
    gain = options.absolute.split("=")
    if len(gain) != 2:
        print("--absolute=A=E")
        exit(1)
    a, e = gain
    e = np.array(e, dtype=float)
    if a.startswith("(") and a.endswith(")"):
        a = np.array(a[1:-1].split(","), dtype=float)
        scale_factor = e / a
    else:
        a = np.array(a, dtype=float)
        scale_factor = e / a
    scale_factor = np.array([scale_factor])
    print("Scaling by ", scale_factor)

if options.gains:
    if not options.absolute:
        scale_factor = np.array([0.7])
    f = np.loadtxt("aida_gains.txt", comments="#", dtype=np.dtype("i,1U,i,f"))
    final = np.zeros((n_dssds, n_x + n_y))
    for e in f:
        dssd, xy, strip, en = e
        if xy == "Y": strip = strip + n_x
        final[dssd, strip] = en
    print("Average gain is {0}".format(np.mean(final)))
    scale_write_gains(final)
    exit(0)


Spn_u = np.empty((n_dssds, n_x, n_y))
Spn_v = np.empty((n_dssds, n_x, n_y))

if options.use:
    print("Loading Spns from spns.npy")
    (Spn_u, Spn_v) = np.load("spns.npy")
else:
    data = np.loadtxt("aida_calibration.txt", comments="#")
    print("Loaded calibration data... {0} entries found".format(len(data)))
    print("Checking how it looks")
    N = np.zeros((n_dssds, n_x, n_y))
    Ns = np.zeros((n_dssds, n_x + n_y))

    for e in data:
        dssd, xstrip, xamp, ystrip, yamp = e
        dssd = int(dssd)
        xstrip = int(xstrip)
        ystrip = int(ystrip)
        N[dssd - 1, xstrip, ystrip] = N[dssd - 1, xstrip, ystrip] + 1
        Ns[dssd -1, xstrip] += 1
        Ns[dssd -1, ystrip + n_x] += 1

    print("{0}/{1} pixels have data".format(np.count_nonzero(N), len(N.flatten())))
    print("{0}/{1} strips have data".format(np.count_nonzero(Ns), len(Ns.flatten())))

    def L(Ap, An, Sp):
        return 1/(0.01**2 + (np.log(Ap/An) - np.log(Sp))**2)

    Spns = np.empty((n_dssds, n_x, n_y), dtype='O')
    Priors = np.empty((n_dssds, n_x, n_y), dtype='O')

    it = np.nditer(Spns, flags=['multi_index', 'refs_ok'])
    pb = p.ShadyBar('Creating arrays', max=len(Spns.flatten()))
    for e in it:
        Spns[it.multi_index] = np.linspace(0.6, 1.4, 5000)
        Priors[it.multi_index] = np.ones(5000)
        Priors[it.multi_index] /= sum(Priors[it.multi_index])
        pb.next()
    pb.finish()

    for e in p.ShadyBar('Analysing calibration data').iter(data):
        dssd, xstrip, xamp, ystrip, yamp = e
        dssd = int(dssd)
        xstrip = int(xstrip)
        ystrip = int(ystrip)
        Priors[dssd -1, xstrip, ystrip] *= L(xamp, yamp, Spns[dssd -1, xstrip, ystrip])
        Priors[dssd -1, xstrip, ystrip] /= sum(Priors[dssd - 1, xstrip, ystrip])


    it = np.nditer(Spn_u, flags=['multi_index', 'refs_ok'])
    pb = p.ShadyBar('Calculating measured pixel gains', max=len(Spns.flatten()))
    for e in it:
        idx = it.multi_index
        Spn_u[idx] = np.sum(Priors[idx] * Spns[idx])
        Spn_v[idx] = np.sum((Spn_u[idx] - Spns[idx])**2 * Priors[idx])
        pb.next()
    pb.finish()

if options.save:
    np.save("spns", (Spn_u, Spn_v))
    print("Saved pixel gains to spns.npy")

spinner = None #Spinner("Fitting strip gains ")

def chi2(xin, d):
    spinner.next()
    x = np.insert(xin, 0, 1)
    #x.shape= (n_dssds, n_x + n_y)
    r = np.zeros((n_x,n_y))
    pns = x[:n_x]
    nns = x[n_x:]
    vals = np.outer(nns, 1/pns)
    r = (((Spn_u[d] - vals)/Spn_v[d]))
    return r.flatten()

final = np.empty((n_dssds, n_x + n_y))

for i in range(n_dssds):
    spinner = Spinner("Fitting strip gains for DSSD {0} ".format(i))
    x0 = np.ones(n_x + n_y - 1)
    gains = opt.least_squares(chi2, x0, verbose=1, loss='cauchy', bounds=(0, np.inf), x_scale='jac', args=(i,))
    spinner.finish()
    final[i, 0] = 1
    final[i, 1:] = gains.x

#final.shape = (n_dssds, n_x + n_y)

print("Average gain is {0}".format(np.mean(final)))
final /= np.mean(final)
scale_write_gains(final)

