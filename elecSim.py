# created by Guang Yang April 2019
#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
import random
import math
from optparse import OptionParser
from array import array

lar_active_vols = [ "vol3DST" , "volCube" , "volcube" ]
neutron_mass = 939.5654133
speedOfLight = 30 #cm per ns


def loop( events, tgeo, tout, nfiles, okruns, GeoNow, doDetSim ):
    if len(okruns) == 0:
        print "There are no runs in this TTree...skipping!"
        return

    print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])    






if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)
    parser.add_option('--geo', help='geometry number', default="Geo4")
    parser.add_option('--doDetSim', help='if do detSim', default=False)


    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_p3pi = array('f',3*[0.0])
    tout.Branch('p3pi',t_p3pi,'p3pi[3]/F')

    loaded = False
    #if os.path.exists("EDepSimEvents/EDepSimEvents.so"):
    #    print "Found EDepSimEvents.so"
    #    ROOT.gSystem.Load("EDepSimEvents/EDepSimEvents.so")
    #    loaded = True

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"
        
    print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    nfiles = 0
    okruns = []
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/%02d/full3DST.%s.%d.edepsim.root" % (args.topdir, run/1000, neutrino, run)
        print fname

        # see if it is an OK file
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            print "File is crap: %s" % fname
            continue
        nfiles += 1
        okruns.append( run )

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    print "OK runs: ", sorted(okruns)
    print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns), args.geo, args.doDetSim )

    fout.cd()
    tout.Write()













