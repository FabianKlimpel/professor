# -*- python -*-

from professor2.histos import *
from professor2.paramsio import *


def read_histos(path, rebin=None):
    "Load histograms from file, into a dict of path -> yoda.Histo[DataBin]"
    histos = {}
    if path.endswith(".root"):
        ## Try ROOT
        try:
            import ROOT
            raise Exception("Can't yet read ROOT format files!")
        except:
            print "Can't yet read ROOT format files... TODO!"
    else:
        ## Try YODA (handles .yoda, .aida, and .flat)
        try:
            import yoda
            s2s = []
            for ao in yoda.read(path, asdict=False):
                if type(ao) == yoda.core.Histo1D and rebin is not None:
                    ao.rebin(rebin)
                s2s.append(ao.mkScatter())
            #s2s = [ao.mkScatter() for ao in yoda.read(path, asdict=False)]
            for s2 in s2s:
                bins = [DataBin(p.xMin, p.xMax, p.y, p.yErrAvg) for p in s2.points]
                #bins = [DataBin(p.xMin, p.xMax, p.y*(p.xMax-p.xMin), p.yErrAvg) for p in s2.points]
                # bins = [DataBin(p.xMin, p.xMax, p.y, p.yErrs) for p in s2.points]
                histos[s2.path] = Histo(bins, s2.path)
        except Exception, e:
            print "Can't load histos from file '%s': %s" % (path, e)
    return histos


def read_rundata(dirs, pfname="params.dat", debug=False, rebin=None): #, formats="yoda,root,aida,flat"):
    params, histos = {}, {}
    import os, glob
    for num, d in enumerate(dirs):
        run = os.path.basename(d)
        files = glob.glob(os.path.join(d, "*"))
        #if debug:
            #print "Reading from %s"%run
        #else:
            #print "\r%.1f per cent read"%((float(num+1)/len(dirs))*100),
        for f in files:
            ## Params file
            if os.path.basename(f) == pfname:
                params[run] = read_paramsfile(f)
            ## Histo file
            else:
                ## Read as a path -> Histo dict
                hs = read_histos(f, rebin=rebin)
                ## Restructure into the path -> run -> Histo return dict
                for path, hist in hs.iteritems():
                    histos.setdefault(path, {})[run] = hist
        ## Check that a params file was found and read in this dir
        if run not in params.keys():
            raise Exception("No params file '%s' found in run dir '%s'" % (pfname, d))
    return params, histos


# TODO: remove this alias
load_rundata = read_rundata