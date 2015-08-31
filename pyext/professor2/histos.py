# -*- python -*-

class Histo(object):
    "A simple histogram -- just a Bin container with an optional path name"

    def __init__(self, bins=None, path=None):
        self.bins = bins if bins else []
        self.path = path

    @property
    def nbins(self):
        return len(self.bins)

    # TODO: NO!!! Only YODA should write YODA format... or we're back into consistency hell. And anyway look at the mess required to make this work
    def toYODA(self, ppoint=None, manpath=None):
        # TODO: follow var naming convention
        if self.path is None and not manpath is None:
            P=manpath
        elif self.path is not None and manpath is None:
            P=self.path
        elif self.path is not None and manpath is not None:
            P=manpath
        else:
            print "No path given!"
            return ""
        s="# BEGIN YODA_SCATTER2D %s\n"%P
        s+="Path=%s\n"%P
        s+="Type=Scatter2D\n"
        s+="# xval   xerr-   xerr+   yval    yerr-   yerr+\n"
        for b in self.bins:
            if ppoint is not None:
                s+="%e\t%e\t%e\t%e\t%e\t%e\n"%(b.xmid, b.xmid-b.xmin, b.xmax-b.xmid, b.val(ppoint), b.errs(ppoint)[0], b.errs(ppoint)[1])
            else:
                s+="%e\t%e\t%e\t%e\t%e\t%e\n"%(b.xmid, b.xmid-b.xmin, b.xmax-b.xmid, b.val, b.err, b.err)
        s+="# END YODA_SCATTER2D\n"
        return s



class Bin(object):
    "A base class for binned data, handling the common x-edge stuff"

    def __init__(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax

    @property
    def xmid(self):
        return (self.xmin + self.xmax) / 2.0

    @property
    def xedges(self):
        return (self.xmin, self.xmax)

    def __cmp__(self, other):
        return cmp(self.xmin, other.xmin)


class DataBin(Bin):
    "A bin containing a data value and its error(s)"

    def __init__(self, xmin, xmax, val=None, errs=None):
        Bin.__init__(self, xmin, xmax)
        self.val = val
        self._errs = errs

    # TODO: return numerical 0 if _errs is None?
    @property
    def err(self):
        "Get a scalar error value, by averaging if necessary"
        if self._errs is None:
            return 0.0
        elif hasattr(self._errs, "__len__"):
            assert len(self._errs) == 2
            return sum(self._errs) / 2.0
        return self._errs
    @err.setter
    def err(self, e):
        "Set a scalar error value"
        assert not hasattr(self._errs, "__len__")
        self._errs = e

    @property
    def errs(self):
        "Get a pair of error values, by construction if necessary"
        if self._errs is None:
            return (0.0, 0.0)
        elif hasattr(self._errs, "__len__"):
            assert len(self._errs) == 2
            return self._errs
        return (self._errs, self._errs)
    @errs.setter
    def errs(self, e):
        "Set a pair of error values"
        if e is None:
            self._errs = None
        elif hasattr(e, "__len__"):
            assert len(e) == 2
            self._errs = e
        else:
            self._errs = [e,e]


class IpolBin(Bin):
    """
    A bin containing a value interpolation and its error(s)

    TODO:
     * Implement optional limits on the range of return values... for val and errs?
     * Provide ierr and ierrs getter/setter pairs cf. err/errs on DataBin? They can't be averaged, so not sure it makes sense...
     * Allow ipol'd error handling, with wrapped relative error parameterisation as an option?
    """

    def __init__(self, xmin, xmax, ival=None, ierrs=None):
        Bin.__init__(self, xmin, xmax)
        self.ival = ival
        self.ierrs = ierrs

    def val(self, params):
        return self.ival.value(params)

    def err(self, params):
        if self.ierrs is None:
            return 0.0
        elif hasattr(self.ierrs, "__len__"):
            assert len(self.ierrs) == 2
            return (self.ierrs[0].value(params) + self.ierrs[1].value(params))/2.0
        else:
            return self.ierrs.value(params)

    def errs(self, params):
        if self.ierrs is None:
            return (0.0, 0.0)
        elif hasattr(self.ierrs, "__len__"):
            assert len(self.ierrs) == 2
            return (self.ierrs[0].value(params), self.ierrs[1].value(params))
        else:
            e = self.ierrs.value(params)
            return (e, e)