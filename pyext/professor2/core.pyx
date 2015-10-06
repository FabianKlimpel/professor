#cython: embedsignature=True

cimport professor as c
cimport cython.operator.dereference as deref


cdef class ParamPoints:
    cdef c.ParamPoints* _ptr

    def __cinit__(self, pvec):
        self._ptr = new c.ParamPoints(pvec)

    def __del__(self):
        del self._ptr


cdef class Ipol:
    """An interpolation of a scalar function built from a list of values across
    a set of parameter point anchors.

    The main workhorse object in Professor. The interpolation coefficients are
    calculated lazily, i.e. when first used.
    """
    cdef c.Ipol* _ptr

    def __cinit__(self, *args):
        if len(args) == 1 and type(args[0]) is str:
            self._ptr = new c.Ipol(args[0])
        else:
            pp = ParamPoints(args[0])
            vals = list(args[1])
            order = int(args[2])
            name = ""
            threshold = 1e-10 #< ???
            if len(args) == 4:
                try:
                    threshold = float(args[3])
                except:
                    name = str(args[3])
            if len(args) == 5:
                    name = str(args[3])
                    threshold = float(args[4])

            self._ptr = new c.Ipol(deref(pp._ptr), vals, order, name, threshold)

    def __del__(self):
        del self._ptr

    @property
    def coeffs(self):
        return self._ptr.coeffs()

    @property
    def dim(self):
        return self._ptr.order()

    @property
    def order(self):
        return self._ptr.order()

    @property
    def name(self):
        return self._ptr.name()

    def value(self, *params, vmin=None, vmax=None):
        """Calculate the value of this interpolation at the given params point,
        forcing return within the range vmin..vmax"""
        if len(params) == 1 and hasattr(params[0], "__iter__"):
            params = params[0]
        v = self._ptr.value(params)
        if vmin is not None and v < vmin:
            return vmin
        if vmax is not None and v > vmax:
            return vmax
        return v

    ## Alias
    val = value

    def toString(self, name=""):
        "Produce a persistent string representing this Ipol object"
        return self._ptr.toString(name)

    def __repr__(self):
        return self.toString(self.name)
