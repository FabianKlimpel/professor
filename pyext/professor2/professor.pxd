from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "Professor/Professor.h" namespace "Professor":
    string version()

cdef extern from "Professor/Ipol.h" namespace "Professor":
    int numCoeffs(int dim, int order)


cdef extern from "Professor/ParamPoints.h" namespace "Professor":
    cdef cppclass ParamPoints:
        ParamPoints(const vector[vector[double]]&) except +

cdef extern from "Professor/Ipol.h" namespace "Professor":
    cdef cppclass Ipol:

        Ipol(const ParamPoints&, const vector[double]&, int, const string&, double, bool) except +
        Ipol(ParamPoints&, const vector[double]&, const vector[double]&, int&, bool, int&, vector[double]&, const string&) except +
        Ipol(const string&) except +
        #Ipol(ParamPoints&, const vector[double]&, const vector[double]&, int&, const string&) except +
        #Ipol(const ParamPoints&, const vector[double]&, int&, int&, int&, const string&, int&) except +

        string name() except +
        int order() except +
        int dim() except +

        const vector[double]& coeffs() except +
        const double& coeff(size_t) except +

        const vector[vector[int]]& structure() except +
        vector[double] longVector(const vector[double]&) except +

        const ParamPoints& params() except +

        double value(const vector[double]&) except +

        double derivative(const vector[double]&) except +
        const vector[double] gradient(const vector[double]&) except +

        string toString() except +
        string toString(const string&) except +

        const vector[double]& minParamVals() except +
        const vector[double]& maxParamVals() except +
        void setParamLimits(const vector[double]&, const vector[double]&) except +
        void setMinParamVals(const vector[double]&) except +
        void setMaxParamVals(const vector[double]&) except +


cdef extern from "Professor/SimplexInterpolation.h" namespace "Professor":
    cdef cppclass SimplexInterpolation:

        SimplexInterpolation( const ParamPoints& p, const vector[double]& ) except +

#        string name() except +
#        int order() except +
#        int dim() except +

#        const ParamPoints& params() except +

        double value(const vector[double]&) except +

        double derivative(const vector[double]&) except +
        const vector[double] gradient(const vector[double]&) except +

        string toString() except +
        string toString(const string&) except +

#        const vector[double]& minParamVals() except +
#        const vector[double]& maxParamVals() except +
#        void setParamLimits(const vector[double]&, const vector[double]&) except +
#        void setMinParamVals(const vector[double]&) except +
#        void setMaxParamVals(const vector[double]&) except +

