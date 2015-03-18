#include "Ipol.h"
#include <eigen3/Eigen/SVD>

using namespace Eigen;

double Ipol::value(vector <double> P) {
  if (_coeffs.size() == 0 && _values.size()>0) {
    cout << "Coeffs empty, values not, calculating" << endl;
    _calcCoeffs();
  }
  vector<double> LV = getLongVector(getDP(P, _center), _coeffs, order());
  double v = 0.0;
  for (int i=0; i< LV.size();i++) {
    v += LV[i]*_coeffs[i];
  }
  return v;
}

void Ipol::_calcCoeffs() {
  assert(_pts->points().size() == _values.size());
  //if (_center.size() <1) {
    //this->setCenter(this->calcCenter()); //TODO add debug message?
  //}
  //if (m_min.size() <1) {
    //this->m_min = this->calcMin();
  //}
  //if (m_max.size() <1) {
    //this->m_max = this->calcMax();
  //}
  int ncoeff = numOfCoefficients(_center.size(), order());
  if (ncoeff > _pts->points().size()) {
    cout << "Error: not enough ("<< ncoeff <<" vs. " <<_pts->points().size()<< ") anchor points, aborting" <<endl;
    abort();
  }
  MatrixXd DP = MatrixXd(_pts->points().size(), ncoeff);
  VectorXd MC = VectorXd(_pts->points().size());

  vector<double> tempLV;
  vector<double> tempDP;
  // Populate the to be inversed matrix
  for (int a=0;a<_pts->points().size();a++) {
    tempLV = getLongVector(getDP(_pts->points()[a], _center), order());
    for (int i=0;i<tempLV.size();i++) {
      DP(a, i) = tempLV[i];
    }
    // The vector of values (corresponding to anchors)
    MC[a] = _values[a];
  }
  VectorXd co = DP.jacobiSvd(ComputeThinU|ComputeThinV).solve(MC);
  vector<double> temp;
  for (int i=0;i<ncoeff;i++) {
    temp.push_back(co[i]);
  }
/*
  tuple<int, vector<double> > pb(order, temp); // TODO: do we want coeffs more multiple orders
  */
  _coeffs = temp;
}

// Tested and working
int Ipol::numOfCoefficients(int dim, int order) {
    int ntok = 1;
    int r = min(order, dim);
    for(int i=0; i<r;++i) {
      ntok=ntok*(dim+order-i)/(i+1);
    }
  return ntok;
}
//int ProfKeeper::binomial(int n, int k) {
    //int ntok = 1;
    //int r = min(k, n-k);
    //for(int i=0; i<r;++i) {
      //ntok=ntok*(n-i)/(i+1);
    //}
    //return ntok;
//}

vector<double> Ipol::getLongVector(vector<double> p, vector<double> coeffs, int order) {
  if (order < 1 || order > 6) {
    std::cout << "ERROR degree " << order << " not implemented, exiting" << std::endl;
    exit(1);
  }
  if (coeffs.size() != numOfCoefficients(p.size(), order)) {
    std::cout << "ERROR invalid number of coefficients: " << coeffs.size() << " supplied, " << numOfCoefficients(p.size(), order) << " required, exiting" << std::endl;
  }
  if (order == 1) return getLongVector1D(p);
  if (order == 2) return getLongVector2D(p);
  if (order == 3) return getLongVector3D(p);
  if (order == 4) return getLongVector4D(p);
  if (order == 5) return getLongVector5D(p);
  if (order == 6) return getLongVector6D(p);

}

vector<double> Ipol::getLongVector(vector<double> p, int order) {
  if (order < 1 || order > 6) {
    std::cout << "ERROR degree " << order << " not implemented, exiting" << std::endl;
    exit(1);
  }
  if (order == 1) return getLongVector1D(p);
  if (order == 2) return getLongVector2D(p);
  if (order == 3) return getLongVector3D(p);
  if (order == 4) return getLongVector4D(p);
  if (order == 5) return getLongVector5D(p);
  if (order == 6) return getLongVector6D(p);
}

vector<double> Ipol::getDP(vector<double> P, vector<double> C) {
  vector<double> dp;
  for (int i=0; i<P.size();i++) {
    dp.push_back(P[i] - C[i]);
  }
  return dp;
}

vector<double> Ipol::getLongVector1D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }

  assert(retvec.size() == numOfCoefficients(nop,1));
  return retvec;
}

vector<double> Ipol::getLongVector2D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      if (i<=j) {
        retvec.push_back(p[i]*p[j]);
      }
    }
  }

  assert(retvec.size() == numOfCoefficients(nop,2));
  return retvec;
}

vector<double> Ipol::getLongVector3D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      if (i<=j) {
        retvec.push_back(p[i]*p[j]);
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        if (i<=j && i<=k && j<=k) {
          retvec.push_back(p[i]*p[j]*p[k]);
        }
      }
    }
  }

  assert(retvec.size() == numOfCoefficients(nop,3));
  return retvec;
}

vector<double> Ipol::getLongVector4D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      if (i<=j) {
        retvec.push_back(p[i]*p[j]);
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        if (i<=j && i<=k && j<=k) {
          retvec.push_back(p[i]*p[j]*p[k]);
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          if (i<=j && i<=k && i<=l &&
                      j<=k && j<=l && 
                              k<=l) {
            retvec.push_back(p[i]*p[j]*p[k]*p[l]);
          }
        }
      }
    }
  }

  assert(retvec.size() == numOfCoefficients(nop,4));
  return retvec;
}
vector<double> Ipol::getLongVector5D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      if (i<=j) {
        retvec.push_back(p[i]*p[j]);
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        if (i<=j && i<=k && j<=k) {
          retvec.push_back(p[i]*p[j]*p[k]);
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          if (i<=j && i<=k && i<=l && j<=k && j<=l && k<=l) {
            retvec.push_back(p[i]*p[j]*p[k]*p[l]);
          }
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          for (int m=0;m<nop;m++) {
            if (
                i<=j && i<=k && i<=l && i<=m &&
                        j<=k && j<=l && j<=m &&
                                k<=l && k<=m &&
                                        l<=m
               ) {
              retvec.push_back(p[i]*p[j]*p[k]*p[l]*p[m]);
            }
          }
        }
      }
    }
  }

  assert(retvec.size() == numOfCoefficients(nop,5));
  return retvec;
}
vector<double> Ipol::getLongVector6D(vector<double> p) {
  int nop = p.size();
  vector<double> retvec;
  retvec.push_back(1.0);    // This is the offset, for alpha
  for (int i=0;i<nop;i++) { // Linear coefficients, for beta
    retvec.push_back(p[i]);
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      if (i<=j) {
        retvec.push_back(p[i]*p[j]);
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        if (i<=j && i<=k && j<=k) {
          retvec.push_back(p[i]*p[j]*p[k]);
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          if (i<=j && i<=k && i<=l && j<=k && j<=l && k<=l) {
            retvec.push_back(p[i]*p[j]*p[k]*p[l]);
          }
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          for (int m=0;m<nop;m++) {
            if (
                i<=j && i<=k && i<=l && i<=m &&
                        j<=k && j<=l && j<=m &&
                                k<=l && k<=m &&
                                        l<=m
               ) {
              retvec.push_back(p[i]*p[j]*p[k]*p[l]*p[m]);
            }
          }
        }
      }
    }
  }
  for (int i=0;i<nop;i++) {
    for (int j=0;j<nop;j++) {
      for (int k=0;k<nop;k++) {
        for (int l=0;l<nop;l++) {
          for (int m=0;m<nop;m++) {
            for (int n=0;n<nop;n++) {
              if (
                  i<=j && i<=k && i<=l && i<=m && i<=n &&
                          j<=k && j<=l && j<=m && j<=n &&
                                  k<=l && k<=m && k<=n &&
                                          l<=m && l<=n &&
                                                  m<=n
                 ) {
                retvec.push_back(p[i]*p[j]*p[k]*p[l]*p[m]*p[n]);
              }
            }
          }
        }
      }
    }
  }

  assert(retvec.size() == numOfCoefficients(nop,6));
  return retvec;
}