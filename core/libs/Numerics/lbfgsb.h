/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/
/*************************************************************************
NEOS, November 1994. (Latest revision June 1996.)
Optimization Technology Center.
Argonne National Laboratory and Northwestern University.

Written by Ciyou Zhu in collaboration with
R.H. Byrd, P. Lu-Chen and J. Nocedal.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.
      
This software is freely available, but we  expect  that  all  publications
describing  work using this software, or all commercial products using it,
quote at least one of the references given below:
    * R. H. Byrd, P. Lu and J. Nocedal.  A Limited  Memory  Algorithm  for
      Bound Constrained Optimization, (1995), SIAM Journal  on  Scientific
      and Statistical Computing , 16, 5, pp. 1190-1208.
    * C. Zhu, R.H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B,
      FORTRAN routines for  large  scale  bound  constrained  optimization
      (1997), ACM Transactions on Mathematical Software,  Vol 23,  Num. 4,
      pp. 550 - 560.
*************************************************************************/

#ifndef _lbfgsb_h
#define _lbfgsb_h

#include "ap.h"

/*-----------------------------------------------
This routines must be defined by you:

void funcgrad(const ap::real_1d_array& x, double& f, ap::real_1d_array& g);
-----------------------------------------------*/

namespace ap
{
/// Class that creates a callback-like connection to a CMTK functional class.
class FunctionAndGradient
{
public:
  /// Evaluate function value and gradient.
  virtual void Evaluate( const ap::real_1d_array& x, double& f, ap::real_1d_array& g ) = 0;
};
} // namespace ap

/*************************************************************************
The  subroutine  minimizes  the  function  F(x) of N arguments with simple
constraints using a quasi-Newton method (LBFGS scheme) which is  optimized
to use a minimum amount of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm (instead  of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

This subroutine uses the FuncGrad subroutine which calculates the value of
the function F and gradient G in point X. The programmer should define the
FuncGrad subroutine by himself.  It should be noted  that  the  subroutine
doesn't need to waste  time for memory allocation of array G, because  the
memory is allocated in calling the  subroutine.  Setting  a  dimension  of
array G each time when calling a subroutine will excessively slow down  an
algorithm.

The programmer could also redefine the LBFGSNewIteration subroutine  which
is called on each new step. The current point X, the function value F  and
the gradient G are passed  into  this  subroutine.  It  is  reasonable  to
redefine the subroutine for better debugging, for  example,  to  visualize
the solution process.

Input parameters:
    N       -   problem dimension. N>0
    M       -   number of  corrections  in  the  BFGS  scheme  of  Hessian
                approximation  update.  Recommended value:  3<=M<=7.   The
                smaller value causes worse convergence,  the  bigger  will
                not  cause  a  considerably  better  convergence, but will
                cause a fall in the performance. M<=N.
    X       -   initial solution approximation.
                Array whose index ranges from 1 to N.
    EpsG    -   positive number which defines a precision of  search.  The
                subroutine finishes its work if the condition ||G|| < EpsG
                is satisfied, where ||.|| means Euclidian norm, G - gradient
                projection onto a feasible set, X - current approximation.
    EpsF    -   positive number which defines a precision of  search.  The
                subroutine  finishes  its  work if on iteration number k+1
                the condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
                is satisfied.
    EpsX    -   positive number which defines a precision of  search.  The
                subroutine  finishes  its  work if on iteration number k+1
                the condition |X(k+1)-X(k)| <= EpsX is satisfied.
    MaxIts  -   maximum number of iterations.
                If MaxIts=0, the number of iterations is unlimited.
    NBD     -   constraint type. If NBD(i) is equal to:
                * 0, X(i) has no constraints,
                * 1, X(i) has only lower boundary,
                * 2, X(i) has both lower and upper boundaries,
                * 3, X(i) has only upper boundary,
                Array whose index ranges from 1 to N.
    L       -   lower boundaries of X(i) variables.
                Array whose index ranges from 1 to N.
    U       -   upper boundaries of X(i) variables.
                Array whose index ranges from 1 to N.

Output parameters:
    X       -   solution approximation.
Array whose index ranges from 1 to N.
    Info    -   a return code:
                    * -2 unknown internal error,
                    * -1 wrong parameters were specified,
                    * 0 interrupted by user,
                    * 1 relative function decreasing is less or equal to EpsF,
                    * 2 step is less or equal to EpsX,
                    * 4 gradient norm is less or equal to EpsG,
                    * 5 number of iterations exceeds MaxIts.

FuncGrad routine description. User-defined.
Input parameters:
    X   -   array whose index ranges from 1 to N.
Output parameters:
    F   -   function value at X.
    G   -   function gradient.
            Array whose index ranges from 1 to N.
The memory for array G has already been allocated in the calling subroutine,
and it isn't necessary to allocate it in the FuncGrad subroutine.

    NEOS, November 1994. (Latest revision June 1996.)
    Optimization Technology Center.
    Argonne National Laboratory and Northwestern University.

    Written by Ciyou Zhu in collaboration with
    R.H. Byrd, P. Lu-Chen and J. Nocedal.
*************************************************************************/

namespace
ap
{
void
lbfgsbminimize
(
  FunctionAndGradient *const functionAndGradient,
  const int& n,
  const int& m,
  ap::real_1d_array& x,
  const double& epsg,
  const double& epsf,
  const double& epsx,
  const int& maxits,
  const ap::integer_1d_array& nbd,
  const ap::real_1d_array& l,
  const ap::real_1d_array& u,
  int& info );
void lbfgsbactive(const int& n,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     ap::real_1d_array& x,
     ap::integer_1d_array& iwhere,
     bool& prjctd,
     bool& cnstnd,
     bool& boxed);
void lbfgsbbmv(const int& m,
     const ap::real_2d_array& sy,
     ap::real_2d_array& wt,
     const int& col,
     const ap::real_1d_array& v,
     ap::real_1d_array& p,
     int& info,
     ap::real_1d_array& workvec);
void lbfgsbcauchy(const int& n,
     const ap::real_1d_array& x,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     const ap::real_1d_array& g,
     ap::integer_1d_array& iorder,
     ap::integer_1d_array& iwhere,
     ap::real_1d_array& t,
     ap::real_1d_array& d,
     ap::real_1d_array& xcp,
     const int& m,
     const ap::real_2d_array& wy,
     const ap::real_2d_array& ws,
     const ap::real_2d_array& sy,
     ap::real_2d_array& wt,
     const double& theta,
     const int& col,
     const int& head,
     ap::real_1d_array& p,
     ap::real_1d_array& c,
     ap::real_1d_array& wbp,
     ap::real_1d_array& v,
     int& nint,
     const ap::real_1d_array& sg,
     const ap::real_1d_array& yg,
     const double& sbgnrm,
     int& info,
     ap::real_1d_array& workvec);
void lbfgsbcmprlb(const int& n,
     const int& m,
     const ap::real_1d_array& x,
     const ap::real_1d_array& g,
     const ap::real_2d_array& ws,
     const ap::real_2d_array& wy,
     const ap::real_2d_array& sy,
     ap::real_2d_array& wt,
     const ap::real_1d_array& z,
     ap::real_1d_array& r,
     ap::real_1d_array& wa,
     const ap::integer_1d_array& index,
     const double& theta,
     const int& col,
     const int& head,
     const int& nfree,
     const bool& cnstnd,
     int& info,
     ap::real_1d_array& workvec,
     ap::real_1d_array& workvec2);
void lbfgsberrclb(const int& n,
     const int& m,
     const double& factr,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     int& task,
     int& info,
     int& k);
void lbfgsbformk(const int& n,
     const int& nsub,
     const ap::integer_1d_array& ind,
     const int& nenter,
     const int& ileave,
     const ap::integer_1d_array& indx2,
     const int& iupdat,
     const bool& updatd,
     ap::real_2d_array& wn,
     ap::real_2d_array& wn1,
     const int& m,
     const ap::real_2d_array& ws,
     const ap::real_2d_array& wy,
     const ap::real_2d_array& sy,
     const double& theta,
     const int& col,
     const int& head,
     int& info,
     ap::real_1d_array& workvec,
     ap::real_2d_array& workmat);
void lbfgsbformt(const int& m,
     ap::real_2d_array& wt,
     const ap::real_2d_array& sy,
     const ap::real_2d_array& ss,
     const int& col,
     const double& theta,
     int& info);
void lbfgsbfreev(const int& n,
     int& nfree,
     ap::integer_1d_array& index,
     int& nenter,
     int& ileave,
     ap::integer_1d_array& indx2,
     const ap::integer_1d_array& iwhere,
     bool& wrk,
     const bool& updatd,
     const bool& cnstnd,
     const int& iter);
void lbfgsbhpsolb(const int& n,
     ap::real_1d_array& t,
     ap::integer_1d_array& iorder,
     const int& iheap);
void lbfgsblnsrlb(const int& n,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     ap::real_1d_array& x,
     const double& f,
     double& fold,
     double& gd,
     double& gdold,
     const ap::real_1d_array& g,
     const ap::real_1d_array& d,
     ap::real_1d_array& r,
     ap::real_1d_array& t,
     const ap::real_1d_array& z,
     double& stp,
     double& dnrm,
     double& dtd,
     double& xstep,
     double& stpmx,
     const int& iter,
     int& ifun,
     int& iback,
     int& nfgv,
     int& info,
     int& task,
     const bool& boxed,
     const bool& cnstnd,
     int& csave,
     ap::integer_1d_array& isave,
     ap::real_1d_array& dsave);
void lbfgsbmatupd(const int& n,
     const int& m,
     ap::real_2d_array& ws,
     ap::real_2d_array& wy,
     ap::real_2d_array& sy,
     ap::real_2d_array& ss,
     const ap::real_1d_array& d,
     const ap::real_1d_array& r,
     int& itail,
     const int& iupdat,
     int& col,
     int& head,
     double& theta,
     const double& rr,
     const double& dr,
     const double& stp,
     const double& dtd);
void lbfgsbprojgr(const int& n,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     const ap::real_1d_array& x,
     const ap::real_1d_array& g,
     double& sbgnrm);
void lbfgsbsubsm(const int& n,
     const int& m,
     const int& nsub,
     const ap::integer_1d_array& ind,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     ap::real_1d_array& x,
     ap::real_1d_array& d,
     const ap::real_2d_array& ws,
     const ap::real_2d_array& wy,
     const double& theta,
     const int& col,
     const int& head,
     int& iword,
     ap::real_1d_array& wv,
     ap::real_2d_array& wn,
     int& info);
void lbfgsbdcsrch(const double& f,
     const double& g,
     double& stp,
     const double& ftol,
     const double& gtol,
     const double& xtol,
     const double& stpmin,
     const double& stpmax,
     int& task,
     ap::integer_1d_array& isave,
     ap::real_1d_array& dsave,
     int& addinfo);
void lbfgsbdcstep(double& stx,
     double& fx,
     double& dx,
     double& sty,
     double& fy,
     double& dy,
     double& stp,
     const double& fp,
     const double& dp,
     bool& brackt,
     const double& stpmin,
     const double& stpmax);
bool additionallbfgsbstoppingcriterion(int iter,
     const ap::real_1d_array& x,
     double f,
     const ap::real_1d_array& g);
bool lbfgsbdpofa(ap::real_2d_array& a, const int& n);
void lbfgsbdtrsl(ap::real_2d_array& t,
     const int& n,
     ap::real_1d_array& b,
     const int& job,
     int& info);
void lbfgsbnewiteration(const ap::real_1d_array& x,
     double f,
     const ap::real_1d_array& g);
} // namespace ap

#endif
