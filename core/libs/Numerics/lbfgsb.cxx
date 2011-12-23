/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
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

#include "lbfgsb.h"

namespace
ap
{

void lbfgsbactive
( const int& n,
  const ap::real_1d_array& l,
  const ap::real_1d_array& u,
  const ap::integer_1d_array& nbd,
  ap::real_1d_array& x,
  ap::integer_1d_array& iwhere,
  bool& prjctd,
  bool& cnstnd,
  bool& boxed)
{
    int nbdd;
    int i;

    nbdd = 0;
    prjctd = false;
    cnstnd = false;
    boxed = true;
    for(i = 1; i <= n; i++)
    {
        if( nbd(i)>0 )
        {
            if( nbd(i)<=2&&x(i)<=l(i) )
            {
                if( x(i)<l(i) )
                {
                    prjctd = true;
                    x(i) = l(i);
                }
                nbdd = nbdd+1;
            }
            else
            {
                if( nbd(i)>=2&&x(i)>=u(i) )
                {
                    if( x(i)>u(i) )
                    {
                        prjctd = true;
                        x(i) = u(i);
                    }
                    nbdd = nbdd+1;
                }
            }
        }
    }
    for(i = 1; i <= n; i++)
    {
        if( nbd(i)!=2 )
        {
            boxed = false;
        }
        if( nbd(i)==0 )
        {
            iwhere(i) = -1;
        }
        else
        {
            cnstnd = true;
            if( nbd(i)==2&&u(i)-l(i)<=0 )
            {
                iwhere(i) = 3;
            }
            else
            {
                iwhere(i) = 0;
            }
        }
    }
}


void lbfgsbbmv(const int&,
     const ap::real_2d_array& sy,
     ap::real_2d_array& wt,
     const int& col,
     const ap::real_1d_array& v,
     ap::real_1d_array& p,
     int& info,
     ap::real_1d_array& workvec)
{
    int i;
    int k;
    int i2;
    ap::real_value_type s;

    if( col==0 )
    {
        return;
    }
    p(col+1) = v(col+1);
    for(i = 2; i <= col; i++)
    {
        i2 = col+i;
        s = 0.0;
        for(k = 1; k <= i-1; k++)
        {
            s = s+sy(i,k)*v(k)/sy(k,k);
        }
        p(i2) = v(i2)+s;
    }
    ap::vmove(workvec.getvector(1, col), p.getvector(col+1, col+col));
    lbfgsbdtrsl(wt, col, workvec, 11, info);
    ap::vmove(p.getvector(col+1, col+col), workvec.getvector(1, col));
    if( info!=0 )
    {
        return;
    }
    for(i = 1; i <= col; i++)
    {
        p(i) = v(i)/sqrt(sy(i,i));
    }
    ap::vmove(workvec.getvector(1, col), p.getvector(col+1, col+col));
    lbfgsbdtrsl(wt, col, workvec, 1, info);
    ap::vmove(p.getvector(col+1, col+col), workvec.getvector(1, col));
    if( info!=0 )
    {
        return;
    }
    for(i = 1; i <= col; i++)
    {
        p(i) = -p(i)/sqrt(sy(i,i));
    }
    for(i = 1; i <= col; i++)
    {
        s = 0;
        for(k = i+1; k <= col; k++)
        {
            s = s+sy(k,i)*p(col+k)/sy(i,i);
        }
        p(i) = p(i)+s;
    }
}


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
     const ap::real_value_type& theta,
     const int& col,
     const int& head,
     ap::real_1d_array& p,
     ap::real_1d_array& c,
     ap::real_1d_array& wbp,
     ap::real_1d_array& v,
     int& nint,
     const ap::real_1d_array&,
     const ap::real_1d_array&,
     const ap::real_value_type& sbgnrm,
     int& info,
     ap::real_1d_array& workvec)
{
    bool xlower;
    bool xupper;
    bool bnded;
    int i;
    int j;
    int col2;
    int nfree;
    int nbreak;
    int pointr;
    int ibp;
    int nleft;
    int ibkmin;
    int iter;
    ap::real_value_type f1;
    ap::real_value_type f2;
    ap::real_value_type dt;
    ap::real_value_type dtm;
    ap::real_value_type tsum;
    ap::real_value_type dibp;
    ap::real_value_type zibp;
    ap::real_value_type dibp2;
    ap::real_value_type bkmin;
    ap::real_value_type tu = 0;
    ap::real_value_type tl = 0;
    ap::real_value_type wmc;
    ap::real_value_type wmp;
    ap::real_value_type wmw;
    ap::real_value_type tj;
    ap::real_value_type tj0;
    ap::real_value_type neggi;
    ap::real_value_type f2org;
    ap::real_value_type tmpv;

    if( sbgnrm<=0 )
    {
        ap::vmove(xcp.getvector(1, n), x.getvector(1, n));
        return;
    }
    bnded = true;
    nfree = n+1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0;
    col2 = 2*col;
    f1 = 0;
    for(i = 1; i <= col2; i++)
    {
        p(i) = 0;
    }
    for(i = 1; i <= n; i++)
    {
        neggi = -g(i);
        if( iwhere(i)!=3&&iwhere(i)!=-1 )
        {
            tl = 0;
            tu = 0;
            if( nbd(i)<=2 )
            {
                tl = x(i)-l(i);
            }
            if( nbd(i)>=2 )
            {
                tu = u(i)-x(i);
            }
            xlower = nbd(i)<=2&&tl<=0;
            xupper = nbd(i)>=2&&tu<=0;
            iwhere(i) = 0;
            if( xlower )
            {
                if( neggi<=0 )
                {
                    iwhere(i) = 1;
                }
            }
            else
            {
                if( xupper )
                {
                    if( neggi>=0 )
                    {
                        iwhere(i) = 2;
                    }
                }
                else
                {
                    if( fabs(neggi)<=0 )
                    {
                        iwhere(i) = -3;
                    }
                }
            }
        }
        pointr = head;
        if( iwhere(i)!=0&&iwhere(i)!=-1 )
        {
            d(i) = 0;
        }
        else
        {
            d(i) = neggi;
            f1 = f1-neggi*neggi;
            for(j = 1; j <= col; j++)
            {
                p(j) = p(j)+wy(i,pointr)*neggi;
                p(col+j) = p(col+j)+ws(i,pointr)*neggi;
                pointr = pointr%m+1;
            }
            if( nbd(i)<=2&&nbd(i)!=0&&neggi<0 )
            {
                nbreak = nbreak+1;
                iorder(nbreak) = i;
                t(nbreak) = tl/(-neggi);
                if( nbreak==1||t(nbreak)<bkmin )
                {
                    bkmin = t(nbreak);
                    ibkmin = nbreak;
                }
            }
            else
            {
                if( nbd(i)>=2&&neggi>0 )
                {
                    nbreak = nbreak+1;
                    iorder(nbreak) = i;
                    t(nbreak) = tu/neggi;
                    if( nbreak==1||t(nbreak)<bkmin )
                    {
                        bkmin = t(nbreak);
                        ibkmin = nbreak;
                    }
                }
                else
                {
                    nfree = nfree-1;
                    iorder(nfree) = i;
                    if( fabs(neggi)>0 )
                    {
                        bnded = false;
                    }
                }
            }
        }
    }
    if( theta!=1 )
    {
        ap::vmul(p.getvector(col+1, col+col), theta);
    }
    ap::vmove(xcp.getvector(1, n), x.getvector(1, n));
    if( nbreak==0&&nfree==n+1 )
    {
        return;
    }
    for(j = 1; j <= col2; j++)
    {
        c(j) = 0;
    }
    f2 = -theta*f1;
    f2org = f2;
    if( col>0 )
    {
        lbfgsbbmv(m, sy, wt, col, p, v, info, workvec);
        if( info!=0 )
        {
            return;
        }
        tmpv = ap::vdotproduct(v.getvector(1, col2), p.getvector(1, col2));
        f2 = f2-tmpv;
    }
    dtm = -f1/f2;
    tsum = 0;
    nint = 1;
    if( nbreak!=0 )
    {
        nleft = nbreak;
        iter = 1;
        tj = 0;
        while(true)
        {
            tj0 = tj;
            if( iter==1 )
            {
                tj = bkmin;
                ibp = iorder(ibkmin);
            }
            else
            {
                if( iter==2 )
                {
                    if( ibkmin!=nbreak )
                    {
                        t(ibkmin) = t(nbreak);
                        iorder(ibkmin) = iorder(nbreak);
                    }
                }
                lbfgsbhpsolb(nleft, t, iorder, iter-2);
                tj = t(nleft);
                ibp = iorder(nleft);
            }
            dt = tj-tj0;
            if( dtm<dt )
            {
                break;
            }
            tsum = tsum+dt;
            nleft = nleft-1;
            iter = iter+1;
            dibp = d(ibp);
            d(ibp) = 0;
            if( dibp>0 )
            {
                zibp = u(ibp)-x(ibp);
                xcp(ibp) = u(ibp);
                iwhere(ibp) = 2;
            }
            else
            {
                zibp = l(ibp)-x(ibp);
                xcp(ibp) = l(ibp);
                iwhere(ibp) = 1;
            }
            if( nleft==0&&nbreak==n )
            {
                dtm = dt;
                if( col>0 )
                {
                    ap::vadd(c.getvector(1, col2), p.getvector(1, col2), dtm);
                }
                return;
            }
            nint = nint+1;
            dibp2 = ap::sqr(dibp);
            f1 = f1+dt*f2+dibp2-theta*dibp*zibp;
            f2 = f2-theta*dibp2;
            if( col>0 )
            {
                ap::vadd(c.getvector(1, col2), p.getvector(1, col2), dt);
                pointr = head;
                for(j = 1; j <= col; j++)
                {
                    wbp(j) = wy(ibp,pointr);
                    wbp(col+j) = theta*ws(ibp,pointr);
                    pointr = pointr%m+1;
                }
                lbfgsbbmv(m, sy, wt, col, wbp, v, info, workvec);
                if( info!=0 )
                {
                    return;
                }
                wmc = ap::vdotproduct(c.getvector(1, col2), v.getvector(1, col2));
                wmp = ap::vdotproduct(p.getvector(1, col2), v.getvector(1, col2));
                wmw = ap::vdotproduct(wbp.getvector(1, col2), v.getvector(1, col2));
                ap::vsub(p.getvector(1, col2), wbp.getvector(1, col2), dibp);
                f1 = f1+dibp*wmc;
                f2 = f2+2.0*dibp*wmp-dibp2*wmw;
            }
            f2 = ap::maxreal(ap::machineepsilon*f2org, f2);
            if( nleft>0 )
            {
                dtm = -f1/f2;
                continue;
            }
            else
            {
                if( bnded )
                {
                    f1 = 0;
                    f2 = 0;
                    dtm = 0;
                }
                else
                {
                    dtm = -f1/f2;
                }
            }
            break;
        }
    }
    if( dtm<=0 )
    {
        dtm = 0;
    }
    tsum = tsum+dtm;
    ap::vadd(xcp.getvector(1, n), d.getvector(1, n), tsum);
    if( col>0 )
    {
        ap::vadd(c.getvector(1, col2), p.getvector(1, col2), dtm);
    }
}


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
     const ap::real_value_type& theta,
     const int& col,
     const int& head,
     const int& nfree,
     const bool& cnstnd,
     int& info,
     ap::real_1d_array& workvec,
     ap::real_1d_array& workvec2)
{
    int i;
    int j;
    int k;
    int pointr;
    ap::real_value_type a1;
    ap::real_value_type a2;

    if( !cnstnd&&col>0 )
    {
        for(i = 1; i <= n; i++)
        {
            r(i) = -g(i);
        }
    }
    else
    {
        for(i = 1; i <= nfree; i++)
        {
            k = index(i);
            r(i) = -theta*(z(k)-x(k))-g(k);
        }
        ap::vmove(workvec2.getvector(1, 2*m), wa.getvector(2*m+1, 4*m));
        lbfgsbbmv(m, sy, wt, col, workvec2, wa, info, workvec);
        ap::vmove(wa.getvector(2*m+1, 4*m), workvec2.getvector(1, 2*m));
        if( info!=0 )
        {
            info = -8;
            return;
        }
        pointr = head;
        for(j = 1; j <= col; j++)
        {
            a1 = wa(j);
            a2 = theta*wa(col+j);
            for(i = 1; i <= nfree; i++)
            {
                k = index(i);
                r(i) = r(i)+wy(k,pointr)*a1+ws(k,pointr)*a2;
            }
            pointr = pointr%m+1;
        }
    }
}


void lbfgsberrclb(const int& n,
     const int& m,
     const ap::real_value_type& factr,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     int& task,
     int& info,
     int& k)
{
    int i;

    if( n<=0 )
    {
        task = 2;
    }
    if( m<=0 )
    {
        task = 2;
    }
    if( m>n )
    {
        task = 2;
    }
    if( factr<0 )
    {
        task = 2;
    }
    for(i = 1; i <= n; i++)
    {
        if( nbd(i)<0||nbd(i)>3 )
        {
            task = 2;
            info = -6;
            k = i;
        }
        if( nbd(i)==2 )
        {
            if( l(i)>u(i) )
            {
                task = 2;
                info = -7;
                k = i;
            }
        }
    }
}


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
     const ap::real_value_type& theta,
     const int& col,
     const int& head,
     int& info,
     ap::real_1d_array& workvec,
     ap::real_2d_array& workmat)
{
    int ipntr;
    int jpntr;
    int iy;
    int iis;
    int jy;
    int js;
    int is1;
    int js1;
    int k1;
    int i;
    int k;
    int col2;
    int pbegin;
    int pend;
    int dbegin;
    int dend;
    int upcl;
    ap::real_value_type temp1;
    ap::real_value_type temp2;
    ap::real_value_type temp3;
    ap::real_value_type temp4;
    ap::real_value_type v;
    int j;

    if( updatd )
    {
        if( iupdat>m )
        {
            for(jy = 1; jy <= m-1; jy++)
            {
                js = m+jy;
                ap::vmove(wn1.getcolumn(jy, jy, m-1), wn1.getcolumn(jy+1, jy+1, m));
                ap::vmove(wn1.getcolumn(js, js, js+m-jy-1), wn1.getcolumn(js+1, js+1, js+m-jy));
                ap::vmove(wn1.getcolumn(jy, m+1, m+m-1), wn1.getcolumn(jy+1, m+2, m+m));
            }
        }
        pbegin = 1;
        pend = nsub;
        dbegin = nsub+1;
        dend = n;
        iy = col;
        iis = m+col;
        ipntr = head+col-1;
        if( ipntr>m )
        {
            ipntr = ipntr-m;
        }
        jpntr = head;
        for(jy = 1; jy <= col; jy++)
        {
            js = m+jy;
            temp1 = 0;
            temp2 = 0;
            temp3 = 0;
            for(k = pbegin; k <= pend; k++)
            {
                k1 = ind(k);
                temp1 = temp1+wy(k1,ipntr)*wy(k1,jpntr);
            }
            for(k = dbegin; k <= dend; k++)
            {
                k1 = ind(k);
                temp2 = temp2+ws(k1,ipntr)*ws(k1,jpntr);
                temp3 = temp3+ws(k1,ipntr)*wy(k1,jpntr);
            }
            wn1(iy,jy) = temp1;
            wn1(iis,js) = temp2;
            wn1(iis,jy) = temp3;
            jpntr = jpntr%m+1;
        }
        jy = col;
        jpntr = head+col-1;
        if( jpntr>m )
        {
            jpntr = jpntr-m;
        }
        ipntr = head;
        for(i = 1; i <= col; i++)
        {
            iis = m+i;
            temp3 = 0;
            for(k = pbegin; k <= pend; k++)
            {
                k1 = ind(k);
                temp3 = temp3+ws(k1,ipntr)*wy(k1,jpntr);
            }
            ipntr = ipntr%m+1;
            wn1(iis,jy) = temp3;
        }
        upcl = col-1;
    }
    else
    {
        upcl = col;
    }
    ipntr = head;
    for(iy = 1; iy <= upcl; iy++)
    {
        iis = m+iy;
        jpntr = head;
        for(jy = 1; jy <= iy; jy++)
        {
            js = m+jy;
            temp1 = 0;
            temp2 = 0;
            temp3 = 0;
            temp4 = 0;
            for(k = 1; k <= nenter; k++)
            {
                k1 = indx2(k);
                temp1 = temp1+wy(k1,ipntr)*wy(k1,jpntr);
                temp2 = temp2+ws(k1,ipntr)*ws(k1,jpntr);
            }
            for(k = ileave; k <= n; k++)
            {
                k1 = indx2(k);
                temp3 = temp3+wy(k1,ipntr)*wy(k1,jpntr);
                temp4 = temp4+ws(k1,ipntr)*ws(k1,jpntr);
            }
            wn1(iy,jy) = wn1(iy,jy)+temp1-temp3;
            wn1(iis,js) = wn1(iis,js)-temp2+temp4;
            jpntr = jpntr%m+1;
        }
        ipntr = ipntr%m+1;
    }
    ipntr = head;
    for(iis = m+1; iis <= m+upcl; iis++)
    {
        jpntr = head;
        for(jy = 1; jy <= upcl; jy++)
        {
            temp1 = 0;
            temp3 = 0;
            for(k = 1; k <= nenter; k++)
            {
                k1 = indx2(k);
                temp1 = temp1+ws(k1,ipntr)*wy(k1,jpntr);
            }
            for(k = ileave; k <= n; k++)
            {
                k1 = indx2(k);
                temp3 = temp3+ws(k1,ipntr)*wy(k1,jpntr);
            }
            if( iis<=jy+m )
            {
                wn1(iis,jy) = wn1(iis,jy)+temp1-temp3;
            }
            else
            {
                wn1(iis,jy) = wn1(iis,jy)-temp1+temp3;
            }
            jpntr = jpntr%m+1;
        }
        ipntr = ipntr%m+1;
    }
    for(iy = 1; iy <= col; iy++)
    {
        iis = col+iy;
        is1 = m+iy;
        for(jy = 1; jy <= iy; jy++)
        {
            js = col+jy;
            js1 = m+jy;
            wn(jy,iy) = wn1(iy,jy)/theta;
            wn(js,iis) = wn1(is1,js1)*theta;
        }
        for(jy = 1; jy <= iy-1; jy++)
        {
            wn(jy,iis) = -wn1(is1,jy);
        }
        for(jy = iy; jy <= col; jy++)
        {
            wn(jy,iis) = wn1(is1,jy);
        }
        wn(iy,iy) = wn(iy,iy)+sy(iy,iy);
    }
    info = 0;
    if( !lbfgsbdpofa(wn, col) )
    {
        info = -1;
        return;
    }
    col2 = 2*col;
    for(js = col+1; js <= col2; js++)
    {
        ap::vmove(workvec.getvector(1, col), wn.getcolumn(js, 1, col));
        lbfgsbdtrsl(wn, col, workvec, 11, info);
        ap::vmove(wn.getcolumn(js, 1, col), workvec.getvector(1, col));
    }
    for(iis = col+1; iis <= col2; iis++)
    {
        for(js = iis; js <= col2; js++)
        {
            v = ap::vdotproduct(wn.getcolumn(iis, 1, col), wn.getcolumn(js, 1, col));
            wn(iis,js) = wn(iis,js)+v;
        }
    }
    for(j = 1; j <= col; j++)
    {
        ap::vmove(workmat.getrow(j, 1, col), wn.getrow(col+j, col+1, col+col));
    }
    info = 0;
    if( !lbfgsbdpofa(workmat, col) )
    {
        info = -2;
        return;
    }
    for(j = 1; j <= col; j++)
    {
        ap::vmove(wn.getrow(col+j, col+1, col+col), workmat.getrow(j, 1, col));
    }
}


void lbfgsbformt(const int&,
     ap::real_2d_array& wt,
     const ap::real_2d_array& sy,
     const ap::real_2d_array& ss,
     const int& col,
     const ap::real_value_type& theta,
     int& info)
{
    int i;
    int j;
    int k;
    int k1;
    ap::real_value_type ddum;

    for(j = 1; j <= col; j++)
    {
        wt(1,j) = theta*ss(1,j);
    }
    for(i = 2; i <= col; i++)
    {
        for(j = i; j <= col; j++)
        {
            k1 = ap::minint(i, j)-1;
            ddum = 0;
            for(k = 1; k <= k1; k++)
            {
                ddum = ddum+sy(i,k)*sy(j,k)/sy(k,k);
            }
            wt(i,j) = ddum+theta*ss(i,j);
        }
    }
    info = 0;
    if( !lbfgsbdpofa(wt, col) )
    {
        info = -3;
    }
}


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
     const int& iter)
{
    int iact;
    int i;
    int k;

    nenter = 0;
    ileave = n+1;
    if( iter>0&&cnstnd )
    {
        for(i = 1; i <= nfree; i++)
        {
            k = index(i);
            if( iwhere(k)>0 )
            {
                ileave = ileave-1;
                indx2(ileave) = k;
            }
        }
        for(i = 1+nfree; i <= n; i++)
        {
            k = index(i);
            if( iwhere(k)<=0 )
            {
                nenter = nenter+1;
                indx2(nenter) = k;
            }
        }
    }
    wrk = ileave<n+1||nenter>0||updatd;
    nfree = 0;
    iact = n+1;
    for(i = 1; i <= n; i++)
    {
        if( iwhere(i)<=0 )
        {
            nfree = nfree+1;
            index(nfree) = i;
        }
        else
        {
            iact = iact-1;
            index(iact) = i;
        }
    }
}


void lbfgsbhpsolb(const int& n,
     ap::real_1d_array& t,
     ap::integer_1d_array& iorder,
     const int& iheap)
{
    int i;
    int j;
    int k;
    int indxin;
    int indxou;
    ap::real_value_type ddum;
    ap::real_value_type dout;

    if( iheap==0 )
    {
        for(k = 2; k <= n; k++)
        {
            ddum = t(k);
            indxin = iorder(k);
            i = k;
            while(true)
            {
                if( i>1 )
                {
                    j = i/2;
                    if( ddum<t(j) )
                    {
                        t(i) = t(j);
                        iorder(i) = iorder(j);
                        i = j;
                        continue;
                    }
                }
                break;
            }
            t(i) = ddum;
            iorder(i) = indxin;
        }
    }
    if( n>1 )
    {
        i = 1;
        dout = t(1);
        indxou = iorder(1);
        ddum = t(n);
        indxin = iorder(n);
        while(true)
        {
            j = i+i;
            if( j<=n-1 )
            {
                if( t(j+1)<t(j) )
                {
                    j = j+1;
                }
                if( t(j)<ddum )
                {
                    t(i) = t(j);
                    iorder(i) = iorder(j);
                    i = j;
                    continue;
                }
            }
            break;
        }
        t(i) = ddum;
        iorder(i) = indxin;
        t(n) = dout;
        iorder(n) = indxou;
    }
}


void lbfgsblnsrlb(const int& n,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     ap::real_1d_array& x,
     const ap::real_value_type& f,
     ap::real_value_type& fold,
     ap::real_value_type& gd,
     ap::real_value_type& gdold,
     const ap::real_1d_array& g,
     const ap::real_1d_array& d,
     ap::real_1d_array& r,
     ap::real_1d_array& t,
     const ap::real_1d_array& z,
     ap::real_value_type& stp,
     ap::real_value_type& dnrm,
     ap::real_value_type& dtd,
     ap::real_value_type& xstep,
     ap::real_value_type& stpmx,
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
     ap::real_1d_array& dsave)
{
    int i;
    ap::real_value_type a1;
    ap::real_value_type a2;
    ap::real_value_type v;
    ap::real_value_type ftol;
    ap::real_value_type gtol;
    ap::real_value_type xtol;
    ap::real_value_type big;
    int addinfo;

    addinfo = 0;
    big = 1.0E10;
    ftol = 1.0E-3;
    gtol = 0.9E0;
    xtol = 0.1E0;
    if( task!=1 )
    {
        v = ap::vdotproduct(d.getvector(1, n), d.getvector(1, n));
        dtd = v;
        dnrm = sqrt(dtd);
        stpmx = big;
        if( cnstnd )
        {
            if( iter==0 )
            {
                stpmx = 1;
            }
            else
            {
                for(i = 1; i <= n; i++)
                {
                    a1 = d(i);
                    if( nbd(i)!=0 )
                    {
                        if( a1<0&&nbd(i)<=2 )
                        {
                            a2 = l(i)-x(i);
                            if( a2>=0 )
                            {
                                stpmx = 0;
                            }
                            else
                            {
                                if( a1*stpmx<a2 )
                                {
                                    stpmx = a2/a1;
                                }
                            }
                        }
                        else
                        {
                            if( a1>0&&nbd(i)>=2 )
                            {
                                a2 = u(i)-x(i);
                                if( a2<=0 )
                                {
                                    stpmx = 0;
                                }
                                else
                                {
                                    if( a1*stpmx>a2 )
                                    {
                                        stpmx = a2/a1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if( iter==0&&!boxed )
        {
            stp = ap::minreal(1/dnrm, stpmx);
        }
        else
        {
            stp = 1;
        }
        ap::vmove(t.getvector(1, n), x.getvector(1, n));
        ap::vmove(r.getvector(1, n), g.getvector(1, n));
        fold = f;
        ifun = 0;
        iback = 0;
        csave = 0;
    }
    v = ap::vdotproduct(g.getvector(1, n), d.getvector(1, n));
    gd = v;
    if( ifun==0 )
    {
        gdold = gd;
        if( gd>=0 )
        {
            info = -4;
            return;
        }
    }
    lbfgsbdcsrch(f, gd, stp, ftol, gtol, xtol, ap::real_value_type(0), stpmx, csave, isave, dsave, addinfo);
    xstep = stp*dnrm;
    if( csave!=4&&csave!=3 )
    {
        task = 1;
        ifun = ifun+1;
        nfgv = nfgv+1;
        iback = ifun-1;
        if( stp==1 )
        {
            ap::vmove(x.getvector(1, n), z.getvector(1, n));
        }
        else
        {
            for(i = 1; i <= n; i++)
            {
                x(i) = stp*d(i)+t(i);
            }
        }
    }
    else
    {
        task = 5;
    }
}


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
     ap::real_value_type& theta,
     const ap::real_value_type& rr,
     const ap::real_value_type& dr,
     const ap::real_value_type& stp,
     const ap::real_value_type& dtd)
{
    int j;
    int pointr;
    ap::real_value_type v;

    if( iupdat<=m )
    {
        col = iupdat;
        itail = (head+iupdat-2)%m+1;
    }
    else
    {
        itail = itail%m+1;
        head = head%m+1;
    }
    ap::vmove(ws.getcolumn(itail, 1, n), d.getvector(1, n));
    ap::vmove(wy.getcolumn(itail, 1, n), r.getvector(1, n));
    theta = rr/dr;
    if( iupdat>m )
    {
        for(j = 1; j <= col-1; j++)
        {
            ap::vmove(ss.getcolumn(j, 1, j), ss.getcolumn(j+1, 2, j+1));
            ap::vmove(sy.getcolumn(j, j, col-1), sy.getcolumn(j+1, j+1, col));
        }
    }
    pointr = head;
    for(j = 1; j <= col-1; j++)
    {
        v = ap::vdotproduct(d.getvector(1, n), wy.getcolumn(pointr, 1, n));
        sy(col,j) = v;
        v = ap::vdotproduct(ws.getcolumn(pointr, 1, n), d.getvector(1, n));
        ss(j,col) = v;
        pointr = pointr%m+1;
    }
    if( stp==1 )
    {
        ss(col,col) = dtd;
    }
    else
    {
        ss(col,col) = stp*stp*dtd;
    }
    sy(col,col) = dr;
}


void lbfgsbprojgr(const int& n,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     const ap::integer_1d_array& nbd,
     const ap::real_1d_array& x,
     const ap::real_1d_array& g,
     ap::real_value_type& sbgnrm)
{
    int i;
    ap::real_value_type gi;

    sbgnrm = 0;
    for(i = 1; i <= n; i++)
    {
        gi = g(i);
        if( nbd(i)!=0 )
        {
            if( gi<0 )
            {
                if( nbd(i)>=2 )
                {
                    gi = ap::maxreal(x(i)-u(i), gi);
                }
            }
            else
            {
                if( nbd(i)<=2 )
                {
                    gi = ap::minreal(x(i)-l(i), gi);
                }
            }
        }
        sbgnrm = ap::maxreal(sbgnrm, fabs(gi));
    }
}


void lbfgsbsubsm(const int&,
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
     const ap::real_value_type& theta,
     const int& col,
     const int& head,
     int& iword,
     ap::real_1d_array& wv,
     ap::real_2d_array& wn,
     int& info)
{
    int pointr;
    int col2;
    int ibd = 0;
    int jy;
    int js;
    int i;
    int j;
    int k;
    ap::real_value_type alpha;
    ap::real_value_type dk;
    ap::real_value_type temp1;
    ap::real_value_type temp2;

    if( nsub<=0 )
    {
        return;
    }
    pointr = head;
    for(i = 1; i <= col; i++)
    {
        temp1 = 0;
        temp2 = 0;
        for(j = 1; j <= nsub; j++)
        {
            k = ind(j);
            temp1 = temp1+wy(k,pointr)*d(j);
            temp2 = temp2+ws(k,pointr)*d(j);
        }
        wv(i) = temp1;
        wv(col+i) = theta*temp2;
        pointr = pointr%m+1;
    }
    col2 = 2*col;
    lbfgsbdtrsl(wn, col2, wv, 11, info);
    if( info!=0 )
    {
        return;
    }
    for(i = 1; i <= col; i++)
    {
        wv(i) = -wv(i);
    }
    lbfgsbdtrsl(wn, col2, wv, 1, info);
    if( info!=0 )
    {
        return;
    }
    pointr = head;
    for(jy = 1; jy <= col; jy++)
    {
        js = col+jy;
        for(i = 1; i <= nsub; i++)
        {
            k = ind(i);
            d(i) = d(i)+wy(k,pointr)*wv(jy)/theta+ws(k,pointr)*wv(js);
        }
        pointr = pointr%m+1;
    }
    for(i = 1; i <= nsub; i++)
    {
        d(i) = d(i)/theta;
    }
    alpha = 1;
    temp1 = alpha;
    for(i = 1; i <= nsub; i++)
    {
        k = ind(i);
        dk = d(i);
        if( nbd(k)!=0 )
        {
            if( dk<0&&nbd(k)<=2 )
            {
                temp2 = l(k)-x(k);
                if( temp2>=0 )
                {
                    temp1 = 0;
                }
                else
                {
                    if( dk*alpha<temp2 )
                    {
                        temp1 = temp2/dk;
                    }
                }
            }
            else
            {
                if( dk>0&&nbd(k)>=2 )
                {
                    temp2 = u(k)-x(k);
                    if( temp2<=0 )
                    {
                        temp1 = 0;
                    }
                    else
                    {
                        if( dk*alpha>temp2 )
                        {
                            temp1 = temp2/dk;
                        }
                    }
                }
            }
            if( temp1<alpha )
            {
                alpha = temp1;
                ibd = i;
            }
        }
    }
    if( alpha<1 )
    {
        dk = d(ibd);
        k = ind(ibd);
        if( dk>0 )
        {
            x(k) = u(k);
            d(ibd) = 0;
        }
        else
        {
            if( dk<0 )
            {
                x(k) = l(k);
                d(ibd) = 0;
            }
        }
    }
    for(i = 1; i <= nsub; i++)
    {
        k = ind(i);
        x(k) = x(k)+alpha*d(i);
    }
    if( alpha<1 )
    {
        iword = 1;
    }
    else
    {
        iword = 0;
    }
}


void lbfgsbdcsrch(const ap::real_value_type& f,
     const ap::real_value_type& g,
     ap::real_value_type& stp,
     const ap::real_value_type& ftol,
     const ap::real_value_type& gtol,
     const ap::real_value_type& xtol,
     const ap::real_value_type& stpmin,
     const ap::real_value_type& stpmax,
     int& task,
     ap::integer_1d_array& isave,
     ap::real_1d_array& dsave,
     int& addinfo)
{
    bool brackt;
    int stage;
    ap::real_value_type finit;
    ap::real_value_type ftest;
    ap::real_value_type fm;
    ap::real_value_type fx;
    ap::real_value_type fxm;
    ap::real_value_type fy;
    ap::real_value_type fym;
    ap::real_value_type ginit;
    ap::real_value_type gtest;
    ap::real_value_type gm;
    ap::real_value_type gx;
    ap::real_value_type gxm;
    ap::real_value_type gy;
    ap::real_value_type gym;
    ap::real_value_type stx;
    ap::real_value_type sty;
    ap::real_value_type stmin;
    ap::real_value_type stmax;
    ap::real_value_type width;
    ap::real_value_type width1;
    ap::real_value_type xtrapl;
    ap::real_value_type xtrapu;

    xtrapl = 1.1E0;
    xtrapu = 4.0E0;
    while(true)
    {
        if( task==0 )
        {
            if( stp<stpmin )
            {
                task = 2;
                addinfo = 0;
            }
            if( stp>stpmax )
            {
                task = 2;
                addinfo = 0;
            }
            if( g>=0 )
            {
                task = 2;
                addinfo = 0;
            }
            if( ftol<0 )
            {
                task = 2;
                addinfo = 0;
            }
            if( gtol<0 )
            {
                task = 2;
                addinfo = 0;
            }
            if( xtol<0 )
            {
                task = 2;
                addinfo = 0;
            }
            if( stpmin<0 )
            {
                task = 2;
                addinfo = 0;
            }
            if( stpmax<stpmin )
            {
                task = 2;
                addinfo = 0;
            }
            if( task==2 )
            {
                return;
            }
            brackt = false;
            stage = 1;
            finit = f;
            ginit = g;
            gtest = ftol*ginit;
            width = stpmax-stpmin;
            width1 = width/0.5;
            stx = 0;
            fx = finit;
            gx = ginit;
            sty = 0;
            fy = finit;
            gy = ginit;
            stmin = 0;
            stmax = stp+xtrapu*stp;
            task = 1;
            break;
        }
        else
        {
            if( isave(1)==1 )
            {
                brackt = true;
            }
            else
            {
                brackt = false;
            }
            stage = isave(2);
            ginit = dsave(1);
            gtest = dsave(2);
            gx = dsave(3);
            gy = dsave(4);
            finit = dsave(5);
            fx = dsave(6);
            fy = dsave(7);
            stx = dsave(8);
            sty = dsave(9);
            stmin = dsave(10);
            stmax = dsave(11);
            width = dsave(12);
            width1 = dsave(13);
        }
        ftest = finit+stp*gtest;
        if( stage==1&&f<=ftest&&g>=0 )
        {
            stage = 2;
        }
        if( brackt&&(stp<=stmin||stp>=stmax) )
        {
            task = 3;
            addinfo = 1;
        }
        if( brackt&&stmax-stmin<=xtol*stmax )
        {
            task = 3;
            addinfo = 2;
        }
        if( stp==stpmax&&f<=ftest&&g<=gtest )
        {
            task = 3;
            addinfo = 3;
        }
        if( stp==stpmin&&(f>ftest||g>=gtest) )
        {
            task = 3;
            addinfo = 4;
        }
        if( f<=ftest&&fabs(g)<=gtol*(-ginit) )
        {
            task = 4;
            addinfo = -1;
        }
        if( task==3||task==4 )
        {
            break;
        }
        if( stage==1&&f<=fx&&f>ftest )
        {
            fm = f-stp*gtest;
            fxm = fx-stx*gtest;
            fym = fy-sty*gtest;
            gm = g-gtest;
            gxm = gx-gtest;
            gym = gy-gtest;
            lbfgsbdcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
            fx = fxm+stx*gtest;
            fy = fym+sty*gtest;
            gx = gxm+gtest;
            gy = gym+gtest;
        }
        else
        {
            lbfgsbdcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
        }
        if( brackt )
        {
            if( fabs(sty-stx)>=0.66*width1 )
            {
                stp = stx+0.5*(sty-stx);
            }
            width1 = width;
            width = fabs(sty-stx);
        }
        if( brackt )
        {
            stmin = ap::minreal(stx, sty);
            stmax = ap::maxreal(stx, sty);
        }
        else
        {
            stmin = stp+xtrapl*(stp-stx);
            stmax = stp+xtrapu*(stp-stx);
        }
        stp = ap::maxreal(stp, stpmin);
        stp = ap::minreal(stp, stpmax);
        if( (brackt&&(stp<=stmin||stp>=stmax))||(brackt&&stmax-stmin<=xtol*stmax) )
        {
            stp = stx;
        }
        task = 1;
        break;
    }
    if( brackt )
    {
        isave(1) = 1;
    }
    else
    {
        isave(1) = 0;
    }
    isave(2) = stage;
    dsave(1) = ginit;
    dsave(2) = gtest;
    dsave(3) = gx;
    dsave(4) = gy;
    dsave(5) = finit;
    dsave(6) = fx;
    dsave(7) = fy;
    dsave(8) = stx;
    dsave(9) = sty;
    dsave(10) = stmin;
    dsave(11) = stmax;
    dsave(12) = width;
    dsave(13) = width1;
}


void lbfgsbdcstep(ap::real_value_type& stx,
     ap::real_value_type& fx,
     ap::real_value_type& dx,
     ap::real_value_type& sty,
     ap::real_value_type& fy,
     ap::real_value_type& dy,
     ap::real_value_type& stp,
     const ap::real_value_type& fp,
     const ap::real_value_type& dp,
     bool& brackt,
     const ap::real_value_type& stpmin,
     const ap::real_value_type& stpmax)
{
    ap::real_value_type gamma;
    ap::real_value_type p;
    ap::real_value_type q;
    ap::real_value_type r;
    ap::real_value_type s;
    ap::real_value_type sgnd;
    ap::real_value_type stpc;
    ap::real_value_type stpf;
    ap::real_value_type stpq;
    ap::real_value_type theta;

    sgnd = dp*(dx/fabs(dx));
    if( fp>fx )
    {
        theta = 3*(fx-fp)/(stp-stx)+dx+dp;
        s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
        gamma = s*sqrt(ap::sqr(theta/s)-dx/s*(dp/s));
        if( stp<stx )
        {
            gamma = -gamma;
        }
        p = gamma-dx+theta;
        q = gamma-dx+gamma+dp;
        r = p/q;
        stpc = stx+r*(stp-stx);
        stpq = stx+dx/((fx-fp)/(stp-stx)+dx)/2*(stp-stx);
        if( fabs(stpc-stx)<fabs(stpq-stx) )
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc+(stpq-stpc)/2;
        }
        brackt = true;
    }
    else
    {
        if( sgnd<0 )
        {
            theta = 3*(fx-fp)/(stp-stx)+dx+dp;
            s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
            gamma = s*sqrt(ap::sqr(theta/s)-dx/s*(dp/s));
            if( stp>stx )
            {
                gamma = -gamma;
            }
            p = gamma-dp+theta;
            q = gamma-dp+gamma+dx;
            r = p/q;
            stpc = stp+r*(stx-stp);
            stpq = stp+dp/(dp-dx)*(stx-stp);
            if( fabs(stpc-stp)>fabs(stpq-stp) )
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            brackt = true;
        }
        else
        {
            if( fabs(dp)<fabs(dx) )
            {
                theta = 3*(fx-fp)/(stp-stx)+dx+dp;
                s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dx), fabs(dp)));
                gamma = s*sqrt(ap::maxreal(ap::real_value_type(0), ap::sqr(theta/s)-dx/s*(dp/s)));
                if( stp>stx )
                {
                    gamma = -gamma;
                }
                p = gamma-dp+theta;
                q = gamma+(dx-dp)+gamma;
                r = p/q;
                if( r<0&&gamma!=0 )
                {
                    stpc = stp+r*(stx-stp);
                }
                else
                {
                    if( stp>stx )
                    {
                        stpc = stpmax;
                    }
                    else
                    {
                        stpc = stpmin;
                    }
                }
                stpq = stp+dp/(dp-dx)*(stx-stp);
                if( brackt )
                {
                    if( fabs(stpc-stp)<fabs(stpq-stp) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                    if( stp>stx )
                    {
                        stpf = ap::minreal(stp+0.66*(sty-stp), stpf);
                    }
                    else
                    {
                        stpf = ap::maxreal(stp+0.66*(sty-stp), stpf);
                    }
                }
                else
                {
                    if( fabs(stpc-stp)>fabs(stpq-stp) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                    stpf = ap::minreal(stpmax, stpf);
                    stpf = ap::maxreal(stpmin, stpf);
                }
            }
            else
            {
                if( brackt )
                {
                    theta = 3*(fp-fy)/(sty-stp)+dy+dp;
                    s = ap::maxreal(fabs(theta), ap::maxreal(fabs(dy), fabs(dp)));
                    gamma = s*sqrt(ap::sqr(theta/s)-dy/s*(dp/s));
                    if( stp>sty )
                    {
                        gamma = -gamma;
                    }
                    p = gamma-dp+theta;
                    q = gamma-dp+gamma+dy;
                    r = p/q;
                    stpc = stp+r*(sty-stp);
                    stpf = stpc;
                }
                else
                {
                    if( stp>stx )
                    {
                        stpf = stpmax;
                    }
                    else
                    {
                        stpf = stpmin;
                    }
                }
            }
        }
    }
    if( fp>fx )
    {
        sty = stp;
        fy = fp;
        dy = dp;
    }
    else
    {
        if( sgnd<0 )
        {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    stp = stpf;
}


bool additionallbfgsbstoppingcriterion(int,
     const ap::real_1d_array&,
     ap::real_value_type,
     const ap::real_1d_array&)
{
    bool result;

    result = false;
    return result;
}


bool lbfgsbdpofa(ap::real_2d_array& a, const int& n)
{
    bool result;
    ap::real_value_type t;
    ap::real_value_type s;
    ap::real_value_type v;
    int j;
    int jm1;
    int k;

    for(j = 1; j <= n; j++)
    {
        s = 0.0;
        jm1 = j-1;
        if( jm1>=1 )
        {
            for(k = 1; k <= jm1; k++)
            {
                v = ap::vdotproduct(a.getcolumn(k, 1, k-1), a.getcolumn(j, 1, k-1));
                t = a(k,j)-v;
                t = t/a(k,k);
                a(k,j) = t;
                s = s+t*t;
            }
        }
        s = a(j,j)-s;
        if( s<=0.0 )
        {
            result = false;
            return result;
        }
        a(j,j) = sqrt(s);
    }
    result = true;
    return result;
}


void lbfgsbdtrsl(ap::real_2d_array& t,
     const int& n,
     ap::real_1d_array& b,
     const int& job,
     int& info)
{
    ap::real_value_type temp;
    ap::real_value_type v;
    int cse;
    int j;
    int jj;

    for(j = 1; j <= n; j++)
    {
        if( t(j,j)==0.0 )
        {
            info = j;
            return;
        }
    }
    info = 0;
    cse = 1;
    if( job%10!=0 )
    {
        cse = 2;
    }
    if( job%100/10!=0 )
    {
        cse = cse+2;
    }
    if( cse==1 )
    {
        b(1) = b(1)/t(1,1);
        if( n<2 )
        {
            return;
        }
        for(j = 2; j <= n; j++)
        {
            temp = -b(j-1);
            ap::vadd(b.getvector(j, n), t.getcolumn(j-1, j, n), temp);
            b(j) = b(j)/t(j,j);
        }
        return;
    }
    if( cse==2 )
    {
        b(n) = b(n)/t(n,n);
        if( n<2 )
        {
            return;
        }
        for(jj = 2; jj <= n; jj++)
        {
            j = n-jj+1;
            temp = -b(j+1);
            ap::vadd(b.getvector(1, j), t.getcolumn(j+1, 1, j), temp);
            b(j) = b(j)/t(j,j);
        }
        return;
    }
    if( cse==3 )
    {
        b(n) = b(n)/t(n,n);
        if( n<2 )
        {
            return;
        }
        for(jj = 2; jj <= n; jj++)
        {
            j = n-jj+1;
            v = ap::vdotproduct(t.getcolumn(j, j+1, j+1+jj-1-1), b.getvector(j+1, j+1+jj-1-1));
            b(j) = b(j)-v;
            b(j) = b(j)/t(j,j);
        }
        return;
    }
    if( cse==4 )
    {
        b(1) = b(1)/t(1,1);
        if( n<2 )
        {
            return;
        }
        for(j = 2; j <= n; j++)
        {
            v = ap::vdotproduct(t.getcolumn(j, 1, j-1), b.getvector(1, j-1));
            b(j) = b(j)-v;
            b(j) = b(j)/t(j,j);
        }
        return;
    }
}


void lbfgsbnewiteration(const ap::real_1d_array&,
     ap::real_value_type,
     const ap::real_1d_array&)
{

}

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
void
lbfgsbminimize
(
  FunctionAndGradient *const functionAndGradient,
  const int& n,
  const int& m,
  ap::real_1d_array& x,
  const ap::real_value_type& epsg,
  const ap::real_value_type& epsf,
  const ap::real_value_type& epsx,
  const int& maxits,
  const ap::integer_1d_array& nbd,
  const ap::real_1d_array& l,
  const ap::real_1d_array& u,
  int& info )
{
  ap::real_value_type f;
  ap::real_1d_array g;
  ap::real_1d_array xold;
  ap::real_1d_array xdiff;
  ap::real_2d_array ws;
  ap::real_2d_array wy;
  ap::real_2d_array sy;
  ap::real_2d_array ss;
  ap::real_2d_array yy;
  ap::real_2d_array wt;
  ap::real_2d_array wn;
  ap::real_2d_array snd;
  ap::real_1d_array z;
  ap::real_1d_array r;
  ap::real_1d_array d;
  ap::real_1d_array t;
  ap::real_1d_array wa;
  ap::real_1d_array sg;
  ap::real_1d_array sgo;
  ap::real_1d_array yg;
  ap::real_1d_array ygo;
  ap::integer_1d_array index;
  ap::integer_1d_array iwhere;
  ap::integer_1d_array indx2;
  int csave;
  ap::boolean_1d_array lsave;
  ap::integer_1d_array isave;
  ap::real_1d_array dsave;
  int task = 0;
  bool prjctd;
  bool cnstnd;
  bool boxed;
  bool updatd;
  bool wrk;
  int i;
  int k;
  int nintol;
  int iback;
  int nskip;
  int head;
  int col;
  int iter;
  int itail;
  int iupdat;
  int nint;
  int nfgv;
  int internalinfo;
  int ifun;
  int iword;
  int nfree;
  int ileave;
  int nenter;
  ap::real_value_type theta;
  ap::real_value_type fold;
  ap::real_value_type dr;
  ap::real_value_type rr;
  ap::real_value_type dnrm;
  ap::real_value_type xstep;
  ap::real_value_type sbgnrm;
  ap::real_value_type ddum;
  ap::real_value_type dtd;
  ap::real_value_type gd;
  ap::real_value_type gdold;
  ap::real_value_type stp;
  ap::real_value_type stpmx;
  ap::real_value_type tf;
  ap::real_1d_array workvec;
  ap::real_1d_array workvec2;
  ap::real_1d_array dsave13;
  ap::real_1d_array wa0;
  ap::real_1d_array wa1;
  ap::real_1d_array wa2;
  ap::real_1d_array wa3;
  ap::real_2d_array workmat;
  ap::integer_1d_array isave2;

  workvec.setbounds(1, m);
  workvec2.setbounds(1, 2*m);
  workmat.setbounds(1, m, 1, m);
  isave2.setbounds(1, 2);
  dsave13.setbounds(1, 13);
  wa0.setbounds(1, 2*m);
  wa1.setbounds(1, 2*m);
  wa2.setbounds(1, 2*m);
  wa3.setbounds(1, 2*m);
  g.setbounds(1, n);
  xold.setbounds(1, n);
  xdiff.setbounds(1, n);
  ws.setbounds(1, n, 1, m);
  wy.setbounds(1, n, 1, m);
  sy.setbounds(1, m, 1, m);
  ss.setbounds(1, m, 1, m);
  yy.setbounds(1, m, 1, m);
  wt.setbounds(1, m, 1, m);
  wn.setbounds(1, 2*m, 1, 2*m);
  snd.setbounds(1, 2*m, 1, 2*m);
  z.setbounds(1, n);
  r.setbounds(1, n);
  d.setbounds(1, n);
  t.setbounds(1, n);
  wa.setbounds(1, 8*m);
  sg.setbounds(1, m);
  sgo.setbounds(1, m);
  yg.setbounds(1, m);
  ygo.setbounds(1, m);
  index.setbounds(1, n);
  iwhere.setbounds(1, n);
  indx2.setbounds(1, n);
  lsave.setbounds(1, 4);
  isave.setbounds(1, 23);
  dsave.setbounds(1, 29);
  col = 0;
  head = 1;
  theta = 1;
  iupdat = 0;
  updatd = false;
  iter = 0;
  nfgv = 0;
  nint = 0;
  nintol = 0;
  nskip = 0;
  nfree = n;
  internalinfo = 0;
  ap::lbfgsberrclb(n, m, epsf, l, u, nbd, task, internalinfo, k);
  if( task==2||maxits<0||epsg<0||epsx<0 )
    {
    info = -1;
    return;
    }
  ap::lbfgsbactive(n, l, u, nbd, x, iwhere, prjctd, cnstnd, boxed);
  ap::vmove(xold.getvector(1, n), x.getvector(1, n));
  functionAndGradient->Evaluate(x, f, g);
  nfgv = 1;
  ap::lbfgsbprojgr(n, l, u, nbd, x, g, sbgnrm);
  if( sbgnrm<=epsg )
    {
    info = 4;
    return;
    }
  while(true)
    {
    iword = -1;
    if( !cnstnd&&col>0 )
      {
      ap::vmove(z.getvector(1, n), x.getvector(1, n));
      wrk = updatd;
      nint = 0;
      }
    else
      {
      ap::vmove(wa0.getvector(1, 2*m), wa.getvector(1, 2*m));
      ap::vmove(wa1.getvector(1, 2*m), wa.getvector(2*m+1, 4*m));
      ap::vmove(wa2.getvector(1, 2*m), wa.getvector(4*m+1, 6*m));
      ap::vmove(wa3.getvector(1, 2*m), wa.getvector(6*m+1, 8*m));
      ap::lbfgsbcauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt, theta, col, head, wa0, wa1, wa2, wa3, nint, sg, yg, sbgnrm, internalinfo, workvec);
      ap::vmove(wa.getvector(1, 2*m), wa0.getvector(1, 2*m));
      ap::vmove(wa.getvector(2*m+1, 4*m), wa1.getvector(1, 2*m));
      ap::vmove(wa.getvector(4*m+1, 6*m), wa2.getvector(1, 2*m));
      ap::vmove(wa.getvector(6*m+1, 8*m), wa3.getvector(1, 2*m));
      if( internalinfo!=0 )
	{
	internalinfo = 0;
	col = 0;
	head = 1;
	theta = 1;
	iupdat = 0;
	updatd = false;
	continue;
	}
      nintol = nintol+nint;
      ap::lbfgsbfreev(n, nfree, index, nenter, ileave, indx2, iwhere, wrk, updatd, cnstnd, iter);
      }
    if( nfree!=0&&col!=0 )
      {
      if( wrk )
	{
	ap::lbfgsbformk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd, m, ws, wy, sy, theta, col, head, internalinfo, workvec, workmat);
	}
      if( internalinfo!=0 )
	{
	internalinfo = 0;
	col = 0;
	head = 1;
	theta = 1;
	iupdat = 0;
	updatd = false;
	continue;
	}
      ap::lbfgsbcmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head, nfree, cnstnd, internalinfo, workvec, workvec2);
      if( internalinfo==0 )
	{
	ap::lbfgsbsubsm(n, m, nfree, index, l, u, nbd, z, r, ws, wy, theta, col, head, iword, wa, wn, internalinfo);
	}
      if( internalinfo!=0 )
	{
	internalinfo = 0;
	col = 0;
	head = 1;
	theta = 1;
	iupdat = 0;
	updatd = false;
	continue;
	}
      }
    for(i = 1; i <= n; i++)
      {
      d(i) = z(i)-x(i);
      }
    task = 0;
    while(true)
      {
      ap::lbfgsblnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, stp, dnrm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, internalinfo, task, boxed, cnstnd, csave, isave2, dsave13);
      if( internalinfo!=0||iback>=20||task!=1 )
	{
	break;
	}
      functionAndGradient->Evaluate(x, f, g);
      }
    if( internalinfo!=0 )
      {
      ap::vmove(x.getvector(1, n), t.getvector(1, n));
      ap::vmove(g.getvector(1, n), r.getvector(1, n));
      f = fold;
      if( col==0 )
	{
	if( internalinfo==0 )
	  {
	  internalinfo = -9;
	  nfgv = nfgv-1;
	  ifun = ifun-1;
	  iback = iback-1;
	  }
	task = 2;
	iter = iter+1;
	functionAndGradient->NextIteration( iter );
	info = -2;
	return;
	}
      else
	{
	if( internalinfo==0 )
	  {
	  nfgv = nfgv-1;
	  }
	internalinfo = 0;
	col = 0;
	head = 1;
	theta = 1;
	iupdat = 0;
	updatd = false;
	continue;
	}
      }
    iter = iter+1;
    functionAndGradient->NextIteration( iter );
    ap::lbfgsbnewiteration(x, f, g);
    ap::lbfgsbprojgr(n, l, u, nbd, x, g, sbgnrm);
    if( sbgnrm<=epsg )
      {
      info = 4;
      return;
      }
    ap::vmove(xdiff.getvector(1, n), xold.getvector(1, n));
    ap::vsub(xdiff.getvector(1, n), x.getvector(1, n));
    tf = ap::vdotproduct(xdiff.getvector(1, n), xdiff.getvector(1, n));
    tf = sqrt(tf);
    if( tf<=epsx )
      {
      info = 2;
      return;
      }
    ddum = ap::maxreal(fabs(fold), ap::maxreal(fabs(f), ap::real_value_type(1)));
    if( fold-f<=epsf*ddum )
      {
      info = 1;
      return;
      }
    if( iter>maxits&&maxits>0 )
      {
      info = 5;
      return;
      }
    if( additionallbfgsbstoppingcriterion(iter, x, f, g) )
      {
      info = 0;
      return;
      }
    ap::vmove(xold.getvector(1, n), x.getvector(1, n));
    for(i = 1; i <= n; i++)
      {
      r(i) = g(i)-r(i);
      }
    rr = ap::vdotproduct(r.getvector(1, n), r.getvector(1, n));
    if( stp==1 )
      {
      dr = gd-gdold;
      ddum = -gdold;
      }
    else
      {
      dr = (gd-gdold)*stp;
      ap::vmul(d.getvector(1, n), stp);
      ddum = -gdold*stp;
      }
    if( dr<=ap::machineepsilon*ddum )
      {
      nskip = nskip+1;
      updatd = false;
      }
    else
      {
      updatd = true;
      iupdat = iupdat+1;
      ap::lbfgsbmatupd(n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head, theta, rr, dr, stp, dtd);
      ap::lbfgsbformt(m, wt, sy, ss, col, theta, internalinfo);
      if( internalinfo!=0 )
	{
	internalinfo = 0;
	col = 0;
	head = 1;
	theta = 1;
	iupdat = 0;
	updatd = false;
	continue;
	}
      }
    }
}

} // namespace ap
