#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#define pi 3.1415927
#define eps 0.000000001
#define t 1.0

void d2bands(double qx, double qy, double &e, double &e1, double &e2);

int main(){
  std::ofstream out("2d.dat");
  int nq=64;
  int ne=1000;
  std::vector<double> xx(ne);
  double qx, qy, eq, e1, e2;
  double de=8.0/ne;
  double dq=2*pi/nq;
  
  int im, ip;
  double e, em, ep, dq1, dq2, am, ap;
  for(int i=0; i<ne; i++)
    xx[i]=0;
  for(int i=1; i<= nq; i++){
    qx=-pi+(i-0.5)*dq;
    for(int j=1; j<= nq; j++){
      qy=-pi+(j-0.5)*dq;
      
      d2bands(qx,qy,eq,e1,e2);
      e1=fabs(e1)+eps;
      e2=fabs(e2)+eps;
      eq=eq-(e1+e2)*dq/2.0;
      im=eq/de-2.0;
      if(im<0) im=0;
      ip=im+(e1+e2)*dq/de+4.0;
      if(ip>ne-1) ip=ne-1;
      //std::cout << e1 << "  " << e2 << "  " << e3 << "  " << ep << "  " << im << "  " << ip << "\n";
      for(int xxi=im; xxi<=ip; xxi++){
        em=de*xxi;
        ep=em+de;
        dq1=(em-eq)/e1;
        dq2=(em-eq)/e2;
        am=dq1*dq2;  // the prefix 1/2
        if(dq1>dq) am=am-(dq1-dq)*(dq1-dq)*dq2/dq1;
        if(dq2>dq) am=am-(dq2-dq)*(dq2-dq)*dq1/dq2;
        if(dq1*dq2>(dq1+dq2)*dq) am=2.0*dq*dq;
        if(dq1<0 || dq2<0) am=0;
          
        dq1=(ep-eq)/e1;
        dq2=(ep-eq)/e2;
        ap=dq1*dq2;  // the prefix 1/2
        if(dq1>dq) ap=ap-(dq1-dq)*(dq1-dq)*dq2/dq1;
        if(dq2>dq) ap=ap-(dq2-dq)*(dq2-dq)*dq1/dq2;
        if(dq1*dq2>(dq1+dq2)*dq) ap=2.0*dq*dq;
        if(dq1<0 || dq2<0) ap=0;          
        xx[xxi]=xx[xxi]+fabs(ap-am);
        // std::cout << fabs(ap-am) << std::endl;
      }
    }
  }
  double x=1.0/(2.0*pi)/(2.0*pi)/2.0;
  
  for(int i=0; i<ne; i++){
    xx[i]=xx[i]*x;
    xx[i]=xx[i]/de;
    e=(i+0.5)*de;
    out << e << "  " << xx[i] << std::endl;
  }
}

void d2bands(double qx, double qy, double &e, double &e1, double &e2){
  e=2.0*t+t*(cos(qx)+cos(qy));
  e1=-t*sin(qx);
  e2=-t*sin(qy);
}
