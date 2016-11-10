#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#define pi 3.1415927
#define eps 0.0000001
#define t 1.0

void d3bands(double qx, double qy, double qk, double &e, double &e1, double &e2, double &e3);

int main(){
  std::ofstream out("3d.dat");
  int nq=60;
  int ne=500;
  std::vector<double> xx(ne);
  double qx, qy, qz, eq, e1, e2, e3;
  double de=8.0/ne;
  double dq=2*pi/nq;
  
  int im, ip;
  double e, em, ep, dq1, dq2, dq3, am, ap;
  for(int i=0; i<ne; i++)
    xx[i]=0;
  for(int i=1; i<= nq; i++){
    qx=-pi+(i-0.5)*dq;
    for(int j=1; j<= nq; j++){
      qy=-pi+(j-0.5)*dq;
      for(int k=1; k<= nq; k++){
        qz=-pi+(k-0.5)*dq;
        d3bands(qx,qy,qz,eq,e1,e2,e3);
        e1=fabs(e1)+eps;
        e2=fabs(e2)+eps;
        e3=fabs(e3)+eps;
        eq=eq-(e1+e2+e3)*dq/2.0;
        im=eq/de-2.0;
        if(im<0) im=0;
        ip=im+(e1+e2+e3)*dq/de+4.0;
        if(ip>ne-1) ip=ne-1;
        //std::cout << e1 << "  " << e2 << "  " << e3 << "  " << ep << "  " << im << "  " << ip << "\n";
        for(int xxi=im; xxi<=ip; xxi++){
          em=de*xxi;
          ep=em+de;
          dq1=(em-eq)/e1;
          dq2=(em-eq)/e2;
          dq3=(em-eq)/e3;
          am=dq1*dq2*dq3;  // the prefix 1/6
          if(dq1>dq) am=am-(dq1-dq)*(dq1-dq)*(dq1-dq)*dq2*dq3/dq1/dq1;
          if(dq2>dq) am=am-(dq2-dq)*(dq2-dq)*(dq2-dq)*dq1*dq3/dq2/dq2;
          if(dq3>dq) am=am-(dq3-dq)*(dq3-dq)*(dq3-dq)*dq1*dq2/dq3/dq3;
          if(dq1*dq2>(dq1+dq2)*dq) am=am+pow((dq1*dq2-dq*(dq1+dq2)),3.0)*dq3/(dq1*dq2)/(dq1*dq2);
          if(dq1*dq3>(dq1+dq3)*dq) am=am+pow((dq1*dq3-dq*(dq1+dq3)),3.0)*dq2/(dq1*dq3)/(dq1*dq3);
          if(dq2*dq3>(dq2+dq3)*dq) am=am+pow((dq2*dq3-dq*(dq2+dq3)),3.0)*dq1/(dq2*dq3)/(dq2*dq3);
          if(dq1*dq2*dq3>(dq1*dq2+dq1*dq3+dq2*dq3)*dq) am=6.0*dq*dq*dq;
          if(dq1<0 || dq2<0 || dq3<0) am=0;

          dq1=(ep-eq)/e1;
          dq2=(ep-eq)/e2;
          dq3=(ep-eq)/e3;
          ap=dq1*dq2*dq3;  // the prefix 1/6
          if(dq1>dq) ap=ap-(dq1-dq)*(dq1-dq)*(dq1-dq)*dq2*dq3/dq1/dq1;
          if(dq2>dq) ap=ap-(dq2-dq)*(dq2-dq)*(dq2-dq)*dq1*dq3/dq2/dq2;
          if(dq3>dq) ap=ap-(dq3-dq)*(dq3-dq)*(dq3-dq)*dq1*dq2/dq3/dq3;
          if(dq1*dq2>(dq1+dq2)*dq) ap=ap+pow((dq1*dq2-dq*(dq1+dq2)),3.0)*dq3/(dq1*dq2)/(dq1*dq2);
          if(dq1*dq3>(dq1+dq3)*dq) ap=ap+pow((dq1*dq3-dq*(dq1+dq3)),3.0)*dq2/(dq1*dq3)/(dq1*dq3);
          if(dq2*dq3>(dq2+dq3)*dq) ap=ap+pow((dq2*dq3-dq*(dq2+dq3)),3.0)*dq1/(dq2*dq3)/(dq2*dq3);
          if(dq1*dq2*dq3>(dq1*dq2+dq1*dq3+dq2*dq3)*dq) ap=6.0*dq*dq*dq;          
          if(dq1<0 || dq2<0 || dq3<0) ap=0;
          
          xx[xxi]=xx[xxi]+fabs(ap-am);
        }
      }
    }
  }
  double x=1.0/(2.0*pi)/(2.0*pi)/(2.0*pi)/6.0;
  
  for(int i=0; i<ne; i++){
    xx[i]=xx[i]*x;
    xx[i]=xx[i]/de;
    e=(i+0.5)*de;
    out << e << "  " << xx[i] << std::endl;
  }
}

void d3bands(double qx, double qy, double qz, double &e, double &e1, double &e2, double &e3){
  e=3.0*t+t*(cos(qx)+cos(qy)+cos(qz));
  e1=-t*sin(qx);
  e2=-t*sin(qy);
  e3=-t*sin(qz);
}
 
