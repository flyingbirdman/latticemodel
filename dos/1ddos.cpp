#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#define pi 3.1415927
#define eps 0.000000001
#define t 1.0

void d1bands(double qx, double &e,double &e1);

int main(){
  std::ofstream out("1d.dat");
  int nq=300;
  int ne=1000;
  std::vector<double> xx(ne);
  double qx, eq, e1;
  double de=4.0/ne;
  double dq=2*pi/nq;
  
  int im, ip;
  double e, em, ep, dq1, am, ap;
  for(int i=0; i<ne; i++)
    xx[i]=0;
  for(int i=1; i<= nq; i++){
    qx=-pi+(i-0.5)*dq;
    d1bands(qx,eq,e1);
    e1=fabs(e1)+eps;
    eq=eq-e1*dq/2.0;
    im=eq/de-2.0;
    if(im<0) im=0;
    ip=im+e1*dq/de+4.0;
    if(ip>ne-1) ip=ne-1;
    for(int j=im; j<=ip; j++){
      em=de*j;
      ep=em+de;
      dq1=(em-eq)/e1;
      am=dq1;
      if(dq1>dq) am=dq;
      if(dq1<0) am=0;
      dq1=(ep-eq)/e1;
      ap=dq1;
      if(dq1>dq) ap=dq;
      if(dq1<0) ap=0;
      xx[j]= xx[j]+fabs(ap-am);
    }
  }
  double x=1.0/2.0/pi;
  for(int i=0; i<ne; i++){
    xx[i]=xx[i]*x;
    xx[i]=xx[i]/de;
    e=(i+0.5)*de;
    out << e << "    " << xx[i] << std::endl;
  }
}


void d1bands(double qx, double &e,double &e1){
  e= t+t*cos(qx);  
  e1= -t*sin(qx);
}
