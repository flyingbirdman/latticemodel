#include "1dlattice.h"
#include "nb.cpp"
#include "nrlib.h"
#include "vector_operation.h"
template<class T>
double D1Mesh<T>::kx_i(int index){
  if(index<0 || index>=Lq){
    std::clog << "index:  " << index << " ;[d2mesh] kx <0 || kx>=Lq" << std::endl;
    std::exit(1.0);
  }
  return px1+(px2-px1)*(index+0.5)/Lq;
}

template<class T>
void D1Mesh<T>::cal_grid(){
  DoubleVector d2(pt.d1);
  for(int i=0; i<Lq; i++){
    ki[0]=kx_i(i);    
    esol->cal_ei(ki,pt.en,veci);
    
    ki[0]=ki[0]-eps;
    esol->cal_ei(ki,pt.d1,veci);
        
    ki[0]=ki[0]+2.0*eps;
    esol->cal_ei(ki,d2,veci);
    for(int j=0; j<numOfBand; j++)
      pt.d1[j]=(d2[j]-pt.d1[j])/eps/2.0;    
    grid.push_back(pt);
  }
}

template<class T>
void D1Mesh<T>::d1FourTransform(int ii, std::ostream &out){
  ComplexVector h_k(Lq);
  double delta= (px2-px1);
  for(int i=0; i<Lq; i++)
    h_k[i]=grid[i].en[ii];
  fft_1d(h_k, 1);
  print_fft_1d(h_k, delta, out);
}

template<class T>
void D1Mesh<T>::datToFile(std::ostream &out){
  out << "####px1,px2,py1,py2,Lq,numOfBand ####" << std::endl;
  out << px1 << "   " << px2 << "   "  << Lq << "   " << numOfBand << std::endl;
  out << "############################################" << std::endl;
  for(int i=0; i<grid.size(); i++){
    for(int k=0; k<numOfBand; k++)
        out << grid[i].en[k] << "  " << grid[i].d1[k] << "   ";
    out << std::endl;
  }
}

template<class T>
void D1Mesh<T>::datFromFile(std::istream &in){
  std::string inbuf;
  in >> inbuf;
  in >> px1 >> px2 >> Lq >> numOfBand;
  in >> inbuf;
  pt.en.resize(numOfBand);
  pt.d1.resize(numOfBand);
  for(int i=0; i<Lq; i++){
    for(int k=0; k<numOfBand; k++)
      in >> pt.en[k]  >> pt.d1[k];
    grid.push_back(pt);
  }
}

template<class T>
void D1Dos<T>::dosFunc(){
  double eq, e1;
  
  int ip, im;
  double em, ep, dq1, am, ap;
  for(int i=0; i<ne; i++)
    for(int j=0; j<numOfBand; j++)
      xx[i][j]=0;
  for(int i=0; i<nq; i++){
    for(int ie=0; ie<numOfBand;ie++)
    {
      eq=d1m->grid[i].en[ie];
      e1=fabs(d1m->grid[i].d1[ie])+eps;
      eq=eq-e1*dq/2.0;
      im=(eq-emin)/de-2.0;
      if(im<0) im=0;
      ip=im+e1*dq/de+4.0;
      if(ip>ne-1) ip=ne-1;
      for(int j=im;  j<=ip;  j++){
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
        xx[j][ie+1]= xx[j][ie+1]+fabs(ap-am);
      }
    }
  }
  for(int i=0; i<ne; i++){
    xx[i][0]=(i+0.5)*de;
    for(int j=1; j<=numOfBand; j++){
      xx[i][j]=xx[i][j]/zoneArea;
      xx[i][j]=xx[i][j]/de;
    }
  }
}

template<class T>
void D1Dos<T>::hilbertTransform(){
  double x,y;
  for(int ie=0; ie<numOfBand; ie++){
    for(int i=0; i<ne; i++){
      yy[i][0]=xx[i][0];
      x=0;
      for(int j=0; j<ne; j++){
        // if(i==j && i==0) continue;
        // y=i+j;
        // x=x-xx[j][ie+1]/y;
        // uncomment for correlation function
        if(i==j) continue;
        y=i-j;
        x=x+xx[j][ie+1]/y;
      }
      yy[i][ie+1]=-x/pi;
    }
  }
}

template<class T>
void D1Dos<T>::xyToFile(std::ostream &out){  
  for(int i=0; i<ne; i++){
    out << xx[i][0] << "  ";
    for(int ie=0; ie<numOfBand; ie++)
      out << yy[i][ie+1] << "  " << xx[i][ie+1] << "  ";
    out << std::endl;
  }
}
