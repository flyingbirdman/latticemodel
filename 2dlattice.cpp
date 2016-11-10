#include "2dlattice.h"
#include "nb.cpp"
#include "nrlib.h"
#include "vector_operation.h"
template<class T>
double D2Mesh<T>::kx_i(int index){
  if(index<0 || index>=Lq){
    std::clog << "index:  " << index << " ;[d2mesh] kx <0 || kx>=Lq" << std::endl;
    std::exit(1.0);
  }
  return px1+(px2-px1)*(index+0.5)/Lq;
}

template<class T>
double D2Mesh<T>::ky_i(int index){
  if(index<0 || index>=Lq){
    std::clog << "index:  " << index << " ;[d2mesh] ky<0 || ky>=Lq" << std::endl;
    std::exit(1.0);
  }
  return py1+(py2-py1)*(index+0.5)/Lq;
}

template<class T>
void D2Mesh<T>::cal_grid(){
  for(int i=0; i<Lq; i++){    
    for(int j=0; j<Lq; j++){
      ki[0]=kx_i(i);
      ki[1]=ky_i(j);
      esol->cal_ei(ki,pt.en,veci);

      ki[0]=ki[0]+eps;
      esol->cal_ei(ki,pt.d1,veci);
            
      ki[0]=ki[0]-eps;
      ki[1]=ki[1]+eps;
      esol->cal_ei(ki,pt.d2,veci);
      for(int k=0; k<numOfBand; k++){
        pt.d1[k]=(pt.d1[k]-pt.en[k])/eps;
        pt.d2[k]=(pt.d2[k]-pt.en[k])/eps;
      }
      std::cout << i << "  " << j << "  " << pt.en[0] << "  " << pt.d1[0] << "  " << pt.d2[0] << std::endl;
      grid[i].push_back(pt);
    }
  }
}

template<class T>
void D2Mesh<T>::d2FourTransform(int ii, std::ostream &out){
  ComplexMatrix h_k(Lq);
  double delta1= (px2-px1);
  double delta2= (py2-py1);
  for(int i=0; i<Lq; i++)
    h_k[i].resize(Lq);
  for(int i=0; i<Lq; i++)
    for(int j=0; j<Lq; j++)
      h_k[i][j]=grid[i][j].en[ii];    
  fft_2d(h_k, 1);
  print_fft_2d(h_k, delta1, delta2, out);
}

template<class T>
void D2Mesh<T>::datToFile(std::ostream &out){
  out << "####px1,px2,py1,py2,Lq,numOfBand ####" << std::endl;
  out << px1 << "   " << px2 << "   " << py1 << "   " << py2 << "   " << Lq << "   " << numOfBand << std::endl;
  out << "############################################" << std::endl;
  for(int i=0; i<grid.size(); i++)
    for(int j=0; j<grid[0].size(); j++){
      for(int k=0; k<numOfBand; k++)
        out << grid[i][j].en[k] << "  " << grid[i][j].d1[k] << "   " << grid[i][j].d2[k] << "   ";
      out << std::endl;
    }
}

template<class T>
void D2Mesh<T>::datFromFile(std::istream &in){
  std::string inbuf;
  in >> inbuf;
  in >> px1 >> px2 >> py1 >> py2 >> Lq >> numOfBand;
  in >> inbuf;
  pt.en.resize(numOfBand);
  pt.d1.resize(numOfBand);
  pt.d2.resize(numOfBand);
  grid.resize(Lq);
  for(int i=0; i<Lq; i++)
    for(int j=0; j<Lq; j++){
      for(int k=0; k<numOfBand; k++)
        in >> pt.en[k]  >> pt.d1[k] >> pt.d2[k];
      grid[i].push_back(pt);
    }
}

template<class T>
void D2Dos<T>::get_ZoneArea(){
  DoubleVector q2(2);
  zoneArea=0;
  for(int i=0; i<nq; i++){
    q2[0]= qx1+(i+0.5)*dq;
    for(int j=0; j<nq; j++){
      q2[1]= qy1+(j+0.5)*dq;
      if(inZoneOrNot(q2)) zoneArea += dq*dq;
    }
  }
}

template<class T>
void D2Dos<T>::dosFunc(){
  DoubleVector q2(2);
  double eq, e1, e2, em, ep;
  double dq1, dq2, am, ap;
  int ip, im;
  for(int i=0; i<ne; i++)
    for(int j=0; j<numOfBand; j++)
      xx[i][j]=0;
  for(int i=0; i<nq; i++){
    q2[0]=d2m->kx_i(i);
    for(int j=0; j<nq; j++){
      q2[1]=d2m->ky_i(j);
      if(!inZoneOrNot(q2)) continue;
      for(int ie=0; ie<numOfBand;ie++)
      {
        eq=d2m->grid[i][j].en[ie];
        e1=fabs(d2m->grid[i][j].d1[ie])+eps;
        e2=fabs(d2m->grid[i][j].d2[ie])+eps;
        eq=eq-(e1+e2)*dq/2.0;
        im=(eq-emin)/de-2.0;
        if(im<0) im=0;
        ip=im+(e1+e2)*dq/de+4.0;
        if(ip>ne-1) ip=ne-1;
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
          xx[xxi][ie+1]=xx[xxi][ie+1]+fabs(ap-am);
        }
      }
    }
  }
  for(int i=0; i<ne; i++){
    xx[i][0]=(i+0.5)*de;
    for(int j=1; j<=numOfBand; j++){
      xx[i][j]=xx[i][j]/zoneArea/2.0;
      xx[i][j]=xx[i][j]/de;
    }
  }
}

template<class T>
void D2Dos<T>::hilbertTransform(){
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
void D2Dos<T>::xyToFile(std::ostream &out){  
  for(int i=0; i<ne; i++){
    out << xx[i][0] << "  ";
    for(int ie=0; ie<numOfBand; ie++)
      out << yy[i][ie+1] << "  " << xx[i][ie+1] << "  ";
    out << std::endl;
  }
}
