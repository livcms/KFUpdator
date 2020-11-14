#include "TrackingTools/KalmanUpdators/interface/prototype.h" 
#include <vector> 
#include <algorithm> 
#include <iostream> 


output testfunction(std::vector<double> xin, std::vector<double> Cin, std::vector<double> rin, std::vector<double> Vin, std::vector<double> VMeasin ){
 
  std::vector<double> rMeasin = {xin[3], xin[4]};

  std::vector<double> r2(2);
 
  // r = r-RMeas
  std::transform (rin.begin(), rin.end(), rMeasin.begin(), r2.begin(), std::minus<double>());
  
  //covariance matrix of residuals
  std::vector<double> R(3);
  std::transform(Vin.begin(), Vin.end(), VMeasin.begin(), R.begin(), std::plus<double>());
  
  //invert R              
  double det_R = R[0]*R[2] - R[1]*R[1]; 
                   
  std::vector<double> InvertedR = {R[2]/det_R, -R[1]/det_R,  R[0]/det_R}; 
  
  //Kalman gain 

  std::vector<double> K = {Cin[6]*InvertedR[0] + Cin[10]*InvertedR[1], Cin[6]*InvertedR[1]+Cin[10]*InvertedR[2], 
                           Cin[7]*InvertedR[0] + Cin[11]*InvertedR[1], Cin[7]*InvertedR[1] + Cin[11]*InvertedR[2], 
                           Cin[8]*InvertedR[0] + Cin[12]*InvertedR[1], Cin[8]*InvertedR[1] + Cin[12]*InvertedR[2],
                           Cin[9]*InvertedR[0] + Cin[13]*InvertedR[1], Cin[9]*InvertedR[1] + Cin[13]*InvertedR[2],
                           Cin[13]*InvertedR[0] + Cin[14]*InvertedR[1], Cin[13]*InvertedR[1] + Cin[14]*InvertedR[2]}; 
  
  std::vector<double> K_Mult_r  = {K[0]*r2[0]+K[1]*r2[1], K[2]*r2[0]+K[3]*r2[1], K[4]*r2[0] + K[5]*r2[1], K[6]*r2[0] + K[7]*r2[1],
                                  K[8]*r2[0]+K[9]*r2[1]};
                                                               
  std::vector<double> fsv(5);
  std::transform(xin.begin(), xin.end(), K_Mult_r.begin(), fsv.begin(), std::plus<double>());
 

   std::vector<double> CFull = {Cin[0], Cin[1], Cin[3], Cin[6], Cin[10], Cin[1], Cin[2], Cin[4], Cin[7], Cin[11], Cin[3], Cin[4], Cin[5], Cin[8], Cin[12], Cin[6], Cin[7], Cin[8], Cin[9], Cin[13], Cin[10], Cin[11], Cin[12], Cin[13], Cin[14]};

  std::vector<double> K1 = {K[0],K[2], K[4], K[6], K[8]};
  std::vector<double>  K2  = {K[1],K[3], K[5], K[7], K[9]};
  std::vector<double> res1(5);
  std::vector<double> res2(5);
  std::vector<double> A;
  for(int i = 0; i<5; i++){
 //
  double myconst = CFull[3+5*i];
  double myconst2 = CFull[4+5*i];
  std::transform(K1.begin(), K1.end(), res1.begin(), [&myconst](auto& c){return c*myconst;});
            std::transform(K2.begin(), K2.end(), res2.begin(), [&myconst2](auto& c){return c*myconst2;});
                std::transform(res1.begin(), res1.end(), res2.begin(), res1.begin(), std::plus<double>());
                    A.insert(A.end(), res1.begin(), res1.end());
                     }


   std::vector<double> A1 = {A[15], A[16], A[17], A[18], A[19]};
   std::vector<double> A2 = {A[20], A[21], A[22], A[23], A[24]};
   std::vector<double> res3(5);
   std::vector<double> res4(5);
   std::vector<double> KA;
   for(int i = 0; i<5; i++){
     double Kelem = K[2*i];
     double Kelem2 = K[2*i+1];
     std::transform(A1.begin(), A1.end(), res3.begin(), [&Kelem](auto& c){return c*Kelem;});
     std::transform(A2.begin(), A2.end(), res4.begin(), [&Kelem2](auto& c){return c*Kelem2;});
     std::transform(res3.begin(), res3.end(), res4.begin(), res3.begin(), std::plus<double>());
     KA.insert(KA.end(), res3.begin(), res3.end());
   }
   std::vector<double> VFull = {Vin[0], Vin[1], Vin[1], Vin[2]};
    std::vector<double> VKT;
    std::vector<double> res5 (5);
    std::vector<double> res6 (5);
    std::vector<double> K3 = {K[0], K[2], K[4], K[6], K[8]};
    std::vector<double> K4 = {K[1], K[3], K[5], K[7], K[9]};
    for(int i = 0; i<2; i++){
      double V = VFull[2*i];
      double V1 = VFull[2*i+1];
      std::transform(K3.begin(), K3.end(), res5.begin(), [&V](auto& c){return c*V;});
      std::transform(K4.begin(), K4.end(), res6.begin(), [&V1](auto& c){return c*V1;});
      std::transform(res5.begin(), res5.end(), res6.begin(), res5.begin(), std::plus<double>());
      VKT.insert(VKT.end(), res5.begin(), res5.end());
    }


  std::vector<double> KVKT;
  std::vector<double> VKT1 = {VKT[0], VKT[1], VKT[2], VKT[3], VKT[4]};
  std::vector<double> VKT2 = {VKT[5], VKT[6], VKT[7], VKT[8], VKT[9]};
  std::vector<double> res7 (5);
  std::vector<double> res8 (5);
  for(int i = 0; i<5; i++){
    double Kelem = K[2*i];
    double Kelem2 = K[2*i+1];
    std::transform(VKT1.begin(), VKT1.end(), res7.begin(), [&Kelem](auto& c){return c*Kelem;});
    std::transform(VKT2.begin(), VKT2.end(), res8.begin(), [&Kelem2](auto& c){return c*Kelem2;});
    std::transform(res7.begin(), res7.end(), res8.begin(), res7.begin(), std::plus<double>());
    KVKT.insert(KVKT.end(), res7.begin(), res7.end());
  }

  std::vector<double> fse;
  for (int i =0; i<5; i++){
    for(int j = 0; j<5; j++){
    fse.push_back(CFull[5*i+j] - A[5*i+j] - A[i+5*j] + KA[5*i+j] + KVKT[5*i+j]);
}
}

 output finalresult; 
 finalresult.finalStateVector = fsv; 
 finalresult.finalStateError = {fse[0], fse[5], fse[6], fse[10], fse[11], fse[12], fse[15], fse[16], fse[17], fse[18], fse[20], fse[21], fse[22], fse[23], fse[24]}; 
  return(finalresult); 
} 
 
