(* ::Package:: *)

OnlyLowEnergySPheno = True;


MINPAR={{1,Lambda1IN},
        {2,LamH1IN},
        {3,LamH2IN},
        {4,LamH12IN},
        {5,Lam14IN},
        {6,Lam24IN},
        {7,Lam22IN},
        {8,Lam13IN},
        {9,Lam31IN},
        {10,MS12input},
        {11,MS22input}
        };


ParametersToSolveTadpoles = {mu2};

BoundaryLowScaleInput={
   {v, vSM}, 
 {Ye, YeSM*vSM},
 {Yd, YdSM*vSM},
 {Yu, YuSM*vSM},
 {g1, g1SM},
 {g2, g2SM},
 {g3, g3SM},
  {Lambda1,Lambda1IN},
  {LamH1,LamH1IN},
  {LamH2,LamH2IN},
  {LamH12,LamH12IN},
  {Lam14,Lam14IN},
  {Lam24,Lam24IN},
  {Lam22,Lam22IN},
  {Lam13,Lam13IN},
  {Lam31,Lam31IN},
  {LamS12,LamS12IN},
  {MS12, MS12Input},
  {MS22, MS22Input}
};




ListDecayParticles = {Fu,Fe,Fd,hh};
ListDecayParticles3B = {{Fu,"Fu.f90"},{Fe,"Fe.f90"},{Fd,"Fd.f90"}};

FlagLoopContributions = True;  

DefaultInputValues ={Lambda1IN -> 0.28, LamH1IN -> 0.0, LamH2IN -> 1.0, LamH12IN -> 0.000026, Lam14IN -> 0, Lam24IN -> 0, Lam22IN -> 0, Lam13IN -> 0, Lam31IN -> 0, MS12input -> 250000, MS22input -> 255025};
