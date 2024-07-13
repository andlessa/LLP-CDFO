(* ::Package:: *)

Off[General::spell]

Model`Name = "SSDMv2";
Model`NameLaTeX ="Two Singlet scalar Dark Matter";
Model`Authors = "Andre Lessa (based on SM model by F.Staub and SSDM by Diego Restrepo and arXiv:2404.19057)";
Model`Date = "2024-07-02";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global Symmetries *)

Global[[1]] = {Z[2], Z2};


(* Gauge Groups *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,1};
Gauge[[2]]={WB, SU[2], left,        g2,True,1};
Gauge[[3]]={G,  SU[3], color,       g3,False,1};


(* Matter Fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3,1};  
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1,1};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3,1};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3,1};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1,1};

ScalarFields[[1]] =  {H, 1, {Hp, H0},     1/2, 2,  1,1};
ScalarFields[[2]] =  {S1, 1, ss1,     0, 1,  1, -1};
ScalarFields[[3]] =  {S2, 1, ss2,     0, 1,  1, -1};
RealScalars = {S1,S2};



        
(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput]= {
	{LagHC, {AddHC->True}},
	{LagNoHC,{AddHC->False}}
};


LagNoHC = -(mu2 conj[H].H + Lambda1/2 conj[H].H.conj[H].H + MS12/2 S1.S1 + MS22/2 S2.S2 
           + (LamH1 S1.S1.conj[H].H  + LamH2 S2.S2.conj[H].H + LamH12 S1.S2.conj[H].H)  
(*           - (v^2/2)*(LamH1 S1.S1.conj[H].H  + LamH2 S2.S2.conj[H].H + LamH12 S1.S2.conj[H].H) *)
           + Lam14 S1.S1.S1.S1  + Lam24 S2.S2.S2.S2   + Lam22 S1.S1.S2.S2
           + Lam13 S1.S2.S2.S2  + Lam31 S1.S1.S1.S2);
LagHC =  -(Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q);



			  		  

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};     
        
        
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
{    {H0, {v, 1/Sqrt[2]}, {Ah, I/Sqrt[2]},{hh, 1/Sqrt[2]}}     };
 

DEFINITION[EWSB][MatterSector]=   
    {{{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}};  


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  DL, conj[DR]},
 Fe ->{  EL, conj[ER]},
 Fu ->{  UL, conj[UR]},
 Fv ->{  vL, 0}};

DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};



