(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`ItoSDEs`", "CFMLab`StockStat`", "Statistics`NormalDistribution`"]

SDESolver1::usage = "SDESolver1[b_, c_, s_, t0_, t1_, K_]
computes the Monte-Carlo solution trajectory of a scalar-valued SDE."

SDESolver2::usage = "SDESolver2[b_, c_, s_, t0_, t1_, K_]
computes the Monte-Carlo solution trajectory of a vector-valued SDE."

SDESolver::usage = "SDESolver[b_, c_, s_, t0_, t1_, K_]
computes the Monte-Carlo solution trajectory of an SDE, whether it is scalar or vector
valued."

ItoDrift::usage = "ItoDrift[g_, {b_, c_}, {t_, S_}] computes the drift of the
process according to the Ito's chain rule."

ItoDiffusion::usage = "ItoDiffusion[g_, c_, {t_, S_}] computes the diffusion of the
process according to the Ito's chain rule."

StockExtract::usage = "StockExtract[i_, S\[Tau]_] extracts individual stock-price trajectories
from the vector-price trajectory."

Begin["`Private`"]

SDESolver1[b_, c_, s_, t0_, t1_, K_] := 
  Module[{dt, dB, G}, dt = N[(t1 - t0)/K]; 
    dB = Table[Random[NormalDistribution[0, Sqrt[dt]]], {K}]; 
    G[{t_, S_}, dB_] := {t + dt, S + b[t, S]*dt + c[t, S]*dB}; 
    FoldList[G, {t0, s}, dB]]

SDESolver2[b_, c_, s_, t0_, t1_, K_] := 
   Module[{n, dt, dB, G}, n = Length[Transpose[c[t0,s]]]; 
     dt = (t1 - t0)/K; dB = Table[Random[NormalDistribution[0, 
         Sqrt[dt]]], {K}, {n}]; G[{t_, S_}, dB_] := 
      {t + dt, S + b[t, S]*dt + c[t, S] . dB}; 
     FoldList[G, {t0, s}, dB]]
    
SDESolver[b_, c_, s_, t0_, t1_, K_] := 
  If[Head[b[t0, s]] === List, SDESolver2[b, c, s, t0, t1, K], 
   SDESolver1[b, c, s, t0, t1, K]]
   
StockExtract[i_, S\[Tau]_] := ({#1[[1]], #1[[2,i]]} & ) /@ S\[Tau]


Der[u_,V_,0]:=u

Der[u_,V_,1]:=D[u,#]&/@V

Der[u_,V_,k_]:=Nest[Der[#,V,1]&,u,k]

ItoDrift[g_, {b_, c_}, {t_, S_}] := 
  D[g, t] + Der[g,S,1]. b + 
   (1/2)*Plus @@ Flatten[c . Transpose[c]*Der[g,S,2]]

ItoDiffusion[g_, c_, {t_, S_}] := Der[g,S,1].c

End[]

EndPackage[]