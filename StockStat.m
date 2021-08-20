(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`StockStat`", "Statistics`NormalDistribution`", 
"Statistics`MultinormalDistribution`"]

Av::usage = "Av[List_] computes the average of a List."

Ave::usage = "Ave[List_] computes the average of a List without flattening."

\[CapitalSigma]::usage = "\[CapitalSigma][PriceList_] computes the standard deviation
of the price data, assuming the Log-Normal model."

AEst::usage = "AEst[PriceList_] computes the appreciation rate
of the price data, assuming the Log-Normal model."

ConfidenceInterval::usage = "ConfidenceInterval[a_, \[Alpha]_, \[Sigma]_, T_]
computes the confidence interval."

A\[CapitalSigma]Est::usage = "A\[CapitalSigma]Est[PriceList_, p_] computes the 
standard deviation, appreciation rate, and confidence interval for the appreciation rate
of the price data, assuming the Log-Normal model."

StockCovariance::usage = "StockCovariance[stocksSDE_] computes the covariance matrix
of the price data, assuming the vector Log-Normal model."

Diag::usage = "Diag[Matrix_] extracts the diagonal of square matrix."

GrowthRateEstimate::usage = "GrowthRateEstimate[stocksSDE_] computes the appreciation vector
of the price data, assuming the vector Log-Normal model."

GrowthRatePosterior::usage = "GrowthRatePosterior[stocksSDE_]."

GrowthRateEstimateWithPrior::usage = "
GrowthRateEstimateWithPrior[stocksSDE_, ap_, \[Sigma]p_]."

GrowthRatePosteriorWithPrior::usage = "
GrowthRatePosteriorWithPrior[stocksSDE_, ap_, \[Sigma]p_]."

Begin["`Private`"]

Av[List_] := Plus @@ Flatten[List]/Length[Flatten[List]]

\[CapitalSigma][PriceList_] := Module[{dS, S, dt}, 
   dS = Drop[Transpose[PriceList][[2]], 1] - 
      Drop[Transpose[PriceList][[2]], -1]; 
    S = Drop[Transpose[PriceList][[2]], -1]; 
    dt = Drop[Transpose[PriceList][[1]], 1] - 
      Drop[Transpose[PriceList][[1]], -1]; 
    Sqrt[Av[dS^2/(S^2*dt)]]]

AEst[PriceList_] := Module[{\[Sigma], T}, \[Sigma] = \[CapitalSigma][PriceList]; 
    T = PriceList[[-1,1]] - PriceList[[1,1]]; 
    \[Sigma]^2/2 + (Log[PriceList[[-1,2]]] - Log[PriceList[[1,2]]])/T]

ConfidenceInterval[a_, \[Alpha]_, \[Sigma]_, T_] := 
  {Quantile[NormalDistribution[a, 
     \[Sigma]/Sqrt[T]], (1 - \[Alpha])/2],
Quantile[NormalDistribution[a, \[Sigma]/Sqrt[T]], 
    1 - (1 - \[Alpha])/2]}

A\[CapitalSigma]Est[PriceList_, p_] := Module[{\[Sigma], T, aest}, 
   \[Sigma] = \[CapitalSigma][PriceList]; T = PriceList[[-1,1]] - 
      PriceList[[1,1]]; aest = \[Sigma]^2/2 + 
      (1/T)*(Log[PriceList[[-1,2]]] - Log[PriceList[[1,2]]]); 
    {{aest, ConfidenceInterval[aest, p, \[Sigma], T]}, \[Sigma]}]

Ave[List_] := Plus@@List/Length[List]

StockCovariance[stocksSDE_]:=Module[{},preS = Transpose[stocksSDE][[2]]; 
preT = Transpose[stocksSDE][[1]]; 
S = Drop[preS, -1]; dS = (Transpose[{#1}] & ) /@ 
    (Drop[preS, 1] - S); 
dT = Drop[preT, 1] - Drop[preT, -1]; 
Ave[(#1 . Transpose[#1] & ) /@ (dS/S)/dT]]

Diag[Matrix_] := Table[Matrix[[i,i]], {i, 1, Length[Matrix]}]

GrowthRateEstimate[stocksSDE_] := 
  Log[stocksSDE[[-1,2]]/stocksSDE[[1,2]]]/
    (stocksSDE[[-1,1]] - stocksSDE[[1,1]]) + 
   (1/2)*Diag[StockCovariance[stocksSDE]]
   
GrowthRatePosterior[stocksSDE_] := MultinormalDistribution[
   GrowthRateEstimate[stocksSDE], StockCovariance[stocksSDE]/
    (stocksSDE[[-1,1]] - stocksSDE[[1,1]])]

GrowthRateEstimateWithPrior[stocksSDE_, ap_, \[Sigma]p_] := 
  Inverse[Inverse[\[Sigma]p . Transpose[\[Sigma]p]] + 
     Inverse[StockCovariance[stocksSDE]/(stocksSDE[[-1,1]] - 
        stocksSDE[[1,1]])]] . 
   (Inverse[\[Sigma]p . Transpose[\[Sigma]p]] . ap + 
    Inverse[StockCovariance[stocksSDE]/(stocksSDE[[-1,1]] - 
        stocksSDE[[1,1]])] . GrowthRateEstimate[stocksSDE])

GrowthRatePosteriorWithPrior[stocksSDE_, ap_, \[Sigma]p_] := 
  MultinormalDistribution[GrowthRateEstimateWithPrior[
    stocksSDE, ap, \[Sigma]p], 
   Inverse[Inverse[\[Sigma]p . Transpose[\[Sigma]p]] + 
     Inverse[StockCovariance[stocksSDE]/(stocksSDE[[-1,1]] - 
        stocksSDE[[1,1]])]]]

End[]

EndPackage[]