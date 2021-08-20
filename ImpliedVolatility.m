(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`ImpliedVolatility`", "CFMLab`StockStat`", "CFMLab`BlackScholes`",
"Global`"]


IndividualCallVolatilities::usage = "IndividualCallVolatilities[aod_, r_]
computes impled volatilities corresponding to particular call options."

IndividualPutVolatilities::usage = "IndividualPutVolatilities[aod_, r_]
computes impled volatilities corresponding to particular call options."

ConstantCallVolatility::usage = "ConstantCallVolatility[aod_, r_, (p_)?EvenQ]
computes the constant implied volatility corresponding to all call options,
by minimizing the least p-distance between the theoretical and observed prices."

ConstantPutVolatility::usage = "ConstantPutVolatility[aod_, r_, (p_)?EvenQ]
computes the constant implied volatility corresponding to all put options,
by minimizing the least p-distance between the theoretical and observed prices."

TimeDependentCallVolatility::usage = "TimeDependentCallVolatility[aod_, qq_, r_, (p_)?EvenQ]
computes the time dependent implied volatility corresponding to all call options,
by minimizing the least p-distance between the theoretical and observed prices."

TimeDependentPutVolatility::usage = "TimeDependentPutVolatility[aod_, qq_, r_, (p_)?EvenQ]
computes the time dependent implied volatility corresponding to all put options,
by minimizing the least p-distance between the theoretical and observed prices."

Begin["`Private`"]


IndividualCallVolatilities[aod_, r_] := 
  Function[y, 
    (\[Sigma] /. FindRoot[{aod[[2,y,#1,3]] == VC[aod[[1,4]], 
           aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], r, 
           \[Sigma]]}, {\[Sigma], 0.7}] & ) /@ 
     Range[1, Length[aod[[2,y]]]]] /@ 
   Range[1, Length[aod[[2]]]]

IndividualPutVolatilities[aod_, r_] := 
  Function[y, 
    (\[Sigma] /. FindRoot[{aod[[2,y,#1,3]] == VP[aod[[1,4]], 
           aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], 
           0.05, \[Sigma]]}, {\[Sigma], 0.8}] & ) /@ 
     Range[1, Length[aod[[2,y]]]]] /@ 
   Range[1, Length[aod[[2]]]]

ConstantCallVolatility[aod_, r_] := ConstantCallVolatility[
   aod, r, 2]

ConstantCallVolatility[aod_, r_, (p_)?EvenQ] := 
  Module[{f}, f = (Plus @@ #1 & )[((#1^p & ) /@ #1 & )[
       Flatten[Function[y, (aod[[2,y,#1,3]] - VC[aod[[1,4]], 
              aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], 
              r, \[Sigma]] & ) /@ Range[1, Length[aod[[2,y]]]]] /@ 
         Range[1, Length[aod[[2]]]]]]]; 
    \[Sigma] /. FindMinimum[f, {\[Sigma], 0.5}][[2]]]

ConstantPutVolatility[aod_, r_] := ConstantPutVolatility[aod, 
   r, 2]

ConstantPutVolatility[aod_, r_, (p_)?EvenQ] := 
  Module[{f}, f = (Plus @@ #1 & )[((#1^p & ) /@ #1 & )[
       Flatten[Function[y, (aod[[2,y,#1,3]] - VP[aod[[1,4]], 
              aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], 
              r, \[Sigma]] & ) /@ Range[1, Length[aod[[2,y]]]]] /@ 
         Range[1, Length[aod[[2]]]]]]]; 
    \[Sigma] /. FindMinimum[f, {\[Sigma], 0.5}][[2]]]

TimeDependentCallVolatility[aod_, qq_, r_] := 
  TimeDependentCallVolatility[aod, qq, r, 2]

Step[t_, V_, T_] := V[[1]] + 
   Plus @@ ((#1[[2]] - #1[[1]] & ) /@ Partition[V, 2, 1]*
     (UnitStep[t - #1] & ) /@ T)

P[k_] := ToExpression /@ Take[CharacterRange["a", "z"], -k]

S[k_][t_] := Step[t, P[k], Take[aod[[1,3]], k - 1]]



TimeDependentCallVolatility[aod_, qq_, r_, (p_)?EvenQ] := 
  Module[{ss, w, f}, kk = Min[qq, Length[aod[[2]]]]; 
    w[t_, S_, T_, k_, rr_] = GeneralCall[t, S, T, k, rr & , 
      rr & , S[kk]]; f = (Plus @@ #1 & )[((#1^p & ) /@ #1 & )[
       Flatten[Function[y, (aod[[2,y,#1,3]] - w[aod[[1,4]], 
              aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], 
              r] & ) /@ Range[1, Length[aod[[2,y]]]]] /@ 
         Range[1, Length[aod[[2]]]]]]]; 
    ss[t_] = S[kk][t] /. FindMinimum[Evaluate[f], 
        Evaluate[Sequence @@ ({#1, 0.5} & ) /@ P[kk]]][[2]]; ss]


TimeDependentPutVolatility[aod_, qq_, r_] := 
  TimeDependentPutVolatility[aod, qq, r, 2]



TimeDependentPutVolatility[aod_, qq_, r_, (p_)?EvenQ] := 
  Module[{ss, w, f}, kk = Min[qq, Length[aod[[2]]]]; 
    w[t_, S_, T_, k_, rr_] = GeneralPut[t, S, T, k, rr & , 
      rr & , S[kk]]; f = (Plus @@ #1 & )[((#1^p & ) /@ #1 & )[
       Flatten[Function[y, (aod[[2,y,#1,3]] - w[aod[[1,4]], 
              aod[[1,5]], aod[[2,y,#1,1]], aod[[2,y,#1,2]], 
              r] & ) /@ Range[1, Length[aod[[2,y]]]]] /@ 
         Range[1, Length[aod[[2]]]]]]]; 
    ss[t_] = S[kk][t] /. FindMinimum[Evaluate[f], 
        Evaluate[Sequence @@ ({#1, 0.5} & ) /@ P[kk]]][[2]]; 
    ss]


End[]

EndPackage[]
