(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`BlackScholes`"]

VC::usage = "VC[t, S, T, k, r (or a, b), \[Sigma]] computes
the fair price of an European call option at time t, at stock-price S, for
the expiration time T, for the strike price k, when the interest rate is r,
and when the stock-volatility is \[Sigma]."

VP::usage = "VP[t, S, T, k, r (or a, b), \[Sigma]] computes
the fair price of an European put option at time t, at stock-price S, for
the expiration time T, for the strike price k, when the interest rate is r,
and when the stock-volatility is \[Sigma]."

GeneralCall::usage = "GeneralCall[t, S, T, k, a, b, \[Sigma]] computes
the fair price of an European call option at time t, at stock-price S, for
the expiration time T, for the strike price k, when the parameters a, b,
and the stock-volatility \[Sigma] are all pure functions of time."

GeneralPut::usage = "GeneralPut[t, S, T, k, a, b, \[Sigma]] computes
the fair price of an European put option at time t, at stock-price S, for
the expiration time T, for the strike price k, when the parameters a, b,
and the stock-volatility \[Sigma] are all pure functions of time."


Begin["`Private`"]

VC[t_, S_, T_, k_, r_, \[Sigma]_] := 
  (1/2)*(Erf[(2*Log[S/k] - (t - T)*(\[Sigma]^2 + 2*r))/
       (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])]*S + S + 
    E^(r*(t - T))*k*(Erfc[(2*Log[S/k] - (t - T)*(2*r - \[Sigma]^2))/
        (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])] - 2))

VC[t_, S_, T_, k_, a_, b_, \[Sigma]_] := (1/2)*E^(b*(t - T))*
   (E^(a*(T - t))*S*(Erf[(2*Log[S/k] - (t - T)*(\[Sigma]^2 + 2*a))/
        (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])] + 1) + 
    k*(Erfc[(2*Log[S/k] - (t - T)*(2*a - \[Sigma]^2))/(2*Sqrt[2]*
         Sqrt[(T - t)*\[Sigma]^2])] - 2))

VP[t_, S_, T_, k_, r_, \[Sigma]_] := 
  (1/2)*(E^(r*(t - T))*k*Erfc[(2*Log[S/k] - (t - T)*(2*r - \[Sigma]^2))/
       (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])] - 
    S*Erfc[(2*Log[S/k] - (t - T)*(\[Sigma]^2 + 2*r))/
       (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])])

VP[t_, S_, T_, k_, a_, b_, \[Sigma]_] := (1/2)*E^(b*(t - T))*
   (k*Erfc[(2*Log[S/k] - (t - T)*(2*a - \[Sigma]^2))/
       (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])] - E^(a*(T - t))*S*
     Erfc[(2*Log[S/k] - (t - T)*(\[Sigma]^2 + 2*a))/
       (2*Sqrt[2]*Sqrt[(T - t)*\[Sigma]^2])])

GeneralCall[t_, S_, T_, k_, a_, b_, \[Sigma]_] := 
  ((1/2)*(E^Integrate[a[\[Tau]], {\[Tau], t, T}]*S*
      (Erf[(2*Integrate[a[\[Tau]], {\[Tau], t, T}] + \
Integrate[\[Sigma][\[Tau]]^2, 
           {\[Tau], t, T}] + \
2*Log[S/k])/(2*Sqrt[2]*Sqrt[Integrate[\[Sigma][\
\[Tau]]^2, 
            {\[Tau], t, T}]])] + 1) - 
     k*Erfc[(-2*Integrate[a[\[Tau]], {\[Tau], t, T}] + Integrate[\[Sigma][\
\[Tau]]^2, 
          {\[Tau], t, T}] - 2*Log[S/k])/(2*Sqrt[2]*
         Sqrt[Integrate[\[Sigma][\[Tau]]^2, {\[Tau], t, T}]])]))/
   E^Integrate[b[\[Tau]], {\[Tau], t, T}]

GeneralPut[t_, S_, T_, k_, a_, b_, \[Sigma]_] := 
  ((1/2)*(k*Erfc[-((-2*Integrate[a[\[Tau]], {\[Tau], t, T}] + Integrate[\
\[Sigma][\[Tau]]^2, 
           {\[Tau], t, T}] - \
2*Log[S/k])/(2*Sqrt[2]*Sqrt[Integrate[\[Sigma][\
\[Tau]]^2, 
            {\[Tau], t, T}]]))] - E^Integrate[a[\[Tau]], {\[Tau], t, T}]*S*
      Erfc[(2*Integrate[a[\[Tau]], {\[Tau], t, T}] + \
Integrate[\[Sigma][\[Tau]]^2, 
          {\[Tau], t, T}] + 2*Log[S/k])/(2*Sqrt[2]*
         Sqrt[Integrate[\[Sigma][\[Tau]]^2, {\[Tau], t, T}]])]))/
   E^Integrate[b[\[Tau]], {\[Tau], t, T}]

End[]

EndPackage[]