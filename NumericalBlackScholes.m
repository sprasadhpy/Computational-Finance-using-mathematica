(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`NumericalBlackScholes`", "CFMLab`StockStat`", "CFMLab`BlackScholes`", 
"LinearAlgebra`Tridiagonal`"]


ForwardParabolicSolver::usage = "ForwardParabolicSolver[{a2_, a1_, a0_}, 
{f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}, IOrder_:2]
computes the solution of the general one-dimensional forward parabolic PDE."

BackwardParabolicSolver::usage = "BackwardParabolicSolver[{a2_, a1_, a0_}, 
{f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}, IOrder_:2]
computes the solution of the general one-dimensional backward parabolic PDE."

NumericalCallFairPrice::usage = "NumericalCallFairPrice[{\[Sigma]_, r_, D0_}, 
{k_, T_}, {t0_, Nt_}, {s0_, s1_, Ns_}, IOrder_:2] computes the solution of
the Black-Scholes PDE."

NumericalPutFairPrice::usage = "NumericalPutFairPrice[{\[Sigma]_, r_, D0_}, 
{k_, T_}, {t0_, Nt_}, {s0_, s1_, Ns_}, IOrder_:2] computes the solution of
the Black-Scholes PDE."

DupireCallSolution::usage = "DupireCallSolution[{\[Sigma]_, r_}, {S_, t0_}, 
{t1_, Nt_}, {s0_, s1_, Ns_}, IOrder_:2] computes the solution of
the Dupire PDE."

DupirePutSolution::usage = "DupirePutSolution[{\[Sigma]_, r_}, {S_, t0_}, 
{t1_, Nt_}, {s0_, s1_, Ns_}, IOrder_:2] computes the solution of
the Dupire PDE."

DupireObstaclePut::usage = "DupireObstaclePut[{\[Sigma]_, r_, D0_,S_}, 
{x0_, x1_, Nx_}, {t0_, t1_, Nt_}, IOrder_:{1, 2}] computes the solution of
the obstacle problem for the Dupire PDE."

ForwardNonCylindricalParabolicSolver::usage = "
ForwardNonCylindricalParabolicSolver[g_, {b_, a2_, a1_, a0_}, 
   {f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}, IOrder_:1]
computes the solution of the general one-dimensional forward parabolic PDE
in a non-cylindrical domain."

BackwardNonCylindricalParabolicSolver::usage = "
BackwardNonCylindricalParabolicSolver[g_, {b_, a2_, a1_, a0_}, 
   {f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}]
computes the solution of the general one-dimensional backward parabolic PDE
in a non-cylindrical domain."

Begin["`Private`"]


ForwardParabolicSolver[{a2_, a1_, a0_}, {f_, p_}, {s0_, s1_, Ns_}, 
   {t0_, t1_, Nt_}, IOrder_:2] := Module[{ds, ss, dt, tt, Arguments, PPList, 
    A2List, A1List, A0List, RightHandSide, A1PlusList, A1MinusList, 
    LowerDiagonal, Diagonal, UpperDiagonal, AAA, EllipticSolver, 
    ParabolicSolution}, ds = (s1 - s0)/Ns; ss = N[Range[s0, s1, ds]]; 
    dt = (t1 - t0)/Nt; tt = N[Range[t0, t1, dt]]; 
    Arguments = Outer[{#1, #2} & , tt, ss]; 
    PPList = ((p @@ #1 & ) /@ #1 & ) /@ Arguments; 
    A2List = ((a2 @@ #1 & ) /@ #1 & ) /@ Arguments; 
    A1List = ((a1 @@ #1 & ) /@ #1 & ) /@ Arguments; 
    A0List = ((a0 @@ #1 & ) /@ #1 & ) /@ Arguments; 
    RightHandSide = ((f @@ #1 & ) /@ #1 & ) /@ Arguments; 
    A1PlusList = ((Max[#1, 0] & ) /@ #1 & ) /@ A1List; 
    A1MinusList = ((-Min[#1, 0] & ) /@ #1 & ) /@ A1List; 
    LowerDiagonal = (ReplacePart[#1, 0, -1] & ) /@ (Drop[#1, 1] & ) /@ 
       (A2List/ds^2 - A1PlusList/ds); Diagonal = 
     (ReplacePart[#1, 1, -1] & ) /@ (ReplacePart[#1, 1, 1] & ) /@ 
       (A0List + A1MinusList/ds + A1PlusList/ds - 2*(A2List/ds^2) - 1/dt); 
    UpperDiagonal = (ReplacePart[#1, 0, 1] & ) /@ (Drop[#1, -1] & ) /@ 
       (A2List/ds^2 - A1MinusList/ds); 
    AAA = Drop[Transpose[{LowerDiagonal, Diagonal, UpperDiagonal, 
        RightHandSide, PPList}], 1]; EllipticSolver := 
     TridiagonalSolve[#2[[1]], #2[[2]], #2[[3]], 
       ReplacePart[ReplacePart[-#2[[4]] - (1/dt)*#1, #2[[5,1]], 1], 
        #2[[5,-1]], -1]] & ; ParabolicSolution = FoldList[EllipticSolver, 
      PPList[[1]], AAA]; Interpolation[(Append[#1[[1]], #1[[2]]] & ) /@ 
      Transpose[{Flatten[Arguments, 1], Flatten[ParabolicSolution]}], 
     InterpolationOrder -> IOrder]]

BackwardParabolicSolver[{a2_, a1_, a0_}, {f_, p_}, {s0_, s1_, Ns_}, 
   {t0_, t1_, Nt_}, IOrder_:2] := Module[{ForwardSol}, 
   ForwardSol = ForwardParabolicSolver[{a2[t1 - #1, #2] & , 
       a1[t1 - #1, #2] & , a0[t1 - #1, #2] & }, {-f[t1 - #1, #2] & , 
       p[t1 - #1, #2] & }, {s0, s1, Ns}, {0, t1 - t0, Nt}, IOrder]; 
    ForwardSol[t1 - #1, #2] & ]


NumericalCallFairPrice[{\[Sigma]_, r_, D0_}, {k_, T_}, {t0_, Nt_}, 
   {s0_, s1_, Ns_}, IOrder_:2] := Module[{ds, ss, dt, tt, Arguments,
    \[Sigma]\[Sigma], rr}, ds = (s1 - s0)/Ns; ss = N[Range[s0, s1, ds]]; 
    dt = (T - t0)/Nt; tt = N[Range[t0, T, dt]]; 
    Arguments = Outer[{#1, #2} & , tt, ss]; 
    \[Sigma]\[Sigma] = Av[(\[Sigma] @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    rr = Av[(r @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    BackwardParabolicSolver[{(1/2)*#2^2*\[Sigma][#1, #2]^2 & , 
      (r[#1, #2] - D0[#1, #2])*#2 & , -r[#1, #2] & }, 
     {0 & , If[#1 < T, VC[#1, #2, T,k,rr, \[Sigma]\[Sigma]], Max[0, #2 - k]] & }, 
     {s0, s1, Ns}, {t0, T, Nt}, IOrder]]

NumericalPutFairPrice[{\[Sigma]_, r_, D0_}, {k_, T_}, {t0_, Nt_}, 
   {s0_, s1_, Ns_}, IOrder_:2] := Module[{ds, ss, dt, tt, Arguments,
    \[Sigma]\[Sigma], rr}, ds = (s1 - s0)/Ns; ss = N[Range[s0, s1, ds]]; 
    dt = (T - t0)/Nt; tt = N[Range[t0, T, dt]]; 
    Arguments = Outer[{#1, #2} & , tt, ss];  
    \[Sigma]\[Sigma] = Av[(\[Sigma] @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    rr = Av[(r @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    BackwardParabolicSolver[{(1/2)*#2^2*\[Sigma][#1, #2]^2 & , 
      (r[#1, #2] - D0[#1, #2])*#2 & , -r[#1, #2] & }, 
     {0 & , If[#1 < T, VP[#1, #2, k, T, rr, \[Sigma]\[Sigma]], Max[0, k - #2]] & }, 
     {s0, s1, Ns}, {t0, T, Nt}, IOrder]]

DupireCallSolution[{\[Sigma]_, r_}, {S_, t0_}, {t1_, Nt_}, {s0_, s1_, Ns_}, 
   IOrder_:2] := Module[{ds, ss, dt, tt, Arguments, \[Sigma]\[Sigma], rr}, 
   ds = (s1 - s0)/Ns; ss = N[Range[s0, s1, ds]]; dt = (t1 - t0)/Nt; 
    tt = N[Range[t0, t1, dt]]; Arguments = Outer[{#1, #2} & , tt, ss]; 
    \[Sigma]\[Sigma] = Av[(\[Sigma] @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    rr = Av[(r @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    ForwardParabolicSolver[{(1/2)*#2^2*\[Sigma][#1, #2]^2 & , (-r[#1, #2])*#2 & , 
      0 & }, {0 & , If[#1 > t0, VC[t0, S, #1,#2, rr,\[Sigma]\[Sigma]], 
        Max[0, S - #2]] & }, {s0, s1, Ns}, {t0, t1, Nt}, IOrder]]

DupirePutSolution[{\[Sigma]_, r_}, {S_, t0_}, {t1_, Nt_}, {s0_, s1_, Ns_}, 
   IOrder_:2] := Module[{ds, ss, dt, tt, Arguments, \[Sigma]\[Sigma],rr}, 
   ds = (s1 - s0)/Ns; ss = N[Range[s0, s1, ds]]; dt = (t1 - t0)/Nt; 
    tt = N[Range[t0, t1, dt]]; Arguments = Outer[{#1, #2} & , tt, ss]; 
    \[Sigma]\[Sigma] = Av[(\[Sigma] @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    rr = Av[(r @@ #1 & ) /@ Flatten[Arguments, 1]]; 
    ForwardParabolicSolver[{(1/2)*#2^2*\[Sigma][#1, #2]^2 & , (-r[#1, #2])*#2 & , 
      0 & }, {0 & , If[#1 > t0, VP[t0, S, #1, #2, rr, \[Sigma]\[Sigma]], 
        Max[0,  #2 - S]] & }, {s0, s1, Ns}, {t0, t1, Nt}, IOrder]]






DupireObstaclePut[{\[Sigma]_, r_, D0_,S_}, {x0_, x1_, Nx_}, {t0_, t1_, Nt_}, 
   IOrder_:{1, 2}] := Module[{a2, a1, a0, p, ds, ss, dt, tt, Arguments, 
    PPList, A2List, A1List, A0List, A1PlusList, A1MinusList, Prev, FB0, 
    LowerDiagonal, Diagonal, UpperDiagonal, FoldingList, 
    MaximalBoundaryValueSolver, BoundaryValueSolver, ReplaceAtPlaces, 
    DupireObstacleProblemSolution}, a2[T_, k_] := (1/2)*\[Sigma][T, k]^2*k^2; 
    a1[T_, k_] := (-(r[T, k] - D0[T, k]))*k; a0[T_, k_] := -D0[T, k]; 
    p[T_, k_] := Max[k - S, 0]; ds = (x1 - x0)/Nx; 
    ss = N[Range[x0, x1, ds]]; dt = (t1 - t0)/Nt; 
    tt = N[Range[t0, t1, dt]]; Arguments = Outer[{#1, #2} & , tt, ss]; 
    PPList = ((p @@ #1 & ) /@ #1 & ) /@ Arguments; 
    A2List = ((a2 @@ #1 & ) /@ #1 & ) /@ Drop[Arguments, 1]; 
    A1List = ((a1 @@ #1 & ) /@ #1 & ) /@ Drop[Arguments, 1]; 
    A0List = ((a0 @@ #1 & ) /@ #1 & ) /@ Drop[Arguments, 1]; 
    A1PlusList = ((Max[#1, 0] & ) /@ #1 & ) /@ A1List; 
    A1MinusList = ((-Min[#1, 0] & ) /@ #1 & ) /@ A1List; 
    Prev = (p @@ #1 & ) /@ Arguments[[1]]; 
    FB0 = {Select[ss, #1 >= S & ][[1]], 
      Position[ss, Select[ss, #1 >= S & ][[1]]][[1,1]]}; 
    LowerDiagonal = (ReplacePart[#1, 0, -1] & ) /@ 
      (Drop[#1, 1] & ) /@ (A2List/ds^2 - A1PlusList/ds); 
    Diagonal = (ReplacePart[#1, 1, 1] & ) /@ 
      (ReplacePart[#1, 1, -1] & ) /@ (-2*(A2List/ds^2) + A1PlusList/ds + 
        A1MinusList/ds + A0List - 1/dt); UpperDiagonal = 
     (ReplacePart[#1, 0, 1] & ) /@ (Drop[#1, -1] & ) /@ 
       (A2List/ds^2 - A1MinusList/ds); FoldingList = 
     Transpose[{LowerDiagonal, Diagonal, UpperDiagonal, 
       Drop[PPList, 1]}]; MaximalBoundaryValueSolver[FB0_, 
      LowerDiagonal_, Diagonal_, UpperDiagonal_, Prev_] := 
     Module[{y, x}, y = {}; 
       Do[x = {FB0 + {i*ds, i}, BoundaryValueSolver[FB0 + {i*ds, i}, 
            LowerDiagonal, Diagonal, UpperDiagonal, Prev]}; 
         y = Append[y, x]; If[i >= 1 && y[[-2,2]][[FB0[[2]]]] > 
            y[[-1,2]][[FB0[[2]]]], Break[]], {i, 0, 1000}]; y[[-2]]]; 
    BoundaryValueSolver[FB0_, LowerDiagonal_, Diagonal_, UpperDiagonal_, 
      Prev_] := Module[{ModifiedLowerDiagonal, ModifiedDiagonal, 
       ModifiedUpperDiagonal, RHS, ModifiedRHS, k}, 
      k = FB0[[2]]; RHS = -(Prev/dt); ModifiedLowerDiagonal = 
        ReplaceAtPlaces[LowerDiagonal, Table[0, 
          {Length[LowerDiagonal] + 1 - k}], Range[k - 1, 
          Length[LowerDiagonal] - 1]]; ModifiedDiagonal = 
        ReplaceAtPlaces[Diagonal, Table[1, {Length[Diagonal] - k}], 
         Range[k, Length[Diagonal] - 1]]; ModifiedUpperDiagonal = 
        ReplaceAtPlaces[UpperDiagonal, Table[0, 
          {Length[UpperDiagonal] + 1 - k}], Range[k, 
          Length[UpperDiagonal]]]; ModifiedRHS = ReplaceAtPlaces[RHS, 
         (Prev[[#1]] & ) /@ Range[k, Length[RHS]], 
         Range[k, Length[RHS]]]; TridiagonalSolve[ModifiedLowerDiagonal, 
        ModifiedDiagonal, ModifiedUpperDiagonal, ModifiedRHS]]; 
    ReplaceAtPlaces[expr_, new_, places_] := 
     Fold[ReplacePart[#1, #2[[1]], #2[[2]]] & , expr, 
      Transpose[{new, places}]]; DupireObstacleProblemSolution = 
     FoldList[MaximalBoundaryValueSolver[#1[[1]], #2[[1]], #2[[2]], 
        #2[[3]], #1[[2]]] & , {FB0, PPList[[1]]}, FoldingList]; 
    DupirePut = Interpolation[({#1[[1,1]], #1[[1,2]], #1[[2]]} & ) /@ 
       Transpose[{(Flatten[#1, 1] & )[Arguments], 
         Flatten[(#1[[2]] & ) /@ DupireObstacleProblemSolution]}], 
      InterpolationOrder -> IOrder]; FreeBoundary = 
     Interpolation[Transpose[{tt, (#1[[1,1]] & ) /@ 
         DupireObstacleProblemSolution}], InterpolationOrder -> 1]; 
    {FreeBoundary, DupirePut}]





ForwardNonCylindricalParabolicSolver[g_, {b_, a2_, a1_, a0_}, 
   {f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}, IOrder_:1] := 
  Module[{ds, ss, dt, tt, Arguments, PPList, A2List, A1List, 
    A0List, RightHandSide, A1PlusList, A1MinusList, 
    LowerDiagonal, Diagonal, UpperDiagonal, AAA, EllipticSolver, 
    ParabolicSolution}, ds = (s1 - s0)/Ns; 
    ss = N[Range[s0, s1, ds]]; dt = (t1 - t0)/Nt; 
    tt = N[Range[t0, t1, dt]]; Arguments = Outer[{#1, #2} & , tt, 
      ss]; PPList = ((p @@ #1 & ) /@ #1 & ) /@ Arguments; 
    zxc = Function[X, Select[X, #1[[2]] >= g[#1[[1]]] & ]] /@ 
      Arguments; DirichletCondAt = 
     Drop[(Union[#1, {1}, {Ns + 1}] & ) /@ 
       Table[Flatten[(Position[Arguments[[i]], #1] & ) /@ 
          zxc[[i]], 2], {i, 1, Length[zxc]}], 1]; 
    DirichletCondAt2 = Flatten[Table[({i, #1} & ) /@ 
        DirichletCondAt[[i]], {i, 1, Length[DirichletCondAt]}], 
      1]; BCoefList = ((b @@ #1 & ) /@ #1 & ) /@ 
      Drop[Arguments, 1]; A2List = ((a2 @@ #1 & ) /@ #1 & ) /@ 
      Drop[Arguments, 1]; A1List = ((a1 @@ #1 & ) /@ #1 & ) /@ 
      Drop[Arguments, 1]; A0List = ((a0 @@ #1 & ) /@ #1 & ) /@ 
      Drop[Arguments, 1]; RightHandSide = 
     ((f @@ #1 & ) /@ #1 & ) /@ Drop[Arguments, 1]; 
    A1PlusList = ((Max[#1, 0] & ) /@ #1 & ) /@ A1List; 
    A1MinusList = ((-Min[#1, 0] & ) /@ #1 & ) /@ A1List; 
    LowerDiagonal = (Drop[#1, 1] & ) /@ 
      Fold[ReplacePart[#1, 0, #2] & , -(A1PlusList/ds) + 
        A2List/ds^2, DirichletCondAt2]; 
    Diagonal = Fold[ReplacePart[#1, 1, #2] & , 
      A0List + A1MinusList/ds + A1PlusList/ds - (2*A2List)/ds^2 - 
       BCoefList/dt, DirichletCondAt2]; UpperDiagonal = 
     (Drop[#1, -1] & ) /@ Fold[ReplacePart[#1, 0, #2] & , 
       -(A1PlusList/ds) + A2List/ds^2, DirichletCondAt2]; 
    AAA = Transpose[{LowerDiagonal, Diagonal, UpperDiagonal, 
       RightHandSide, Drop[PPList, 1], BCoefList, 
       DirichletCondAt}]; EllipticSolver := 
     Function[{X, Y}, TridiagonalSolve[Y[[1]], Y[[2]], Y[[3]], 
       Fold[ReplacePart[#1, Y[[5,#2]], #2] & , 
        -Y[[4]] - (Y[[6]]/dt)*X, Y[[7]]]]]; 
    ParabolicSolution = FoldList[EllipticSolver, PPList[[1]], 
      AAA]; Interpolation[(Append[#1[[1]], #1[[2]]] & ) /@ 
      Transpose[{Flatten[Arguments, 1], 
        Flatten[ParabolicSolution]}], InterpolationOrder -> 
      IOrder]]
      
BackwardNonCylindricalParabolicSolver[g_, {b_, a2_, a1_, a0_}, 
   {f_, p_}, {s0_, s1_, Ns_}, {t0_, t1_, Nt_}] := 
  Module[{ForwardSol}, 
   ForwardSol = ForwardNonCylindricalParabolicSolver[
      g[t1 - #1] & , {b[t1 - #1, #2] & , a2[t1 - #1, #2] & , 
       a1[t1 - #1, #2] & , a0[t1 - #1, #2] & }, 
      {-f[t1 - #1, #2] & , p[t1 - #1, #2] & }, {s0, s1, Ns}, 
      {0, t1 - t0, Nt}]; ForwardSol[t1 - #1, #2] & ]

End[]

EndPackage[]