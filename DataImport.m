(* Copyright 2002 by Srdjan Stojanovic  *)


BeginPackage["CFMLab`DataImport`"]

YearLength::usage = "YearLength is a constant, an average length of a year in seconds."

YearClock::usage = "YearClock[y_] is date (in years) since the year 2000."

GetPriceData::usage = "GetPriceData[Ticker_] imports the closing stock price data for
the chosen stock, from the MarketData directory."

GetAllData::usage = "GetAllData[Ticker_,OptionType_,Vol_:1] combines all the
option data and the stock data."

Begin["`Private`"]

SetDirectory[
  ToFileName[{$TopDirectory, "AddOns", "Applications", "CFMLab", 
      "MarketData"}]];

YearLength = (1/100)*FromDate[{2000, 1, 1, 0, 0, 0}];
 
YearClock[y_] := 
   N[(FromDate[y] - FromDate[{2000, 1, 1, 0, 0, 0}])/
     YearLength]; 

GetPriceData[Ticker_] := 
  Get[ToFileName[{$TopDirectory, "AddOns", "Applications", 
     "CFMLab", "MarketData"}, StringJoin[ToString[Ticker], 
     "PriceData.txt"]]]



GetAllData["QQQ", OptionType_, Volume_:300, NumberOfExpirations_:3] := 
  Module[{InputData, AllData, F, AllData2, AllCallData, AllCallData2, 
    AllPutData, AllPutData2, x, y, ExpTimes, z}, 
   InputData = InputData = OpenRead[ToFileName[{$TopDirectory, "AddOns", 
         "Applications", "CFMLab", "MarketData"}, 
        "OptionMetricsData2.txt"]]; 
    AllData = Drop[ReadList[InputData, Table[Word, {9}]], -1]; 
    Close[InputData]; F := Module[{xx, yy, zz}, 
       xx = Partition[Append[Prepend[(#1[[1]] & ) /@ StringPosition[#1, 
              "/"], 0], StringLength[#1] + 1], 2, 1]; 
        yy = ToExpression[Function[y, StringTake[#1, y + {1, -1}]] /@ 
           xx]; zz = ({If[#1[[3]] > 50, #1[[3]] + 1900, #1[[3]] + 2000], 
            #1[[1]], #1[[2]], 16, 0, 0} & )[yy]; YearClock[zz]] & ; 
    AllData2 = ({#1[[1]], #1[[2]], F[#1[[2]]], #1[[3]], F[#1[[3]]], 
        ToExpression[#1[[4]]], #1[[5]], ToExpression[#1[[6]]], 
        ToExpression[#1[[7]]], ToExpression[#1[[8]]], 
        ToExpression[#1[[9]]], N[(1/2)*(ToExpression[#1[[6]]] + 
           ToExpression[#1[[7]]])]} & ) /@ Drop[AllData, 1]; 
    AllCallData = Select[AllData2, #1[[7]] == "C" & ]; 
    AllPutData = Select[AllData2, #1[[7]] == "P" & ]; 
    AllCallData2 = Transpose[(Transpose[AllCallData][[#1]] & ) /@ 
       {5, 6, -1, -3}]; AllPutData2 = Transpose[
      (Transpose[AllPutData][[#1]] & ) /@ {5, 6, -1, -3}]; 
    x = If[ToString[OptionType] === "Call", AllCallData2, AllPutData2]; 
    y = Select[x, #1[[-1]] >= Volume & ]; 
    ExpTimes = Union[Transpose[y][[1]]]; 
    ExpTimes = Take[ExpTimes, Min[NumberOfExpirations, 
       Length[ExpTimes]]]; z = Function[t, Select[y, #1[[1]] == t & ]] /@ 
      ExpTimes; {{"QQQ", ToString[OptionType], ExpTimes, 
      0.8314167853101158, 76.765617}, z}]

End[]

EndPackage[]
