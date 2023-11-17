(* Standard Clausen function *)
Cl[m_, x_] := (PolyLog[m, Exp[I x]] + PolyLog[m, Exp[-I x]])/2 /; !EvenQ[m];

Cl[m_, x_] := (PolyLog[m, Exp[I x]] - PolyLog[m, Exp[-I x]])/(2 I) /; EvenQ[m];
