(* Glaisher-Clausen function *)
Sl[m_, x_] := (-1)^(m/2-1) (2Pi)^m/(2 Factorial[m]) BernoulliB[m, x/(2Pi)] /; EvenQ[m];

Sl[m_, x_] := (-1)^((m+1)/2) (2Pi)^m/(2 Factorial[m]) BernoulliB[m, x/(2Pi)] /; !EvenQ[m];
