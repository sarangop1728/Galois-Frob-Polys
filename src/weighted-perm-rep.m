function v(d,g)
  return Vector([x^g : x in {1..d}]);
end function;

function permMatrix(d,G)
  return Transpose(Matrix([v(d,g): g in G]));
end function;

function MapToSlopes(d,NP)
  graph := &join[{<i,NP[i] - NP[2*d + 1 - i]>, <-i,-NP[i] + NP[2*d + 1 - i]>} : i in [1..d]];
  return map<X2(d) -> Rationals() | graph>;
end function;

intrinsic NewtonMatrix(G::GrpPerm, NP::SeqEnum) -> AlgMatElt
{
Return the Newton hyperplane matrix corresponding to the pair G and the 
Newton polygon NP (a sequence of Newton slopes)
}
  d := #NP div 2;
  f := MapToSlopes(d,NP);
  return Transpose(Matrix([Vector([f(x^g) : x in {1..d}]): g in G]));
end intrinsic;

intrinsic AngleRank(G::GrpPerm, NP::SeqEnum) -> RngIntElt
{
Return the angle rank of the pair G and the 
Newton polygon NP (a sequence of Newton slopes)
}
  return Rank(NewtonMatrix(G, NP));
end intrinsic;
