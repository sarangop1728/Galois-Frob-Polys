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

intrinsic ExistsAdmissibleDecompositions(G::GrpPerm,wt::SeqEnum) -> Bool
  v := wt[1];
  w := wt[2];
  d := &+[#entry : entry in w]/2;
  d := Integers()!d;
  W2d := W2(d);
  Stab := Stabilizer(W2d,w);
  for D in Subgroups(Stab) do
      Ds := D`subgroup;
      validD := true;
      for i in [1 .. d] do
          sum := 0;
          for elt in Orbit(Ds,i) do
              if elt in [1 .. d] then
                  sum := sum + v[elt];
              else
                  sum := sum + v[5 - elt];
              end if;
          end for;
          if not sum in Integers() then
              validD := false;
          end if;
            end for;
      if validD then
          if IsCyclic(Ds) then
              return true;
          else
              for G0 in NormalSubgroups(Ds) do
                    G0s := G0`subgroup;
                    for G1 in NormalSubgroups(G0s) do
                        G1s := G1`subgroup;
                        if IsCyclic(quo<Ds | G0s>) and IsCyclic(quo<G0s | G1s>) and GCD(#quo<G0s | G1s>, #G1s) eq 1 and IsPrimePower(#G1s) then
                           return true;
                        end if;
                    end for;
              end for;
          end if;
      end if;
  end for;
  return false;
end intrinsic;
