function conjugateToCorrectComplexConj(cc)
  /*
  Input: cc = (1,2)(3,5)(4,6) in S6.
  Output: Permutation s such that s*cc*s^(-1) = (1,4)(2,5)(3,6).
  */
  n := Degree(Parent(cc));
  Sn := SymmetricGroup(n);
  d := n div 2;
  s := Identity(Sn);
  for j in [1..d] do
    l := [1..n]; l[d+j] := j^cc; l[j^cc] := d+j;
    t := Sn!l;
    s := t*s;
    cc := t*cc*t^(-1);
  end for;
  return s;
end function;


function complexConjCandidates(H)
  /*
  Returns the elements of H of the same cycle type as cc.
  */
  d := Degree(H) div 2;
  return [g : g in H | CycleStructure(g)[1] eq <2,d> ];
end function;


function isComplexConj(roots,c,q)
  /*
  Determines if the candidate c is complex conjugation, for a given list of
  roots and prime power q.
  */
  n := #roots;
  return {roots[i]*roots[i^c]/q eq 1 : i in [1..n]} eq {true};
end function;

function permutationFromOrderedRoots(g,tau,roots)
  /*
  Given an element of the Galois group G (which acts on the set of roots),
  return the permutation corresponding to the ordering of the roots.
  */
  return [Index(roots, tau(g)(r)): r in roots];
end function;


function PermutationToW2d(perm)
  n := Degree(Parent(perm));
  d := n div 2;
  X := X2(d);
  W := W2(d);
  return W!{@ X[i^perm] : i in [1..n] @};
end function;


intrinsic NewtonGaloisGroup(Poly::RngUPolElt) -> SeqEnum, GrpPerm, GrpPermElt, FldNum
{
  Input: a squarefree q-Weil polynomial P.
  Output:
  1. The list of (roots,valuations) of P in the splitting field L, with respect
  to an indexing of the roots (\bar(j) = 2d+1-j, and nondecreasing valuations
  with respect to some prime pp above p in L).
  2. The Galois group H of P, as a permutation group on the roots of P, with
  respect to the indexing chosen in (1).
  3. The complex conjugation element cc of H.
  4. The splitting field L of P.
}
  Pol<T> := PolynomialRing(Rationals());
  P := Pol!Poly;
  assert IsSquarefree(P);
  n := Degree(P);
  d := n div 2;
  _, q := IsPower(Integers()!Coefficient(P,0), Degree(P) div 2);
  p := PrimeDivisors(q)[1];
  // Splitting field and Galois group.
  L<a> := SplittingField(P);
  G, _, tau := AutomorphismGroup(L);
  roots := [r[1]: r in Roots(P,L)];
  // Order according to valulation.
  pOL := ideal<Integers(L) | p>;
  pp := Factorisation(pOL)[1][1];
  roots_valuations := [<r,Valuation(r,pp)> :  r in roots];
  Sort(~roots_valuations, func< a, b | a[2] - b[2]>);
  roots := [rv[1]: rv in roots_valuations];
  valuations := [rv[2]: rv in roots_valuations];
  // Get the subgroup of Sn corresponding to G, according to this ordering of the
  // roots.
  Sn := SymmetricGroup(n);
  gens := [];
  for g in Generators(G) do
    Append(~gens, permutationFromOrderedRoots(g,tau,roots));
  end for;
  H := sub<Sn | gens>;
  // Get the cc of H.
  cc := Identity(H);
  for c in complexConjCandidates(H) do
    if isComplexConj(roots,c,q) then
      cc := cc*c;
    end if;
  end for;
  s := conjugateToCorrectComplexConj(cc);
  H := Conjugate(H,s);
  cc := H!(s*cc*s^(-1));
  // Reorder the list of roots to coincide with the permutation action of H.
  roots_valuations := [roots_valuations[i^s]: i in [1..n]];
  // Coerce H to W2(d).
  W := W2(d);
  H := sub<W | [PermutationToW2d(h): h in Generators(H)]>;
  cc := H!PermutationToW2d(cc);
return roots_valuations, H, cc, L;
end intrinsic;
