intrinsic X2(d::RngIntElt) -> SetIndx
{
Define the set which we denote X_(2d) in the paper. Note the set is 
ordered, so that we can sort subgroups of W_(2d)
}   
  set := [1..d] cat Reverse([-d..-1]);
  set := {@ x : x in set @};
  return set;
end intrinsic;


intrinsic W2(d::RngIntElt) -> GrpPerm, GrpPermElt
{
Define W_(2d) as a subgroup of the symmetric group on X_(2d). The
second return value is the complex conjugation element.
}    
  S2d := SymmetricGroup(X2(d));
  cc := &*[S2d!(i,-i) : i in [1..d]];
  W2d := Centralizer(S2d, cc);
  return W2d, W2d!cc;
end intrinsic;


intrinsic AdmissibleSubgroupClasses(d::RngIntElt) -> SeqEnum[SetEnum], GrpPerm
{
Returns the conjugacy classes of subgroups of W_(2d) which are 
transitive and contain complex conjugation. Second return is W_(2d).
}
  G, cc := W2(d);
  cc := sub<G | cc>;
  subs := Subgroups(G, cc);
  classes := [H : H in subs | IsTransitive(H`subgroup)];
  return [Class(G, H`subgroup) : H in classes], G;
end intrinsic;


intrinsic WeightConjugacyClass(w::SeqEnum, H::GrpPerm) -> SetEnum
{
Given a weighting w, return the w-conjugacy class
}
  d := &+[#entry : entry in w]/2;
  d := Integers()!d;
  W2d := W2(d);
  G := Stabilizer(W2d,w);
  return {Conjugate(H,g) : g in G};
end intrinsic;

function CoerceToW2d(perm)
  n := Degree(Parent(perm));
  d := n div 2;
  map := [1..d] cat [-d..-1];
  return W2(d)!([map[i^perm] : i in [1..n]]);
end function;
