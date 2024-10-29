alphabet := "abcdefghijklmnopqrstuvwxyz";
alphabet := [alphabet[i] : i in [1..26]];
TNums := TransitiveGroupIdentification;

//////////////////////////////////////////////////
// Start with subgroup sorting
//////////////////////////////////////////////////

function IsLex(l1, l2)
  // Return true if and only if l1 <= l2 in lex order.
  // Assume that l1 and l2 are sequences of positive integers.
  if #l1 lt #l2 then
    l1 cat:= [0 : i in [1..#l2-#l1]];
  elif #l2 lt #l1 then
    l2 cat:= [0 : i in [1..#l1-#l2]];
  end if;

  for i in [1..#l1] do
    if l1[i] lt l2[i] then
      return true;
    elif l2[i] lt l1[i] then
      return false;
    end if;
  end for;

  return true;
end function;


function IsLexLists(L1, L2)
  // Return true if and only if L1 <= L2 in lex order.
  // Here L1 and L2 are lists of lists 
  if #L1 lt #L2 then
    L1 cat:= [[] : i in [1..#L2-#L1]];
  elif #L2 lt #L1 then
    L2 cat:= [[] : i in [1..#L1-#L2]];
  end if;
 
  for i in [1..#L1] do
    if IsLex(L1[i], L2[i]) and (L1[1] ne L2[i]) then
      return true;
    elif IsLex(L2[i], L1[i]) and (L1[1] ne L2[i]) then
      return false;
    end if;
  end for;

  return true;
end function;


intrinsic LexCompareSubgroups(G::GrpPerm, H::GrpPerm) -> BoolElt
{Given 2 subgroups of the same ambient symmetric group, return 
true if G <= H under some well defined ordering from
CanonicalInvariant
}
  GG := StandardGroup(G);
  HH := StandardGroup(H);

  LG := CanonicalInvariant(GG);
  LH := CanonicalInvariant(HH);

  return IsLexLists(LG, LH);  
end intrinsic;


function LexSrt(G, H)
  flag := LexCompareSubgroups(G, H);
  if G eq H then
    return 0;
  elif flag then
    return -1;
  else
    return 1;
  end if;
end function;


function OrderThenTransitiveThenLexSrtList(GG, HH)
  // assume GG and HH come in sorted and each group in GG (resp. HH)
  // is isomorphic. Sort first by order (bigger is better) then transitive
  // label (bigger is better)
  if #GG[1] eq #HH[1] then
    if IsIsomorphic(GG[1], HH[1]) then
      return LexSrt(GG[1], HH[1]);
    else
      x1,y1 := TNums(GG[1]);
      x2,y2 := TNums(HH[1]);

      if y1 ne y2 then
        return -(y1 - y2);
      else
        return -(x1 - x2);
      end if;
    end if;
  else
    return -(#GG[1] - #HH[1]);
  end if;
end function;


intrinsic SortSubgroups(GG::SeqEnum[GrpPerm]) -> SeqEnum[GrpPerm]
{
Sorts GG
}
  return Sort(GG, LexSrt);
end intrinsic;


intrinsic SortClasses(GGs::SeqEnum) -> SeqEnum
{
Sorts GGs
}
  ret := [];
  for GG in GGs do
    HH := SortSubgroups(GG);
    Append(~ret, HH);
  end for;

  Sort(~ret, OrderThenTransitiveThenLexSrtList);
  return ret;
end intrinsic;


//////////////////////////////////////////////////
// Now for labelling
//////////////////////////////////////////////////

function IsDihedral(G)
  // true if and only if G is dihedral
  if #G mod 2 eq 1 then
    return false;
  elif #G eq 4 then
    return false;  
  end if;

  D := DihedralGroup(#G div 2);
  return IsIsomorphic(G, D);
end function;


function IsV4(G)
  if #G eq 4 then
    if not IsCyclic(G) then
      return true;
    end if;
  end if;
  return false;
end function;


function Getd(G)
  return #GSet(G) div 2;
end function;


function IsW(G)
  if #G eq #W2(Getd(G)) then
    return true;
  else
    return false;
  end if;
end function;       


intrinsic MakeAbstractW2dLabel(G::GrpPerm) -> MonStgElt
{
Make the first 3 things in the W2d labelling
e.g., W6.6.t or D6.6.t
}
  d := Getd(G);
  
  ret := "";
  if IsCyclic(G) then
    ret cat:= Sprintf("C%o", #G);
  elif IsV4(G) then
    ret cat:= "V4";
  elif IsDihedral(G) then
    ret cat:= Sprintf("D%o", #G div 2);
  elif IsW(G) then
    ret cat:= Sprintf("W%o", 2*d);
  elif IsAlternating(G) then
    _,_,r := GuessAltsymDegree(G);
    ret cat:= Sprintf("A%o", r);
  elif IsSymmetric(G) then
    _,_,r := GuessAltsymDegree(G);
    ret cat:= Sprintf("S%o", r);
  else
    n, q := TNums(G);
    ret cat:= Sprintf("%oT%o", q, n);
  end if;

  ret cat:= Sprintf(".%o", 2*d);

  if IsTransitive(G) then
    ret cat:= ".t";
  else
    ret cat:= ".nt";
  end if;
  
  return ret;  
end intrinsic;


intrinsic MakeW2dLabelCCs(GGs::SeqEnum[SeqEnum]) -> SeqEnum, SeqEnum[MonStgElt]
{
Given a sequence of conjugacy classes of subgroups of W_(2d) reorder them
canonically and then label each group.
}
  HHs := SortClasses(GGs);
  labs := [];
  for HH in HHs do
    start := MakeAbstractW2dLabel(HH[1]);
    Append(~labs, start);
  end for;

  cclabs := [];
  for lab in {@ x : x in labs @} do
    n := #[x : x in labs | x eq lab];
    these_labs := [lab cat "." cat alphabet[i] : i in [1..n]];
    cclabs cat:= these_labs;
  end for;

  all_labs := [];
  for i in [1..#cclabs] do
    these_labs := [Sprintf("%o.%o", cclabs[i], j) : j in [1..#HHs[i]]];
    Append(~all_labs, these_labs);
  end for;
  
  return HHs, all_labs;
end intrinsic;


intrinsic MakeW2dLabelCCs(GGs::SeqEnum[SetEnum]) -> SeqEnum, SeqEnum[MonStgElt]
{
Given a sequence of conjugacy classes of subgroups of W_(2d) reorder them
canonically and then label each group.
}
  HHs := [[h : h in HH] : HH in GGs];
  return MakeW2dLabelCCs(HHs);
end intrinsic;


intrinsic AdmissibleSubgroupsAndLabels(d::RngIntElt) -> SeqEnum, SeqEnum[MonStgElt], GrpPerm
{
Does the obvious thing   <-- make this better lol
}
  GGs, W := AdmissibleSubgroupClasses(d);
  HHs, labs := MakeW2dLabelCCs(GGs);
  return HHs, labs, W;
end intrinsic;
