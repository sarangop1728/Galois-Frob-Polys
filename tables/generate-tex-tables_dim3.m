AttachSpec("../src/spec");

GG, labs := AdmissibleSubgroupsAndLabels(3);

A := <[0,0,0,1,1,1], [{1,2,3},{-3,-2,-1}]>;
B := <[0,0,1/2,1/2,1,1], [{1,2},{3,-3},{-2,-1}]>;
C := <[0,1/2,1/2,1/2,1/2,1], [{1},{2,3,-3,-2},{-1}]>;
D := <[1/3,1/3,1/3,2/3,2/3,2/3], [{1,2,3},{-3,-2,-1}]>;
E := <[1/2,1/2,1/2,1/2,1/2,1/2], [{1,2,3,-3,-2,-1}]>;

wts := [A,B,C,D,E];

//////////////////////////////////////////////////

function make_table(wt, GG, labs)
  // sort conj classes
  cclasses := {};
  all_subs := &cat GG;
  all_labs := &cat labs;
  
  for G in all_subs do
    Include(~cclasses, WeightConjugacyClass(wt[2], G));
  end for;

  cclasses := [[G : G in cl] : cl in cclasses];
//  cclasses := SortClasses(cclasses);
  
  cclasses_labs := [];
  for cl in cclasses do
    this_cl := [];
    for G in cl do
      i := Index(all_subs, G);
      Append(~this_cl, all_labs[i]);
    end for;
    Append(~cclasses_labs, this_cl);
  end for;

  // make table
  ret := "\\hline \n";
  for cl in cclasses_labs do
    i := Index(all_labs, cl[1]);
    ark_str := Sprintf("\\multirow{%o}{*}{%o}", #cl, AngleRank(all_subs[i], wt[1]) );
    ret cat:= Sprintf("\\wl{%o} & %o &&& \\\\ \n", cl[1], ark_str);
    for G in cl[2..#cl] do
      ret cat:= Sprintf("\\wl{%o} &&&& \\\\ \n", G);
    end for;
    ret cat:= "\\hline \n";
  end for;

  return ret;
end function;

for wt in wts do
  make_table(wt, GG, labs);
  printf "Newton polygon: %o", wt[1];
  "\n\n--------------------------------------------------";
  "--------------------------------------------------\n\n";
end for;
