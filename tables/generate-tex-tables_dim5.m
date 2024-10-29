AttachSpec("../src/spec");

GG, labs := AdmissibleSubgroupsAndLabels(5);
ord := <[0,0,0,0,0,1,1,1,1,1], [{1,2,3,4,5},{-5,-4,-3,-2,-1}]>;

//////////////////////////////////////////////////

function easy_sort(l1,l2)
  n := StringToInteger(Split(l1, ".")[5]);
  m := StringToInteger(Split(l2, ".")[5]);
  return (n - m);
end function;

function srt(a,b)
  return easy_sort(a[2], b[2]);
end function;

function bar(n)
  if n lt 0 then
    return Sprintf("\\bar{%o}", -n);
  else
    return Sprint(n);
  end if;
end function;

function repr_elt(g)
  s := Sprint(g);
  s := &cat Split(s, ",");
  s := &cat Split(s, "(");
  s := Split(s, ")");
  s := [Split(ss, " ") : ss in s];
  s := [[StringToInteger(a) : a in ss] : ss in s];
  s := [[bar(n) : n in ss] : ss in s];
  s := [&cat ss : ss in s];
  s := ["(" cat ss cat ")" : ss in s];
  return &cat s;
end function;
  

function make_table(wt, GG, labs)
  // sort conj classes
  cclasses := {@ @};
  all_subs := &cat GG;
  all_labs := &cat labs;
  
  for i in [1..#all_subs] do
    G := all_subs[i];
    lab := all_labs[i];
    cc := WeightConjugacyClass(wt[2], G);
    cc_labs := [];
    for H in cc do
      j := Index(all_subs, H);
      Append(~cc_labs, <H, all_labs[j]>);
    end for;    
    Sort(~cc_labs, srt);
    Include(~cclasses, cc_labs);
  end for;

  reprs := [x[1] : x in cclasses];
  
  // make table
  gens_table := "\\hline \n";
  table := "\\hline \n";
  for cl in reprs do
    gen_str := &cat [repr_elt(g) cat "," : g in FewGenerators(cl[1])];
    gens_table cat:= Sprintf("\\texttt{%o} & $ %o $ \\\\ \n", cl[2], gen_str[1..#gen_str-1]);
    gens_table cat:= "\\hline \n";
    table cat:= Sprintf("\\texttt{%o} & %o &&& \\\\ \n", cl[2], AngleRank(cl[1], wt[1]));
    table cat:= "\\hline \n";
  end for;

  return table, gens_table;
end function;

table, gens_table := make_table(ord, GG, labs);

SetColumns(1024);
gens_table;
