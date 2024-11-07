AttachSpec("../src/spec");

Pol<T> := PolynomialRing(Rationals());

// Each list corresponds to a table in the paper, and each entry is the Weil
// polynomial corresponding to the LMFDB label displayed there.
dim2_A_s := [T^4 - 2*T^3 + 3*T^2 - 4*T + 4, T^4 - 3*T^3 + 5*T^2 - 6*T + 4, T^4 - 3*T^3 + 5*T^2 - 9*T + 9];
dim2_B_s := [T^4 - T^3 - 2*T + 4];
dim2_C_s := [T^4 - 2*T^3 + 2*T^2 - 4*T + 4, T^4 - 2*T^3 + 4*T^2 - 8*T + 16, T^4 - 4*T^2 + 4];
dim2_A_ns := [T^4 - 3*T^3 + 8*T^2 - 9*T + 9, T^4 + 3*T^2 + 4];
dim2_B_ns := [T^4 - 3*T^3 + 6*T^2 - 6*T + 4, T^4 - 7*T^3 + 20*T^2 - 28*T + 16];
dim2_C_ns := [T^4 - 2*T^3 + 4*T^2 - 4*T + 4, T^4 + 4, T^4 - 6*T^3 + 16*T^2 - 24*T + 16, T^4 - 8*T^2 + 16];
dim3_A := [T^6 - 3*T^5 + 5*T^4 - 7*T^3 + 10*T^2 - 12*T + 8, T^6 - 3*T^5 + 6*T^4 - 9*T^3 + 12*T^2 - 12*T + 8, T^6 - 3*T^3 + 8, T^6 - 2*T^5 + 3*T^3 - 8*T + 8, T^6 - 4*T^5 + 9*T^4 - 15*T^3 + 18*T^2 - 16*T + 8, T^6 - 10*T^5 + 48*T^4 - 151*T^3 + 336*T^2 - 490*T + 343];
dim3_B := [T^6 - T^5 - T^4 + 2*T^3 - 2*T^2 - 4*T + 8, T^6 - 2*T^5 - T^4 + 6*T^3 - 4*T^2 - 32*T + 64, T^6 - 2*T^5 + T^4 + 2*T^2 - 8*T + 8];
dim3_C := [T^6 - T^5 - 4*T + 8, T^6 - T^5 + 2*T^4 + 8*T^2 - 16*T + 64, T^6 - T^5 - 4*T^3 - 16*T + 64];
dim3_D := [T^6 - 2*T^5 + 2*T^4 - 2*T^3 + 4*T^2 - 8*T + 8, T^6 - 3*T^5 + 9*T^4 - 15*T^3 + 27*T^2 - 27*T + 27, T^6 - 2*T^3 + 8, T^6 - 35*T^3 + 343];
dim3_E := [T^6 - 9*T^3 + 27];

paper_tables_dim2 := [dim2_A_s, dim2_B_s, dim2_C_s, dim2_A_ns, dim2_B_ns, dim2_C_ns];
paper_tables_dim3 := [dim3_A, dim3_B, dim3_C, dim3_D, dim3_E];

function realFactors(P)
  _, q := IsPower(Integers()!Coefficient(P,0), Degree(P) div 2);
  Pol<T> := PolynomialRing(Rationals());
  square, sqrtq := IsSquare(q);
  if square then
    real_q_polys := {T^2-q, T-sqrtq, T+sqrtq};
  else
    real_q_polys := {T^2-q};
  end if;
  return [f: f in Factorization(P) | f[1] in real_q_polys];
end function;

// Input a q-Weil polynomial and output the angle rank
slowAngleRank := function(poly)
    P := &*[f[1]: f in Factorization(poly)];
	L := SplittingField(P);
	OL := RingOfIntegers(L);
	R := Roots(P, L);
	idealFactorizations := [];
	setOfPrimes := [];
	for root in R do
		I := ideal<OL | root[1]>;
		factors := Factorization(I);
		for factor in factors do
			if not factor[1] in setOfPrimes then
				Append(~setOfPrimes, factor[1]);
			end if;
		end for;
		Append(~idealFactorizations, factors);
	end for;
	// degree = number of rows
	// set of primes = number of columns
	// no repeated roots in poly
	M := ZeroMatrix(Integers(), Degree(P), #setOfPrimes);
	for i in [1 .. Degree(P)] do
		for j in [1 .. #setOfPrimes] do
			factorization := idealFactorizations[i];
			prime := setOfPrimes[j];
			for k in [1 .. #factorization] do
				if factorization[k][1] eq prime then
					M[i,j] := factorization[k][2];
				end if;
			end for;
		end for;
	end for;
	return Rank(M) - 1;
end function;

// We verify first the claims of the tables in dimension two.
for table in paper_tables_dim2 do
  Sprintf("\nTable A.%o\n", Index(paper_tables_dim2,table));
  Sprint("w-conjugacy class | angle rank\n");
  for poly in table do
    if #realFactors(poly) eq 0 then         
      rv_list, G, _, _ := NewtonGaloisGroup(poly);
      NP := [rv[2]: rv in rv_list];
      Sprintf("%-17o | %o", MakeAbstractW2dLabel(G), AngleRank(G,NP));
    else
      P := &*[f[1] : f in Factorization(P)];
      G, _, Gdata := GaloisGroup(P);
      error if GaloisProof(P,S) eq false, "Incorrect Galois group";
      Sprintf("%-17o | %o", GroupName(G), slowAngleRank(P));
    end if;
  end for;                
end for;            
