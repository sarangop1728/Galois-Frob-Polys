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

Pol<T> := PolynomialRing(Rationals());
// Copied from "examples-in-tables.sage"
paper_labels := ["2.2.ac_d", "2.2.ad_f", "2.3.ad_f", "2.2.ab_a", "2.2.ac_c", "2.4.ac_e", "2.2.a_ae", "2.3.ad_i", "2.2.a_d", "2.2.ad_g", "2.4.ah_u", "2.2.ac_e", "2.2.a_a", "2.4.ag_q", "2.4.a_ai", "3.2.ad_f_ah", "3.2.ad_g_aj", "3.2.a_a_ad", "3.2.ac_a_d", "3.2.ae_j_ap", "3.7.ak_bw_afv", "3.2.ab_ab_c", "3.4.ac_ab_g", "3.2.ac_b_a", "3.2.ab_a_a", "3.4.ab_c_a", "3.4.ab_a_ae", "3.2.ac_c_ac", "3.3.ad_j_ap", "3.2.a_a_ac", "3.7.a_a_abj","3.3.a_a_aj"];
paper_polynomials := [T^4 - 2*T^3 + 3*T^2 - 4*T + 4, T^4 - 3*T^3 + 5*T^2 - 6*T + 4, T^4 - 3*T^3 + 5*T^2 - 9*T + 9, T^4 - T^3 - 2*T + 4,T^4 - 2*T^3 + 2*T^2 - 4*T + 4, T^4 - 2*T^3 + 4*T^2 - 8*T + 16, T^4 - 4*T^2 + 4, T^4 - 3*T^3 + 8*T^2 - 9*T + 9, T^4 + 3*T^2 + 4, T^4 - 3*T^3 + 6*T^2 - 6*T + 4, T^4 - 7*T^3 + 20*T^2 - 28*T + 16, T^4 - 2*T^3 + 4*T^2 - 4*T + 4, T^4 + 4, T^4 - 6*T^3 + 16*T^2 - 24*T + 16, T^4 - 8*T^2 + 16, T^6 - 3*T^5 + 5*T^4 - 7*T^3 + 10*T^2 - 12*T + 8, T^6 - 3*T^5 + 6*T^4 - 9*T^3 + 12*T^2 - 12*T + 8, T^6 - 3*T^3 + 8, T^6 - 2*T^5 + 3*T^3 - 8*T + 8, T^6 - 4*T^5 + 9*T^4 - 15*T^3 + 18*T^2 - 16*T + 8, T^6 - 10*T^5 + 48*T^4 - 151*T^3 + 336*T^2 - 490*T + 343, T^6 - T^5 - T^4 + 2*T^3 - 2*T^2 - 4*T + 8, T^6 - 2*T^5 - T^4 + 6*T^3 - 4*T^2 - 32*T + 64, T^6 - 2*T^5 + T^4 + 2*T^2 - 8*T + 8, T^6 - T^5 - 4*T + 8, T^6 - T^5 + 2*T^4 + 8*T^2 - 16*T + 64, T^6 - T^5 - 4*T^3 - 16*T + 64, T^6 - 2*T^5 + 2*T^4 - 2*T^3 + 4*T^2 - 8*T + 8, T^6 - 3*T^5 + 9*T^4 - 15*T^3 + 27*T^2 - 27*T + 27, T^6 - 2*T^3 + 8, T^6 - 35*T^3 + 343, T^6 - 9*T^3 + 27];
angle_ranks := [];
for P in paper_polynomials do
    i := Index(paper_polynomials,P);
    delta_P := slowAngleRank(P);
    Append(~angle_ranks,delta_P);
    print i, paper_labels[i], delta_P;
end for;
