AttachSpec("src/spec");
import "src/weighted-perm-rep.m": permMatrix;

// Returns Weil polynomial from first g point counts on curve.
function WeilPolyFromPointCount(points, q)
    g := #points;
    R<T> := PolynomialRing(Rationals());
    Q<t> := PowerSeriesRing(Rationals());
    u := &+[points[j]*t^j/j : j in [1..g]];
    L := Exp(u)*(1-t)*(1-q*t);
    L_coefficients := [Coefficients(L)[i] : i in [1..g+1]];
    for i in [1..g] do
        Append(~L_coefficients,q^i*L_coefficients[g-i+1]);
    end for;
    P_coefficients := [L_coefficients[2*g-i+1] : i in [0..2*g]];
    return R!P_coefficients;
end function;

// .............................................................
// .............................................................
// Example: Jacobian of the hyperelliptic curve with affine
// equation y^2 = x^9-1 over finite field with 19 elements.
// .............................................................
// .............................................................

q := 19;
_<x> := PolynomialRing(GF(q));
C := HyperellipticCurve(x^9-1);
F := AlgorithmicFunctionField(FunctionField(C));
g := Genus(F);
Sprintf("\nC: %o, of genus %o.", C, g);

points := [NumberOfPlacesOfDegreeOneECF(C, i) : i in [1..4]];
Sprintf("The first %o point counts of C are: %o.\n",g,points);

R<T> := PolynomialRing(Rationals());
P := R!WeilPolyFromPointCount(points, q);
Sprintf("From the point counts, we calculate that the Frobenius polynomial of the Jacobian is given by \n P(T) = %o\n",P);

P1, P2 := Explode(Factorization(P));
P1 := P1[1];
P2 := P2[1];
Sprintf("P(T) factors are P1(T)*P2(T), where \n P1(T) = %o \n P2(T) = %o\n",P1,P2);

roots_NP, H, cc, L := NewtonGaloisGroup(P);
Sprintf("The polynomial P(T) has Galois group %o. \n And permutation representation with image: %o.\n", GroupName(H),H);

RootsA := [rnp[1]: rnp in roots_NP];
NP := [rnp[2]: rnp in roots_NP];
Sprintf("The Newton polygon of P(T) has slopes: %o.\n",NP);

X := X2(4);
RootsE := [r: r in RootsA | Evaluate(P1,r) eq 0 ];
IndE := [Index(RootsA,r): r in RootsE];
RootsB := [r: r in RootsA | Evaluate(P2,r) eq 0 ];
IndB := [Index(RootsA,r): r in RootsB];
a,b,c, abar, bbar, cbar := Explode(RootsA[IndB]);
pi, pibar := Explode(RootsA[IndE]);
Sprintf("The eigenvalues a,b,c, abar, bbar, cbar of B correspond to the indices %o. \n", [X[i]: i in IndB]);
Sprintf("The eigenvalues pi, pibar of E correspond to the indices %o. \n", [X[i]: i in IndE]);
Sprintf("To calculate the angle rank, we consider the matrix: \n%o, \nwhose entries come from the indexing of the roots and the action of the group H.\n", permMatrix(g,H));
Sprintf("The Newton hyperplane matrix is obtained from this one by replacing each index i with alpha_i/alphabar_i, and taking the valuations of the entries. We obtain:\n%o.\n",NewtonMatrix(H,NP));
Sprintf("We have the multiplicative relation a*b/(pi*c) = %o. This was found by calculating the kernel of the matrix above.\n", a*b/(pi*c));