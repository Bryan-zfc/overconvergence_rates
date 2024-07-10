// The functions dk and epsilon are needed to compute a basis for the Bi as in Lauder's splitting
intrinsic dk(k::RngIntElt) -> RngIntElt
{Compute d_k as in the definition of Lauder's splitting
}
    n := Floor(k/12);
    if k mod 12 eq 2 then
        return n;
    else
        return n+1;
    end if;
end intrinsic;

intrinsic epsilon(x::RngIntElt) -> RngIntElt
{Computes epsilon of x as in Lauder's splitting
}
    if x mod 4 eq 0 then
        return 0;
    else 
        return 1;
    end if;
end intrinsic;


// Create the list of tuples (i,j) such that g_{i,j} is a basis element of B_i for i \leq N

intrinsic Tuples_in_B_i(N::RngIntElt, prime::RngIntElt) -> SeqEnum
{Computes the list of tuples [i, j] where i ranges from 1 to N and j ranges from the lower bound to the upper bound
    determined as in Lauder's splitting.
}
    B := [[0, 0]];
    for i in [1..N] do
        lower_bounds := dk((i-1)*(prime-1));
        upper_bound := dk(i*(prime-1)) - 1;
        for j in [lower_bounds..upper_bound] do
            B := B cat [[i, j]];
        end for;
    end for;
    return B;
end intrinsic;

intrinsic CreateCoeffMatrix(prime::RngIntElt, B::SeqEnum[SeqEnum], Ring::Rng, N::RngIntElt) -> Mtrx
{Computes the matrix whose rows have the Fourier coefficients of g_ij for i,j in B
}
    // Define the power series ring
    R<q> := PowerSeriesRing(Rationals(), #B + 10); // The +10 is just to be extra certain

    // Compute all the required modular forms
    M := ModularForms(Gamma0(1), 12);
    Delta := R!qExpansion(Basis(M)[2], #B + 2);
    E4 := Eisenstein(4, q);
    E6 := Eisenstein(6, q);
    Epmin1 := Eisenstein(prime - 1, q);

    // Create the matrix with rows as the Fourier coefficients of g_{i,j}
    KR := KMatrixSpace(Rationals(), #B, #B);
    M := KR!0;
    for n in [1..#B] do
        i := B[n][1];
        j := B[n][2];
        e := epsilon(i * (prime - 1));
        a := (i * (prime - 1) - 12 * j - 6 * e) / 4;
        f := (Delta^j) * (E4^a) * (E6^e) * (Epmin1^(N - i)); // This is g_{i,j}
        for y in [1..#B] do
            M[n, y] := Coefficient(f, y);
        end for;
    end for;

    // Force the matrix over the specified ring
    MS := RMatrixSpace(Ring, #B, #B);
    M := MS!M;
    return M;
end intrinsic;
