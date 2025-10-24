//Create Vandemonde Matrix over the ring R, with List as inputs, and of size #List x n

intrinsic Vandermonde( R::Rng, List::SeqEnum , n::RngIntElt) -> Mtrx
  //Create the correct matrix space and a zero matrix that we will change into the Vandermonde matrix	
  
  Mat_space:= RMatrixSpace(R,#List,n);
  M := Mat_space!0	 	
  for i := 1 to #List do
    for j := 1 to n do
      M[i,j] := List[i]^(j-1);
    end for;
  end for; 
  return M;
end intrinsic;


//Create the functions d_k and epsilon (as in the paper).

intrinsic d_k( k::RngIntElt ) -> RngIntElt
	n := Floor(k/12);
	if k mod 12 eq 2 then
		return n;
	else 
		return n+1;
	end if;
end intrinsic;

intrinsic epsilon( x::RngIntElt ) -> RngIntElt
	if x mod 4 eq 0 then
		return 0;
	else 
		return 1;
	end if;
end intrinsic;

//For a given prime and N, compute all the tuples (i,j) such that g_(i,j) is in B_{i}, up to j=N-1
intrinsic Create_basis_tuples( prime::RngIntElt , N::RngIntElt) -> SeqEnum

	Basis_tuples := [[0,0]];

	for i in [1..N] do
		lower_bound := d_k( (i-1) * (prime-1) );
		upper_bound := d_k( i* (prime-1) ) - 1;
		
		for j in [lower_bound..upper_bound] do
			Basis_tuples := Basis_tuples cat [ [i,j] ];
		end for;

	end for;
	return Basis_tuples;
end intrinsic;


//Create the matrix spaces and the coefficient matrix of the g_{i,j}'s which form a basis of the B, for a list B containing tuples (i,j) and over the ring R
intrinsic Create_coeff_matrix( prime::RngIntElt , B::SeqEnum[SeqEnum] , R::Rng) -> Mtrx

	Q<q>:= PowerSeriesRing(Rationals(),#B+1);
	Mod_forms_space := ModularForms(Gamma0(1),12);
	Delta := Q!qExpansion(Basis(Mod)[2], #B+1);
	E_4 := Q!Eisenstein(4,q);
	E_6 := Q!Eisenstein(6,q);
	E_pmin1 := Q!Eisenstein(prime-1,q);
	Matrix_Space_over_Q:=KMatrixSpace(Rationals(),#B,#B);
	Matrix_Space_Over_R := RMatrixSpace(R,#B,#B);
	M:=Matrix_Space_over_Q!0;
	for n in [1..#B] do
		i := B[n][1];
		j := B[n][2];
		e := epsilon(i*(prime-1));
		a := (i*(prime-1)-12*j-6*e)/4;
		g_ij :=  (Delta^j)*(E4^a)*(E6^e)*(Epmin1^(N-i)); //this is g_{i,j}
		for k in [1..#B] do
			M[n,k] := Coefficient(g_ij,k-11);
		end for;
	end for;
	M := Matrix_Space_Over_R!M;
	return M;
end intrinsic;

//Create E*/VE* for  prime, up to q^N, of classical weight
intrinsic Create_Eisenstein_ratio( prime::RngIntElt , N::RngIntElt , weight::RngIntElt RngUPolElt) -> RngUPolElt
        R<q> := PowerSeriesRing(Rationals(), N);
        E := Eisenstein(weight, q);
        E_star := E - prime^(weight-1) * Evaluate(E, q^prime);
        V_E_star := Evaluate(E_star,q^(prime));
        return E_star/V_E_star;
end intrinsic;


//Compute the coefficients (which are modular forms!) in the Katz Expansion of the modular function "f" for a given prime, up till B_n  and over the ring Z[[q]]/ (q^N, p^C) (where N is related to n), see the paper. 

intrinsic Katz_exp_coef( prime::RngIntElt , C::RngIntElt, f::RngUpolElt , n::RngIntElt)
	
	N := d_k(n * (prime -1) );
	Z_p := quo<Integers() | prime^C>;

	#Find the values (i,j) such that the g_{i,j} form a basis for Bi

	Basis_tuples := Create_basis_tuples(prime, N);
	basis_size := #Basis_tuples;

	//Create the modular functions needed 
	Q<q>:= PowerSeriesRing(Rationals(),basis_size);

	Mod_forms_space := ModularForms( Gamma0(1), 12);

	Delta := Q!qExpansion( Basis(Mod_forms_space)[2], N + 1 + basis_size);
	E_4 := Q!Eisenstein(4,q);
	E_6 := Q!Eisenstein(6,q);
	E_pmin1 := Q!Eisenstein(prime-1, q);
	f_tilde :=  E_pmin1^(N)*Q!f;

	//Create the matrix spaces and the coefficient matrix of the g_{i,j}'s which form a basis of the B_i and solve the system of equations	
	
	Coeff_matrix := Create_coeff_matrix := function(prime, Basis_tuples, Q)
	Coeff_vector_of_f := Vector([Z_p!Coefficient(f_tilde,i-1) : i in [1..basis_size] ] );
	Solution := Solution(Coeff_matrix, Coeff_vector_of_f);

	S<q>:= PowerSeriesRing(Z_p,basis_size);
	g := S!0;

	for i in [1..basis_size] do
		if Basis_tuples[i][1] eq j then 
			n := Bsis_tuples[i][2];
			e := epsilon(j*(prime-1));
			a := (j*(prime-1)-12*n-6*e)/4;
			f :=  (S!Delta^n)*(S!E4^a)*(S!E6^e); 
			g := g + Solution[i]*(S!f);
		end if;
	end for;
	return g;
end intrinsic;
