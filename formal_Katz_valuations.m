Attach(aux_functions.m)

\\ The valuations of b_(i,j) appearing in the formal Katz expansion foe the given prime, for all (i,j) where 0 ≤ j ≤ i ≤ r. The input lambda is the number of weights (increase to get sharper results, see Chapter 2 in "COMPUTATIONS ON OVERCONVERGENCE RATES RELATED TO THE EISENSTEIN FAMILY").


Valuations_Formal_Katz := function( prime,lambda, r)
   
  \\Create the quotient ring of the $p$-adic integers, with accuracy lambda + 1
	Quot_ring := quo<Integers() | prime^(lambda+1)>;

	\\Create the list of weights in the classical case and in the overconvergent case
	List_of_weights_classical := [ n * (prime-1) : n in [1..lambda] ];
	List_of_padic_weights:= [ Quot_ring!( (prime+1)^w-1 ) : w in List_of_weights_classical];
	
	\\Create the Vandermonde Matrix over Quot_ring and with the list of p-adic weights as input
	Vander_matrix:= Vandermonde(Quot_ring, List_of_padic_weights, lambda); 

	\\Compute the valuations appearing in the kernel of the Vandermonder Matrix. This is needed as our matrix is generally singular and thus our system of equations will not have a 	unique solution (but we know it has at least one solution). 

	Sequence_Basis := ElementToSequence( Basis( Kernel( Transpose(Vander_matrix) ) ) );
	Valuations_Kernel := [lambda + 1 : l in [1..lambda]];

	for i in [1..lambda] do
		for j in [1..#Sequence_Basis] do
			Valuations_Kernel[i] := Minimum(Valuations_Kernel[i], Valuation(Integers()!Sequence_Basis[j][i],prime));
		end for;
	end for;


	\\Compute the coefficients in the Katz Expansions of the Eisenstein series (as in step 2 of algorithm 2 in the paper) for all the weights
	Coefficients_Katz_Expansions := [ Katz_exp_coef := function(prime, lambda+1,  Create_Eisenstein_ratio(prime, b, weight) , r) : weigth in List_of_weights_classical];

	\\Solve the linear equations and compute the valuations (up to an element in the kernel). It will return "-1" as valuation if it cannot be determined (in which case one can increase lambda, i.e. the number of weights).
	Answer := [ [ r, j-1, -1]: j in [1..lambda] ];
  
	for j in [1..r] do
		Coeff_vector := [ Coefficient( Coefficients_Katz_Expansions[i],j)   : i in [1..lambda]];
		Transposed_coeff_vector_over_Quot_ring := Transpose( RMatrixSpace(Quot_ring , 1, lambda) ! Vector(Coeff_vector) );
		Solution_vector:= Solution( Transpose(Vander_matrix), Vector(Transposed_coeff_vector_over_Quot_ring));

		for k in [1..lambda] do
			if Valuation(Integers() ! Solution_vector[k],prime) ne Infinity() 
			   and Valuation(Integers() ! Solution_vector[k],prime) lt Valuations_Kernel[k]
			   and Valuation(Integers()! Solution_vector[k],prime) lt lambda + 1 then
				val_up_to_kernel := Valuation(Integers()! Solution_vector[k] , prime);
				Answer[k][3] := Rationals!val_up_to_kernel;
			end if;
		end for;

	return Answer;
end function;
