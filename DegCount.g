# Input: degree. Output: the number of all the irreps of dimension <= deg.

DegTotal := function(deg)
	local count, p, ld, pmax, ldmax, pset, ldset;
	
	count := 1;

	# p odd, ld = 1.	
	pmax := 2*deg + 1;
	pset := Filtered([3..pmax], x -> IsPrime(x));
	
	for p in pset do
		if p + 1 <= deg then
			count := count + (p-3)/2;
		fi;

		if p - 1 <= deg then 
			count := count + (p-1)/2;
		fi;

		if (p+1)/2 <= deg then
			count := count + 2;
		fi;

		if (p-1)/2 <= deg then
			count := count + 2;
		fi;

		if p <= deg then
			count := count + 1;
		fi;
	od;

	# p odd, ld >= 2.
	pmax := Maximum([RootInt(2*deg+1, 2), 3]);
	ldmax := Maximum([LogInt(deg, 3) + 2, 2]);
	pset := Filtered([3..pmax], x -> IsPrime(x));
	ldset := [2..ldmax];
	
	for p in pset do
		for ld in ldset do
			if (p+1)*p^(ld-1) <= deg then
				count := count + p^(ld-2)*(p-1)^2/2;
			fi;

			if (p-1)*p^(ld-1) <= deg then
				count := count + p^(ld-2)*(p^2-1)/2;
			fi;

			if p^(ld-2)*(p^2-1)/2 <= deg then
				count := count + Sum([1..ld-1], si -> 4*(p-1)*p^(ld-si-1));
				count := count + 4;
			fi;
		od;
	od;

	# p = 2, ld = 1.
	if 1 <= deg then 
		count := count + 1;
	fi;

	if 2 <= deg then 
		count := count + 1;
	fi;

	# p = 2, ld = 2.
	if 3 <= deg then 
		count := count + 4;
	fi;

	if 1 <= deg then 
		count := count + 2;
	fi;

	if 2 <= deg then 
		count := count + 1;
	fi;

	# p = 2, ld = 3.
	
	if 6 <= deg then 
		count := count + 6;
	fi;

	if 4 <= deg then 
		count := count + 2;
	fi;

	if 2 <= deg then 
		count := count + 4;
	fi;

	if 3 <= deg then 
		count := count + 8;
	fi;

	# p = 2, ld = 4.
	if 24 <= deg then 
		count := count + 2;
	fi;

	if 8 <= deg then 
		count := count + 6;
	fi;

	if 6 <= deg then 
		count := count + 20;
	fi;

	if 3 <= deg then 
		count := count + 16;
	fi;

	if 12 <= deg then 
		count := count + 2;
	fi;

	# p = 2, ld = 5.
	if 48 <= deg then 
		count := count + 4;
	fi;

	if 16 <= deg then 
		count := count + 12;
	fi;

	if 12 <= deg then 
		count := count + 40;
	fi;

	if 24 <= deg then 
		count := count + 4;
	fi;

	if 6 <= deg then 
		count := count + 32;
	fi;

	# p = 2, ld > 5.
	ldmax := Maximum([LogInt(deg, 2) + 4, 6]);
	ldset := [6..ldmax];

	for ld in ldset do
		if 3*2^(ld-1) <= deg then
			count := count + 2^(ld-3);
		fi;
		
		if 2^(ld-1) <= deg then
			count := count + 3*2^(ld-3);
		fi;

		if 3*2^(ld-2) <= deg then
			count := count + 2^(ld-3);
		fi;

		if 3*2^(ld-3) <= deg then
			count := count + 10*2^(ld-3);
		fi;

		if 3*2^(ld-4) <= deg then
			count := count + 32 + Sum([3..ld-3], si -> 16*2^(ld-si-2));
		fi;
	od;

	return count;
end;

# Input: degree. Output: the number of all the irreps of dimension = deg.

DegExact := function(deg)
	if deg = 1 then
		return DegTotal(deg);
	else
		return DegTotal(deg) - DegTotal(deg-1);
	fi;
end;





	

































			




















