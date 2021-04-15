RepRp1:=function (p,a)
	local S,T, n, i, j;
	n:=(p+1)/2;
	S:=IdentityMat(n);
	for i in [0..n-1] do
		for j in [i..n-1] do
			if i=0 and j > 0 then
				S[i+1][j+1]:=Sqrt(2);
			elif i>0 and j >0 then
				S[i+1][j+1]:=E(p)^(2*a*i*j) + E(p)^(-2*a*i*j);
			fi;
		od;
	od;
	for i in [1..n] do
		for j in [1..i-1] do
			S[i][j]:=S[j][i];
		od;
	od;
      S:=Sqrt(Jacobi(-1,p))*Jacobi(a,p)*S/Sqrt(p);
	T:=DiagonalMat(List([0..n-1], i->E(p)^(a*i^2)));
	return [S, T, n];
end;

RepRp2:=function (p,a)
	local S,T, n, i, j;
	n:=(p-1)/2;
	S:=IdentityMat(n);
	for i in [1..n] do
		for j in [i..n] do
			S[i][j]:=E(p)^(2*a*i*j) - E(p)^(-2*a*i*j);
		od;
	od;
	for i in [1..n] do
		for j in [1..i-1] do
			S[i][j]:=S[j][i];
		od;
	od;
      S:=Sqrt(Jacobi(-1,p))*Jacobi(a,p)*S/Sqrt(p);
	T:=DiagonalMat(List([1..n], i->E(p)^(a*i^2)));
	Print(S^4=S^0, (S*T)^3=S^0);
	return [S, T, n];
end;
