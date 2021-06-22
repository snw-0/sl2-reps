# Orbit of M/A. Note that A does not act freely.

orbx := function(x)
	local a, ans;
	ans := [];
	for a in A do
		Add(ans, Prod(x, a));
	od;

	return ans;
end;

O := AsSet(List(M, x -> AsSet(orbx(x))));

for x in O do
	Print(x, "\n");
od;
