# Sqrt of roots of unity.

SR := function(x)
	local l, i;
	l := Order(x);
	for i in [1..l] do
		if E(l)^i = x then
			return E(2*l)^i;
		break;
		fi;
	od;
end;