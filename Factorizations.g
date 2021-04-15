_FactorizationsStep := function(factors)
	local c, m, o, rest, r;

	if Length(factors) <= 1 then
		return [factors];
	fi;

	o := [];

	for c in Combinations([1 .. Length(factors)]) do
		if Length(c) = 0 then
			continue;
		fi;

		m := Product(factors{c});

		rest := [1 .. Length(factors)];
		SubtractSet(rest, c);
		rest := factors{rest};

		for r in _FactorizationsStep(rest) do
			Add(r, m);
			Add(o, r);
		od;
	od;

	return o;
end;

Factorizations := function(n)
	local o, f;

	o := [];

	for f in _FactorizationsStep(Factors(n)) do
		Sort(f);
		f := Reversed(f);
		Add(f, 1);
		if not f in o then
			Add(o,f);
		fi;
	od;

	StableSort(o);
	return Reversed(o);
end;
