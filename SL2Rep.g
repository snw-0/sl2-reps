# For now, since we haven't installed the codes as a package, we need to change the directory to read the files. Need to change by hand.

#ChangeDirectoryCurrent("path-to-package-folder");

# For me, it is
# ChangeDirectoryCurrent("/Users/MacBookPro/Dropbox/N_Wilson/SL2\ Irrep\ code");

#Read("./RepN.g");
#Read("./RepD.g");
#Read("./RepRs_regular.g");
#Read("./RepRu.g");

# arg is a list of additional parameters. Length(arg) varies for different types.
# Form of arg:
# Type N: [[i,j]]
# Type D: [[i,j]]
# Type R: [sigma, r, t, [i,j]]
# Type Ru: [[i,j]]

Rep := function(type, p, ld, arg)
	if type = "N" then
		if Length(arg) = 1 then
			Print("Correct number of parameters. Computing the representation...\n");
			return RepN(p, ld, arg[1]);
		else
			Print("Wrong number of parameters... Cannot proceed!\n");
		fi;
	elif type = "D" then
		if Length(arg) = 1 then
			Print("Correct number of parameters. Computing the representation...\n");
			return RepD(p, ld, arg[1]);
		else
			Print("Wrong number of parameters... Cannot proceed!\n");
		fi;
	elif type = "R" then
		if Length(arg) = 4 then
			Print("Correct number of parameters. Computing the representation...\n");
			return RepR(p, ld, arg[1], arg[2], arg[3], arg[4]);
		else
			Print("Wrong number of parameters... Cannot proceed!\n");
		fi;
	elif type = "Ru" then
		if Length(arg) = 1 then
			Print("Correct number of parameters. Computing the representation...\n");
			return RepRUnary(p, ld, arg[1]);
		else
			Print("Wrong number of parameters... Cannot proceed!\n");
		fi;
	else
		Print("Wrong type... Cannot proceed!\n");
	fi;
end;

