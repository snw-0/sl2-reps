gap> START_TEST( "Exceptional irrep test" );
gap> TestPrint := function(rho); Print(rho.name, " : ", rho.degree, " ", rho.level, " ", Trace(rho.S), " ", Trace(rho.T), "\n"); end;;
gap> test := SL2IrrepsExceptional();; for rho in test do TestPrint(rho); od;
Xi_6 tensor R_2^0(1,3)_1 : 3 4 1 -1
Xi_9 tensor R_3^0(1,3,nu)_1 : 6 8 0 0
Xi_6 tensor R_4^2(1,3,[0,1]) : 6 16 0 -1-E(4)
Xi_6 tensor R_4^2(3,3,[0,1]) : 6 16 0 -1+E(4)
N_3([1,0])+ tensor R_4^0(1,7,[1,1])+ : 12 16 0 0
N_3([1,3])+ tensor R_4^0(1,7,[1,1])+ : 12 16 0 0
Xi_9 tensor R_5^2(1,1,[0,0])_1 : 12 32 0 0
Xi_9 tensor R_5^2(1,1,[0,1])_1 : 12 32 0 0
Xi_9 tensor R_5^2(3,1,[0,0])_1 : 12 32 0 0
Xi_9 tensor R_5^2(3,1,[0,1])_1 : 12 32 0 0
Xi_6 tensor R_6^4(1,1,nu)_1 : 12 64 0 -1-E(4)
Xi_6 tensor R_6^4(1,3,nu)_1 : 12 64 0 -1-E(4)
Xi_6 tensor R_6^4(3,1,nu)_1 : 12 64 0 -1+E(4)
Xi_6 tensor R_6^4(3,3,nu)_1 : 12 64 0 -1+E(4)
Xi_6 tensor R_6^4(5,1,nu)_1 : 12 64 0 -1-E(4)
Xi_6 tensor R_6^4(5,3,nu)_1 : 12 64 0 -1-E(4)
Xi_6 tensor R_6^4(7,1,nu)_1 : 12 64 0 -1+E(4)
Xi_6 tensor R_6^4(7,3,nu)_1 : 12 64 0 -1+E(4)
gap> STOP_TEST( "except.tst" );
