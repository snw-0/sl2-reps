#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type N.

#! @Chapter Methods for Constructing Representations

#! @Section Representations of Type N
#!
#! TODO: Description of type $N$.

#! @Arguments p,lambda
#! @Returns a record describing $M$
#! @Description
#!  Constructs information about the underlying module $M$ of type $N$.
DeclareGlobalFunction( "SL2Reps_ModuleN" );

#! @Arguments p,lambda,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of type $N$ with level $p^\lambda$
#!  corresponding the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepN" );
