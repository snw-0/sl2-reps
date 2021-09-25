#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type D.

#! @Chapter Methods for Constructing Representations

#! @Section Representations of Type D
#!
#! TODO: Description of type $D$.

#! @Arguments p,lambda
#! @Returns a record describing $M$
#! @Description
#!  Constructs information about the underlying module $M$ of type $D$.
DeclareGlobalFunction( "SL2Reps_ModuleD" );

#! @Arguments p,lambda,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of type $D$ with level $p^\lambda$
#!  corresponding the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepD" );
