#
# SL2Reps: Constructs representations of SL2(Z).
#
# Representations of type D.

#! @Chapter Methods for constructing representations

#! @Section Representations of type D
#!
#! TODO: Description of type $D$.

#! @Arguments p,ld
#! @Returns a record describing $(M,Q)$
#! @Description
#!  Constructs information about the underlying quadratic module $(M,Q)$ of type $D$.
DeclareGlobalFunction( "SL2Reps_ModuleD" );

#! @Arguments p,ld,chi_index
#! @Returns a list of lists of the form $[S,T]$
#! @Description
#!  Constructs the irreducible representation(s) of type $D$ with level $p^\ld$
#!  corresponding to the character $\chi$ indexed by `chi_index`.
DeclareGlobalFunction( "SL2Reps_RepD" );
