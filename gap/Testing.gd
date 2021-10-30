#
# SL2Reps: Constructs representations of SL2(Z).
#
# Testing functions.

#! @Chapter Testing
#! @ChapterTitle Methods for testing

#! @Section Testing

#! @Arguments p,ld
#! @Returns the group $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$ with conjugacy classes set to the format we use.
DeclareGlobalFunction( "SL2Reps_SL2Conj" );

#! @Arguments S,T,p,ld
#! @Returns a list representing a character of $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$
#! @Description Converts the modular data $(S,T)$, which must have level dividing $p^\lambda$,
#! into a character of $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$, presented in a
#! form matching the conjugacy classes used in `SL2Reps_SL2Conj`.
DeclareGlobalFunction( "SL2Reps_ChiST" );

#! @Arguments p, lambda
#! @Returns a boolean
#! @Description Constructs and tests all irreps of level dividing $p^\lambda$ by checking their
#! positions in `Irr(G)`.
DeclareGlobalFunction( "SL2Reps_IrrepPositionTest" );
