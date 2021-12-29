#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Testing functions.
#


#! @Chapter Testing
#! @ChapterTitle Methods for testing
#! By the Chinese Remainder Theorem, it suffices to test irreps of prime power level,
#! so those are the irreps handled by the functions in this section.

#! @Section Testing

#! @Arguments p,lambda
#! @Returns the group $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$ with conjugacy classes set to the format we use.
DeclareGlobalFunction( "SL2WithConjClasses" );

#! @Arguments S,T,p,lambda
#! @Returns a list representing a character of $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$.
#! @Description Converts the modular data $(S,T)$, which must have level dividing $p^\lambda$,
#! into a character of $\mathrm{SL}_2(\mathbb{Z}/p^\lambda\mathbb{Z})$, presented in a
#! form matching the conjugacy classes used in `SL2WithConjClasses`.
DeclareGlobalFunction( "SL2ChiST" );

#! @Arguments p,lambda
#! @Returns a boolean.
#! @Description Constructs and tests all non-trivial irreps of level dividing $p^\lambda$ by checking their
#! positions in `Irr(G)` (see <URL> <Link>https://www.gap-system.org/Manuals/doc/ref/chap71.html#X873B3CC57E9A5492</Link><LinkText>Section 71.8-2 of the GAP Manual</LinkText></URL>).
#! Note that this function will print information on the irreps involved if `InfoSL2Reps` is set to level 1 or higher; see Section <Ref Sect="Chapter_Introduction_Section_Usage"/>.
DeclareGlobalFunction( "SL2TestPositions" );

#! @Arguments p,lambda
#! @Returns a boolean.
#! @Description Constructs and tests all irreps of level $p^\lambda$, confirming that the
#! $S$-matrix is symmetric and unitary and the $T$ matrix is diagonal.
#! Note that this function will print information on the irreps involved if `InfoSL2Reps` is set to level 1 or higher; see Section <Ref Sect="Chapter_Introduction_Section_Usage"/>.
DeclareGlobalFunction( "SL2TestSymmetry" );
