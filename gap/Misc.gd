#
# SL2Reps: Constructing symmetric representations of SL(2,Z).
#
# Miscellaneous and helper functions.
#


# InfoClass
DeclareInfoClass( "InfoSL2Reps" );
SetInfoLevel(InfoSL2Reps, 0);

# Private helper functions
DeclareGlobalFunction( "_SL2SqrtOfRootOfUnity" );
DeclareGlobalFunction( "_SL2QuadNonRes" );
DeclareGlobalFunction( "_SL2Factorizations" );
DeclareGlobalFunction( "_SL2RecordIrrep" );
DeclareGlobalFunction( "_SL2ConcatNames");
DeclareGlobalFunction( "_SL2ConjClassesOdd" );
DeclareGlobalFunction( "_SL2ConjClassesEven" );
DeclareGlobalFunction( "_SL2ConjClasses" );
