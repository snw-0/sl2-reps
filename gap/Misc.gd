#
# SL2Reps: Constructs representations of SL2(Z).
#

DeclareInfoClass( "InfoSL2Reps" );

# Private helper functions
DeclareGlobalFunction( "_SL2SqrtOfRootOfUnity" );
DeclareGlobalFunction( "_SL2QuadNonRes" );
DeclareGlobalFunction( "_SL2Factorizations" );
DeclareGlobalFunction( "_SL2RecordIrrep" );
DeclareGlobalFunction( "_SL2ConcatNames");
DeclareGlobalFunction( "_SL2ConjClassesOdd" );
DeclareGlobalFunction( "_SL2ConjClassesEven" );
DeclareGlobalFunction( "_SL2ConjClasses" );

# DeclareGlobalFunction( "_SL2Reps_CharNorm" );
# DeclareGlobalFunction( "_SL2Reps_RepChi" );
# DeclareGlobalFunction( "_SL2Reps_ClassMap" );
# DeclareGlobalFunction( "_SL2Reps_ChiTest" );
