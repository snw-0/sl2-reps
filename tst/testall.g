#
# SL2Reps: Constructs representations of SL2(Z).
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "SL2Reps" );

TestDirectory(DirectoriesPackageLibrary("SL2Reps", "tst"),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace")));

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
