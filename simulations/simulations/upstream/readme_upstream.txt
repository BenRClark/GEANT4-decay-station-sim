
This is the directory for analyzing upstream clover crystals.

You can ignore the Clover_Analyzer.C and Clover_Analyzer_Downstream.C files. The new analysis code is called upstream_Clover_Analyzer.C

Use upstream_Clover_Analyzer.C to analyze clover detectors that have picked up energy peaks below 123keV. The .SRM file in this directory is probably correct, but make sure that it has all the appropriate energies that you wish to fit before executing this program.

move "Clover_%d_Results.txt" files to the "results_upstream" directory once generated. Can then execute "FitEfficiency.C" to generate an efficiency curve.

I am defining upstream to mean any detector which picks up all energy peaks we are interested in. i.e. those that detect energies less than 123keV.
