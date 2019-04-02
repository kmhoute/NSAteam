# NSAteam Directory
** Compiled/Compressed Codes: **
  Signal and Beampattern Formation
    ULA_Analysis is a function; includes FullArrayAnalysis and -WithFilter 
    NSA_Analysis is a function; includes Nested -ArrayAnalysis and -WithFilter
                    --> can possibly incorporate ProdMinMUSIC
                    --> does/can it use real and simulated data?
    BP_Formation is a CLASS; contains functions for Nested and Linear
                  ***to call the class functions: BP_Formation.Nested(#,#,#,#) or BP_Formation.Linear(#,#,#,#)
                    --> includes ProductMinBeampattern and BeampatternLinear array

  Physical Array
    Rearrange_DataFiles is a function; combines Rearrange_Data and rearrange_all
                                       --> Rearrange_Data is still necessary with current code setup 
                                       --> rearrange all can be deleted
    ifourierTrans  (idk what this is)
    -* Rearrange_Data  --> which one do we actually use?

  Comparison Tests
    NestedProcessorComp (real); uses majority of NSA_Analysis(data=real, filter=yes, processor=both) 
    NestedBeampatternComp (simulation) -- *** not fully working; something wrong with peak finding function ***
    psl_check analyzes psl through math derivation of the 3/MN pt; 
                --> need to abjust -13dB to -13.6dB and check/adjust function usage
    ScatterTest (real) *** not really sure what this did ****
  
  Estimates&Analysis (4/01/19)
    RUNdirEst.m is a CLASS (04/01/19); combines directionEstimatesVersion2.m and RUNdirectionEstimates.m
                            --> includes ways to vary SNR, snapshots, and array parameters
                            --> could use work to reduce/simplify code
                              * plotting needs editing on line width and plot size
                              
  ***add from camille: testingallMandP, changingL_SLH, and merged; from Hossam: genDATA_noise
  *** can add an extras folder or make an extras class and include the random test functions



ALL CODE:
 * '--' in front means it can be eliminated
DrKay:
  -- FullArrayAnalysis (real and simulation) 
  -- NestedArrayAnalysis (real and simulation)
  -- ProductMinBeampattern (simulation)
  ProdMinMUSIC
  
RevisedFromDrKay
  -- FullArrayAnalysisWithFilter (real and simulation)
  -- NestedWithFiltering (real and simulation)
  -- NestedBeampatternFunc (sim)

Comparison
  NestedProcessorComp (real); uses majority of NSA_Analysis(data=real, filter=yes, processor=both) 
  -- ProdMinBPComp (simulation) -- just simplified function arguments; incorporated in BP_Formation
  NestedBeampatternComp (simulation) -- ***not fully working; something wrong with peak finding function***
  psl_check analyzes psl through math derivation of the 3/MN pt; 
              --> need to abjust -13dB to -13.6dB and check/adjust function usage
  ScatterTest (real) ***not really sure what this did ****

Other
  ifourierTrans
  -- rearrange_all   --> which one do we actually use?
  -* Rearrange_Data  --> which one do we actually use?
  
Compile
  ULA_Analysis is a function; includes FullArrayAnalysis and -WithFilter 
  NSA_Analysis is a function; includes Nested -ArrayAnalysis and -WithFilter
                  --> can possibly incorporate ProdMinMUSIC
                  --> does/can it use real and simulated data?
  Rearrange_DataFiles is a function; combines Rearrange_Data and rearrange_all
                                     --> Rearrange_Data is still necessary with current code setup 
                                     --> rearrange all can be deleted
  BP_Formation is a CLASS; contains functions for Nested and Linear
                ***to call the class functions: BP_Formation.Nested(#,#,#,#) or BP_Formation.Linear(#,#,#,#)
                  --> includes ProductMinBeampattern and BeampatternLinear array
