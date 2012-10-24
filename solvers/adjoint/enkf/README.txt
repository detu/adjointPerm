

General:

The EnKF package contains routines to perform history matching experiments using the MRST model construction 
and simulation functionality. The included algorithms are sequential or iterated global and local EnKF and 
EnRML. Measurements can be well rates, water cut, bottom hole pressure and grid cell saturation. Multiple 
measurement times can be specified per update time. The models can be simulated forward after the final update.
The static parameters permeability and porosity, and dynamic variables pressure and saturation can be updated.
Also generic structural model parameters can be specified and updated (see example 2). Two routines are included
for visualization of the results. The possible input options are further explained in the example inputSettings 
files. The code has been tested with Matlab version R2011a.


Instructions:

1. Place the folder enkf-2012a somewhere in the MRST tree.
2. The main routine is mrstEnKF.m.
3. The input for EnKF experiments should be specified in inputSettings.m.
4. To run one of the examples, rename inputSettings#.m to inputSettings.m, where # is the example number.
5. For example 3, the file NORNE.GRDECL should be placed in the folder specified in testGrid3.m. Example 1 uses 
the file with realizations specified in testGrid1.m.
6. Results can be visualized with the routines plotGrid.m and plotProduction.m.
