% Low-level input routines for ECLIPSE deck.
%
% Files
%   readGRDECL      - Read Eclipse file assuming corner-point grid format.
%   readPARAMS      - Read Eclipse or FrontSim data file.
%   readpvt         - Read Black Oil PVT data from ECLIPSE parameter file into table.
%   readVector      - Input vector of floating point numbers from GRDECL file.
%   cacheDir        - Define platform-dependent directory for temporary files.
%   readCache       - Read cached variables associated with given file into caller's workspace.
%   readEclipseFile - Input ECLIPSE or FrontSim Deck specification
%   readWellKW      - Read well definitions from an ECLIPSE Deck specification.
%   writeCache      - Write callers workspace to matfile in directory ./.cache/
%   writeGRDECL     - Write a GRDECL structure out to permanent file on disk.

%{
#COPYRIGHT#
%}
