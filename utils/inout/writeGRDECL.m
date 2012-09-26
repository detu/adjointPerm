function success = writeGRDECL(grdecl, filename)
%Write a GRDECL structure out to permanent file on disk.
%
% SYNOPSIS:
%   writeGRDECL(grdecl, file)
%
% PARAMETERS:
%   grdecl - A corner-point GRDECL structure which may be passed to
%            function 'processGRDECL' in order to construct a grid
%            structure.
%
%   file   - Name of output file.  String.  This name is passed directly to
%            function 'fopen' using mode 'wt'.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   fopen, readGRDECL, processGRDECL.

%{
#COPYRIGHT#
%}

% $Date: 2009-09-07 12:36:50 +0200 (ma, 07 sep 2009) $
% %Revision:$

   success = true;

   format = @(s,n) [repmat([s, ' '], [1, n]), '\n'];

   [fid, msg] = fopen(filename, 'wt');
   if fid < 0, error(msg); end

   fprintf(fid, 'SPECGRID\n%d %d %d 1 F\n/\n\n', grdecl.cartDims);
   fprintf(fid, 'COORD\n');
   fprintf(fid, format('%12f', 6), grdecl.COORD);
   fprintf(fid, '/\n\n');

   fprintf(fid, 'ZCORN\n');
   fprintf(fid, format('%12f', 8), grdecl.ZCORN);
   fprintf(fid, '/\n\n');

   fprintf(fid, 'ACTNUM\n');
   if isfield(grdecl, 'ACTNUM'),
      fprintf(fid, format('%d', 20), grdecl.ACTNUM);
   else
      % No field 'ACTNUM', assume all cells active.
      fprintf(fid, format('%d', 20), ones([prod(grdecl.cartDims(:)), 1]));
   end
   fprintf(fid, '/\n');

   fclose(fid);
end
