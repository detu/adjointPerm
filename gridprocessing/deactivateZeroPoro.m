function grdecl = deactivateZeroPoro(grdecl, varargin)
%Mark cells of zero porosity as inactive.
%
% SYNOPSIS:
%   grdecl = deactivateZeroPoro(grdecl)
%   grdecl = deactivateZeroPoro(grdecl, 'PORO', phi)
%   grdecl = deactivateZeroPoro(grdecl, rock)
%
% PARAMETERS:
%   grdecl     - Eclipse case specification as defined by function
%                readGRDECL.
%
%   'PORO'/phi - Explicit definition of cell porosity.  OPTIONAL.
%                Argument 'phi' must be a vector of length
%                PROD(grdecl.cartDims) (i.e. total number of cells in the
%                case specification), each value giving the porosity in the
%                corresponding cell.
%
%   rock       - Rock data structure assumed to contain field rock.poro
%                giving cell porosity for each external cell in the grid.
%                OPTIONAL.
%
% RETURNS:
%   grdecl - Updated Eclipse case specification.  The field grdecl.ACTNUM
%            is either created or updated to mark zero-porosity cells as
%            inactive.
%
% EXAMPLE:
%   case = [datadir, filesep, 'case.grdecl']
%   grdecl = readGRDECL(case)
%   nc = prod(grdecl.cartDims);
%   grdecl = deactivateZeroPoro(grdecl, 'PORO', ...
%                               [zeros(100,1); repmat(0.3, nc-100, 1)])
%   G = buildMatchingGrid(grdecl)
%
% NOTE:
%   Function deactivateZeroPoro must be called prior to calling function
%   buildMatchingGrid.
%
% SEE ALSO:
%   readGRDECL, buildMatchingGrid.

%{
#COPYRIGHT#
%}

% $Id: deactivateZeroPoro.m 1949 2009-03-31 10:16:08Z bska $

if ~isempty(grdecl) && isstruct(grdecl),
   if isfield(grdecl, 'ACTNUM'),
      ACTNUM = grdecl.ACTNUM;
   else
      ACTNUM = ones([prod(grdecl.cartDims), 1]);
   end
   
   PORO = parse_poro(grdecl, varargin{:});

   grdecl.ACTNUM = ACTNUM & (PORO > 0);
end


function PORO = parse_poro(grdecl, varargin)
opt = struct(varargin{:});
if isempty(opt),
   opt = struct('poros', ones([prod(grdecl.cartDims), 1]));
end

if isfield(grdecl, 'PORO'),
   % deactivateZeroPoro(grdecl)
   PORO = grdecl.PORO;
elseif isfield(opt, 'PORO'),
   % deactivateZeroPoro(grdecl, 'PORO', phi)
   PORO = opt.PORO;
elseif isfield(opt, 'poros'),
   % deactivateZeroPoro(grdecl, rock) or
   % deactivateZeroPoro(grdecl) when 'grdecl' has no 'PORO' field.
   PORO = opt.poro;
end
