classdef CompressibleFluidBase < handle
   %Base class for compressible fluids

   % $Date: 2009-10-01 16:38:31 +0200 (to, 01 okt 2009) $
   % $Revision: 2926 $

   properties(GetAccess=public, SetAccess=protected, Abstract=true)
      surfaceDensity
      info
      viscosity
      miscible
   end

   methods (Abstract)
      relperm(self, s)
      pvt    (self, p, z)
   end
end
