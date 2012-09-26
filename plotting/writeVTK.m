function writeVTK(vertex, face, colour, filename)
%Write PATCH (vertex,face,colour) data to file on disk in VTK format.
%
% SYNOPSIS:
%   writeVTK(vertex, face, colour, filename)
%
% PARAMETERS:
%   vertex - Unique point list.
%   face   - Vector of face indices to vertex
%   colour - rgb colour or scalar field values.
%
%   filename -
%            Name of output file.  String.  This name is passed directly to
%            function 'fopen' using mode 'wt'.
%
% RETURNS:
%   Nothing.
%
% SEE ALSO:
%   fopen, plotFaces.

%{
#COPYRIGHT#
%}

% $Id: writeVTK.m 2023 2009-04-17 18:04:21Z bska $

   [fid, msg] = fopen(filename, 'wt');
   if fid < 0, error(msg); end

   % Write heading.
   num_verts   = size(vertex, 1);
   num_faces   = size(face  , 1);
   num_colours = size(colour, 1);

   fprintf(fid, '# vtk DataFile Version 5.0\n');
   fprintf(fid, 'Faces\n');
   fprintf(fid, 'ASCII\n');
   fprintf(fid, 'DATASET POLYDATA\n\n');

   assert (any(size(vertex,2) == [2, 3]));
   fprintf(fid, 'POINTS %d float\n', num_verts);
   if size(vertex,2) == 2,
      vertex = [vertex, zeros([num_verts,1])];
   end

   % Write vertices (inverted Z coordinate to make output look like M).
   vertex(:,3) = -vertex(:,3);
   fprintf(fid, '%16.4f %16.4f %16.4f\n', vertex .');
   fprintf(fid, '\n');

   % Write faces.
   tot_np = length(find(~isnan(face)));
   fprintf(fid, 'POLYGONS %d %d\n', num_faces, num_faces + tot_np);
   for f = 1 : num_faces,
      row = find(~isnan(face(f,:)));
      nrp = length(row);
      fprintf(fid, '%d ', nrp);
      for p = 1 : nrp,
         fprintf(fid, ' %d', face(f,row(p))-1);
      end
      fprintf(fid, '\n');
   end
   fprintf(fid, '\n');

   % Write dataset attributes
   if num_colours == 1, colour = colour(ones([num_faces, 1]), :); end

   fprintf(fid, 'CELL_DATA %d\n', num_faces);

   if size(colour,2) == 1,
      % Write scalars.
      fprintf(fid, 'SCALARS scalars float\n');
      fprintf(fid, 'LOOKUP_TABLE default\n');
      fprintf(fid, '%8.4f\n', colour);
   else
      % Write colours.
      assert (size(colour,2) == 3);
      fprintf(fid, 'COLOR_SCALARS rgb 3\n');
      fprintf(fid, '%8.4f %8.4 %8.4f\n', colour.');
   end

   fclose(fid);
end
