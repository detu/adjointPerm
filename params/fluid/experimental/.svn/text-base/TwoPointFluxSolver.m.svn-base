classdef TwoPointFluxSolver < PressureSolver
    %

    properties(Access=public)
        ht  % n x 2 matrix of one-sided transmissibilities

        grid
        boundary
    end

    methods(Access=public)
        %
        function self = TwoPointFluxSolver(g, rock, b)
            self.grid     = g; % grid should be an instance of a class derived from handle!
            self.boundary = b; % likewise!
            self.ht       = self.computeHalfTransmissibilities(rock);
        end



        function [resSol, wellSol, S, rhs] = solve(self, resSol, wellSol, fluidprop)
            [S, rhs] = self.assemble(resSol, wellSol, fluidprop);
            resSol.cellPressure = S\rhs;
            resSol.facePressure = [];
            resSol.cellFlux = [];
            t = self.getTransmissibilities(fluidprop.totalMobility(resSol.sw, resSol.cellPressure));
            resSol.faceFlux     = self.computeFaceFlux(resSol, wellSol, t);
        end
    end

    methods(Access=protected)

        function [A, RHS] = assemble(self, resSol, wellSol, fluidprop) %#ok
            ncells = numel(resSol.cellPressure);
            N   = double(self.grid.faces.neighbors);
            ind = find(prod(N, 2));

            lt    = fluidprop.totalMobility(resSol.sw, resSol.cellPressure);
            %trans = self.getTransmissibilities(lt);

            t  = self.ht(ind,:) .* lt (N(ind,:));
            t2 = self.getTransmissibilities(fluidprop.totalMobility(resSol.sw, resSol.cellPressure));

            S  = sparse(reshape(N(ind,:),          [], 1), ...
                        reshape(fliplr(N(ind, :)), [], 1), ...
                        reshape(1./t,              [], 1));

            [ii, jj, ss] = find(S);
            A = sparse(ii, jj, 1 ./ ss);

            d = sum (A, 2);
            A = A - spdiags(d, 0, ncells, ncells);

            RHS = zeros(ncells,1);

            % Dirichlet
            [f,v]  = self.boundary.getDirichletBoundary();
            if any(v)
                c      = sum(N(f, :), 2);
                RHS(c) = -v.*self.grid.faces.areas(f).*t2(c);
                ind = sub2ind(size(A), c, c);
                A(ind) = A(ind)-t2(c);
                %[-v.*self.grid.faces.areas(f).*t2(c), t2(c)]
            end
            % Sett diagonalen til T(f, nonzero)

            % Neumann
            [f,v]  = self.boundary.getNeumannBoundary();
            if any(v)
                c      = sum(N(f, :), 2);
                RHS(c) = v;
                A(sub2ind(size(A), c, c)) = t2(c);
            end
        end

        function faceFlux = computeFaceFlux(self, resSol, wellSol, t) %#ok

            faceFlux = zeros([self.grid.faces.num, 1]);
            bdry     = prod(double(self.grid.faces.neighbors), 2) == 0;

            % Internal faces
            faceFlux(~bdry) = -t(~bdry).*diff(resSol.cellPressure(self.grid.faces.neighbors(~bdry,:)), 1, 2);

            % Dirichlet
            [f,v]=self.boundary.getDirichletBoundary();
            cells = sum(double(self.grid.faces.neighbors(f, :)), 2);
            faceFlux(f) = -t(f).* ( v - resSol.cellPressure(cells));

            % Neumann
            [f,v] = self.boundary.getNeumannBoundary();
            faceFlux(f) = v;

        end

        function facePressure = computeFacePressure(self, resSol, wellSol, t) %#ok
            assert(false);
        end

        function t = getTransmissibilities(self, lt)
            t = zeros([self.grid.faces.num, 1]);
            bdry     = prod(double(self.grid.faces.neighbors), 2) == 0;
            t( bdry) = sum(self.ht(bdry,:), 2);
            t(~bdry) = TwoPointFluxSolver.harmonicAverage(self.ht(~bdry,1), self.ht(~bdry, 2));

        end



        function ht = computeHalfTransmissibilities(self, rock)
            % Compute half-transmissibilities

            FCON    = 1;%0.00852701472;
            permCol = [1, 2, 3, ...
                       2, 4, 5, ...
                       3, 5, 6];

            dimProd = double(self.grid.cells.numFaces);
            cf      = double(self.grid.cellFaces);

            ht  = zeros([self.grid.faces.num, 2]);

            areas   = self.grid.faces.areas    (cf(:,1)  );
            %outNormals = g.faces.normals  (cf(:,1),:) .* cf(:, [2,2,2]);

            % assumes normal sign is +1 if normal points from
            % g.faces.neighbors(*,1) to g.faces.neighbors(*,2).
            %Neighbors  =  g.faces.neighbors(sub2ind(size(g.faces.neighbors), cf(:,1), 1 + (cf(:,2) > 0)));
            side       =  1 + (cf(:,2) < 0);

            faceCent   = self.grid.faces.centroids(cf(:,1),:);
            perm = FCON * TwoPointFluxSolver.getPerm(self.grid, rock);

            ix = 0;
            for k = 1 : self.grid.cells.num
                nF  = dimProd(k);
                pR  = (ix  + 1 : ix  + nF)';

                a = areas(pR);
                K = reshape(perm(k, permCol), [3, 3]);
                C = faceCent(pR,:) - self.grid.cells.centroids(k(ones([nF,1])), :);
                %N = outNormals(pR,:);

                ck = C*K;
                t  = sqrt(sum(ck.*ck, 2)).*a;

                ht(sub2ind(size(ht) , double(self.grid.cellFaces(pR,1)), side(pR))) = t;

                ix = ix + nF;
            end
        end

    end

    % Helper functions
    methods(Static=true)
               
        function t = harmonicAverage(t1, t2)
           t = 2.*t1.*t2./(t1+t2);
        end
        
        function perm = getPerm(grid, rock)
            if isfield(grid.cells, 'indexMap')
                im = grid.cells.indexMap;
            else
                im = (1 : grid.cells.num) .';
            end

            if isempty(rock),
                error(id('Rock:Empty'), ...
                    'Empty ''rock'' structure is not supported.')
            elseif size(rock.perm, 2) == 1,                  % isotropic
                perm = rock.perm(im)*[1, 0, 0, 1, 0, 1];
            elseif size(rock.perm, 2) == 3,                  % diagonal
                perm = [rock.perm(im,1), zeros([numel(im), 2]), ...
                    rock.perm(im,2), zeros([numel(im), 1]), ...
                    rock.perm(im,3)];
            else                                             % full
                perm = rock.perm(im,:);
            end
        end
    end

end
