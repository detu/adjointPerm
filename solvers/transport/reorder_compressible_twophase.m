function z = reorder_compressible_twophase(z, tf, q, gm, pv, fluid)
  muo = 1.0;
  muw = 1.0;

  f = fopen ('tmp.txt', 'w');
  fprintf(f, ['/fluid/Phase0/Viscosity=%f\n', ...
              '/fluid/Phase1/Viscosity=%f\n', ...
              '/transport/sparse_linear_solver/solver_library=trilinos\n', ...
              '/transport/sparse_linear_solver/solver_type=GMRES\n', ...
              '/transport/sparse_linear_solver/use_multilevel_preconditioner=false\n',...
              '/transport/max_iterations=100\n', ...
              '/transport/tolerance=1e-6\n',...
              '/transport/relaxation_factor=0.8\n'], fluid.mu(end, 1), fluid.mu(end, 2));

    fclose(f);


    z1       = z(:,1);
    z2       = z(:,2);
    b1       = fluid.B(:,1);
    b2       = fluid.B(:,2);
    sources  = full(q);%*[1,0]; % /B_source
    pressure = ones(numel(z1), 1);
%    pm       = 1;
%    pd       = 0.01;
%    visc     = [0.1 0.1];

    [i,j,v]  = find(gm);
    Bf       = (fluid.B(i,:)+fluid.B(j,:))./2;
    [r1, r2] = matlab_reorder_compressible(z1, z2, tf, sources, -1.0 * gm, Bf, pv, b1, b2,...
                                           pressure, 'tmp.txt');
%    [r1, r2] = matlab_reorder_compressible(z1, z2, tf, sources, -1.0 * gm, fluid.B, pv, 'tmp.txt',...
%                                            pressure, pm, pd);
    z(:,1)   = r1(:);
    z(:,2)   = r2(:);

    if any(isnan(r1(:))) || any(isnan(r2(:)))
        error('NaN in return value from matlab_reorder_compressible');
    end
end
