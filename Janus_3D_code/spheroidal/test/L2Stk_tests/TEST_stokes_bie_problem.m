classdef TEST_stokes_bie_problem < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
        end
    end

    methods (Test)
        function testInteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 3/sqrt(2);
            target_distances = 0.3; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 3/sqrt(2);
            target_distances = 3; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = false;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorNeumannProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 3/sqrt(2);
            target_distances = 0.4; % Distance from the surface to evaluate the potential
            plt = true;
            neumann = true;
            interior = true;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorNeumannProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 3/sqrt(2);
            target_distances = 1; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end
    end
end