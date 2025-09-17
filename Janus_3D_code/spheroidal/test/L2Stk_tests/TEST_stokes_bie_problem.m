classdef TEST_stokes_bie_problem < matlab.unittest.TestCase
    properties
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            addpath(genpath('../../../spheroidal'));
            addpath(genpath('../../../support'))
        end
    end

    methods (Test)
        function testInteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 10;
            target_distances = 0.2; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorDirichletProblemTwoSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 2;
            u0 = [5 5];
            target_distances = 0.2*ones(1, ns);
            plt = true;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorDirichletProblemThreeSpheroids(testCase)
            p = 16;
            eta = 30;
            ns = 3;
            u0 = [1.5 2 1.7];
            target_distances = 0.2*ones(1, ns);
            plt = true;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 2/sqrt(3);
            target_distances = 1; % Distance from the surface to evaluate the potential
            plt = true;
            neumann = false;
            interior = false;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemMultipleSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 2;
            u0 = [21/sqrt(41) 1.8];
            target_distances = 1*ones(1, ns);
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
            u0 = 2/sqrt(3);
            target_distances = 1e-2; % Distance from the surface to evaluate the potential
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
            % u0 = 11/sqrt(21); % AR = 1.1
            u0 = 2/sqrt(3); % AR = 2
            % u0 = 4/sqrt(15); % AR = 4
            target_distances = 0.5; % Distance from the surface to evaluate the potential
            plt = true;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorNeumannProblemTwoSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 15;
            ns = 2;
            u0 = [5 5];
            target_distances = 1*ones(1, ns);
            plt = false;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end
    end
end