classdef TEST_stokes_bie_problem < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;

        u0_prolate = 4/sqrt(3);
        u0_oblate = 8/sqrt(3);
        a_prolate;
        a_oblate;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);
        end
    end

    methods (Test)
        function testInteriorDirichletProblemOneSpheroid(testCase)
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 13/sqrt(69);
            target_distances = 0.4; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, truesoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorDirichletProblemMultipleSpheroids(testCase)
        end

        function testExteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 13/sqrt(69);
            target_distances = 3; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = false;
            
            [soln, truesoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemMultipleSpheroid(testCase)
        end
    end
end