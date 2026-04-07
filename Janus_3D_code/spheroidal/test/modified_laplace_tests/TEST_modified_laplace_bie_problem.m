classdef TEST_modified_laplace_bie_problem < matlab.unittest.TestCase
    properties
        p = 16
        lambda = 1
    end

    methods (TestClassSetup)
        function setup(testCase)
            test_dir = fileparts(mfilename('fullpath'));
            repo_root = fullfile(test_dir, '..', '..', '..');
            addpath(genpath(fullfile(repo_root, 'spheroidal')));
            addpath(genpath(fullfile(repo_root, 'support')));
        end
    end

    methods (Test)
        function testBIEConvergenceWithOrder(testCase)
            p_list = [8, 12, 16];
            u0 = 2/sqrt(3);
            fit_err = zeros(size(p_list));
            sol_err = zeros(size(p_list));
            condS = zeros(size(p_list));

            for i = 1:numel(p_list)
                [~, ~, ~, condK, info] = modified_laplace_bie_problem( ...
                    p_list(i), testCase.lambda, u0, 5e-1, false, "exterior_dirichlet");
                sol_err(i) = info.rel_soln_err;
                condS(i) = condK;
            end

            testCase.verifyTrue(all(sol_err(2:end) < sol_err(1:end-1)), ...
                'Solution error should improve as p increases.');
            testCase.verifyLessThan(sol_err(end), 1e-5, ...
                'Highest-order solution error is too large.');
        end

        function testBIESLPTargetReconstructionAtMultipleDistances(testCase)
            distances = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 1e-2, 1e-1, 1];
            u0 = 2/sqrt(3);
            rel_err = zeros(size(distances));

            for i = 1:numel(distances)
                [~, ~, ~, ~, info] = modified_laplace_bie_problem( ...
                    testCase.p, testCase.lambda, u0, distances(i), false, "exterior_dirichlet");
                rel_err(i) = info.rel_soln_err;
            end

            testCase.verifyLessThan(rel_err, 8e-5, ...
                'Near-singular target evaluation is not sufficiently accurate for some target distances.');
            testCase.verifyTrue(all(diff(rel_err) < 0), ...
                'Target error should decrease as targets move farther from the surface.');
        end

        function testExteriorNeumannSLPProblem(testCase)
            target_distance = 1e-1;
            u0 = 2/sqrt(3);
            [~, ~, ~, condK, info] = modified_laplace_bie_problem( ...
                testCase.p, testCase.lambda, u0, target_distance, false, "exterior_neumann");
            testCase.verifyLessThan(info.rel_soln_err, 5e-5, ...
                'Exterior Neumann solution is too inaccurate.');
        end
    end
end
