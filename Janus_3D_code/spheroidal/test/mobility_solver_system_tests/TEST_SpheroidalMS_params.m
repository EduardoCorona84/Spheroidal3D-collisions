%{
Note SpheroidalMS_params is not an actual object within the codebase.
However, the code files relating to it are quite mysterious (wrt Rbs_mobility).
The test functions below are to illustrate (for me) the functionality of it.
%}

classdef TEST_SpheroidalMS_params < matlab.unittest.TestCase
    properties
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
        end
    end

    methods (Test)
        function testSpheroidalMSInitializeParameters(testCase)
            %{
            The idea is that the user should supply the system of spheres/spheroids.
            %}
            p = 16;
            C = [0 5 3];
            equ_radii = 1;
            polar_radii = 1/2;
            eps = 1e-3;
            mdist = 1e-2;
            out = true;
            dense = true;

            Fparams = struct();
            Fparams.parbd.p = p;
            Fparams.parbd.Ct = C;
            Fparams.parbd.equ_radii = equ_radii;
            Fparams.parbd.polar_radii = polar_radii;
            Fparams.parbd.eps = eps;
            Fparams.parbd.mdist = mdist;
            Fparams.parbd.out = out;
            Fparams.denseMV = dense;
            Fparams.tdisc = "euler";

            Fparams_updated = SpheroidalMS_initparams(Fparams);
            parbd = Fparams_updated.parbd;

            % Xrp vs. Xp
            %{
            Rotation seems to be done outside...
            Xrp represents spheres/spheroids in its local frame
            Xp represents UNROTATED spheres/spheroids in its global frame

            Not sure where Xp is used however...
            %}
            Xp = parbd.Xp;
            Xrp = parbd.Xrp;
            scatter3(Xrp(:,1), Xrp(:,2), Xrp(:,3));
            hold on;
            scatter3(Xp(:,1), Xp(:,2), Xp(:,3));
            hold off;
        end
    end
end
