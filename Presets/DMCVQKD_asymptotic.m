%Preset input for discrete-modulated continuous variable QKD
%In this preset we are using the asymptotic solver.

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = DMCVQKD_asymptotic()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handle protocolDescription(parameters) and expectations(parameters); description is a struct with multiple fields
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('DMCVheterodyneDescription');
    channelModel=str2func('DMCVheterodyneChannel');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["amp_ps","phase_ps","cutoffN","recon","eta","noise","alphaValue","f"]; %DMCVQKD

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    parameters.scan.eta = 10.^(-0.2*(50:5:50)/10); %channel transmittance
    

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    
    parameters.fixed.amp_ps = 0; %amplitude post-selection
    parameters.fixed.phase_ps = 0; %phase post-selection
    parameters.fixed.cutoffN = 12; %photon number cutoff
    parameters.fixed.recon = 1; %use reverse reconciliation, 1 for on
    parameters.fixed.noise = 0.005; %excess noise ksi
    parameters.fixed.alphaValue = 0.5; %signal amplitude
    parameters.fixed.f = 0.95; %error-correction efficiency - note that for CVQKD, by convention it is called beta and <1


    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
%     parameters.optimize.alphaValue = [0.2,0.5,0.7];
    
    
end

%set the running options for the solvers
function solverOptions=setOptions()


    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'sdpt3';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 2; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 10; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 1; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 3;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
end