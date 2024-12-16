%Preset input for measurement-device-independent BB84 with asymptotic solver
%There is no squashing model, and channel only contains misalignment.

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = MDIBB84WCP_decoy_pD_w()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handle protocolDescription(parameters) and expectations(parameters); description is a struct with multiple fields
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('MDIBB84Description_w');
    channelModel=str2func('MDIBB84WCPChannel_pD_w');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

%     parameters.names = ["eta","pz","pd","f","thetaA","thetaB","mu1","mu2","mu3","etad","q"]; %MDIBB84
    parameters.names = ["eta","pz","pd","f","thetaA","thetaB","mu1","mu2","mu3"]; %MDIBB84
    
    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
%     parameters.scan.thetaA = 0:pi/64:pi/4; %misalignment angle between Alice-Charlie
    parameters.scan.eta = 10.^(-0.2*(0:5:120)/10); %channel transmittance between Alice-Charlie (half of total distance)

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    parameters.fixed.pz = 0.5; %basis choice probability (for Z basis); only Z basis is used for coding
    parameters.fixed.pd = 1e-6; %dark count probability
    parameters.fixed.f = 1; %error-correction efficiency
    parameters.fixed.thetaA = 0.183; %misalignment angle between Alice-Charlie sin^2(theta) 
    parameters.fixed.thetaB = 0; %misalignment angle between Bob-Charlie  sin^2(theta)
    parameters.fixed.mu1 = 0.4; %first intensity (signal)
    parameters.fixed.mu2 = 0.1; %second intensity (decoy 1)
    parameters.fixed.mu3 = 0.0007; %third intensity (decoy 2)

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
%     parameters.optimize.mu1 = [0.25,0.25,0.6]; %first intensity (signal)
    
    
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
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
    solverOptions.solver1.name = 'asymptotic_inequality';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-8; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = ''; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
end