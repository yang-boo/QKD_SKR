%Preset input for prepare-and-measure BB84 with asymptotic solver
%There is squashing model for Bob's POVM, and channel model contains loss, misalignment, and dark count.
%In this preset we are using the asymptotic solver.

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = pmBB84_asymptotic_w()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handles of description/channel/error-correction files (functions)
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('pmBB84Description_2023using_w');
    channelModel=str2func('pmBB84Channel_2023using');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["ed","pz","pd","eta","etad","f","fullstat"]; %BB84 Lossy

    %%%%%%%%%%%%%%%% 1.parameter to scan over 参数扫描 %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    %必须命名至少一个要扫描的参数（如果只对固定值感兴趣，可以是单点数组）
    
    parameters.scan.eta = 10.^( -0.2*(0:5:0) /10); %channel transmittance 信道透射率
    
    %%%%%%%%%%%%%%%% 2.fixed parameters 固定参数 %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    %可选；常量值可以是数值或数组/矩阵
    

    parameters.fixed.ed = 0.033;
    %misalignment, defined as single-photon error, i.e. sin^2(theta) 
    % 偏差，定义为单光子误差，即 sin^2(theta)
    parameters.fixed.pz = 0.5; 
    %basis choice probability (for Z basis) 基选择概率（针对 Z 基）
    parameters.fixed.pd = 1e-5; %1e-6; 
    %dark count probability 暗计数率
    parameters.fixed.etad = 0.2; %0.045; 
    %detector efficiency 探测器效率
    parameters.fixed.fullstat = 1; 
    %using full statistics or using QBER/Gain observables only
    %使用完整的统计数据或仅使用 QBER/Gain 可观测数据
    parameters.fixed.f = 1; 
    %error correction efficiency 纠错效率


    %%%%%%%%%%%%%%%% 3.optimizable parameters 可优化参数 %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    %可选；声明可优化参数会自动调用本地搜索优化器，
    %格式必须为 [lowerBound, initialValue, upperBound] （下限值、初始值、上限值）
    
%     parameters.optimize.pz = [0.1,0.5,0.9];
%     parameters.optimize.f = [1.0,1.2,2.0];
end

%set the running options for the solvers
%设置求解器的运行选项
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration 每个主迭代的输出
    %2.output at each solver 1 FW iteration 每个求解器 1 FW 迭代的输出
    %3.show all SDP solver output 显示所有 SDP 求解器输出
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting 参数优化器设置 %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; 
    %choose between 'coordinateDescent' and 'bruteForce'
    %在 "coordinateDescent "和 "bruteForce "之间做出选择
%     solverOptions.optimizer.name = 'localSearch_Adam';
    solverOptions.optimizer.linearResolution = 3; 
    %resolution in each dimension (for brute force search and coordinate descent)
    %每个维度的分辨率（用于暴力搜索和坐标下降）
    solverOptions.optimizer.maxIterations = 1;  
    %max number of iterations (only for coordinate descent)
    %最大迭代次数（仅适用于坐标下降）
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; 
    %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    %选择内置的 "fminbnd "算法和自定义的 "iterative"算法进行线性搜索（仅适用于坐标下降算法）
    solverOptions.optimizer.iterativeDepth = 2; 
    %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    %选择迭代级别的深度；函数计数 = 深度 * 线性分辨率（仅适用于坐标下降和使用 "iterative "时）
    solverOptions.optimizer.maxSteps = 5; 
    %max number of steps (only for gradient descent and ADAM)
    %最大步数（仅适用于梯度下降和 ADAM）
    solverOptions.optimizer.optimizerVerboseLevel = 2; 
    %0:only output optimized result; 1:output a progress bar; 2:output at each function call
    %0:只输出优化结果；1:输出进度条；2:每次调用函数时输出

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; 
    %1e-6 for asymptotic, 2.5e-3 for finite;
    %1e-6（渐近），2.5e-3（有限）；
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; 
    %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    %最小化 norm(rho0-rho) 或 -lambda_min(rho)，有限大小使用方法 1，渐近 v1 使用方法 2
    
    %default options默认选项
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; 
    %true for testing gap, false for testing f1-f0
    %测试 gap 时为 true，测试 f1-f0 时为 false
    solverOptions.solver1.removeLinearDependence = 'rref'; 
    %'qr' for QR decomposition, 'rref' for row echelon form, 'empty' for not removing linear dependence
    %'qr'表示 QR 分解，'rref'表示行梯形，'空'表示不去除线性相关性

    %%%%%%%%%%%%%%%%% step2Solver setting step2Solver 设置 %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end