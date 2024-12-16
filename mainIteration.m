%% FUNCTION NAME: mainIteration
% The main iteration which processes incoming parameter and description
% files, and calls the solver module many times to calculate the key rate.主要迭代处理传入的参数和说明文件，并多次调用求解器模块计算密钥率。
% The main iteration scans through all sample data points, and optionally 
% performs parameter optimization.主迭代扫描所有样本数据点，并执行参数优化。
% can optionally uncomment the "parfor" on line 51 and comment out "for" 
% to accelerate with multithread/multiple computers 可以选择取消注释第 51 行的 "parfor"，并注释掉 "for"，以加快多线程/多台计算机的运行速度

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions)
    
    helper=helperFunctions; %load helper function file
    [parameters.names,parameters.order]=helper.getOrder(parameters); %generate sorting order or parameters (depending on the given name list)

    %process the scannable (variable) parameters
    parameters_scan = struct2cell(parameters.scan);
    dimensions=helper.getCellDimensions(parameters_scan); %read the dimensions of each parameter to be scanned
    N=prod(dimensions); %total dimension of data to scan要扫描的数据的总维度

    %process the fixed parameters
    if(~isfield(parameters,'fixed'))
        %no parameters to optimize
        parameters.fixed = struct;
    end
    
    %process the optimizable parameters
    if(~isfield(parameters,'optimize'))
        %no parameters to optimize
        parameters.optimize = struct;
        p_lower=[];
        p_start=[];
        p_upper=[];
        isOptimizing = false;
    else
        optimize_list=cell2mat(struct2cell(parameters.optimize));
        p_lower=optimize_list(:,1)';
        p_start=optimize_list(:,2)';
        p_upper=optimize_list(:,3)';
        isOptimizing = true;
    end

    % %main iteration that test each point in parameters.scan
    % %can optionally uncomment the "parfor" line and comment out "for" to accelerate with multithread/multiple computers
    % %(if parallel computing toolbox is available)
    % %if using parfor, it is advised to set verbose levels to none (here, or in preset file) like below:
%     测试parameters.scan中每个点的主迭代
%     可以选择性地取消注释“parfor”行并注释掉“for”以使用多线程/多台计算机加速
%     （如果有并行计算工具箱）
%     如果使用parfor，建议将详细级别设置为none（此处或预设文件中），如下所示：

%     solverOptions.globalSetting.verboseLevel = 0; %outputs only iteration number
%     solverOptions.optimizer.optimizerVerboseLevel = -1; %outputs nothing when optimizing

%     parfor i=1:N
      for i=1:N
        if(solverOptions.globalSetting.verboseLevel >= 1)
            fprintf('main iteration: %d\n',i);
        end
        
        %%%%%%%%%%%%%% process single-point parameter %%%%%%%%%%%%%%
        indices=helper.expandIndices(i, dimensions); %get indices of parameters.scan 获取parameters.scan的索引
        p_scan = helper.selectCellRow(indices,parameters_scan);
        p_fixed = struct2cell(parameters.fixed)';

        %%%%%%%%%%%%%% optimize parameter (optional) %%%%%%%%%%%%%%
        
        %optimization of parameters (if parameters.optimize is non-empty)
        if(isOptimizing)
            %wrap the key rate function as a single-input function of "p" (numerical array of optimizable parameters)
            rateFunction = @(p) getKeyRate_wrapper(parameters.names,helper.reorder([p_scan,p_fixed,p],parameters.order),protocolDescription,channelModel,leakageEC,solverOptions); %convert getKeyRate as a single-value function f(p) for optimization
            if(solverOptions.globalSetting.verboseLevel >= 1)
                fprintf('begin optimization\n');
            end
            if(strcmp(solverOptions.optimizer.name,'bruteForce'))
                p_optimal = helper.bruteForceSearch(rateFunction,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'coordinateDescent'))
                p_optimal = helper.coordinateDescent(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'gradientDescent'))
                p_optimal = helper.gradientDescent(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'localSearch_Adam'))
                p_optimal = helper.localSearch_Adam(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            end
            if(solverOptions.globalSetting.verboseLevel >= 1)
                fprintf('finished optimization\n');
            end
        else
            p_optimal = [];
        end

        %%%%%%%%%%%%%% evaluate descriptions %%%%%%%%%%%%%%
        
        %generation of single-row parameter list (a cell array)
        %生成单行参数列表（单元格数组）
        p_scan = num2cell(p_scan); %convert p_scan to cell array
        if(~isempty(p_optimal))
            p_optimal = num2cell(p_optimal); %convert p_optimal to cell array
        end
        p_full=[p_scan,p_fixed,p_optimal]; %concatenate with p_fixed
        p_full = helper.reorder(p_full,parameters.order);
        
        %evaluate the protocol description, channel model, and leakage
        thisProtocolDescription=protocolDescription(parameters.names,p_full);
        thisChannelModel=channelModel(thisProtocolDescription,parameters.names,p_full);
        thisLeakageEC=leakageEC(thisChannelModel,parameters.names,p_full);
        

        %%%%%%%%%%%%%% perform calculation %%%%%%%%%%%%%%
        
        %calculate key rate by calling the solver module通过调用求解器模块计算密钥率
        %note that the full parameter list is also passed into the function for the solver to optionally directly access (e.g. security parameter eps, data size N, etc.).
        %注意，完整的参数列表也会传入函数，供求解器选择直接访问（如安全参数 eps、数据大小 N 等）。
        [results(i).lowerBound,results(i).upperBound,results(i).FWBound,results(i).debugInfo] = getKeyRate(thisProtocolDescription,thisChannelModel,thisLeakageEC,solverOptions,p_full,parameters.names);
        
        %also save the current parameter set for reference (1-D cell array)
        %还保存当前参数集以供参考（1-D cell array）
        results(i).debugInfo.current_parameters = p_full;
        results(i).debugInfo.names = parameters.names;
    end
    
    
end

%helper function used for optimization algorithm 用于优化算法的辅助函数
function rate = getKeyRate_wrapper(names,p_full,protocolDescription,channelModel,leakageEC,solverOptions)

    thisProtocolDescription=protocolDescription(names,p_full);
    thisChannelModel=channelModel(thisProtocolDescription,names,p_full);
    thisLeakageEC=leakageEC(thisChannelModel,names,p_full);
    
    solverOptions.solver1.maxiter = 3; %reduce precision for faster speed when optimizing
    solverOptions.globalSetting.verboseLevel = 0;
    
    [lowerBound,~,~,~] = getKeyRate(thisProtocolDescription,thisChannelModel,thisLeakageEC,solverOptions,p_full,names);
    
    rate = lowerBound;
end