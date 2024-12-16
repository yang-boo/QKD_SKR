%% FUNCTION NAME: generalEC
% This function returns a function handle of error-correction leakage, for general formulation of protocols
% 此函数返回纠错泄漏的函数句柄，用于协议的一般表述
% The return function depends on channelModel only, which contains QBER or probability distribution (depending on fine/coarse grain model) and probability of sifting.
% 返回函数只取决于 channelModel，其中包含 QBER 或概率分布（取决于细粒/粗粒模型）以及筛选概率。
%%

function [leakageEC]  = generalEC_BB84(channelModel,names,p)
    %this list varNames should be a subset of the full parameter list declared in the preset file
    varNames = ["ed","pz","pd","eta","etad","f","mu1","mu2","mu3","active","fullstat"];

    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    %from here on the parameters specified in varNames can be used like any other MATLAB variables
    % 函数 findVariables 和 addVariables 会根据 varNames 自动搜索输入 (names,p) 中的参数值，并将其转换为 MATLAB 变量。
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pSift can be an array (e.g. [pz^2, (1-pz)^2])
    leakageEC= 0;
    errorRate = channelModel.errorRate;
    siftingFactor = channelModel.pSift;
    dimension = length(siftingFactor);
    
    if(~isfield(channelModel,'isCVQKD') || channelModel.isCVQKD == 0)
%         DVQKD
        ita_b = 0.5; % Bob端物理系统的光路透过率
        Y0 = pd;    % 暗计数响应率，忽略了double click的概率
        ita = etad*eta*ita_b;     %一定距离下，总的透过率
        e0 = 0.5;       % 暗计数的误码率，通常在不认为 Eve 能攻击探测器的时候均约为0.5

        Qs = Y0 + 1 - exp(-ita*mu1); %信号光的响应率
        Es = (e0*Y0+etad*(1-exp(-ita*mu1)))/Qs; %信号光的误码率
        if(isfield(channelModel,'QU'))
            Qd = Y0 + 1 - exp(-ita*mu2); %诱骗态光的响应率

            Qv = Y0 + 1 - exp(-ita*mu3); %真空态光的响应率

            Y0L = max((mu2*Qv*exp(mu3)-mu3*Qd*exp(mu2))/(mu2-mu3),0); %光子数为0的脉冲的计数率的下限
            Y1L = mu1/(mu1*mu2-mu1*mu3-mu2^2+mu3^2)*(Qd*exp(mu2)-Qv*exp(mu3)-(mu2^2-mu3^2)/mu1^2*(Qs*exp(mu1)-Y0L));%光子数为1的脉冲的计数率

            Y1U = (Qd*exp(mu2)-Qv*exp(mu3))/(mu2-mu3); % 我们推到的光子数为1的脉冲的计数率的上限
            Q1U = Y1U*exp(-mu1)*mu1; % 我们推到的光子数为1的脉冲的总计数率的上限
            Y0U = Qv*exp(mu3)-Y1L*mu3; %光子数为0的脉冲的计数率的上限
            Q0U = Y0U*exp(-mu1);

            leak1 = Es*2;
            leak = f * binaryEntropy(Es);                        
            leak2 = leak - leak1;
            leakU = (Q0U+Q1U)/Qs*leak1 + leak2;
            leakageEC = Qs*leakU;             
        end
    else       
        %CVQKD
        % DR: EC= (IXY+HX_Y)*pSift - f * IXY*pSift, where f<=1.  
        % RR: EC= (IXY+HY_X)*pSift - f * IXY*pSift, where f<=1. 
        if(~isfield(channelModel,'probDist'))
            fprintf('need a probability distribution for error correction!\n');
            return;
        end
        probDist = channelModel.probDist;
        if(isfield(channelModel,'recon') & channelModel.recon == 1)
            %use reverse reconciliation
            for i = 1:dimension
                if(dimension == 1)
                    [HX_Y, HY_X, IXY] = calculateEC(probDist);
                else
                    [HX_Y, HY_X, IXY] = calculateEC(probDist{i});
                end
                leakageEC = leakageEC +(IXY+HY_X)*siftingFactor(i) - f*IXY*siftingFactor(i);
            end
        else
            %use direct reconciliation
            for i = 1:dimension
                if(dimension == 1)
                    [HX_Y, HY_X, IXY] = calculateEC(probDist);
                else
                    [HX_Y, HY_X, IXY] = calculateEC(probDist{i});
                end
                leakageEC = leakageEC +(IXY+HX_Y)*siftingFactor(i) - f*IXY*siftingFactor(i); 
            end
        end
    end    
end

function [HX_Y, HY_X, IXY] = calculateEC(prob_dist)   	
    %a direct call with single return value will return H(X|Y)
    %直接调用时只有一个返回值，将返回 H(X|Y)

    px_dist = sum(prob_dist,2);
    py_dist = sum(prob_dist)';
    
    HXY = - prob_dist(prob_dist > 0)' * log2(prob_dist(prob_dist > 0));
    HX = - px_dist(px_dist > 0)' * log2(px_dist(px_dist > 0));
    HY = - py_dist(py_dist > 0)' * log2(py_dist(py_dist > 0));
	
    HY_X = HXY - HX;
    HX_Y = HXY - HY;
    IXY  = HY - HY_X;
end