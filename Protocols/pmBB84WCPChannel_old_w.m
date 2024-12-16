% %% FUNCTION NAME: pmBB84WCPChannel
% Realistic channel model for prepare-and-measure BB84 using WCP source and decoy states. 
% The transmittance eta, misalignment ed, and dark count pd are considered.
% The expectations correspond to a squashing model with five POVM outcomes (including photon loss).
%
% The flags active and fullstat can switch between active/passive
% detection, and coarse-grain/fine-grain statistics.
%
% Decoy state analysis is performed inside this channel model function.
%
% Additional information including masks (denoting which statistics are
% bounded by decoy state analysis and therefore uncertain) 
% and pSignal (denoting single photon probability) are included.
%%

function channelModel = pmBB84WCPChannel_old_w(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pz","pd","eta","etad","mu1","mu2","mu3","active","fullstat"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    decoys = [mu1,mu2,mu3];
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    nDecoys = length(decoys);
    dimPB = 5;
    px = 1 - pz;
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    signalStates = {[1;0], [0;1], ketPlus, ketMinus};
    probList = [pz/2; pz/2; (1-pz)/2; (1-pz)/2];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    %对飞行量子比特系统进行部分跟踪，以获得局部 rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    expectations = [];
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    
    %simulating the channel
    %模拟信道
    for i=1:length(decoys)
        rawExpectations(:,:,i)=coherentSourceChannel(active,decoys(i),eta,etad,pd,ed,px);
    end
    rawExpectations_signal=coherentSourceChannel(active,decoys(1),eta,etad,pd,ed,px);
    
    %convert 16-D pattern to 5-D POVM
    %Patterns:
    %column: [H,V,+,-,None]
    %row: Alice POVM
    MappingH=zeros(16,1);
    MappingH(9)=1; %1000
    MappingH(13)=0.5; %1100
    MappingV=zeros(16,1);
    MappingV(5)=1; %0100
    MappingV(13)=0.5; %1100
    
    MappingD=zeros(16,1);
    MappingD(3)=1; %0010
    MappingD(4)=0.5; %0011
    
    MappingA=zeros(16,1);
    MappingA(2)=1; %0001
    MappingA(4)=0.5; %0011
    
    MappingN=zeros(16,1);
    MappingN(1)=1; %0000
    MappingN(6)=1; %0101
    MappingN(7)=1; %0110
    MappingN(8)=1; %0111
    MappingN(10)=1; %1001
    MappingN(11)=1; %1010
    MappingN(12)=1; %1011
    MappingN(14)=1; %1101
    MappingN(15)=1; %1110
    MappingN(16)=1; %1111
    
    Mapping=[MappingH,MappingV,MappingD,MappingA,MappingN];

%     %perform decoy analysis on fine-grained raw analysis
%     %对细粒度原始分析进行诱饵分析
%     L1=size(rawExpectations,1);
%     L2=size(rawExpectations,2);
%     Y1_expectations=zeros(L1,L2,2);
%     for i=1:L1
%         for j=1:L2
%             decoy_expectations=squeeze(rawExpectations(i,j,:))';
%             [Y1_expectations(i,j,1),Y1_expectations(i,j,2)]=decoyAnalysis(decoys,decoy_expectations);
%         end
%     end
%     
%     Y1_expectations(:,:,1)
%     Y1_expectations(:,:,2)
%     Y1L = Y1_expectations(:,:,1)*Mapping;
%     Y1U = Y1_expectations(:,:,2)*Mapping;
    
    %first bin into squashed model, then perform decoy analysis
    %首先将数据输入压缩模型，然后进行诱饵分析
    for i=1:length(decoys)
        bipartiteExpectationsWCP(:,:,i)=rawExpectations(:,:,i)*Mapping;
    end

    L1=size(bipartiteExpectationsWCP,1);
    L2=size(bipartiteExpectationsWCP,2);
    Y1L=zeros(L1,L2);
    Y1U=zeros(L1,L2);
    Y0U=zeros(L1,L2);
    for i=1:L1
        for j=1:L2
            decoy_expectations=squeeze(bipartiteExpectationsWCP(i,j,:))';
            [Y1L(i,j),Y1U(i,j),Y0U(i,j)]=decoyAnalysis(decoys,decoy_expectations);
        end
    end
    
    
    %normalize by signal state probabilities
    %按信号态概率进行归一化
    Y1L = diag(probList)*Y1L;
    Y1U = diag(probList)*Y1U;
    
    Y1L_1D = zeros(4*5,1);
    Y1U_1D = zeros(4*5,1);
    for i = 1:4
        for j = 1:5
            Y1L_1D(5*(i-1)+(j-1)+1) = Y1L(i,j);
            Y1U_1D(5*(i-1)+(j-1)+1) = Y1U(i,j);
        end
    end
    
    if(fullstat==1)
        addExpectations(Y1L_1D,'mask',1);
        addExpectations(Y1U_1D,'mask',2);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        coarseGrainExpectations_Y1L = [Y1L_1D(select(1,2));Y1L_1D(select(2,1));Y1L_1D(select(1,1));Y1L_1D(select(2,2));...
                Y1L_1D(select(3,4));Y1L_1D(select(4,3));Y1L_1D(select(3,3));Y1L_1D(select(4,4))];
        coarseGrainExpectations_Y1U = [Y1U_1D(select(1,2));Y1U_1D(select(2,1));Y1U_1D(select(1,1));Y1U_1D(select(2,2));...
                Y1U_1D(select(3,4));Y1U_1D(select(4,3));Y1U_1D(select(3,3));Y1U_1D(select(4,4))];
        %normalization
        normL = 1-sum(coarseGrainExpectations_Y1U);
        normU = 1-sum(coarseGrainExpectations_Y1L);
        
        addExpectations(coarseGrainExpectations_Y1L,'mask',1);
        addExpectations(normL,'mask',1);
        addExpectations(coarseGrainExpectations_Y1U,'mask',2);
        addExpectations(normU,'mask',2);
    end
    
    
    %QBER and Gain statistics for error correction
    signal_simulation = diag(probList)*rawExpectations_signal*Mapping;
    gainz = signal_simulation(1,1) + signal_simulation(1,2) + signal_simulation(2,1) + signal_simulation(2,2);
    gainx = signal_simulation(3,3) + signal_simulation(3,4) + signal_simulation(4,3) + signal_simulation(4,4); 
    
    errorz = signal_simulation(1,2) + signal_simulation(2,1);
    errorx = signal_simulation(3,4) + signal_simulation(4,3); 

    %signal state proportion
    %信号状态比例
    P1 = decoys(1)*exp(-decoys(1));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorz/gainz,errorx/gainx];  
    channelModel.pSift = [gainz,gainx]; 
    channelModel.pSignal = P1;
    channelModel.flag = 0; %0:eR; 1:pD
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%helper function that performs decoy state analysis
%辅助函数，用于执行诱饵态分析
function [Y1L,Y1U,Y0U] = decoyAnalysis(decoys,decoy_expectations)

    cvx_solver mosek

    n_photon=10;
    n_decoy=size(decoys,2);
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);
    decoy_tolerance=1e-12;
    Obj=zeros(1,n_photon+1);
    Obj(2)=1;
    
    Obj_0=zeros(1,n_photon+1);
    Obj_0(1)=1;
    
    %solve for Y1 upper bound
    try
        cvx_begin quiet
            variable Y(n_photon+1)
            maximize Obj*Y
            for k=1:n_photon+1
               Y(k)<=1;
               Y(k)>=0;
            end
            for i = 1:n_decoy
                C=zeros(1,n_photon+1);
                Ptotal=0;
                for k=1:n_photon+1
                    P=Poisson(decoys(i),k-1);
                    C(k)=P;
                    Ptotal=Ptotal+P;
                end
                Ptotal;
                %if(decoy_expectations(i)>0)
                    %ignore the zero observables
                    C*Y<=decoy_expectations(i)+decoy_tolerance;
                    C*Y>=decoy_expectations(i)-decoy_tolerance-(1-Ptotal);
                %end
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        end
    catch
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
    end
    
    Y1U=Y(2);
    
    %solve for Y1 lower bound
    try
        cvx_begin quiet
            variable Y(n_photon+1)
            minimize Obj*Y
            for k=1:n_photon+1
               Y(k)<=1;
               Y(k)>=0;
            end
            for i = 1:n_decoy
                C=zeros(1,n_photon+1);
                Ptotal=0;
                for k=1:n_photon+1
                    P=Poisson(decoys(i),k-1);
                    C(k)=P;
                    Ptotal=Ptotal+P;
                end
                Ptotal;
                %if(decoy_expectations(i)>0)
                    %ignore the zero observables
                    C*Y<=decoy_expectations(i)+decoy_tolerance;
                    C*Y>=decoy_expectations(i)-decoy_tolerance-(1-Ptotal);
                %end
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        end
    catch 
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
    end
    
    Y1L=Y(2);
    
        %solve for Y0 upper bound
    try
        cvx_begin quiet
            variable Y(n_photon+1)
            maximize Obj_0*Y
            for k=1:n_photon+1
               Y(k)<=1;
               Y(k)>=0;
            end
            for i = 1:n_decoy
                C=zeros(1,n_photon+1);
                Ptotal=0;
                for k=1:n_photon+1
                    P=Poisson(decoys(i),k-1);
                    C(k)=P;
                    Ptotal=Ptotal+P;
                end
                Ptotal;
                %if(decoy_expectations(i)>0)
                    %ignore the zero observables
                    C*Y<=decoy_expectations(i)+decoy_tolerance;
                    C*Y>=decoy_expectations(i)-decoy_tolerance-(1-Ptotal);
                %end
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        end
    catch
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
    end
    
    Y0U=Y(1);
end

%helper function that simulates the channel for WCP sources
%辅助函数，用于模拟 WCP 信号源的信道
function expectations=coherentSourceChannel(active,mu,eta,etad,pd,ed,px)
    expectations=zeros(4,16);
    pz=1-px;
    t=eta*etad; %total transmittance
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);

    theta=asin(sqrt(ed));
    PL=sin(pi/4-theta)^2;
    PU=cos(pi/4-theta)^2;
    
    mapping_passive=[[pz*(1-ed),pz*ed,px*PU,px*PL];[pz*ed,pz*(1-ed),px*PL,px*PU];[pz*PL,pz*PU,px*(1-ed),px*ed];[pz*PU,pz*PL,px*ed,px*(1-ed)]];
    mapping_active=[[1-ed,ed,PU,PL];[ed,1-ed,PL,PU];[PL,PU,1-ed,ed];[PU,PL,ed,1-ed]];
    
    for input=1:4
        %iterating over each input state
        
        for output=1:16
            %iterating over each pattern
            a=index1to4(output-1); %detector event, 4 elements corresponding to [H,V,D,A]
            
            Ppattern=1;
            if(active==0)
                %passive basis choice
                for k=1:4
                    %iterating over each detector
                    Pclick=1-Poisson(mu*t*mapping_passive(input,k),0); 
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    if(a(k)==1)
                        Ppattern=Ppattern*Pclick;
                    elseif(a(k)==0)
                        Ppattern=Ppattern*(1-Pclick);
                    end
                end
            else
                %active basis choice
                PpatternZ=1;
                prob_activeZ=[1,1,0,0];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeZ(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternZ=PpatternZ*Pclick;
                    elseif(a(k)==0)
                        PpatternZ=PpatternZ*(1-Pclick);
                    end
                end
                PpatternX=1;
                prob_activeX=[0,0,1,1];
                for k=1:4
                    %iterating over each detector
                    Pclick=(1-Poisson(mu*t*mapping_active(input,k),0));
                    Pclick=1-(1-Pclick)*(1-pd); %effect of dark count
                    Pclick=prob_activeX(k)*Pclick; %effect of basis choice (active)
                    if(a(k)==1)
                        PpatternX=PpatternX*Pclick;
                    elseif(a(k)==0)
                        PpatternX=PpatternX*(1-Pclick);
                    end
                end
                Ppattern=px*PpatternX+pz*PpatternZ;
            end
            expectations(input,output)=Ppattern;
        end
    end
end

%takes in an index of 1 and convert to array of 4
%以 1 为索引，转换为 4 数组
function a=index1to4(index)
    a(1) = floor(index/8);
    index = mod(index,8);
    a(2) = floor(index/4);
    index = mod(index,4);
    a(3) = floor(index/2);
    index = mod(index,2);
    a(4) = index;
end