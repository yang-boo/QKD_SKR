%% FUNCTION NAME: pmBB84Description
% Simple description for prepare-and-measure BB84. 
% The expectations correspond to non-squashing model with four POVM outcomes.
% 函数名：pmBB84Description
% 准备和测量 BB84 的简单描述。
% 期望值对应于有四个 POVM 结果的非挤压模型。
%%

function protocolDescription = pmBB84Description_2023using(names,p)

    %user should supply the list of parameters used in this description/channel file 用户应提供该description/channel文件中使用的参数列表
    %this list varNames should be a subset of the full parameter list declared in the preset file 该列表 varNames 应是预置文件中声明的完整参数列表的子集
    %parameters specified in varNames can be used below like any other MATLAB variables 在 varNames 中指定的参数可以像其他 MATLAB 变量一样在下面使用
    varNames=["pz"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    %函数 findVariables 和 addVariables 会根据 varNames 自动搜索输入 (names,p) 中的参数值，并将其转换为 MATLAB 变量。
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    %可通过调用 addObservables(x) 或 addObservables(x,'mask',maskValue) 自动填入该描述文件的输出。
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin 用户提供的说明开始 %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimA = 2;
    dimB = 2;
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];

    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    % kraus算子对 G 映射进行后处理。寄存器的排序为 R、A、B，二维公告寄存器（Alice 和 Bob 的公告寄存器经过筛选后组合在一起）
    krausOpZ = kron(kron(kron(zket(2,1),sqrt(pz) * diag([1,0]))+kron(zket(2,2),sqrt(pz) * diag([0,1])),sqrt(pz) * diag([1,1])),[1;0]); % for Z basis
    krausOpX = kron(kron(kron(zket(2,1),sqrt((1-pz)/2) * (ketPlus*ketPlus'))+kron(zket(2,2),sqrt((1-pz)/2) * (ketMinus*ketMinus')),sqrt((1-pz)) * (diag([1,1]))),[0;1]); % for X basis
    
    krausOp = {krausOpZ, krausOpX};
 
    % components for the pinching Z map 挤压 Z 映射的成分
    keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2)); 
    keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2));
    keyMap = {keyProj1, keyProj2};

    % Constraints 约束
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)));
    end
    
    % Z and X constraints Z 和 X 约束条件
%     addObservables(kron(diag([1,0]),diag([0,1])) + kron(diag([0,1]), diag([1,0])));
%     addObservables(kron(ketPlus * ketPlus',ketMinus * ketMinus') + kron(ketMinus * ketMinus', ketPlus * ketPlus'));
    
    addObservables(kron(diag([1,0]),diag([0,1])) + kron(diag([0,1]), diag([1,0])));
    addObservables(kron(ketPlus * ketPlus',ketMinus * ketMinus') + kron(ketMinus * ketMinus', ketPlus * ketPlus'));
    
%     % Cross terms
%     addObservables(kron(diag([1,-1,0,0]), ketPlus * ketPlus' - ketMinus * ketMinus'));
%     addObservables(kron(diag([0,0,1,-1]), diag([1,-1])));

    % Normalization 
    addObservables(eye(dimA*dimB));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end 用户提供的说明结束 %%%%%%%%%%%%%%%%%%%%%%%%%
    
    protocolDescription.observables = observables;
    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.dimensions = [dimA,dimB];

end