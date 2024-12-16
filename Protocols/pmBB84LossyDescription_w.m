%% FUNCTION NAME: pmBB84LossyDescription
% Lossy description for prepare-and-measure BB84.
% The observables correspond to a squashing model with five POVM outcomes
% (including photon loss).
%%

function protocolDescription = pmBB84LossyDescription_w(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz","fullstat"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimA = 4;
    dimB = 3;
    dimPB = 5;
    
    pz0=blkdiag([0],diag([1,0]));
    pz1=blkdiag([0],diag([0,1]));
    px0=blkdiag([0],[1,1;1,1]);
    px1=blkdiag([0],[1,-1;-1,1]);
    
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    krausOpZ0 = kron(kron(kron(kron(zket(2,1),diag([1,0,0,0])),sqrt(pz)*pz0)+kron(kron(zket(2,2),diag([0,1,0,0])),sqrt(pz)*pz1),[1;0]),[1;0]); % for Z basis
    krausOpZ1 = kron(kron(kron(kron(zket(2,1),diag([1,0,0,0])),sqrt(pz)*pz1)+kron(kron(zket(2,2),diag([0,1,0,0])),sqrt(pz)*pz0),[1;0]),[0;1]); % for Z basis
    krausOpX0 = kron(kron(kron(kron(zket(2,1),diag([0,0,1,0])),sqrt((1-pz)/2)*px0)+kron(kron(zket(2,2),diag([0,0,0,1])),sqrt((1-pz)/2)*px1),[0;1]),[1;0]); % for X basis
    krausOpX1 = kron(kron(kron(kron(zket(2,1),diag([0,0,1,0])),sqrt((1-pz)/2)*px1)+kron(kron(zket(2,2),diag([0,0,0,1])),sqrt((1-pz)/2)*px0),[0;1]),[0;1]); % for X basis
    krausOp = {krausOpZ0, krausOpZ1, krausOpX0, krausOpX1};
 
    % components for the pinching Z map
    keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2*2)); 
    keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2*2));
    keyMap = {keyProj1, keyProj2};

    % Constraints
    %observables = {};

    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)),'mask',0);
    end
    
    % Normalization
    addObservables(eye(dimA*dimB),'mask',0);
    
    %Constructing the constraints
    %basicBobPOVMs = {1/2*[1,0;0,0],1/2*[0,0;0,1],1/2*1/2*[1,1;1,1],1/2*1/2*[1,-1;-1,1]};
    basicBobPOVMs = {pz*pz0,pz*pz1,(1-pz)/2*px0,(1-pz)/2*px1,diag([1,0,0])};
%     basicBobPOVMs = {pz*[1,0;0,0],pz*[0,0;0,1],(1-pz)*1/2*[1,1;1,1],(1-pz)*1/2*[1,-1;-1,1]}; %no squashing model
    basicAlicePOVMs = {diag([1,0,0,0]),diag([0,1,0,0]),diag([0,0,1,0]),diag([0,0,0,1])};
    
    %Full set of bipartite POVMS
    bipartitePOVMs = cell(dimA*dimPB,1);
    for i = 1:dimA
        for j = 1:dimPB
            bipartitePOVMs{dimPB*(i-1)+(j-1)+1} = kron(basicAlicePOVMs{i},basicBobPOVMs{j});
        end
    end    
    
    if(fullstat==1)
        addObservables(bipartitePOVMs,'mask',1);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        %observables = [observables; bipartitePOVMs{select(1,2)} ; bipartitePOVMs{select(2,1)}; bipartitePOVMs{select(3,4)} ; bipartitePOVMs{select(4,3)}];
        newObservables = {bipartitePOVMs{select(1,2)} ; bipartitePOVMs{select(2,1)}; bipartitePOVMs{select(1,1)} ; bipartitePOVMs{select(2,2)}; ...
            bipartitePOVMs{select(3,4)} ; bipartitePOVMs{select(4,3)}; bipartitePOVMs{select(3,3)} ; bipartitePOVMs{select(4,4)}};
        addObservables(newObservables,'mask',1);
        
        %normalization
        temp = bipartitePOVMs{select(1,2)}+bipartitePOVMs{select(2,1)}+bipartitePOVMs{select(1,1)}+bipartitePOVMs{select(2,2)}+ ...
            bipartitePOVMs{select(3,4)}+bipartitePOVMs{select(4,3)}+ bipartitePOVMs{select(3,3)} + bipartitePOVMs{select(4,4)};
        
        addObservables(eye(dimA*dimB)-temp,'mask',1);
    end
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end