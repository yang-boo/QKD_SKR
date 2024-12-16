%% FUNCTION NAME: krausFunc
%
% This function works as the Kraus operator description of a quantum
% channel. 
% 这个函数就像量子通道的Kraus算子描述一样。
%%

function rhoPrime = krausFunc(rho,krausOperators,transpose)

%     disp(rho);
%     if ishermitian(rho)
%         disp('矩阵rho是Hermitian矩阵');
%     else
%         disp('矩阵rho不是Hermitian矩阵');
%     end
    
    if nargin == 2 || isempty(transpose)
        rhoPrime = 0;
        
        for i = 1:numel(krausOperators)
            rhoPrime = rhoPrime + krausOperators{i}*rho*krausOperators{i}';
            
%            disp(isequal(krausOperators{i}*rho*krausOperators{i}',(krausOperators{i}*rho*krausOperators{i}')'));
            
%             if ishermitian(krausOperators{i}*rho*krausOperators{i}')
%                 disp('矩阵是Hermitian矩阵');
%             else
%                 disp('矩阵不是Hermitian矩阵');
%             end
            
        end
        
%         if ishermitian(rhoPrime)
%             disp('矩阵rhoPrime是Hermitian矩阵');
%         else
%             disp('矩阵rhoPrime不是Hermitian矩阵');
%         end
        
    else
        rhoPrime = 0;
        for i = 1:numel(krausOperators)
            rhoPrime = rhoPrime + krausOperators{i}'*rho*krausOperators{i};
            
%             if ishermitian(rhoPrime)
%                 disp('矩阵rhoPrime是Hermitian矩阵');
%             else
%                 disp('矩阵rhoPrime不是Hermitian矩阵');
%             end
            
        end
    end

end