%% FUNCTION NAME: primalDf
% This function calculates the gradient of primal problem objective
% function.该函数计算原始问题目标函数的梯度。
%%
function dfval = primalDf(rho,keyMap,krausOperators)

    if nargin == 2 || isempty(krausOperators)
        % if there is no post-selection map 如果没有后选择映射
        
        [rho,~]=perturbation_channel(rho);
        
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*rho*keyMap{j};
        end
        [zRho,~]=perturbation_channel(zRho);
        
        dfval = logm(rho)-logm(zRho);
    else
        % if there is a post-selection map.如果有后选择映射
        
%         disp(rho);
        gRho = krausFunc(rho,krausOperators);
        
%         if ishermitian(gRho)
%             disp('矩阵是Hermitian矩阵');
%         else
%             disp('矩阵不是Hermitian矩阵');
%         end
        
        %check validity of gRho and perform perturbation if not valid 检查 gRho 的有效性，如果无效则执行扰动
        [gRho,~]=perturbation_channel(gRho);
        
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*gRho*keyMap{j}';
        end
        
        %check validity of zRho and perform perturbation if not valid 检查 zRho 的有效性，如果无效则执行扰动
        [zRho,~]=perturbation_channel(zRho);
        
        dfval = krausFunc(logm(gRho)-logm(zRho),krausOperators,'transpose');
    end

end