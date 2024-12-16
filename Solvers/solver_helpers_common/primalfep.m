%% FUNCTION NAME: primalfep
%  Calculate $f_{\epsilon}(\rho)$ function, where $\epsilon$ value is
%  carefully chosen. 
%  计算 $f_{\epsilon}(\rho)$ 函数，其中 $\epsilon$ 的值是精心选择的。
%%

function [fval,realEpsilon] = primalfep(rho,keyMap,krausOperators,options)
%     
%     defaultoptions.epsilon = 0; % 0<=epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
%     defaultoptions.perturbation = 1e-16; % a small value added to the minimum eigenvalue;
%     if nargin == 4
%         if ~isfield(options,'epsilon')
%             options.epsilon = defaultoptions.epsilon;
%         end
%         if ~isfield(options,'perturbation')
%             options.perturbation = defaultoptions.perturbation;
%         end
%     else 
%         options = defaultoptions;
%     end


    if nargin == 3 || isempty(krausOperators)
        % in the case that there is no post-selection (so no Kraus operator).
        dim = size(rho,1);
        [rho,epsilon1] = perturbation_channel(rho);
   
    
        zRho = 0;
        for jMapElm = 1 : numel(keyMap)
            zRho = zRho + keyMap{jMapElm} * rho * keyMap{jMapElm};
        end
        
        [zRho,epsilon2] = perturbation_channel(zRho);
        realEpsilon = max(epsilon1, epsilon2);
        fval = real(trace(rho * (logm(rho) - logm(zRho))));
    
    else
        % in the case that there is a post-selection map.
        [rho,epsilon1] = perturbation_channel(rho);
        gRho = krausFunc(rho,krausOperators);
        [gRho,epsilon2] = perturbation_channel(gRho);
        
      
   
    
        zRho = 0;
        
        for jMapElm = 1 : numel(keyMap)
            zRho = zRho + keyMap{jMapElm} * gRho * keyMap{jMapElm};
        end
        
       [zRho,epsilon3] = perturbation_channel(zRho);
       realEpsilon = max([epsilon1,epsilon2,epsilon3]);
       fval = real(trace(gRho * (logm(gRho) - logm(zRho))));
    end

end