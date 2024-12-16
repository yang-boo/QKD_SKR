%helper function that looks for variables based on varNames,
%searching from provided list of names and values
% 辅助函数，根据 varNames 从提供的名称和值列表中查找变量
function varValues = findVariables(varNames,names,values)
    
    varValues = {};
    counter = 0;
    for i=1:length(varNames)
        for j=1:length(names)
            if(strcmp(varNames(i),names(j))) %判断两个字符串是否相等
                varValues = [varValues,values{j}];
                counter = counter + 1;
                break; %repeated entries will be ignored 重复输入将被忽略
            end
        end
    end
    
    if(counter~=length(varNames))
        fprintf('description/channel parameters not a subset of input parameter list!\n');
        varValues = {}; %error output
    end
end