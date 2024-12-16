%% FUNCTION NAME: readResults
% This file helps to parse and plot mat results from the solver. 
% The N-dimensional data is expanded according to dimensions specified
% Only 1D,2D,3D data can be plotted.
%%
%the third optional argument is the plotting style
%available options:
%1.'linear': plot x and y linearly
%2.'linear-log': plot x versus log10(y)
%3.'km-log': plot -log10(x)*10/0.2 (x is assumed to be transmittance and converted to km) versus log(y)
%4.'none': do not plot
%% 函数名称：读取结果
% 该文件有助于解析和绘制来自求解器的矩阵结果。
% 根据指定的维度展开 N 维数据
% 只能绘制一维、二维和三维数据。
%%
% 第三个可选参数是绘图样式
% 可用选项：
%1.'linear'：线性绘制 x 和 y
%2.'linear-log'：绘制 x 与 log10(y)的关系图
%3.'km-log'：绘制 -log10(x)*10/0.2（假定 x 为透射率并转换为 km）与 log(y) 的关系图
%4.'none'：不绘制
function MDIplotResults(results,parameters,style)
    if(~strcmp(style,'none'))
        parameters_scan = struct2cell(parameters.scan);
        dimensions = getCellDimensions(parameters_scan);
        varnames = fieldnames(parameters.scan);
        dim = length(parameters_scan);
        rate1D = [results.lowerBound]; %can also plot upperBound or FWBound

        if(dim>1)
            rate = reshape(zeros(1,length(rate1D)),dimensions);
        end

        for i=1:length(rate1D)
            indices=num2cell(expandIndices(i, dimensions));
            rate(indices{:})=rate1D(i); %convert the cell array to comma-separated list for indexing将单元格数组转换为逗号分隔的列表以进行索引
        end

        rate(rate<0)=0; %truncate the R<0 data points

        handleF=figure(1);
        set(handleF,'Position',[800,600,800,600]);
        
        if(dim==1)
            %1-D plot
            X = parameters_scan{1};
            if(strcmp(style,'linear-log'))
                XAXISL=xlabel(varnames{1});
                YAXISL=ylabel('log_{10}(R)');
            elseif(strcmp(style,'km-log'))
                plot(-log10(X)*50.0*2,log10(rate))
                XAXISL=xlabel('L_{AB}(km)');
                YAXISL=ylabel('log_{10}(R)');
            elseif(strcmp(style,'km'))
                plot(-log10(X)*50.0*2,rate)
                XAXISL=xlabel('L_{AB}(km)');
                YAXISL=ylabel('R');
            else
                plot(X,rate)
                XAXISL=xlabel(varnames{1});
                YAXISL=ylabel('rate');
            end
            XAXISL.FontSize=16;
            YAXISL.FontSize=16;
        elseif(dim==2)
            X = parameters_scan{1};
            Y = parameters_scan{2};
            %2-D contour
            contour(cell2mat(parameters_scan(1)),cell2mat(parameters_scan(2)),rate)
            XAXISL=xlabel(varnames{1});
            YAXISL=ylabel(varnames{2});
            ZAXISL=zlabel('rate');
            XAXISL.FontSize=16;
            YAXISL.FontSize=16;
            ZAXISL.FontSize=16;
        elseif(dim==3)
            X = parameters_scan{1};
            Y = parameters_scan{2};
            Z = parameters_scan{3};
            %3-D contour
            isosurface(X,Y,Z,rate)
            XAXISL=xlabel(varnames{1});
            XAXISL=ylabel(varnames{2});
            XAXISL=zlabel(varnames{3});
            XAXISL.FontSize=16;
            YAXISL.FontSize=16;
            ZAXISL.FontSize=16;
        else
            fprintf('dimension too high for plotting\n')
        end    
    
    end
end

function indices=expandIndices(i, dimensions)
    i = i - 1; % convert to 0-based index
    dim = length(dimensions);
    indices=[];
    residue = i;
    for j=1:dim
        multiple = 1;
        for k=j+1:dim
            multiple = multiple * dimensions(k);
        end
        index = floor(residue/multiple);
        residue = rem(residue,multiple);
        indices=[indices,index];
    end
    indices = indices + 1; % convert to 1-based index
end

function dimensions = getCellDimensions(array)
    L = length(array);
    dimensions = [];
    if L~=0
        for k=1:L
           element = cell2mat(array(k));
           dimensions = [dimensions,length(element)];
        end
    else
        dimensions=[0];
    end
end