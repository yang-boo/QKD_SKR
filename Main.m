%% FUNCTION NAME: Main
% Main entry point function for key rate calculation.
% Please find manual for detailed explanations on the functions and the software structure.
%
% The user can start the program by choosing a preset to run or modifying a custum preset.
% Each preset contains three functions 
% 'setDescription', 'setParameters', and 'setOptions', 
% where respectively the protocol/channel description files, protocol parameters, 
% and software/solver options are inputted.
% 密钥率计算的主要入口点功能。
% 有关功能和软件结构的详细说明，请参见手册。
%
%用户可以通过选择预置运行程序或修改预置来启动程序。
% 每个预置包含三个功能 
%'setDescription', 'setParameters', and 'setOptions'、 
%其中分别包含协议/通道描述文件、协议参数和软件/求解器程序选项、 
%%

format long
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%% Setting MATLAB Library Path 设置 MATLAB 库路径%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %automatically set path for each run 自动设置每次运行的路径
% %please modify installation directories accordingly 请相应修改安装目录
% %DO NOT AUTO SET PATH IF YOUR MATLAB HAS OTHER LIBRARIY DEPENDENCIES IN THE PATH! 如果您的 Matlab 路径中有其他依赖库，请不要自动设置路径！
% %(OR IT WILL CLEAR ALL OTHER DEPENDENCY PATHS)(否则会清除所有其他依赖路径）。
% restoredefaultpath;
% addpath(genpath('.')) %add current software directory and all subfolders
% %Windows path example
% addpath(genpath('C:\cvx2')) %cvx
% addpath(genpath('C:\Program Files\Mosek\9.0\toolbox\R2015a')) %external mosek
% %Mac/Linux path example
% addpath(genpath('~/cvx'))
% addpath(genpath('~/mosek'))

%%%%%%%%%%%%%%%%%%%%% Setting User Input 设置用户输入%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose a preset test case, or use custom setting (where user can fill in protocol/solver/parameter/options independently)
%选择预设测试例子，或使用自定义设置（用户可独立填写协议/求解器/参数/选项）
%available options: 可用选项： 
%1.'pmBB84_asymptotic'
%2.'pmBB84_finite'
%3.'pmBB84WCP_decoy'
%4.'MDIBB84_asymptotic'
%5.'MDIBB84_finite'
%6.'MDIBB84WCP_decoy' (can use parallel toolbox; decoy analysis might take a few minutes) (可使用并行工具箱；诱饵分析可能需要几分钟时间）
%7.'DMCVQKD_asymptotic'
%8.(archived) 'pmBB84Simple_asymptotic'
%9.(archived) 'MDIBB84Simple_asymptotic'
%10.(archived) 'MDIBB84Simple_finite'
%11. any custom setting (can be composed based upon templatePreset) 任何自定义设置（可根据templatePreset组成）

preset='pmBB84WCP_decoy_old_w';
% preset='MDIBB84WCP_decoy_eR_w';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

%%%%%%%%%%%%%%%%%%%%% Run Main Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call main iteration function
results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions);

%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BB84----save the results to file
path_file = fullfile('E:\yangbo\MATLAB Files\Output\','output_BB84_old_0_140_day0725.mat');
save(path_file,'results','parameters');

% % MDI----save the results to file
% path_file = fullfile('E:\yangbo\MATLAB Files\Output\','output_MDI_eR_0_100.mat');
% save(path_file,'results','parameters');

% %can also load a previous session's result to plot it
% %(can comment out main iteration above to skip computation)
% load('output.mat','results','parameters_scan')

%can uncomment this line to output debugging info
% results.debugInfo

%automatically parse and plot the results (optional)
%the third optional argument is the plotting style
%available options for 1D data:
%1.'linear': plot x and y linearly
%2.'linear-log': plot x versus log10(y)
%3.'km-log': plot -log10(x)*10/0.2 (x is assumed to be transmittance and converted to km) versus log(y)
%4.'km': plot km and y linearly
%5.'none': do not plot

% BB84
% plotResults(results,parameters,'none');

% MDI
% MDIplotResults(results,parameters,'none');

format long
clear all;
close all;
clc;

preset='pmBB84WCP_decoy_w';
% preset='MDIBB84WCP_decoy_w';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

%%%%%%%%%%%%%%%%%%%%% Run Main Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call main iteration function
results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions);

%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BB84----save the results to file
path_file = fullfile('E:\yangbo\MATLAB Files\Output\','output_BB84_0_140_day0823.mat');
save(path_file,'results','parameters');
