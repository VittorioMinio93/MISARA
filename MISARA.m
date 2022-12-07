% Copyright (c) 2022, Vittorio Minio 
% All rights reserved.
%
% 
% 1) TERMS OF USE
%   If you use MISARA or any function(s) of it, you need to acknowledge 
%   MISARA by citing the following DOI:
%
%   10.5281/zenodo.7410076
%  
% 
% 2) LICENSE:
%   This program is part of MISARA
%   MISARA is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%      MISARA is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%       You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%% Home function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main module of the package to configure and manage datasets, parameters
% and modules.
% December 2022 
% Vittorio Minio

%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params: user parameters
% h: UI graphics parameters
% h2: UI control parameters
% main_fig: Home figure 
%%%%%%%%%%%%h.variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% menu_main: menu panel in the Home figure   
% menu (*): menu section in the Home figure

%% Add folders to search path
clc;clear;close all
addpath(genpath(['.', filesep, 'parameters']))
addpath(genpath(['.', filesep, 'main_fig']))
addpath(genpath(['.', filesep, 'predefind_mats']))
addpath(genpath(['.', filesep, 'Array_map']))
addpath(genpath(['.', filesep, 'Spectrogram']))
addpath(genpath(['.', filesep, 'Spectral_coherence']))
addpath(genpath(['.', filesep, 'Polarization']))
addpath(genpath(['.', filesep, 'STA_LTA']))
addpath(genpath(['.', filesep, 'Salped']))
addpath(genpath(['.', filesep, 'ZLC']))
addpath(genpath(['.', filesep, 'Semblance']))
addpath(genpath(['.', filesep, 'Radial Semblance']))
addpath(genpath(['.', filesep, 'Music']))
addpath(genpath(['.', filesep, 'func']))
addpath(genpath(['.', filesep, 'Instrumental correction']))
addpath(genpath(['.', filesep, 'CreateData']))
addpath(genpath(['.', filesep, 'Doc']))
addpath(genpath(['.', filesep, 'Doc/Video Tutorial']))
%% Declaration global variables
global params h2
%% Control and load user/default parameters in the Home module
if exist('user_values.mat','file')==2
    load user_values.mat
else
    load default_values.mat    
end

params.version= 'MISARA 1.0.1';%Version of MISARA
%% Home figure
mainfig=figure('name',[params.version, ' - Home'],...
        'Menubar','none', 'NumberTitle','off','units','normalized');%%Create the Home figure  
%%Create the Help menu 
m = uimenu(mainfig,'Text','&Help'); 
%%Create the User manual menu 
mitem1 = uimenu(m,'Text','&User Manual');
mitem1.MenuSelectedFcn = @UserManual;
%%Create the Video tutorial menu 
mitem2 = uimenu(m,'Text','&Video Tutorial');
mitem_a = uimenu(mitem2,'Text','&Case study 1');
mitem_b = uimenu(mitem2,'Text','&Case study 2');
mitem_c = uimenu(mitem2,'Text','&Case study 3');
%%Video tutorial for Case study 1
g(1) = uimenu(mitem_a,'Text','&Formatting data');
g(2) = uimenu(mitem_a,'Text','&Setting parameters');
g(3) = uimenu(mitem_a,'Text','&Array response function and geometry');
g(4) = uimenu(mitem_a,'Text','&Volcanic tremor features');
g(5) = uimenu(mitem_a,'Text','&Volcanic tremor source position');
g(1).MenuSelectedFcn = @VideoTutorial_a;
g(2).MenuSelectedFcn = @VideoTutorial_a;
g(3).MenuSelectedFcn = @VideoTutorial_a;
g(4).MenuSelectedFcn = @VideoTutorial_a;
g(5).MenuSelectedFcn = @VideoTutorial_a;

%%Video tutorial for Case study 2
f(1) = uimenu(mitem_b,'Text','&Seismic signal');
f(2) = uimenu(mitem_b,'Text','&VLP detection');
f(3) = uimenu(mitem_b,'Text','&VLP source location');
f(4)= uimenu(mitem_b,'Text','&LP detection');
f(5) = uimenu(mitem_b,'Text','&LP source location');
f(1).MenuSelectedFcn = @VideoTutorial_b;
f(2).MenuSelectedFcn = @VideoTutorial_b;
f(3).MenuSelectedFcn = @VideoTutorial_b;
f(4).MenuSelectedFcn = @VideoTutorial_b;
f(5).MenuSelectedFcn = @VideoTutorial_b;

%%Video tutorial for Case study 3
k(1) = uimenu(mitem_c,'Text','&Infrasound signal features');
k(2) = uimenu(mitem_c,'Text','&Infrasound signal source position');
k(1).MenuSelectedFcn = @VideoTutorial_c;
k(2).MenuSelectedFcn = @VideoTutorial_c;



config_main;%%Return the the modules location in the Home figure
config_general;%%Return the General Management figure
config_map;%%Return the Instrumental figure
config_analysis;%%Return the Preliminary analysis figure
config_array_analysis;%%Return the Array analysis figure
%% Menu panel
h.menu_main = uibuttongroup(mainfig,'visible','on','units','normalized','Position',[.015 .65 .2 .30],...
    'BackgroundColor','w','HighlightColor',[1 1 1]*.3,'BorderWidth',1,'BorderType','beveledin');%%Create the menu panel in the Home figure 


h.menu(1) = uicontrol( h.menu_main, 'Style','Radio','String','General Management',...
    'BackgroundColor','w', 'units','normalized','position',[.06 .8 .9 .2],...
    'HandleVisibility','off','Userdata',h.panel(1));%%Create the General Management section of the menu panel

h.menu(2) = uicontrol( h.menu_main, 'Style','Radio','String','Data Pre-processing ',...
    'BackgroundColor','w', 'units','normalized','position',[.06 .65 .9 .2],...
    'HandleVisibility','off', 'Userdata',h.panel(2));%%Create the Instrumental section of the menu panel

h.menu(3) = uicontrol( h.menu_main, 'Style','Radio','String','Signal features analysis',...
    'BackgroundColor','w', 'units','normalized','position',[.06 .5 .9 .2],'HandleVisibility','off',...
    'Userdata',h.panel(3));%%Create the Preliminary analysis section of the menu panel

h.menu(4) = uicontrol( h.menu_main, 'Style','Radio','String','Array analysis',...
    'BackgroundColor','w', 'units','normalized','position',[.06 .35 .9 .2],'HandleVisibility','off',...
    'Userdata',h.panel(4));%%Create the Array analysis section of the menu panel

set(h.menu_main,'selectedobject',h.menu(1));%%Set the starting section of the menu panel

set(h.menu_main,'SelectionChangeFcn',@menu_cbk);%%Change the section of the menu panel
%% Create Dataset button
uicontrol(mainfig,'style','pushbutton', 'string','Create Dataset','units','normalized',...
    'BackgroundColor','w','position',[.03 .68 .15 .07], 'callback','CreateDataStruct');

%% Nested functions
%%Open file pdf of MISARA User Manual 
function UserManual(~,~)
open('Doc\MISARA 1.0.1 manual.pdf')
end

%%%%%%%%%%%%%%%Open videos about Case study 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VideoTutorial_a(g,~)
switch g.Text
    case '&Setting parameters'
        file_video=[cd '/Doc/Video_tutorial/Case_study1_setting.mp4'];
        system(file_video)
    case '&Formatting data'
        file_video=[cd '/Doc/Video_tutorial/Case_study1_formatting.mp4'];
        system(file_video)
    case '&Array response function and geometry'
        file_video=[cd '/Doc/Video_tutorial/Case_study1_arrayresponse.mp4'];
        system(file_video)
    case '&Volcanic tremor features'
        file_video=[cd '/Doc/Video_tutorial/Case_study1_VolcaniTremorFeatures.mp4'];
        system(file_video)
    case '&Volcanic tremor source position'
        file_video=[cd '/Doc/Video_tutorial/Case_study1_VolcaniTremorSource.mp4'];
        system(file_video)
  

end
end
%%%%%%%%%%%%%%%%%Open videos about Case study 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VideoTutorial_b(f,~)
switch f.Text
    case '&Seismic signal'
        file_video=[cd '/Doc/Video_tutorial/Case_study2_SignalViewer.mp4'];
        system(file_video)
    case '&LP detection'
        file_video=[cd '/Doc/Video_tutorial/Case_study2_LPDetection.mp4'];
        system(file_video)
    case '&LP source location'
        file_video=[cd '/Doc/Video_tutorial/Case_study2_LPSource.mp4'];
        system(file_video)
    case '&VLP detection'
        file_video=[cd '/Doc/Video_tutorial/Case_study2_VLPDetection.mp4'];
        system(file_video)
    case '&VLP source location'
        file_video=[cd '/Doc/Video_tutorial/Case_study2_VLPSource.mp4'];
        system(file_video)
end  
end

%%%%%%%%%%%%%%%%%%%5Open videos about Case study 3%%%%%%%%%%%%%%%%%%%%%%%%%
function VideoTutorial_c(k,~)
switch k.Text
    case '&Infrasound signal features'
        file_video=[cd '/Doc/Video_tutorial/Case_study3_InfrasoundSignalFeatures.mp4'];
        system(file_video)
    case '&Infrasound signal source position'
        file_video=[cd '/Doc/Video_tutorial/Case_study3_InfrasoundSourceSource.mp4'];
        system(file_video)
end
end