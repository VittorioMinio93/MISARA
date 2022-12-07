%%%%%%%%%%%%% Main window configuration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio

%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main_fig: Home figure 
% panel_color: colors of the panel 
% fig_color: colors of the figure

%% Set Colors and load icons
panel_color= [.97 .97 .97];
fig_color= [.94 .94 .94];
load icons.mat

% Set the dimensions of  the main window   
ps= get(mainfig,'position');
set(mainfig,'position',[.5-1.5*ps(3)/2 .5-1.5*ps(4)/2 1.5*ps(3) 1.5*ps(4)]);

% Add a separation line
uibuttongroup(mainfig,'visible','on','BackgroundColor',panel_color,...
    'units','normalized','Position',[0 .11 1 .01]);

% Add the 11 Modules buttons
% Beam pattern, Signal Viewer, Spectrogram, Spectral coherence,
% Polarization, STA-LTA, Salped, ZLCC, MUSIC, Semblance, Radial Semblance
uicontrol(mainfig,'style','pushbutton', 'string','Beam Pattern','units','normalized',...
                    'position',[.03 .02 .15 .07], 'callback','arraymap');
uicontrol(mainfig,'style','pushbutton', 'string','Spectrogram','units','normalized',...
                    'position',[.19 .02 .15 .07], 'callback','spec_analysis');                                
uicontrol(mainfig,'style','pushbutton', 'string','Spectral coherence','units','normalized',...
                    'position',[.35 .02 .15 .07], 'callback','spectral_coh'); 
uicontrol(mainfig,'style','pushbutton', 'string','Polarization','units','normalized',...
                    'position',[.51 .02 .15 .07], 'callback','polariz');
uicontrol(mainfig,'style','pushbutton', 'string','STA/LTA','units','normalized',...
                    'position',[.67 .02 .15 .07], 'callback','stalta');
uicontrol(mainfig,'style','pushbutton', 'string','Salped','units','normalized',...
                    'position',[.83 .02 .15 .07], 'callback','salpedDetection');
uicontrol(mainfig,'style','pushbutton', 'string','ZLCC','units','normalized',...
                    'position',[.03 .55 .15 .07], 'callback','zlc');
uicontrol(mainfig,'style','pushbutton', 'string','Semblance','units','normalized',...
                    'position',[.03 .35 .15 .07], 'callback','SemLoc');
uicontrol(mainfig,'style','pushbutton', 'string','Radial Semblance','units','normalized',...
                    'position',[.03 .25 .15 .07], 'callback','RadSemLoc');
uicontrol(mainfig,'style','pushbutton', 'string','MUSIC','units','normalized',...
                    'position',[.03 .45 .15 .07], 'callback','musicArray');
uicontrol(mainfig,'style','pushbutton', 'string','Signal Viewer','units','normalized',...
                    'position',[.03 .15 .15 .07], 'callback','instr_corr');               