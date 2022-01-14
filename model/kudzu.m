function varargout = kudzu(varargin)
% KUDZU MATLAB code for kudzu.fig
%      KUDZU, by itself, creates a new KUDZU or raises the existing
%      singleton*.
%
%      H = KUDZU returns the handle to a new KUDZU or the handle to
%      the existing singleton*.
%
%      KUDZU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KUDZU.M with the given input arguments.
%
%      KUDZU('Property','Value',...) creates a new KUDZU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kudzu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kudzu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kudzu

% Last Modified by GUIDE v2.5 21-Jan-2011 09:47:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kudzu_OpeningFcn, ...
                   'gui_OutputFcn',  @kudzu_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before kudzu is made visible.
function kudzu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kudzu (see VARARGIN)

% Choose default command line output for kudzu
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kudzu wait for user response (see UIRESUME)
% uiwait(handles.figure);

set(handles.figure, 'UserData', struct());
InitializeParameters(handles);
InitializeGrid(handles);

% --- Outputs from this function are returned to the command line.
function varargout = kudzu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function growMethod_CreateFcn(hObject, eventdata, handles)
    fid = fopen('kudzuGrow.m');
    while(isempty(strfind(fgetl(fid), 'growAB')))
    end
    M = [];
    V = 0;
    while(feof(fid)==0)
        L = fgetl(fid);
        if(strfind(L, 'case'))
            strfind(L, '%')
            M = [M ; cellstr(L((strfind(L,'%')+1):end))];
            V = V +1;
        end
    end
    
    set(hObject, 'String', cellstr(M));
    set(hObject, 'Value', 1);


function controlMeethod_CreateFcn(hObject, eventdata, handles)



% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
    InitializeGrid(handles);

% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
    if(get(handles.doProfile, 'Value'))
        profile on;
    end

    InitializeParameters(handles);
    Grow(handles);

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
    if(get(handles.doProfile, 'Value'))
        profile viewer;
    end
    
    set(hObject, 'UserData',  1);
    guidata(hObject, handles);


function gridSize_Callback(hObject, eventdata, handles)
    InitializeGrid(handles);
    
function colonyCount_Callback(hObject, eventdata, handles)
    InitializeGrid(handles);

function colonySize_Callback(hObject, eventdata, handles)
    InitializeGrid(handles);

function conversionRate_Callback(hObject, eventdata, handles)
function leafDeathRate_Callback(hObject, eventdata, handles)
function rootDeathRate_Callback(hObject, eventdata, handles)
function photoRate_Callback(hObject, eventdata, handles)
function RKSteps_Callback(hObject, eventdata, handles)
function updateDelta_Callback(hObject, eventdata, handles)
function delayTime_Callback(hObject, eventdata, handles)
function doProfile_Callback(hObject, eventdata, handles)
function seasonLength_Callback(hObject, eventdata, handles)
function growMethod_Callback(hObject, eventdata, handles)
function controlMeethod_Callback(hObject, eventdata, handles)


function delayTime_CreateFcn(hObject, eventdata, handles)
function updateDelta_CreateFcn(hObject, eventdata, handles)
function RKSteps_CreateFcn(hObject, eventdata, handles)
function conversionRate_CreateFcn(hObject, eventdata, handles)
function leafDeathRate_CreateFcn(hObject, eventdata, handles)
function rootDeathRate_CreateFcn(hObject, eventdata, handles)
function photoRate_CreateFcn(hObject, eventdata, handles)
function seasonLength_CreateFcn(hObject, eventdata, handles)
    
function colonyCount_CreateFcn(hObject, eventdata, handles)
function colonySize_CreateFcn(hObject, eventdata, handles)
function gridSize_CreateFcn(hObject, eventdata, handles)


function InitializeParameters(handles)
    F = get(handles.figure, 'UserData');
   
    F.growth.p   = eval(get(handles.photoRate, 'String'));
    F.growth.m   = eval(get(handles.conversionRate, 'String'));
    F.growth.muA = eval(get(handles.leafDeathRate, 'String'));
    F.growth.muB = eval(get(handles.rootDeathRate, 'String'));
    F.growth.RK_steps = floor(eval(get(handles.RKSteps, 'String')));
    F.growth.updateDelta = eval(get(handles.updateDelta, 'String'));
    F.growth.seasonLength = eval(get(handles.seasonLength, 'String'));
    
    F.growth.growMethod = get(handles.growMethod, 'Value');
    
    set(handles.figure, 'UserData', F);
     


function InitializeGrid(handles)
    N = floor(eval(get(handles.gridSize, 'String')));
    colonyCount = floor(eval(get(handles.colonyCount, 'String')));
    colonySize = floor(eval(get(handles.colonySize, 'String')));

    %Build Initial colonies
    handles.grid = image(zeros(N), 'parent', handles.gridAxis);
    guidata(handles.figure, handles);
    set(handles.gridAxis, 'XTick', [])
    set(handles.gridAxis, 'YTick', [])

    A = zeros(N);
    c = ceil(rand(1, colonyCount) * N*N);
    A(c) = .1;
    A = conv2(A, ones(colonySize), 'same');
    B = A;
    
    %Update/Store data in the Grid control
    set(handles.grid, 'UserData', struct('A', A, 'B', B));
    InitializeGridImage(handles);
        
    
function InitializeGridImage(handles)
    %Clear Grid
    G = get(handles.grid, 'UserData');
    G.cLevels = 8;
    set(handles.grid, 'UserData', G);
    
    %Clear History
    F = get(handles.figure, 'UserData');
    F.history = [];
    set(handles.figure, 'UserData', F); 
    
    UpdateFigure(handles);
    
    %Build Palette
    c = G.cLevels;
    x = [1 linspace(.7,  0, c-1); 1 linspace(.8, .6, c-1); 1 linspace(.5,  0, c-1)]'; %green above ground
    y = [1 linspace(.8, .6, c-1); 1 linspace(.5, .4, c-1); 1 linspace(.5,  0, c-1)]'; %brown below ground
    P = x(reshape(repmat(1:c, c, 1), c*c, 1), :) ...
      + y(reshape(repmat(1:c, 1, c), c*c, 1), :);
    colormap(P/2);
    

function UpdateFigure(handles)
    %Build Composite Image
    G = get(handles.grid, 'UserData');
    R = (ceil(G.A * G.cLevels) * G.cLevels) + ceil(G.B * G.cLevels);
    set(handles.grid, 'CData', R);

    %Build the history plot
    F = get(handles.figure, 'UserData');
    if(~isfield(F, 'history'))
        F.history = [];
    end
    F.history = [F.history; sum(sum(G.A)), sum(sum(G.B))];
    set(handles.figure, 'UserData', F);
    s = size(F.history);
    T = (0:s(1)-1) * F.growth.updateDelta;
    plot(T,F.history, 'Parent', handles.history);
    drawnow;
    
        
function Grow(handles)
    G = get(handles.grid, 'UserData');
    F = get(handles.figure, 'UserData');

    delay = str2double(get(handles.delayTime, 'String'));
    
    set(handles.stop, 'UserData', 0)
    
    for t=0:F.growth.updateDelta:(10-F.growth.updateDelta)
        start = tic;
        if(1 == get(handles.stop, 'UserData'))
            break
        end
        [G.A G.B] = kudzuGrow(G.A, G.B, t, F.growth);
        set(handles.grid, 'UserData', G);
        UpdateFigure(handles);
        
        while(toc(start) < delay)
        end
        
    end
    
    set(handles.stop, 'UserData', 0)
