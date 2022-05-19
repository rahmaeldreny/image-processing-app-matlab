function varargout = motion(varargin)
% MOTION MATLAB code for motion.fig
%      MOTION, by itself, creates a new MOTION or raises the existing
%      singleton*.
%
%      H = MOTION returns the handle to a new MOTION or the handle to
%      the existing singleton*.
%
%      MOTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOTION.M with the given input arguments.
%
%      MOTION('Property','Value',...) creates a new MOTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before motion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to motion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help motion

% Last Modified by GUIDE v2.5 14-May-2022 19:06:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @motion_OpeningFcn, ...
                   'gui_OutputFcn',  @motion_OutputFcn, ...
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


% --- Executes just before motion is made visible.
function motion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to motion (see VARARGIN)

% Choose default command line output for motion
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes motion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = motion_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.JPG';'*.bmp';'*.tif'});
         global  sub;
         global fullFileName;
         sub = fullfile(folder, baseFileName);
            img=imread(sub);
img1=imread(fullFileName);
img2=imread(sub);
[rows cols matrciesNo]=size(img1);
[rows cols matrciesNo]=size(img2);
img1=double(img1);
img2=double(img2);
sub=img1-img2;
figure;
imshow(sub,[]);
