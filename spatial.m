
function varargout = spatial(varargin)

% SPATIAL MATLAB code for spatial.fig
%      SPATIAL, by itself, creates a new SPATIAL or raises the existing
%      singleton*.
%
%      H = SPATIAL returns the handle to a new SPATIAL or the handle to
%      the existing singleton*.
%
%      SPATIAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPATIAL.M with the given input arguments.
%
%      SPATIAL('Property','Value',...) creates a new SPATIAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spatial_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spatial_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spatial

% Last Modified by GUIDE v2.5 18-May-2022 20:20:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spatial_OpeningFcn, ...
                   'gui_OutputFcn',  @spatial_OutputFcn, ...
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


% --- Executes just before spatial is made visible.
function spatial_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spatial (see VARARGIN)

% Choose default command line output for spatial
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spatial wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = spatial_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%histo
img=getimage(handles.axes2)

[rows cols ]=size(img);
counts=zeros(1,256);
for i=1:rows
 for   j=1:cols
    greyLevel=img(i,j);
    counts(greyLevel+1)=counts(greyLevel+1)+1;
end
end
figure;
greyLevels=0:255;
bar(greyLevels,counts,'barWidth',1,'faceColor','b');
xlabel('gray level','fontsize',20);
ylabel('pixel count','fontsize',20);
title('histogram','fontsize',20);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%negative       

img=getimage(handles.axes2);

[rows cols matreciesNo]=size(img);
%g_img=rgb2gray(img);
l=2^8;
neg=(l-1)-img;
axes(handles.axes2);
imshow(neg);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%log

c=inputdlg('value of c=');
c=str2num(c{1});
img=getimage(handles.axes2);
[rows cols matreciesNo]=size(img);
% g_img=rgb2gray(img);
doubleImg=im2double(img);
s=(c*log(1+doubleImg))*256;
sl=uint8(s);
axes(handles.axes2);
imshow(sl);



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% powerLow
img=getimage(handles.axes2);
[rows cols matreciesNo]=size(img);
c=inputdlg('value of c=');
c=str2double(c);
y=inputdlg('value of y=');
y=str2double(y);


doubleImg=double(img);
new=c(1)*power(doubleImg,y(1));
new=uint8(new);
disp(c);
disp(y);
axes(handles.axes2);
imshow(new);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%contrast
img=getimage(handles.axes2);
%gray_img=rgb2gray(img);
s_min=inputdlg('value of S_min=');
s_min=str2num(s_min{1});
s_max=inputdlg('value of S_max=');
s_max=str2num(s_max{1});
[rows cols matreciesNo]=size(img);

doubleImg=double(img);

r_min=min(doubleImg(:));
r_max=max(doubleImg(:));
s=((s_max-s_min)/(r_max-r_min))*(doubleImg-r_min)+s_min;

sI=uint8(s);
axes(handles.axes2);
imshow(sI);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%threshold
img=getimage(handles.axes2);
%gray_img=rgb2gray(img);
k=inputdlg('value of k=');
k=str2num(k{1});
[rows cols matreciesNo]=size(img);
doubleImg=double(img);
for i=1:rows
    for j=1:cols
        if doubleImg(i,j)>=k
            doubleImg(i,j)=255;
        else
            doubleImg(i,j)=0;
        end
    end
end
s=uint8(doubleImg);
axes(handles.axes2);
imshow(s);

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%first approach
img=getimage(handles.axes2);
 %gray_img=rgb2gray(img);
min=inputdlg('value of min range=');
min=str2num(min{1});
max=inputdlg('value of max range=');
max=str2num(max{1});
[rows cols matreciesNo]=size(img);

doubleImg=double(img);
%desired range in white
for i=1:rows
    for j=1:cols
        if doubleImg(i,j)>=min&&doubleImg(i,j)<=max
            doubleImg(i,j)=255;
        else
            doubleImg(i,j)=0;
        end
    end
end
s=uint8(doubleImg);
axes(handles.axes2);

imshow(s);

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%bitslicing
img=getimage(handles.axes2);
%gray_img=rgb2gray(img);
[rows cols matreciesNo]=size(img);
newImg=zeros(rows,cols,0); 
figure;
for k=1:8
    for i=1:1:rows
        for j=1:1:cols
       newImg(i,j,k)=bitget(img(i,j),k);
        end
    end

subplot(2,5,k+1);
imshow(newImg(:,:,k));
title ('img for bit number ');
end



% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%add
img=getimage(handles.axes2);
s=inputdlg('value of scalinig=');
s=str2num(s{1});
doubleImg=double(img);
add=img+s;
axes(handles.axes2);
imshow(add);



% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%motion
open('motion.fig');


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% open('scalew.fig');
img=getimage(handles.axes2);
s=inputdlg('value of scalinig=');
s=str2num(s{1});
doubleImg=double(img);
s_img=img*s;
s_img=uint8(s_img);
axes(handles.axes2);
imshow(s_img);





% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [baseFileName, folder] = uigetfile({'*.jpg';'*.png';'*.JPG';'*.bmp';'*.tif'});
         global  fullFileName;
         fullFileName = fullfile(folder, baseFileName);
            img=imread(fullFileName);
            [rows cols matreciesNo]=size(img);
           
axes(handles.axes2);
    imshow(img);


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%2*nd approach of grey slicing
img=getimage(handles.axes2);
% gray_img=rgb2gray(img);
min=inputdlg('value of min range=');
min=str2num(min{1});
max=inputdlg('value of max range=');
max=str2num(max{1});
[rows cols matreciesNo]=size(img);

doubleImg=double(img);
%desired range in white
for i=1:rows
    for j=1:cols
        if doubleImg(i,j)>=min&&doubleImg(i,j)<=max
             doubleImg(i,j)=255;
        else
            doubleImg(i,j);
        end
    end
end
s=uint8(doubleImg);
axes(handles.axes2);
imshow(s);



% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
%histo aquization
GIm=getimage(handles.axes2);
%GIm = rgb2gray(img);
numofpixels=size(GIm,1)*size(GIm,2);

HIm=uint8(zeros(size(GIm,1),size(GIm,2)));

freq=zeros(256,1);

probc=zeros(256,1);

output=zeros(256,1);

for i=1:size(GIm,1)

    for j=1:size(GIm,2)

        value=GIm(i,j);

        freq(value+1)=freq(value+1)+1;

    end

end


sum=0;
no_bins=255;

for i=1:size(freq)

   sum=sum+freq(i);

   probc(i)=sum/numofpixels;

   output(i)=round(probc(i)*no_bins);

end

for i=1:size(GIm,1)

    for j=1:size(GIm,2)

            HIm(i,j)=output(GIm(i,j)+1);

    end

end
axes(handles.axes2);

imshow(HIm);


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%subtraction
img=getimage(handles.axes2);
s=inputdlg('value of scalinig=');
s=str2num(s{1});
doubleImg=double(img);
sub=img-s;
axes(handles.axes2);
imshow(sub);


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%inverselog
% global fullFileName
c=inputdlg('value of c=');
c=str2num(c{1});
img=getimage(handles.axes2)
[rows cols matreciesNo]=size(img);
% g_img=rgb2gray(img);
doubleImg=im2double(img);
s=power(2,c*log(doubleImg))*256;
sl=uint8(s);
axes(handles.axes2);
    imshow(sl);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
filter=get(handles.listbox1,'value');
switch(filter)
    case 2
        %avgfilter
        img=getimage(handles.axes2);
m=inputdlg('value of m=');
m=str2num(m{1});
n=inputdlg('value of n=');
n=str2num(n{1});
%gimg=rgb2gray(img);
gimg=double(img);
[rows,cols]=size(gimg);
 filter_coefficient=1/(m*n);
 filter=zeros(m,n);
 for i=1:m
     for j=1:n
        filter(i,j)=filter_coefficient;
     end 
 end
    a=(m-1)/2;
    b=(n-1)/2;
    outimage=gimg;
    for i=a+1:rows-a
     for j=b+1:cols-b
      current_pixel_with_neighbour=gimg(i-a:i+a,j-b:j+b);
      new_pixel_with_neighbour= current_pixel_with_neighbour.*filter;
      newpixel=sum(new_pixel_with_neighbour(:),'double');
      outimage(i,j)=newpixel;
     end 
    end
    outimage=uint8(outimage);
  axes(handles.axes2);
imshow(outimage);
%%%%%%%%%%%weighted
    case 3
        img=getimage(handles.axes2);
img=double(img);
[rows,cols]=size(img);
mask=[1/16 2/16 1/16;2/16 4/16 2/16;1/16 2/16 1/16];
out=img;
for i=2:rows-1
    for j=2:cols-1
        temp=mask.*img(i-1:i+1,j-1:j+1);
        value=sum(temp(:));
        out(i,j)=value;
    end
end
out=uint8(out);
  axes(handles.axes2);
imshow(out);
%%%%%%%%%%%median
    case 4
        img=getimage(handles.axes2);
%gray_img=rgb2gray(gray_image);
[rows,cols]=size(img);
out=img;
for i=2:rows-1
 for j=2:cols-1
     temp = [img(i-1, j-1) img(i-1, j) img(i-1, j + 1) img(i, j-1) img(i, j) img(i, j + 1) img(i + 1, j-1) img(i + 1, j) img(i + 1, j + 1)];
     temp = sort(temp);
     out(i, j)= temp(5);
end
end
  axes(handles.axes2);
imshow(out);
%%%%%%%%%%%1'st
    case 5
        
img=getimage(handles.axes2);
[rows ,cols ] = size(img);
% gray_img = rgb2gray(img);
img = double(img);
filtered_image = zeros(size(img));
  
% Sobel Operator Mask
Mx = [-1 0 1; -2 0 2; -1 0 1];
My = [-1 -2 -1; 0 0 0; 1 2 1];
for i = 1:size(img, 1) - 2
    for j = 1:size(img, 2) - 2
  
        Gx = sum(sum(Mx.*img(i:i+2, j:j+2)));
        Gy = sum(sum(My.*img(i:i+2, j:j+2)));
                 
        filtered_image(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
    end
end
  axes(handles.axes2);
imshow(filtered_image);
%%%%%%%%%%%%laplace
    case 6
        
img=getimage(handles.axes2);
%img=rgb2gray(img);
img = double(img);
[rows,cols]=size(img);
mask = [0,1,0;1,-4,1;0,1,0];
out = img;
for i=2:rows-1
     for j=2:cols-1
     temp = mask.*img(i-1:i+1,j-1:j+1);
     value = sum(temp(:));
     out(i, j)= value;
    end
end
out = uint8(out);
  axes(handles.axes2);
imshow(out);
%%%%%%%%%%%%%%%%%%%%%compsite laplace
    case 7
        img=getimage(handles.axes2);
img = double(img);
[rows,cols]=size(img);
mask = [0,1,0;1,-5,1;0,1,0];
out = img;
for i=2:rows-1
     for j=2:cols-1
     temp = mask.*img(i-1:i+1,j-1:j+1);
     value = sum(temp(:));
     out(i, j)= value;
    end
end
out = uint8(out);
  axes(handles.axes2);
imshow(out);
end 
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%reset
global  fullFileName;
img=imread(fullFileName);
 [rows cols matreciesNo]=size(img);
         
  axes(handles.axes2);
imshow(img);


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%size up
img=getimage(handles.axes2);
[rows cols matrixno] = size(img);
disp(rows);
disp(cols);
f=inputdlg('value of sampling factor=');
f=str2num(f{1});
     outimg(:,:,matrixno)=resamplingUp(img(:,:,matrixno),f);

axes(handles.axes2);
    imshow(outimg);
    [r,c]=size(outimg);
    disp(r);
    disp(c);

% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%size down
img=getimage(handles.axes2);
[rows cols matricesNo] = size(img);
disp(rows);
disp(cols);
f=inputdlg('value of sampling factor=');
f=str2num(f{1});

    outimg(:,:,matricesNo) = subSampling(img(:,:,matricesNo),f);


axes(handles.axes2);
    imshow(outimg);
    [r,c]=size(outimg);
    disp(r);
    disp(c);

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
value=get(handles.listbox2,'value');
switch(value)
    case 2
        %ilpf

img=getimage(handles.axes2);
[m,n]=size(img);
ft_img=fft2(double(img));
d0=15;

u=0:(m-1);
idx=find(u>m/2);
u(idx)=u(idx)-m;
v=0:(n-1);
idy=find(v>n/2);
v(idy)=v(idy)-n;

[V,U]=meshgrid(v,u);
d=sqrt(U.^2+V.^2);
h=double(d<=d0);
g=h.*ft_img;
out_img=real(ifft2(double(g)));
axes(handles.axes2);
 imshow(out_img,[]);
    case 3
        %%%%%%%%%%blpf
 img=getimage(handles.axes2);
[M, N] = size(img);
FT_img = fft2(double(img));
D0 = 15; 
n=2*2;
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
  

D = sqrt(U.^2+V.^2);

D = D./ D0;

H = 1./((1+D).^n);
  

G = H.*FT_img;
out_img = real(ifft2(double(G)));

axes(handles.axes2);
 imshow(out_img,[]);
    case 4
        %%%%%%%%%glpf
img=getimage(handles.axes2);
[M,N] = size(img);
FT_img = fft2 (double(img));
D0=15;
D0 = (D0^2)*2; % designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V,U] = meshgrid(v,u);
D = sqrt(U.^2+V.^2);
D = -D.^2;
H = exp(D/D0);
G = H.*FT_img;
out_img = real(ifft2(double(G)));
axes(handles.axes2);
 imshow(out_img,[]);
    case 5
 %ihpf
 
img=getimage(handles.axes2);
[M,N] = size(img);
FT_img = fft2 (double(img));
D0=10;
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
  D = sqrt(U.^2+V.^2);
H = double(D > D0);
G = H.*FT_img;
out_img = real(ifft2(double(G)));
axes(handles.axes2);
 imshow(out_img,[]);
    case 6
        %bhpf

img=getimage(handles.axes2);
[M,N] = size(img);
FT_img = fft2 (double(img));
n = 2;
D0=10;
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;
[V, U] = meshgrid(v, u);
  D = sqrt(U.^2+V.^2);
H = 1./(1 + (D0./D).^(2*n));G = H.*FT_img;
out_img = real(ifft2(double(G)));
axes(handles.axes2);
 imshow(out_img,[]);
 case 7
 %ghpf
 
img=getimage(handles.axes2);
[m n]=size(img);
f_transform=fft2(img);
f_shift=fftshift(f_transform);
p=m/2;
q=n/2;
d0=70;
for i=1:m
for j=1:n
distance=sqrt((i-p)^2+(j-q)^2);
low_filter(i,j)=1-exp(-(distance)^2/(2*(d0^2)));
end
end
filter_apply=f_shift.*low_filter;
image_orignal=ifftshift(filter_apply);
image_filter_apply=abs(ifft2(image_orignal));
axes(handles.axes2);
 imshow(image_filter_apply,[]);
end
% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            img=getimage(handles.axes2);
            [rows cols matreciesNo]=size(img);
            if matreciesNo>1
            img=rgb2gray(img);
            end
axes(handles.axes2);
    imshow(img);


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
value=get(handles.listbox3,'value');
global fullFileName;
switch(value)
    case 2
        %and-one
        
        img=getimage(handles.axes2);
        
        s=inputdlg('value of start point=');
         s=str2num(s{1});
            e=inputdlg('value of end point=');
            e=str2num(e{1});
         outimg = logic(img,s,e);
        axes(handles.axes2);
        imshow(outimg);
         case 3
        %and-zero
        img=getimage(handles.axes2);
           s=inputdlg('value of start point=');
         s=str2num(s{1});
            e=inputdlg('value of end point=');
            e=str2num(e{1});
         outimg = and_zero(img,s,e);
        axes(handles.axes2);
        imshow(outimg);
         case 4
        %or-one
        img=getimage(handles.axes2);
           s=inputdlg('value of start point=');
         s=str2num(s{1});
            e=inputdlg('value of end point=');
            e=str2num(e{1});
        outimg=or_one(img,s,e);
        axes(handles.axes2);
        imshow(outimg);
        case 5
        %or-zero
        img=getimage(handles.axes2);
           s=inputdlg('value of start point=');
         s=str2num(s{1});
            e=inputdlg('value of end point=');
            e=str2num(e{1});
        outimg=or_zero(img,s,e);
        axes(handles.axes2);
        imshow(outimg);
end


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
