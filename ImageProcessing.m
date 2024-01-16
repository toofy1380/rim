function varargout = ImageProcessing(varargin)
% IMAGEPROCESSING MATLAB code for ImageProcessing.fig
%      IMAGEPROCESSING, by itself, creates a new IMAGEPROCESSING or raises the existing
%      singleton*.
%
%      H = IMAGEPROCESSING returns the handle to a new IMAGEPROCESSING or the handle to
%      the existing singleton*.
%
%      IMAGEPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEPROCESSING.M with the given input arguments.
%
%      IMAGEPROCESSING('Property','Value',...) creates a new IMAGEPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageProcessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageProcessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageProcessing

% Last Modified by GUIDE v2.5 22-Aug-2023 12:24:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageProcessing_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageProcessing_OutputFcn, ...
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


% --- Executes just before ImageProcessing is made visible.
function ImageProcessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageProcessing (see VARARGIN)

% Choose default command line output for ImageProcessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageProcessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImageProcessing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function Browse_Callback(hObject, eventdata, handles)
% hObject    handle to Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show a file dialog box to allow the user to select an image
[filename, pathname] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.tiff', 'Image Files (*.jpg, *.jpeg, *.png, *.bmp, *.tif, *.tiff)'}, 'Select an image');

% If the user selected an image, update the GUI with the selected image
if isequal(filename,0) || isequal(pathname,0)
    % User cancelled the dialog box
else
    % Load the selected image
    image_path = fullfile(pathname, filename);
    handles.original_image = imread(image_path);
    handles.current_image = handles.original_image; % use original image as current image
    imshow(handles.current_image, 'Parent', handles.axes1);
end
%%%%%%%%%%%
% If the user selected an image, update the GUI with the selected image
if isequal(filename,0) || isequal(pathname,0)
    % User cancelled the dialog box
else
    % Load the selected image
    image_path = fullfile(pathname, filename);
    handles.original_image = imread(image_path);
    handles.current_image = handles.original_image; % use original image as current image
    imshow(handles.current_image, 'Parent', handles.axes1);
end

% Update the handles structure
guidata(hObject, handles);


% --- Executes on button press in Show Original.
function ShowOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to ShowOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Display the original captured frame
axes(handles.axes1);
imshow(handles.original_image);  % Display the original captured frame

% Update the handles structure to use the original image as the current image
handles.current_image = handles.original_image;
guidata(hObject, handles);


% --- Executes on button press in GreyScale.
function GreyScale_Callback(hObject, eventdata, handles)
% hObject    handle to GreyScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Convert the current image to grayscale
gray_image = rgb2gray(handles.current_image);

% Update the GUI with the grayscale image
imshow(gray_image, 'Parent', handles.axes1);

% Update the handles structure
handles.current_image = gray_image;
guidata(hObject, handles);


% --- Executes on button press in BlackAndWhite.
function BlackAndWhite_Callback(hObject, eventdata, handles)
% hObject    handle to Binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Convert the image to grayscale
gray_image = rgb2gray(handles.current_image);

% Determine the threshold value using Otsu's method
threshold = graythresh(gray_image);

% Convert the image to binary using the threshold
binary_image = imbinarize(gray_image, threshold);

% Update the GUI with the binary image
imshow(binary_image, 'Parent', handles.axes1);

% Update the handles structure
handles.current_image = binary_image;
guidata(hObject, handles);

% --- Executes on button press in HorizontalRotation.
function HorizontalRotation_Callback(hObject, eventdata, handles)
% hObject    handle to HorizontalRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Flip the image horizontally
rotated_image = fliplr(handles.current_image);

% Update the GUI with the rotated image
imshow(rotated_image, 'Parent', handles.axes1);

% Update the handles structure
handles.current_image = rotated_image;
guidata(hObject, handles);

% --- Executes on button press in Histogram.
function Histogram_Callback(hObject, eventdata, handles)
% hObject    handle to Histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Calculate and plot the histogram of the current image
histogram(handles.current_image, 'Parent', handles.axes1);

% Update the handles structure
guidata(hObject, handles);


% --- Executes on button press in VerticalRotation.
function VerticalRotation_Callback(hObject, eventdata, handles)
% hObject    handle to VerticalRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Flip the image vertically
rotated_image = flipdim(handles.current_image, 1);

% Update the GUI with the rotated image
imshow(rotated_image, 'Parent', handles.axes1);

% Update the handles structure
handles.current_image = rotated_image;
guidata(hObject, handles);




% --- Executes on button press in LinearConformalTransform.
function LinearConformalTransform_Callback(hObject, eventdata, handles)
% hObject    handle to LinearConformalTransform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Define the transformation matrix
tform = [1 -0.2 0; 0.2 1 0; 0 0 1]; % scale the image horizontally by 1.2 and vertically by 0.8

% Apply the transformation to the image
transformed_image = imwarp(image, affine2d(tform));

% Update the GUI with the transformed image
imshow(transformed_image, 'Parent', handles.axes1);

% Update the handles structure
handles.current_image = transformed_image;
guidata(hObject, handles);





% --- Executes on button press in RotateImage.
function RotateImage_Callback(hObject, eventdata, handles)
% hObject    handle to RotateImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the rotation angle
prompt = {'Enter the rotation angle (in degrees):'};
dlgtitle = 'Rotate Image';
dims = [1 35];
definput = {'30'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
angle = str2double(answer{1});

% Rotate the image using imrotate
rotated_image = imrotate(image, angle, 'nearest', 'crop');

% Update the GUI with the rotated image
axes(handles.axes1);
imshow(rotated_image);

% Update the handles structure
handles.current_image = rotated_image;
guidata(hObject, handles);


% --- Executes on button press in TranslateImage.
function TranslateImage_Callback(hObject, eventdata, handles)
% hObject    handle to TranslateImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the translation distance
prompt = {'Enter the horizontal translation distance (in pixels):', ...
          'Enter the vertical translation distance (in pixels):'};
dlgtitle = 'Translate Image';
dims = [1 35];
definput = {'50', '50'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
dx = str2double(answer{1});
dy = str2double(answer{2});

% Translate the image using imtranslate
translated_image = imtranslate(image, [dx, dy]);

% Update the GUI with the translated image
axes(handles.axes1);
imshow(translated_image);

% Update the handles structure
handles.current_image = translated_image;
guidata(hObject, handles);


% --- Executes on button press in ShearImageHorizontally.
function ShearImageHorizontally_Callback(hObject, eventdata, handles)
% hObject    handle to ShearImageHorizontally (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the shear angle
prompt = {'Enter the shear angle (in degrees):'};
dlgtitle = 'Shear Image Horizontally';
dims = [1 35];
definput = {'30'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
angle = str2double(answer{1});

% Shear the image horizontally using imwarp
tform = affine2d([1 0 0; tand(angle) 1 0; 0 0 1]);
sheared_image = imwarp(image, tform);

% Update the GUI with the sheared image
axes(handles.axes1);
imshow(sheared_image);

% Update the handles structure
handles.current_image = sheared_image;
guidata(hObject, handles)


% --- Executes on button press in ShearImageVertically.
function ShearImageVertically_Callback(hObject, eventdata, handles)
% hObject    handle to ShearImageVertically (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the shear angle
prompt = {'Enter the shear angle (in degrees):'};
dlgtitle = 'Shear Image Vertically';
dims = [1 35];
definput = {'30'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
angle = str2double(answer{1});

% Shear the image vertically using imwarp
tform = affine2d([1 tand(angle) 0; 0 1 0; 0 0 1]);
sheared_image = imwarp(image, tform);

% Update the GUI with the sheared image
axes(handles.axes1);
imshow(sheared_image);

% Update the handles structure
handles.current_image = sheared_image;
guidata(hObject, handles);

% --- Executes on button press in GaussianFilter.
function GaussianFilter_Callback(hObject, eventdata, handles)
% hObject    handle to GaussianFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt the user for the standard deviation of the Gaussian filter
prompt = {'Enter the standard deviation of the Gaussian filter:'};
dlgtitle = 'Gaussian Filter';
dims = [1 35];
definput = {'1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
sigma = str2double(answer{1});

% Apply the Gaussian filter to the current image
filtered_image = imgaussfilt(handles.current_image, sigma);

% Show the filtered image in the GUI
imshow(filtered_image, 'Parent', handles.axes1);

% Update the handles structure to store the filtered image as the new current image
handles.current_image = filtered_image;
guidata(hObject, handles);


function MaxFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MaxFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the filter size
prompt = {'Enter the filter size:'};
dlgtitle = 'Max Filter';
dims = [1 35];
definput = {'3'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
filter_size = str2double(answer{1});

% Apply the max filter using imdilate
max_filtered_image = imdilate(image, ones(filter_size));

% Update the GUI with the max filtered image
axes(handles.axes1);
imshow(max_filtered_image);

% Update the handles structure
handles.current_image = max_filtered_image;
guidata(hObject, handles);



function MinFilter_Callback(hObject, eventdata, handles)
% hObject    handle to the MinFilter button
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current image
image = handles.current_image;

% Prompt the user for the size of the min filter
prompt = {'Enter the size of the min filter (odd integer):'};
dlgtitle = 'Min Filter';
dims = [1 35];
definput = {'3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    % User cancelled the dialog box
    return;
end
filter_size = str2double(answer{1});

% Apply the min filter using imerode
se = strel('square', filter_size);
filtered_image = imerode(image, se);

% Update the GUI with the filtered image
axes(handles.axes1);
imshow(filtered_image);

% Update the handles structure
handles.current_image = filtered_image;
guidata(hObject, handles);


%Face Detection code section


function testing_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
axes(handles.axes1);
imshow('blank.jpg');
axis off;
guidata(hObject, handles);



function varargout = testing_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
% Initialize current_image as an empty matrix
handles.current_image = [];

% start function
function startpc_Callback(hObject, eventdata, handles)
% in code baraye peydakardan webcam id system
webcam_info = imaqhwinfo('winvideo');
webcam_id = webcam_info.DeviceIDs;
handles.vid = videoinput('winvideo' , 1, 'YUY2_640X480');
guidata(hObject, handles);

%face detection function
function facedetect_Callback(hObject, eventdata, handles)
triggerconfig(handles.vid ,'manual');
set(handles.vid, 'TriggerRepeat',inf);
set(handles.vid, 'FramesPerTrigger',1);
handles.vid.ReturnedColorspace = 'rgb';
 handles.vid.Timeout = 5;
start(handles.vid);
while(1)

facedetector = vision.CascadeObjectDetector;                                                 
trigger(handles.vid); 
handles.im = getdata(handles.vid, 1);
bbox = step(facedetector, handles.im);
hello = insertObjectAnnotation(handles.im,'rectangle',bbox,'Face');
imshow(hello);
end
guidata(hObject, handles);


% --- Executes on button press in Stop button for face detection.
function stopc_Callback(hObject, eventdata, handles)
handles.output = hObject;

% Stop the webcam
stop(handles.vid);

% Capture the current frame
captured_frame = getdata(handles.vid, 1);

% Display the captured frame in its original color
axes(handles.axes1);
imshow(captured_frame);  % Display the raw (color) captured frame

% Update the handles structure with the captured frame
handles.current_image = captured_frame;

% Store the original captured frame in the handles structure
handles.original_image = captured_frame;

% Perform additional image processing functions here if needed

% Update the handles structure
guidata(hObject, handles);


