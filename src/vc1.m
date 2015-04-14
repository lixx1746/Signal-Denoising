function varargout = vc1(varargin)
% VC1 MATLAB code for vc1.fig
%      VC1, by itself, creates a new VC1 or raises the existing
%      singleton*.
%
%      H = VC1 returns the handle to a new VC1 or the handle to
%      the existing singleton*.
%
%      VC1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VC1.M with the given input arguments.
%
%      VC1('Property','Value',...) creates a new VC1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vc1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vc1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vc1

% Last Modified by GUIDE v2.5 13-Dec-2011 22:24:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vc1_OpeningFcn, ...
                   'gui_OutputFcn',  @vc1_OutputFcn, ...
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


% --- Executes just before vc1 is made visible.
function vc1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vc1 (see VARARGIN)

% Choose default command line output for vc1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vc1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vc1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function snr_Callback(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%guidata(hObject,handles);
%snr=get(handles.snr,'Value');
global sigdata
global noise_data
type=get(handles.datatype,'Value');
switch type
    case 1
        sigdata=MakeSignal('Blocks', 128);
    case 2
        sigdata=MakeSignal('HeaviSine', 128);
end
snr=get(hObject,'Value');

        t=linspace(0,1,128);
        snrdb=10*log10(snr);
        set(handles.edit4,'string',num2str(snrdb));
        noise_data=awgn(sigdata,snrdb);
        axes(handles.axes1);
        plot(t,sigdata);
        axes(handles.axes2);
        plot(t,noise_data);
        guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in datatype.
function datatype_Callback(hObject, eventdata, handles)
% hObject    handle to datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns datatype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from datatype
global sigdata
global noise_data
%global noise_data1
val = get(hObject,'Value');
switch val
    case 1
        sigdata=MakeSignal('Blocks', 128);
        snr=get(hObject,'Value');
        t=linspace(0,1,128);
        snrdb=10*log10(snr);
        noise_data=awgn(sigdata,snrdb);
        axes(handles.axes1);
        plot(t,sigdata);
        axes(handles.axes2);
        plot(t,noise_data);

    case 2
        sigdata=MakeSignal('HeaviSine', 128);
        snr=get(hObject,'Value');
        t=linspace(0,1,128);
        snrdb=10*log10(snr);
        noise_data=awgn(sigdata,snrdb);
        axes(handles.axes1);
        plot(t,sigdata);
        axes(handles.axes2);
        plot(t,noise_data);

guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function datatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

        
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1
global sigdata
global noise_data
global wc
global wc2
global wceo
global wcoef
%global noise_data1
%snr=get(handles.snr,'Value');
h=get(handles.axes1,'YLim');
n=128;
R_pre=zeros(1,n/2);
        t=linspace(0,1,128);
        %noise_data=awgn(sigdata, snr,0);
        
        qmf = MakeONFilter('Symmlet',8);
        wc = FWT_PO(noise_data,0,qmf);
        %out2=IWT_PO(wc, 0, qmf);
        %err=out2-noise_data;
        var_est=median(abs(wc(2:128)))/0.6745;
        out1 = ThreshWave(noise_data,'S',0,var_est,sqrt(2*log(128)),0,qmf);
        wc2= FWT_PO(out1,0,qmf);
        b=length(find(abs(wc2)>0.001));
        axes(handles.axes3);
        plot(t,out1);
        axis([0 1 h]);
        [out2,wcoef]  = WaveShrink(noise_data,'SURE',0,qmf);
        c=length(find(abs(wcoef)>0.001));
        axes(handles.axes5);
        plot(t,out2);
        axis([0 1 h]);
        wc1(1)=abs(wc(1));
        for i=0:6
        wc1((2^i+1):(2^(i+1)))=abs(wc((2^i+1):(2^(i+1)))).*(0.5^i);
        end
        
         [wcnew,index]=sort(wc1);
         wcnew=fliplr(wcnew);
         srmindex=fliplr(index);
        
         for m= 1:(n/2)
           
         wceo=zeros(1,n);
         wceo(srmindex(1:m))=wc(srmindex(1:m));
         out = IWT_PO(wceo,0,qmf);
         p=(m+1)/n;
         R_emp1=1/n.*sum((noise_data-out).^2);
         a=(1-(p-p*log(p)+log(n)/(2*n)).^0.5);
         r(m)=1/max((1-(p-p*log(p)+log(n)/(2*n)).^0.5),0);
         R_pre1=(R_emp1*r(m));
         R_pre(m)=(round(R_pre1.*1))./1;
         end
        
         minerror=min(R_pre);
         a=min(find(R_pre==min(R_pre)));
         wceo=zeros(1,n);
         wceo(srmindex(1:a))=wc(srmindex(1:a));
         optimal_out = IWT_PO(wceo,0,qmf);
         %wceo1=FWT_PO(optimal_out,0,qmf);
         a1=length(find(abs(wceo)>0.001));
         axes(handles.axes4);
         plot(t,optimal_out);
         axis([0 1 h]);
         set(handles.edit1,'string', num2str((sqrt(sum((sigdata-out1).^2)/128)/std(sigdata))));
          set(handles.edit6,'string', num2str((sqrt(sum((sigdata-out2).^2)/128)/std(sigdata))));
         set(handles.edit2,'string', num2str((sqrt(sum((sigdata- optimal_out).^2)/128)/std(sigdata))));
         set(handles.edit5,'string',num2str(a1));
         set(handles.edit7,'string',num2str(b));
         set(handles.edit8,'string',num2str(c));
        guidata(hObject,handles);


% --- Executes on button press in draw.
function draw_Callback(hObject, eventdata, handles)
% hObject    handle to draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global wc
global wc2
global wceo
global wcoef
% Hint: get(hObject,'Value') returns toggle state of draw
figure(1);
subplot(1,2,1);
PlotWaveCoeff(wc,0,0);
title('noise coe');
subplot(1,2,2);
PlotWaveCoeff(wc2,0,0);
title('VISU coe');
figure(2);
subplot(1,2,1);
PlotWaveCoeff(wc,0,0);
title('noise coe');
subplot(1,2,2);
PlotWaveCoeff(wceo,0,0);
title('VC-based coe');
figure(3);
subplot(1,2,1);
PlotWaveCoeff(wc,0,0);
title('noise coe');
subplot(1,2,2);
PlotWaveCoeff(wcoef,0,0);
title('SURE coe');



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cp.
function cp_Callback(hObject, eventdata, handles)
% hObject    handle to cp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cp
val = get(handles.datatype,'Value');
snr=get(handles.snr,'Value');
n1=128;
%global noise_data1
for k=1:300
    R_pre=zeros(1,n1/2);
    if val==1
       sigdata=MakeSignal('Blocks',n1);
       snrdb=10*log10(snr);
       noise_data=awgn(sigdata,snrdb);  
       qmf = MakeONFilter('Symmlet',8);
       wc = FWT_PO(noise_data,0,qmf);
       var_est=median(abs(wc(2:n1)))/0.6745;
       out1 = ThreshWave(noise_data,'S',0,var_est,sqrt(2*log(n1)),0,qmf);
       wcoef1= FWT_PO(out1,0,qmf);
      
       [out2,wcoef2]  = WaveShrink(noise_data,'SURE',0,qmf);
       [out3,wcoef3]  = WaveShrink(noise_data,'Hybrid',0,qmf);
       [out4,wcoef4]  = WaveShrink(noise_data,'MinMax',0,qmf);
       wc1(1)=abs(wc(1));
        for i=0:(log2(n1)-1)

        wc1((2^i+1):(2^(i+1)))=abs(wc((2^i+1):(2^(i+1)))).*(0.5^i);
        end
        
         [wcnew,index]=sort(wc1);
         srmindex=fliplr(index);
        for m= 1:n1/2
         wceo=zeros(1,n1);
         wceo(srmindex(1:m))=wc(srmindex(1:m));
         out = IWT_PO(wceo,0,qmf);
         p=(m+1)/n1;
         R_emp1=1/n1.*sum((noise_data-out).^2);
         r=1/max((1-(p-p*log(p)+log(n1)/(2*n1)).^0.5),0);
         R_pre1=(R_emp1*r);
         R_pre(m)=(round(R_pre1.*1))./1
        end
         minerror=min(R_pre);
         a=min(find(R_pre==min(R_pre)));
         wceo=zeros(1,n1);
         wceo(srmindex(1:a))=wc(srmindex(1:a));
         out5 = IWT_PO(wceo,0,qmf);
         nrms_visu(k)=sqrt(sum((sigdata-out1).^2)/n1)/std(sigdata);
         nrms_sure(k)=sqrt(sum((sigdata-out2).^2)/n1)/std(sigdata);
         nrms_hybrid(k)=sqrt(sum((sigdata-out3).^2)/n1)/std(sigdata);
         nrms_min(k)=sqrt(sum((sigdata-out4).^2)/n1)/std(sigdata);
         nrms_vc(k)=sqrt(sum((sigdata-out5).^2)/n1)/std(sigdata);
         dof_visu(k)=length(find(abs(wcoef1)>0.001));
         dof_sure(k)=length(find(abs(wcoef2)>0.001));
         dof_hybrid(k)=length(find(abs(wcoef3)>0.001));
         dof_min(k)=length(find(abs(wcoef4)>0.001));
         dof_vc(k)=length(find(abs(wceo)>0.000000000000000001));
         
    else
        sigdata=MakeSignal('HeaviSine', n1); 
        snrdb=10*log10(snr);
        noise_data=awgn(sigdata,snrdb);
         qmf = MakeONFilter('Symmlet',8);
       wc = FWT_PO(noise_data,0,qmf);
       var_est=median(abs(wc(2:n1)))/0.6745;
       out1 = ThreshWave(noise_data,'S',0,var_est,sqrt(2*log(n1)),0,qmf);
       wcoef1= FWT_PO(out1,0,qmf);
       [out2,wcoef2]  = WaveShrink(noise_data,'SURE',0,qmf);
       [out3,wcoef3]  = WaveShrink(noise_data,'Hybrid',0,qmf);
       [out4,wcoef4]  = WaveShrink(noise_data,'MinMax',0,qmf);
       wc1(1)=abs(wc(1));
        for i=0:(log2(n1)-1)

        wc1((2^i+1):(2^(i+1)))=abs(wc((2^i+1):(2^(i+1)))).*(0.5^i);
        end
         [wcnew,index]=sort(wc1);
         srmindex=fliplr(index);
        for m= 1:n1/2
         wceo=zeros(1,n1);
         wceo(srmindex(1:m))=wc(srmindex(1:m));
         out = IWT_PO(wceo,0,qmf);
         p=(m+1)/n1;
         R_emp1=1/n1.*sum((noise_data-out).^2);
         r=1/max((1-(p-p*log(p)+log(n1)/(2*n1)).^0.5),0);
         R_pre1=(R_emp1*r);
         R_pre(m)=(round(R_pre1.*1))./1;
        end
         minerror=min(R_pre);
         a=min(find(R_pre==min(R_pre)));
         wceo=zeros(1,n1);
         wceo(srmindex(1:a))=wc(srmindex(1:a));
         out5 = IWT_PO(wceo,0,qmf);
         nrms_visu(k)=sqrt(sum((sigdata-out1).^2)/n1)/std(sigdata);
         nrms_sure(k)=sqrt(sum((sigdata-out2).^2)/n1)/std(sigdata);
         nrms_hybrid(k)=sqrt(sum((sigdata-out3).^2)/n1)/std(sigdata);
         nrms_min(k)=sqrt(sum((sigdata-out4).^2)/n1)/std(sigdata);
         nrms_vc(k)=sqrt(sum((sigdata-out5).^2)/n1)/std(sigdata);
         dof_visu(k)=length(find(abs(wcoef1)>0.001));
         dof_sure(k)=length(find(abs(wcoef2)>0.001));
         dof_hybrid(k)=length(find(abs(wcoef3)>0.001));
         dof_min(k)=length(find(abs(wcoef4)>0.001));
         dof_vc(k)=length(find(abs(wceo)>0.0010000000000000));
    end
   
   %boxplot(nrms_visu,0,nrms_sure,0,nrms_hybrid,0,nrms_min,0,nrms_vc,0);
  % boxplot('visu',dof_visu,'sure',dof_sure,'hybrid', dof_hybrid,'min',dof_min,'vc',dof_vc);
    
end
figure(4)
d=[nrms_visu';nrms_sure';nrms_hybrid';nrms_min';nrms_vc'];
dof=[dof_visu';dof_sure';dof_hybrid';dof_min';dof_vc'];
name=[repmat('visu  ',300,1);repmat('sure  ',300,1);repmat('hybird',300,1);repmat('min   ',300,1);repmat('vc    ',300,1)];
subplot(1,2,1);
boxplot(d,name);
title('nrms of different method');
subplot(1,2,2);
boxplot(dof,name)
title('dof of different method');

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
