function varargout = HHsimITgui(varargin)
% HHSIMITGUI MATLAB code for HHsimITgui.fig
%      HHSIMITGUI, by itself, creates a new HHSIMITGUI or raises the existing
%      singleton*.
%
%      H = HHSIMITGUI returns the handle to a new HHSIMITGUI or the handle to
%      the existing singleton*.
%
%      HHSIMITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HHSIMITGUI.M with the given input arguments.
%
%      HHSIMITGUI('Property','Value',...) creates a new HHSIMITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HHsimITgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HHsimITgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HHsimITgui

% Last Modified by GUIDE v2.5 24-Feb-2014 14:31:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HHsimITgui_OpeningFcn, ...
    'gui_OutputFcn',  @HHsimITgui_OutputFcn, ...
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


% --- Executes just before HHsimITgui is made visible.
function HHsimITgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HHsimITgui (see VARARGIN)

% Choose default command line output for HHsimITgui
handles.output = hObject;

% ==== Set initial state ===========
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HHsimITgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HHsimITgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function C_in_Na_Callback(hObject, eventdata, handles)
% hObject    handle to C_in_Na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_in_Na as text
%        str2double(get(hObject,'String')) returns contents of C_in_Na as a double
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));


% --- Executes during object creation, after setting all properties.
function C_in_Na_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_in_Na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_out_Na_Callback(hObject, eventdata, handles)
% hObject    handle to C_out_Na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_out_Na as text
%        str2double(get(hObject,'String')) returns contents of C_out_Na as a double
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));


% --- Executes during object creation, after setting all properties.
function C_out_Na_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_out_Na (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_out_K_Callback(hObject, eventdata, handles)
% hObject    handle to C_out_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_out_K as text
%        str2double(get(hObject,'String')) returns contents of C_out_K as a double
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));



% --- Executes during object creation, after setting all properties.
function C_out_K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_out_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function C_in_K_Callback(hObject, eventdata, handles)
% hObject    handle to C_in_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C_in_K as text
%        str2double(get(hObject,'String')) returns contents of C_in_K as a double
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
% Vrest =  calcSimSS(t, V,E_L, E_Na, E_K, gbar_Na,gbar_K, deltaT);
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));


% --- Executes during object creation, after setting all properties.
function C_in_K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_in_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runSimPush.
function runSimPush_Callback(hObject, eventdata, handles)
% hObject    handle to runSimPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ======== get all values ==========

%===simulation time===
simulationTime = str2double(get(handles.simDur, 'String')); %in milliseconds
deltaT=.001;
Fs = 1/deltaT;
t=0:deltaT:simulationTime;
vc = get(handles.VCbutton, 'Value');

% drugs application
pctTTX = 0.01*str2double(get(handles.pctTTX, 'String'));
pctTEA = 0.01*str2double(get(handles.pctTEA, 'String'));
TTX = 1-pctTTX;
TEA = 1 - pctTEA;

%===constant parameters===%
%All of these can be found in Table 3
C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));
% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin

E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
E_L=10.6;
gbar_K=36*TEA;
gbar_Na=120*TTX;
g_L=.3;
C=1;

Vrest = K*((gbar_K*C_out_K+gbar_Na*C_out_Na)/(gbar_K*C_in_K+gbar_Na*C_in_Na));
%===set the initial states===%
V=0; %Baseline voltage
Vset = 58.449;
Vrest = str2double(get(handles.VrestVal, 'String'))
E_K = E_K_orig- Vrest;
E_Na = E_Na_orig - Vrest;
V = Vrest+Vset;

alpha_n = (.01 * ( (10-V) / (exp((10-V)/10)-1) ));
beta_n = .125*exp(-V/80);
alpha_m = (.1*( (25-V) / (exp((25-V)/10)-1) ));
beta_m = 4*exp(-V/18);
alpha_h = .07*exp(-V/20);
beta_h = 1/(exp((30-V)/10)+1);

n(1) = alpha_n/(alpha_n+beta_n); %Equation 9
m(1) = alpha_m/(alpha_m+beta_m); %Equation 18
h(1) = alpha_h/(alpha_h+beta_h); %Equation 18

% Define Vrest
g_Na = (gbar_Na*(m(1).^3).*h(1)); % .01
g_K = (gbar_K*(n(1).^4)); % .039
g_L = .3;

% Vrest = I_ion/g_ion;
% Vrest = ((g_Na*E_Na_orig)+(g_K*E_K_orig))/(g_Na+g_K);

% E_K = -12;
% E_Na = 115;


I_Na = zeros(1,length(t)); %Equations 3 and 14
I_K = zeros(1,length(t)); %Equations 4 and 6
I_L = zeros(1,length(t)); %Equation 5
I_ion = zeros(1,length(t));

I = zeros(1,numel(t));
I1 = zeros(1,numel(t));
I2 = zeros(1,numel(t));
V = zeros(1, numel(t));
%% ======== define stimuli =====================
if vc
    
    voltageLevels1 = str2double(get(handles.stim1Amp, 'String'))+abs(Vrest); %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
    
    stimS1 = str2double(get(handles.stim1S, 'String')); % in ms
    stimE1 = stimS1 + str2double(get(handles.stim1Dur, 'String'));
    if stimE1 ~= stimS1
        if stimS1 == 0, stimS1 = deltaT, end
        V(round(stimS1*Fs):round(stimE1*Fs)) = voltageLevels1;
    end
    
    voltageLevels2 = str2double(get(handles.stim2Amp, 'String'))+abs(Vrest); %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
    
    stimS2 = str2double(get(handles.stim2S, 'String')); % in ms
    stimE2 = stimS2 + str2double(get(handles.stim2Dur, 'String'));
    if stimE2 ~= stimS2
        if stimS2 == 0, stimS2 =deltaT, end
        V(round(stimS2*Fs):round(stimE2*Fs)) = voltageLevels2;
    end
else
    %Set externally applied current across time
    %Here, first 500 timesteps are at current of 50, next 1500 timesteps at
    %current of zero (resets resting potential of neuron), and the rest of
    %timesteps are at constant current
    currentLevels1 = str2double(get(handles.stim1Amp, 'String')); %Change this to see effect of different currents on voltage (Suggested values: 3, 20, 50, 1000)
    stimS1 = str2double(get(handles.stim1S, 'String')); % in ms
    stimE1 = stimS1 + str2double(get(handles.stim1Dur, 'String'));
    if stimE1 ~= stimS1
        if stimS1 == 0, stimS1 =deltaT, end
        I1(round(stimS1*Fs):round(stimE1*Fs)) = currentLevels1;
    end
    
    currentLevels2 = str2double(get(handles.stim2Amp, 'String'));
    stimS2 = str2double(get(handles.stim2S, 'String')); % in ms
    stimE2 = stimS2 + str2double(get(handles.stim2Dur, 'String'));
    if stimE2 ~= stimS2
        if stimS2 == 0, stimS2 =deltaT, end
        I2(round(stimS2*Fs):round(stimE2*Fs)) = currentLevels2;
    end
    I = I1+I2;
    %plot the stimulus
    axes(handles.axes5)
    stairs(t,I,'r','linewidth', 2)
    set(handles.axes5,'Ycol',[1 1 1])
    set(handles.axes5,'Xcol',[1 1 1])
    xlabel('time (ms)')
    set(handles.axes5,'color',[0 0 0])
    
    %     grid on
    %Comment out the above line and uncomment the line below for constant current, and observe effects on voltage timecourse
    %I(1:numel(t)) = currentLevels;
end

for i=1:numel(t)-1 %Compute coefficients, currents, and derivates at each time step
    
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    alpha_n(i) = (.01 * ( (10-V(i)) / (exp((10-V(i))/10)-1) ));
    beta_n(i) = .125*exp(-V(i)/80);
    alpha_m(i) = (.1*( (25-V(i)) / (exp((25-V(i))/10)-1) ));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = .07*exp(-V(i)/20);
    beta_h(i) = 1/(exp((30-V(i))/10)+1);
    
    
    %---calculate the currents---%
    I_Na(i) = (m(i)^3) * gbar_Na * h(i) * (V(i)-E_Na); %Equations 3 and 14
    I_K(i) = (n(i)^4) * gbar_K * (V(i)-E_K); %Equations 4 and 6
    I_L(i) = g_L *(V(i)-E_L); %Equation 5
    I_ion(i) = I(i) - I_K(i) - I_Na(i) - I_L(i);
    
    
    
    %---calculate the derivatives using Euler first order approximation---%
    if ~vc
        V(i+1) = V(i) + deltaT*I_ion(i)/C;
    end
    %     I_ion(i) = (C*(-V(i)+V(i+1)))/deltaT;
    n(i+1) = (n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i))); %Equation 7
    m(i+1) = (m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i))); %Equation 15
    h(i+1) = (h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i))); %Equation 16
    
end


V = V+Vrest; %Set resting potential 

%===plot Voltage===%
axes(handles.axes1)
plot(t,V,'LineWidth',2,'color',[1 1 1])
hold on
set(handles.axes1,'color',[1 1 1])
set(handles.axes1,'Ycol',[1 1 1])
set(handles.axes1,'Xcol',[1 1 1])
legend('\color[rgb]{1, 1, 1} voltage')
legend boxoff
ylabel('Voltage (mv)')
set(handles.axes1,'Color',[0 0 0])

% grid on

% title('Voltage over Time in Simulated Neuron')


%===plot Conductance===%
axes(handles.axes2)
p1 = plot(t,gbar_K*n.^4,'b','LineWidth',2);
hold on
p2 = plot(t,gbar_Na*(m.^3).*h,'r','LineWidth',2);
set(handles.axes2,'color',[1 1 1])
set(handles.axes2,'Ycol',[1 1 1])
set(handles.axes2,'Xcol',[1 1 1])
legend([p1, p2], '\color{blue} g K', '\color{red} g Na')
legend boxoff
ylabel('Conductance (mS)')
set(handles.axes2,'Color',[0 0 0])

% grid on
% title('Conductance for Potassium and Sodium Ions in Simulated Neuron')

% ======== plot n, m and h ================
axes(handles.axes3)
p1 = plot(t,n,'b','LineWidth',2);
hold on
p2 = plot(t,m,'r','LineWidth',2);
p3 = plot(t,h,'g','LineWidth',2);
% title('m h and n')
set(handles.axes3,'Ycol',[1 1 1])
set(handles.axes3,'Xcol',[1 1 1])
set(handles.axes3,'color',[1 1 1])
legend('\color{blue} n','\color{red} m','\color{green} h')
legend boxoff
set(handles.axes3,'Color',[0 0 0])

% grid on

% ========== plot individual currents ================
axes(handles.axes4)
p1 = plot(t,I_Na,'r','LineWidth',2);
hold on
p2 = plot(t,I_K,'b','LineWidth',2);
Itot = I_Na+I_K+I_L-I;
p3 = plot(t,Itot,'g','LineWidth',2);
% title('Currents')
set(handles.axes4,'color',[1 1 1])
set(handles.axes4,'Ycol',[1 1 1])
set(handles.axes4,'Xcol',[1 1 1])
legend('\color{red} I Na', '\color{blue} I K', '\color{green} Total')
legend boxoff
ylabel('Current (mA)')
set(handles.axes4,'Color',[0 0 0])

% grid on

% update show values handles
handles.tShowValt = t;
handles.VshowValt = V;
handles.g_NaShowValt = gbar_Na*(m.^3).*h;
handles.g_KShowValt = gbar_K*n.^4;
handles.I_NaShowValt = I_Na;
handles.I_KShowValt = I_K;
handles.I_totalShowValt = Itot;
handles.mShowValt = m;
handles.nShowValt = n;
handles.hShowValt = h;
% Update handles structure
guidata(hObject, handles);
% ====================== END user code =======================






function stim2S_Callback(hObject, eventdata, handles)
% hObject    handle to stim2S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim2S as text
%        str2double(get(hObject,'String')) returns contents of stim2S as a double


% --- Executes during object creation, after setting all properties.
function stim2S_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim2S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim2Dur_Callback(hObject, eventdata, handles)
% hObject    handle to stim2Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim2Dur as text
%        str2double(get(hObject,'String')) returns contents of stim2Dur as a double


% --- Executes during object creation, after setting all properties.
function stim2Dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim2Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim2Amp_Callback(hObject, eventdata, handles)
% hObject    handle to stim2Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim2Amp as text
%        str2double(get(hObject,'String')) returns contents of stim2Amp as a double


% --- Executes during object creation, after setting all properties.
function stim2Amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim2Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim1S_Callback(hObject, eventdata, handles)
% hObject    handle to stim1S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim1S as text
%        str2double(get(hObject,'String')) returns contents of stim1S as a double


% --- Executes during object creation, after setting all properties.
function stim1S_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim1S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim1Dur_Callback(hObject, eventdata, handles)
% hObject    handle to stim1Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim1Dur as text
%        str2double(get(hObject,'String')) returns contents of stim1Dur as a double


% --- Executes during object creation, after setting all properties.
function stim1Dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim1Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stim1Amp_Callback(hObject, eventdata, handles)
% hObject    handle to stim1Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim1Amp as text
%        str2double(get(hObject,'String')) returns contents of stim1Amp as a double


% --- Executes during object creation, after setting all properties.
function stim1Amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim1Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pctTTX_Callback(hObject, eventdata, handles)
% hObject    handle to pctTTX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pctTTX as text
%        str2double(get(hObject,'String')) returns contents of pctTTX as a double


% --- Executes during object creation, after setting all properties.
function pctTTX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pctTTX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pctTEA_Callback(hObject, eventdata, handles)
% hObject    handle to pctTEA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pctTEA as text
%        str2double(get(hObject,'String')) returns contents of pctTEA as a double


% --- Executes during object creation, after setting all properties.
function pctTEA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pctTEA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function simDur_Callback(hObject, eventdata, handles)
% hObject    handle to simDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simDur as text
%        str2double(get(hObject,'String')) returns contents of simDur as a double


% --- Executes during object creation, after setting all properties.
function simDur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearPush.
function clearPush_Callback(hObject, eventdata, handles)
% hObject    handle to clearPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset')
set(handles.axes1,'color',[0 0 0])
cla(handles.axes2,'reset')
set(handles.axes2,'color',[0 0 0])
cla(handles.axes3,'reset')
set(handles.axes3,'color',[0 0 0])
cla(handles.axes4,'reset')
set(handles.axes4,'color',[0 0 0])
cla(handles.axes5,'reset')
set(handles.axes5,'color',[0 0 0])


% --- Executes on button press in pushShowVal.
function pushShowVal_Callback(hObject, eventdata, handles)
% hObject    handle to pushShowVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tShowVal = str2double(get(handles.tShowVal,'String'));
VshowVal = handles.VshowValt(handles.tShowValt == tShowVal);
set(handles.VshowVal,'String',VshowVal);
g_NaShowVal = handles.g_NaShowValt(handles.tShowValt == tShowVal);
set(handles.g_NaShowVal,'String',g_NaShowVal);
g_KShowVal = handles.g_KShowValt(handles.tShowValt == tShowVal);
set(handles.g_KShowVal,'String',g_KShowVal);
I_NaShowVal = handles.I_NaShowValt(handles.tShowValt == tShowVal);
set(handles.I_NaShowVal,'String',I_NaShowVal);
I_KShowVal = handles.I_KShowValt(handles.tShowValt == tShowVal);
set(handles.I_KShowVal,'String',I_KShowVal);
I_totalShowVal = handles.I_totalShowValt(handles.tShowValt == tShowVal);
set(handles.I_totalShowVal,'String',I_totalShowVal);
mShowVal = handles.mShowValt(handles.tShowValt == tShowVal);
set(handles.mShowVal,'String',mShowVal);
nShowVal = handles.nShowValt(handles.tShowValt == tShowVal);
set(handles.nShowVal,'String',nShowVal);
hShowVal = handles.hShowValt(handles.tShowValt == tShowVal);
set(handles.hShowVal,'String',hShowVal);


function tShowVal_Callback(hObject, eventdata, handles)
% hObject    handle to tShowVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tShowVal as text
%        str2double(get(hObject,'String')) returns contents of tShowVal as a double


% --- Executes during object creation, after setting all properties.
function tShowVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tShowVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in buttonSimType.
function buttonSimType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in buttonSimType
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipanel6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in defValPush.
function defValPush_Callback(hObject, eventdata, handles)
% hObject    handle to defValPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.C_out_K,'String',num2str(22));
set(handles.C_in_K,'String',num2str(400));
set(handles.C_out_Na,'String',num2str(460));
set(handles.C_in_Na,'String',num2str(50));

C_out_Na = str2double(get(handles.C_out_Na, 'String'));
C_in_Na = str2double(get(handles.C_in_Na, 'String'));
C_out_K = str2double(get(handles.C_out_K, 'String'));
C_in_K = str2double(get(handles.C_in_K, 'String'));

% R = 8314.4621; %J/kmol*K
% F = 9.64853399e4; % Faraday constant
K = 0.086173; % R/zF
T = 273.16 + 16.91; % in kelvin
E_K_orig = K*T*(log(C_out_K/C_in_K));
E_Na_orig = K*T*(log(C_out_Na/C_in_Na));
p_K = 1;
p_Na = 0.038;
% Vrest = K*T*log(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));
% Vrest =  calcSimSS(t, V,E_L, E_Na, E_K, gbar_Na,gbar_K, deltaT);
Vrest = 58*log10(((p_K*C_out_K)+(p_Na*C_out_Na))/((p_K*C_in_K)+(p_Na*C_in_Na)));

set(handles.VrestVal,'String', num2str(Vrest));
set(handles.V_K,'String',num2str(E_K_orig));
set(handles.V_Na,'String',num2str(E_Na_orig));

vc = get(handles.VCbutton, 'Value');
if ~vc
    set(handles.stim1S,'String','0');
    set(handles.stim1Dur,'String','0.2');
    set(handles.stim1Amp,'String','100');
    set(handles.stim2S,'String','0'); 
   set(handles.stim2Dur,'String','0');
   set(handles.stim2Amp,'String','0');
else
    set(handles.stim1S,'String','0.5');
    set(handles.stim1Dur,'String','5');
    set(handles.stim1Amp,'String','0');
    set(handles.stim2S,'String','5.5');
    set(handles.stim2Dur,'String','0.5');
    set(handles.stim2Amp,'String','-60');
end



