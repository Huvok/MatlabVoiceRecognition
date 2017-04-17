% TASK Voice Recognition

%                   %Author: Hugo García
%                   %Co-Author(s): 
%                   %Date: 4/17/2017
%                   %PURPOSE: Voice recognition for biometrics security.

%==========================================================================
function [] = voicerecognition()
clc;

%--------------------------------------------------------------------------
%                   %Initial variables

%                   %Minimum distance for speech recognition.
INT_MIN_DISTANCE = 1;

%                   %Number of users to be stored in the database.
INT_MAX_USERS = 8;

% Preallocating array
str = {8}; 
fstr = {8}; 
nbtr = {8};
ste = {8}; 
fste = {8}; 
nbte = {8};
data = {INT_MAX_USERS, 4};
code = {8};

% For making Choice
ch=0;
poss=5;
while ch~=poss

    ch=menu('Speaker Recognition System','1: Human speaker recognition',...
        '2: Technical data of samples',...
        '3: Power Spectrum',...
        '4: Test with other speech files','5: Exit');
    disp('                                                                  ');
    %----------------------------------------------------------------------


%% 2: Technical data of samples
    if ch==2
        ch23=0;
        while ch23~=2
            ch23=menu('Select Option','Plot','Exit');

            if ch23==1

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');

                    ch32=0;
                    while ch32 ~=2
                        ch32=menu('Database Information','Database','Exit');

                        if ch32==1
                            charFilename1=strcat('Sound Database has : #',num2str(sound_number),'words. Enter a database number : #');
                            prompt = {charFilename1};
                            dlg_title = 'Database Information';
                            num_lines = 1;
                            def = {'1'};
                            options.Resize='on';
                            options.WindowStyle='normal';
                            options.Interpreter='tex';
                            an = inputdlg(prompt,dlg_title,num_lines,def);
                            an=cell2mat(an);
                            a = str2double(an);

                            if (isempty(an))

                            else

                                if (a <= sound_number)
                                    charFilename1=strcat('u',num2str(an));
                                    [s, fs]=wavread(charFilename1);
                                    t = 0:1/fs:(length(s) - 1)/fs;
                                    plot(t, s), axis([0, (length(s) - 1)/fs -0.4 0.5]);
                                    charFilename1=sprintf('Plot of signal %s',charFilename1);
                                    title(charFilename1);
                                    xlabel('Time [s]');
                                    ylabel('Amplitude (normalized)')

                                else
                                    warndlg('Invalid Word ','Warning');
                                end
                            end

                        end
                    end
                    close all;

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end

        end
    end
    %----------------------------------------------------------------------

%% 3: linear and logarithmic power spectrum plot
    if ch==3
        M = 100;
        N = 256;

        ch23=0;
        while ch23~=2
            ch23=menu('Select Option','Plot','Exit');

            if ch23==1

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');

                    ch32=0;
                    while ch32 ~=2
                        ch32=menu('Database Information','Database','Exit');

                        if ch32==1
                            charFilename1=strcat('Sound Database has : #',num2str(sound_number),'words. Enter a database number : #');
                            prompt = {charFilename1};
                            dlg_title = 'Database Information';
                            num_lines = 1;
                            def = {'1'};
                            options.Resize='on';
                            options.WindowStyle='normal';
                            options.Interpreter='tex';
                            an = inputdlg(prompt,dlg_title,num_lines,def);
                            an=cell2mat(an);
                            a = str2double(an);

                            if (isempty(an))

                            else

                                if (a <= sound_number)
                                charFilename1=strcat('u',num2str(an));
                                [s fs]=wavread(charFilename1);
                                frames = blockFrames(s, fs, M, N);
                                t = N / 2;
                                tm = length(s) / fs;
                                subplot(121);
                                imagesc([0 tm], [0 fs/2], abs(frames(1:t, :)).^2), axis xy;
                                title('Power Spectrum (M = 100, N = 256)');
                                xlabel('Time [s]');
                                ylabel('Frequency [Hz]');
                                colorbar;
                                %Question 3 (logarithmic)
                                subplot(122);
                                imagesc([0 tm], [0 fs/2], 20 * log10(abs(frames(1:t, :)).^2)), axis xy;
                                title('Logarithmic Power Spectrum (M = 100, N = 256)');
                                xlabel('Time [s]');
                                ylabel('Frequency [Hz]');
                                colorbar;

                                else
                                    warndlg('Invalid Word ','Warning');
                                end
                            end
                        end
                    end
                    close all;

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end

        end
    end
    %----------------------------------------------------------------------

%% 10: Voice Recognition by letting user enter into database and then compare
    if ch==4
        chos=0;
        possibility=5;
        while chos~=possibility,
            chos=menu('Speaker Recognition System','Add a new sound from microphone',...
                'Speaker recognition from microphone',...
                'Database Info','Delete database','Exit');

            %----------------------------------------------------------------------

%% 10.1 Add a new sound from microphone

            if chos==1

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');                   
                    classe = input('Insert a class number (sound ID) that will be used for recognition:');
                    if isempty(classe)
                        classe = sound_number+1;
                        disp( num2str(classe) );
                    end
                    message=('The following parameters will be used during recording:');
                    disp(message);
                    message=strcat('Sampling frequency',num2str(samplingfrequency));
                    disp(message);
                    message=strcat('Bits per sample',num2str(samplingbits));
                    disp(message);
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                    micrecorder = audiorecorder(samplingfrequency,samplingbits,1);
                    disp('Now, speak into microphone...');
                    record(micrecorder,durata);
                    
                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y1 = getaudiodata(micrecorder);
                    y = getaudiodata(micrecorder, 'uint8');

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    sound_number = sound_number+1;
                    data{sound_number,1} = y;
                    data{sound_number,2} = classe;
                    data{sound_number,3} = 'Microphone';
                    data{sound_number,4} = 'Microphone';
                    charFilename1=strcat('u',num2str(sound_number));
                    wavwrite(y1,samplingfrequency,samplingbits,charFilename1)
                    save('sound_database.dat','data','sound_number','-append');
                    msgbox('Sound added to database','Database result','help');
                    disp('Sound added to database');

                else
                    classe = input('Insert a class number (sound ID) that will be used for recognition:');
                    if isempty(classe)
                        classe = 1;
                        disp( num2str(classe) );
                    end
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                    samplingfrequency = input('Insert the sampling frequency (22050 recommended):');
                    if isempty(samplingfrequency )
                        samplingfrequency = 22050;
                        disp( num2str(samplingfrequency) );
                    end
                    samplingbits = input('Insert the number of bits per sample (8 recommended):');
                    if isempty(samplingbits )
                        samplingbits = 8;
                        disp( num2str(samplingbits) );
                    end
                    micrecorder = audiorecorder(samplingfrequency,samplingbits,1);
                    disp('Now, speak into microphone...');
                    record(micrecorder,durata);

                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y1 = getaudiodata(micrecorder);
                    y = getaudiodata(micrecorder, 'uint8');

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    sound_number = 1;
                    data{sound_number,1} = y;
                    data{sound_number,2} = classe;
                    data{sound_number,3} = 'Microphone';
                    data{sound_number,4} = 'Microphone';
                    charFilename1=strcat('u',num2str(sound_number));
                    wavwrite(y1,samplingfrequency,samplingbits,charFilename1)
                    save('sound_database.dat','data','sound_number','samplingfrequency','samplingbits');
                    msgbox('Sound added to database','Database result','help');
                    disp('Sound added to database');
                end
            end

            %----------------------------------------------------------------------

%% 10.2 Voice Recognition from microphone
            if chos==2

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');
                    Fs = samplingfrequency;
                    durata = input('Insert the duration of the recording (in seconds):');
                    if isempty(durata)
                        durata = 3;
                        disp( num2str(durata) );
                    end
                    micrecorder = audiorecorder(samplingfrequency,samplingbits,1);
                    disp('Now, speak into microphone...');
                    record(micrecorder,durata);

                    while (isrecording(micrecorder)==1)
                        disp('Recording...');
                        pause(0.5);
                    end
                    disp('Recording stopped.');
                    y = getaudiodata(micrecorder);
                    charFilename1='v';
                    wavwrite(y,samplingfrequency,samplingbits,charFilename1);
                    y = getaudiodata(micrecorder, 'uint8');
                    % if the input sound is not mono

                    if size(y,2)==2
                        y=y(:,1);
                    end
                    y = double(y);
                    %----- code for speaker recognition -------
                    disp('MFCC cofficients computation and VQ codebook training in progress...');
                    disp(' ');
                    % Number of centroids required
                    k =16;

                    for ii=1:sound_number
                        % Compute MFCC cofficients for each sound present in database
                        v = mfcc(data{ii,1}, Fs);
                        % Train VQ codebook
                        code{ii} = vqlbg(v, k);
                        disp('...');
                    end
                    disp('Completed.');
                    % Compute MFCC coefficients for input sound
                    v = mfcc(y,Fs);
                    % Current distance and sound ID initialization
                    distmin = Inf;
                    k1 = 0;

                    for ii=1:sound_number
                        d = disteu(v, code{ii});
                        dist = sum(min(d,[],2)) / size(d,1);
                        message=strcat('For User #',num2str(ii),' Dist : ',num2str(dist));
                        disp(message);
             
                        if dist < distmin
                            distmin = dist;
                            k1 = ii;
                        end
                    end

                    if distmin < INT_MIN_DISTANCE
                        min_index = k1;
                        speech_id = data{min_index,2};
                        %-----------------------------------------
                        disp('Matching sound:');
                        message=strcat('File:',data{min_index,3});
                        disp(message);
                        message=strcat('Location:',data{min_index,4});
                        disp(message);
                        message = strcat('Recognized speaker ID: ',num2str(speech_id));
                        disp(message);
                        msgbox(message,'Matching result','help');

                        ch3=0;
                        while ch3~=3
                            ch3=menu('Matched result verification:','Recognized Sound','Recorded sound','Exit');

                            if ch3==1
                                charFilename1=strcat('u',num2str(speech_id));
                                [s fs nb]=wavread(charFilename1);
                                p=audioplayer(s,fs,nb);
                                play(p);
                            end

                            if ch3==2
                                [s fs nb]=wavread('v');
                                p=audioplayer(s,fs,nb);
                                play(p);
                            end
                        end

                    else
                        warndlg('Wrong User . No matching Result.',' Warning ')
                    end
                else
                    warndlg('Database is empty. No matching is possible.',' Warning ')
                end
            end
            %----------------------------------------------------------------------

%% 10.3 Database Info
            if chos==3

                if (exist('sound_database.dat','file')==2)
                    load('sound_database.dat','-mat');
                    message=strcat('Database has #',num2str(sound_number),'words:');
                    disp(message);
                    disp(' ');

                    for ii=1:sound_number
                        message=strcat('Location:',data{ii,3});
                        disp(message);
                        message=strcat('File:',data{ii,4});
                        disp(message);
                        message=strcat('Sound ID:',num2str(data{ii,2}));
                        disp(message);
                        disp('-');
                    end

                    ch32=0;
                    while ch32 ~=2
                        ch32=menu('Database Information','Database','Exit');

                        if ch32==1
                            charFilename1=strcat('Sound Database has : #',num2str(sound_number),'words. Enter a database number : #');
                            prompt = {charFilename1};
                            dlg_title = 'Database Information';
                            num_lines = 1;
                            def = {'1'};
                            options.Resize='on';
                            options.WindowStyle='normal';
                            options.Interpreter='tex';
                            an = inputdlg(prompt,dlg_title,num_lines,def);
                            an=cell2mat(an);
                            a=str2double(an);
                            
                            if (isempty(an))

                            else

                                if (a <= sound_number)
                                    charFilename1=strcat('u',num2str(an));
                                    [s fs nb]=wavread(charFilename1);
                                    p=audioplayer(s,fs,nb);
                                    play(p);
                                else
                                    warndlg('Invalid Word ','Warning');
                                end
                            end
                        end
                    end

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end
            %----------------------------------------------------------------------

%% 10.4 Delete database
            if chos==4
                %clc;
                close all;

                if (exist('sound_database.dat','file')==2)
                    button = questdlg('Do you really want to remove the Database?');

                    if strcmp(button,'Yes')
                        load('sound_database.dat','-mat');

                        for ii=1:sound_number
                            charFilename1=strcat('u',num2str(ii),'.wav');
                            delete(charFilename1);
                        end

                        if (exist('v.wav','file')==2)
                            delete('v.wav');
                        end

                        delete('sound_database.dat');
                        msgbox('Database was succesfully removed from the current directory.','Database removed','help');
                    end

                else
                    warndlg('Database is empty.',' Warning ')
                end
            end
        end
    end

end

close all;
end
%--------------------------------------------------------------------------

%% blockFrames Function
% blockFrames: Puts the signal into frames
%
% Inputs: s contains the signal to analize
% fs is the sampling rate of the signal
% m is the distance between the beginnings of two frames
% n is the number of samples per frame
%
% Output: M3 is a matrix containing all the frames

function M3 = blockFrames(s, fs, m, n)
l = length(s);
nbFrame = floor((l - n) / m) + 1;
for i = 1:n
    for j = 1:nbFrame
        M(i, j) = s(((j - 1) * m) + i); %#ok<AGROW>
    end
end
h = hamming(n);
M2 = diag(h) * M;
for i = 1:nbFrame
    M3(:, i) = fft(M2(:, i)); %#ok<AGROW>
end
end
%--------------------------------------------------------------------------

%% MFCC Function
% MFCC
%
% Inputs: s contains the signal to analize
% fs is the sampling rate of the signal
%
% Output: r contains the transformed signal

function r = mfcc(s, fs)
m = 100;
n = 256;
frame=blockFrames(s, fs, m, n);
m = melfb(20, n, fs);
n2 = 1 + floor(n / 2);
z = m * abs(frame(1:n2, :)).^2;
r = dct(log(z));
end
%--------------------------------------------------------------------------

%% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm
% VQLBG Vector quantization using the Linde-Buzo-Gray algorithm
%
% Inputs: d contains training data vectors (one per column)
% k is number of centroids required
%
% Output: r contains the result VQ codebook (k columns, one for each  centroids)

function r = vqlbg(d,k)
e = .01;
r = mean(d, 2);
dpr = 10000;
for i = 1:log2(k)
    r = [r*(1+e), r*(1-e)];
    while (1 == 1)
        z = disteu(d, r);
        [m,ind] = min(z, [], 2);
        t = 0;
        for j = 1:2^i
            r(:, j) = mean(d(:, find(ind == j)), 2); %#ok<FNDSB>
            x = disteu(d(:, find(ind == j)), r(:, j)); %#ok<FNDSB>
            for q = 1:length(x)
                t = t + x(q);
            end
        end
        if (((dpr - t)/t) < e)
            break;
        else
            dpr = t;
        end
    end
end
end
%--------------------------------------------------------------------------