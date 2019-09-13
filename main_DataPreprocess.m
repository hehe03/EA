clc; clearvars;

% % --------------------- Data 1: MI Data 1------------------------------
% Dataset 1: http://www.bbci.de/competition/iv/desc_1.html
% 7 subjects, 100 trials in each class; 59 EEG channels

dataFolder='./RawData/data1/';
files=dir([dataFolder 'BCICIV_ca*.mat']);
ref=[];
for s=1:length(files)
    s
    load([dataFolder files(s).name]);
    EEG=.1*double(cnt);
    b=fir1(50,2*[8 30]/nfo.fs);
    EEG=filter(b,1,EEG);
    y=mrk.y'; %(-1 for class one or 1 for class two)
    nTrials=length(y);
    X=nan(size(EEG,2),300,nTrials);
    for i=1:nTrials
        X(:,:,i)=EEG(mrk.pos(i)+0.5*nfo.fs:mrk.pos(i)+3.5*nfo.fs-1,:)'; % [0.5-3.5] seconds epoch, channels*Times
    end
    save(['./Data/data1/A' num2str(s) '.mat'],'X','y');
    %---------------- Resting-----------------
    tmp=nan(size(EEG,2),100,nTrials);
    for i=1:nTrials
        tmp(:,:,i)=EEG(mrk.pos(i)+4.25*nfo.fs:mrk.pos(i)+5.25*nfo.fs-1,:)';
    end
    ref=cat(3,ref,tmp);
end
save('./Data/data1/Resting.mat','ref');


%% --------------------- Data 2: MI Data 2a ------------------------------
ref=[];
dataFolder='./RawData/data2/';
files=dir([dataFolder '*T.gdf']);
for s=1:length(files)
    s
    [EEG, h] = sload([dataFolder files(s).name]); % need to enable bioSig toolbox
    EEG(:,end-2:end)=[]; % last three channels are EOG
    for i=1:size(EEG,2)
        EEG(isnan(EEG(:,i)),i)=nanmean(EEG(:,i));
    end
    b=fir1(50,2*[8 30]/h.SampleRate); 
    EEG=filter(b,1,EEG); 
    ids1=h.EVENT.POS(h.EVENT.TYP==769); % left hand
    ids2=h.EVENT.POS(h.EVENT.TYP==770); % right hand
    y=[zeros(length(ids1),1); ones(length(ids2),1)];
    ids=[ids1; ids2];
    X=[];
    for i=length(ids):-1:1
        X(:,:,i)=EEG(ids(i)+.5*h.SampleRate:ids(i)+3.5*h.SampleRate-1,:)';
    end
    [~,index]=sort(ids);
    y=y(index); X=X(:,:,index);
    save(['./Data/data2/A' num2str(s) '.mat'],'X','y');
    %---------------- Resting-----------------
    tmp=[];
    for i=length(ids):-1:1
        tmp(:,:,i)=EEG(ids(i)+round(4.25*h.SampleRate):ids(i)+round(5.25*h.SampleRate)-1,:)';
    end
    tmp=tmp(:,:,index);
    ref=cat(3,ref,tmp);
end
save('./Data/data2/Resting.mat','ref');
