clc; clearvars; close all; warning off all;
rng('default');
%% Compare raw and aligned features on MI
%% Leave-one subject-out
%% need to enable covariancetoolbox

mAcc=cell(1,2); mTime=cell(1,2);
for ds=1:2
    %% make data
    dataFolder=['./Data/data' num2str(ds) '/'];
    files=dir([dataFolder 'A*.mat']);
    Ref=load([dataFolder 'Resting.mat']); % break time for all subjects
    XRaw=[]; yAll=[]; XAlignE=[]; XAlignR=[];
    for s=1:length(files)
        s
        load([dataFolder files(s).name]);
        XRaw=cat(3,XRaw,X);
        yAll=cat(1,yAll,y);
        nTrials=length(y);
        Bt=Ref.ref(:,:,(s-1)*nTrials+1:s*nTrials);
        RtE=mean(covariances(X),3); % reference state, Euclidean space
        RtR=riemann_mean(covariances(Bt)); % reference state, Riemmanian space
        sqrtRtE=RtE^(-1/2); sqrtRtR=RtR^(-1/2);
        XR=nan(size(X,1),size(X,2),nTrials);
        XE=nan(size(X,1),size(X,2),nTrials);
        for j=1:nTrials
            XR(:,:,j)=sqrtRtR*X(:,:,j);
            XE(:,:,j)=sqrtRtE*X(:,:,j);
        end
        XAlignE=cat(3,XAlignE,XE); XAlignR=cat(3,XAlignR,XR);
    end
    
    Accs=cell(1,length(files));
    times=cell(1,length(files));
    for t=1:length(files)    %  target user
        t
        yt=yAll((t-1)*nTrials+1:t*nTrials);
        ys=yAll([1:(t-1)*nTrials t*nTrials+1:end]);
        XtRaw=XRaw(:,:,(t-1)*nTrials+1:t*nTrials);
        XsRaw=XRaw(:,:,[1:(t-1)*nTrials t*nTrials+1:end]);
        XtAlignE=XAlignE(:,:,(t-1)*nTrials+1:t*nTrials);
        XsAlignE=XAlignE(:,:,[1:(t-1)*nTrials t*nTrials+1:end]);
        XtAlignR=XAlignR(:,:,(t-1)*nTrials+1:t*nTrials);
        XsAlignR=XAlignR(:,:,[1:(t-1)*nTrials t*nTrials+1:end]);
        
        %% mdRm
        % raw covariance matrices
        tic
        covTest=covariances(XtRaw);
        covTrain=covariances(XsRaw);
        yPred = mdm(covTest,covTrain,ys);
        Accs{t}(1)=100*mean(yt==yPred);
        times{t}(1)=toc;
        
        %align covariance matrices
        tic
        covTest=covariances(XtAlignR);
        covTrain=covariances(XsAlignR);
        yPred = mdm(covTest,covTrain,ys);
        Accs{t}(2)=100*mean(yt==yPred);
        times{t}(2)=toc;
        
        
        %% CSP+LDA
        %raw trials
        tic
        [fTrain,fTest]=CSPfeature(XsRaw,ys,XtRaw);
        LDA = fitcdiscr(fTrain,ys);
        yPred=predict(LDA,fTest);
        Accs{t}(3)=100*mean(yt==yPred);
        times{t}(3)=toc;
        
        % align trials
        tic
        [fTrain,fTest]=CSPfeature(XsAlignE,ys,XtAlignE);
        LDA = fitcdiscr(fTrain,ys);
        yPred=predict(LDA,fTest);
        Accs{t}(4)=100*mean(yt==yPred);
        times{t}(4)=toc;
    end
    mAcc{ds}=[]; mTime{ds}=[];
    for t=1:length(files)
        mAcc{ds}=cat(1,mAcc{ds},Accs{t});
        mTime{ds}=cat(1,mTime{ds},times{t});
    end
    mAcc{ds}=cat(1,mAcc{ds},mean(mAcc{ds}));
    mTime{ds}=cat(1,mTime{ds},mean(mTime{ds}));
    mAcc{ds}
    mTime{ds}
end
save('MIall.mat','mAcc','mTime');