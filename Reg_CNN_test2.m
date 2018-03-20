% CNN
clc; close all; clear ;

% 실험 정보
N_mark = 28;
N_sub = 21;
Idx_sub = 1 : N_sub;
Idx_trl = 1 : 15;
Idx_sub4train = 1 : 19;
Idx_sub4test = find(countmember(Idx_sub,Idx_sub4train)==0);
% Idx_sub4train(10) = [];
Idx_trl4train = 1 : 5;
Idx_trl4test = find(countmember(Idx_trl,Idx_trl4train)==0);
Label_mark = {'central down lip';'central nose';'central upper lip';'head 1';'head 2';'head 3';'head 4';'jaw';'left central lip';'left cheek';'left dimple';'left down eye';'left down lip';'left eyebrow inside';'left eyebrow outside';'left nose';'left upper eye';'left upper lip';'right central lip';'right cheek';'right dimple';'right down eye';'right down lip';'right eyebrow inside';'right eyebrow outside';'right nose';'right upper eye';'right upper lip'};
delay = 1;
% neuronsHiddenLayer = [30 30];
% % 마커 코 부분 추출
% fpath = fullfile(cd,'DB_v2','DB_markset','mark_nose');
% load(fpath);
    
% EMG/marker delay 조절
N_delay = 5;
Idx_sub_4testing = 3;
Idx_use_mark_type = 1:3;
Idx_use_emg_feat = 1:8;

EMG_feat_for_net1 = 1:4;

% validation 조정
val_subject_indepe = 0;
use_saved_network = 0;




% 마커DB 불러오기 
i_mark = 12;
Path = fullfile('G:\CHS\git\DB_v2\DB_markset_10Hz_basecorr_znorm');
load(fullfile(Path,sprintf('mark_%d',i_mark)));
marker_set = marker_set';
win_sizes = cellfun('length',marker_set);
Marker = cell2mat(marker_set(:));
clear marker_set;

% EMG DB 폴더에 저장
% Path = fullfile(cd,'DB_v2','emg_win_10Hz_Time Alignment','comb_1');
% % Path = fullfile(cd,'DB_v2','emg_raw');
% for i_sub = 1 : 21
%     for i_trl = 1 : 15
%         fname = sprintf('sub_%03d_trl_%03d',i_sub,i_trl);
%         load(fullfile(cd,'DB_v2','emg_win_10Hz_Time Alignment',...
%             'comb_1',fname));
%         load(fullfile(cd,'DB_v2','trg_win',fname));
%         len_win = trg_w(27);
%         for i_win = 1 : len_win
%             fname = sprintf('sub_%03d_trl_%03d_win_%03d',...
%                 i_sub,i_trl,i_win);
%             emg_raw = emg_win{i_win};
%             save(fullfile(cd,'DB_v2','emg_raw',fname),'emg_raw');
%         end
%             disp([i_sub,i_trl]);
%     end
% end
% feat = feat';
% EMG = cell2mat(feat(:));
% clear feat;

% DB를 마구 조작하려면 (이미지는 sub->trl->win 순서대로 정해져 있기
% 떄문에 그 순서대로 sub, trl, win index를 정해주는 것이 좋다.
idx_sub = []; idx_trl = [];
for i_sub = 1 : 21
for i_trl = 1 : 15
    idx_sub = [idx_sub; ...
        i_sub*ones(win_sizes(i_trl,i_sub),1)];
    idx_trl = [idx_trl; ...
        i_trl*ones(win_sizes(i_trl,i_sub),1)];
end
end

% subject independent

% idx_train  = idx_sub==1 .* ...
%     logical(countmember(Idx_sub,Idx_sub4train));
% idx_train = find(idx_train==1);
% 
% idx_test  = idx_sub==1 .* ...
%     logical(countmember(Idx_sub,Idx_sub4test));
% idx_test = find(idx_test==1);

% subject dependent
idx_train  = idx_sub==1 .* ...
    logical(countmember(idx_trl,Idx_trl4train));
idx_train = find(idx_train==1);

idx_test  = idx_sub==1 .* ...
    logical(countmember(idx_trl,Idx_trl4test));
idx_test = find(idx_test==1);

% train/test index가 잘 뽑혔느지 확인
figure;
plot(idx_sub); hold on;
plot(idx_trl);
plot(idx_train,5*ones(length(idx_train),1),'LineWidth',3)
plot(idx_test,5*ones(length(idx_test),1),'LineWidth',3)
hold off;

% Prepare DB of whole Images
EMGDatasetPath = fullfile('G:\CHS\git\DB_v2\emg_raw');

% EMGData = datastore(EMGDatasetPath,'FileExtensions', '.mat');
f_path = dir(EMGDatasetPath);
f_path = struct2cell(f_path);
f_path = f_path(1,:)';
f_path(1:2) = [];

% train DB
Num_train = length(idx_train);
TrainImg = double(zeros(224,224,3,Num_train));
Trainlab = zeros(Num_train,3);
count = 0;
while(1)
    count = count + 1;
    i= idx_train(count);
    load(fullfile(EMGDatasetPath,f_path{i}))
    TrainImg(:,:,:,count) = imresize(cat(3,emg_raw,emg_raw,emg_raw),[224 224]);
    Trainlab(count,:) = Marker(i,1:3);
    if count == Num_train
        break;
    end
end

% test DB
Num_test = length(idx_test);
TestImg = double(zeros(224,224,3,Num_train));
TestLab = zeros(Num_test,3);
count = 0;
while(1)
    count = count + 1;
    i= idx_test(count);
    load(fullfile(EMGDatasetPath,f_path{i}))
    TestImg(:,:,:,count) = imresize(cat(3,emg_raw,emg_raw,emg_raw),[224 224]);
    TestLab(count,:) = Marker(i,1:3);
    if count == Num_test
        break;
    end
end

% get google net
net = googlenet
lgraph = layerGraph(net);
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
% plot(lgraph)
lgraph = removeLayers(lgraph, {'loss3-classifier','prob','output'});

newLayers = [
%     fullyConnectedLayer(500,'Name','FULLY_500')
    fullyConnectedLayer(100,'Name','FULLY_100')
    fullyConnectedLayer(3,'Name','FULLY_3')
    regressionLayer('Name','Regression Layer')];

lgraph = addLayers(lgraph,newLayers);
lgraph = connectLayers(lgraph,'pool5-drop_7x7_s1','FULLY_100');

figure('Units','normalized','Position',[0.1 0.1 0.8 0.8]);
plot(lgraph)


% training option
options = trainingOptions('sgdm', ...
    'Plots','training-progress',...
    'InitialLearnRate',1e-10...
    );

% train Network
rng('default') % For reproducibility
net = trainNetwork(TrainImg,Trainlab,lgraph,options);

predictedTest= predict(net,TestImg);
predictionError = TestLab - predictedTest;

for i_mktype = 1 : 3
    figure(i_mktype)
plot(TestLab(:,i_mktype));
hold on; 
plot(predictedTest(:,i_mktype))
end

% classify(net,I)

% get network 
% layers = [ ...
%     layersTransfer
% 
%     fullyConnectedLayer(1000)
%     fullyConnectedLayer(3)
%     regressionLayer];
% inputlayer = imageInputLayer([408 4 3],'Name','data');

% layers = [ ...
%     imageInputLayer([408 4],'Name','data')
%     convolution2dLayer([3 1],32)
%     reluLayer
%     convolution2dLayer([3 1],32)
%     reluLayer
%     maxPooling2dLayer([2 1])
%     convolution2dLayer([3 1],64)
%     maxPooling2dLayer([2 1])
%     fullyConnectedLayer(3)
%     regressionLayer];
    
%     
%      convolution2dLayer([4 4],96,'Padding','same')
%      reluLayer
%      convolution2dLayer([4 4],96)
%      reluLayer
% %     convolution2dLayer(12,25,'Padding','same')
%     fullyConnectedLayer(3)
%     regressionLayer]



% % visualization of weight vectors
% w1 = net.Layers(2).Weights;
% 
% imshow(w1(:,:,:,1))
% % Scale and resize the weights for visualization
% w1 = mat2gray(w1);
% w1 = imresize(w1,5);
% 
% % Display a montage of network weights. There are 96 individual sets of
% % weights in the first layer.
% figure
% montage(w1)
% title('First convolutional layer weights')

% extracted features from layer
for idx2show = 4;
featureLayer = 'conv5';
trainingFeatures = activations(net, TrainImg(:,:,:,idx2show), featureLayer ...
    ,'OutputAs','channels');
figure
    imshow(TrainImg(:,:,:,idx2show))
for i=1
    temp = imresize(trainingFeatures(:,:,i),5);
    figure
    imshow(temp)
end

featureLayer = 'pool5';
trainingFeatures = activations(net, TrainImg(:,:,:,idx2show), featureLayer ...
    ,'OutputAs','channels');
figure
    imshow(TrainImg(:,:,:,idx2show))
for i=1
    temp = imresize(trainingFeatures(:,:,i),5);
    figure
    imshow(temp)
end

end
% montage(trainingFeatures(:,:,1:3))