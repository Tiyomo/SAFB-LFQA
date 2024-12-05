Feature_table = readtable('NBUfeat.xlsx');   
Dmos = table2array(Feature_table(:, end)); % ��ǩ
Feature = table2array(Feature_table(:, 1:end-1)); % ����

c = 8.46; gamma = 0.013 ;   

numitr = 1000; 
lower = -1; 
upper = 1; 

[Feature_norm, MAX, MIN] = normalization(Feature, lower, upper);

SROCC = zeros(1,numitr);
KROCC = zeros(1,numitr);
RMSE = zeros(1,numitr);
PLCC = zeros(1,numitr);

for trialid = 1:numitr    
    train_ratio = 0.8; % ѵ��������
    num_data = size(Feature,1); % ��������
    num_train = round(train_ratio*num_data); % ѵ��������
    rand_index = randperm(num_data); % �����������
    train_ind = rand_index(1:num_train); % ѵ��������
    test_ind = rand_index(num_train+1:end); % ���Լ�����
    traindata = Feature_norm(train_ind, :);
    trainlabel = Dmos(train_ind, :);
    testdata = Feature_norm(test_ind, :);
    testlabel = Dmos(test_ind, :);
    
    model = libsvmtrain(trainlabel, traindata, ['-s 3 -t 2 -c ',num2str(c),' -g ',num2str(gamma)]);

    [quality, ~, ~] = libsvmpredict(testlabel, testdata, model);

    [srocc,krocc,plcc,rmse] = verify_performance(testlabel,quality);

    SROCC(trialid) = srocc;
    KROCC(trialid) = krocc;
    RMSE(trialid) = rmse;
    PLCC(trialid) = plcc; 
end

fprintf('SROCC: %.4f\n', median(SROCC));
fprintf('RMSE: %.4f\n', median(RMSE));
fprintf('PLCC: %.4f\n', median(PLCC));
fprintf('KROCC: %.4f\n', median(KROCC));