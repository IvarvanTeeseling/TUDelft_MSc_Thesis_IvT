function DIC = dicDataLoad(dicF, dicFr)

% Initiate a waitbar
f = waitbar(1/length(dicFr), 'Uploading MATLAB data. Please wait...');

% Read, extract and store the DIC .mat files selected for processing
for i = 1:length(dicFr)
    
    % Load DIC data
    dic(i) = load(fullfile(dicF(dicFr(i)).folder, dicF(dicFr(i)).name));
    
    % Pre-allocate memory
    if i == 1
        m = size(dic(i).X);
        n = length(dicF);
        DIC.X = zeros([m n]);
        DIC.Y = zeros([m n]);
        DIC.W = zeros([m n]);
        DIC.ey = zeros([m n]);
        DIC.ex = zeros([m n]);
        DIC.e1 = zeros([m n]);
        DIC.e2 = zeros([m n]);
        DIC.s = zeros([m n]);
    end
    
    % Store data in matrix form for easier handeling
    DIC.X(:,:,i) = dic(i).X;
    DIC.Y(:,:,i) = dic(i).Y;
    DIC.W(:,:,i) = dic(i).W;
    DIC.ey(:,:,i) = dic(i).eyy;
    DIC.ex(:,:,i) = dic(i).exx;
    DIC.e1(:,:,i) = dic(i).e1;
    DIC.e2(:,:,i) = dic(i).e2;
    DIC.sig(:,:,i) = dic(i).sigma;
    
    % Update waitbar
    waitbar(i/length(dicFr), f);
end

% Close waitbar
close(f);

end