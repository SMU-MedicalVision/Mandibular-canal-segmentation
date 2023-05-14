function [leftSubPano, rightSubPano, leftSubCanalMask, rightSubCanalMask, sampledVol, sampledCanalMask] = generate_sub(vol, mandible_mask, foramen_cords, canal_mask)

foramen_cords = [foramen_cords(2, :); size(vol, 1) - foramen_cords(1, :); foramen_cords(3:4, :)];
% 参数设置
sampleWidth = 64;
sampleDepth = 96; 
margin = 16;

bboxMask = false(size(vol));
col_center = round(0.5*foramen_cords(2,2) + 0.5*foramen_cords(2,3));%两颏孔中点

[M, N, D] = size(vol);
r1 = min(round(foramen_cords(1, 1)) - margin, M);
d1 = min(round(foramen_cords(3, 1)) + margin, D);
bboxMask(r1:end, col_center+1:end, 1:d1) = 1;

r1 = min(round(foramen_cords(1, 4)) - margin, M);
d1 = min(round(foramen_cords(3, 4)) + margin, D);
bboxMask(r1:end, 1:col_center, 1:d1) = 1;
% bboxMask_zero = false(size(bboxMask));

mandible_mask = mandible_mask & bboxMask;

[Xs, Ys, Zs] = ndgrid(1:M, 1:N, 1:D);%N 维空间中的矩形网格

bwProj = sum(sum(mandible_mask, 1), 3);
bwProj(bwProj < 1) = 1;

xc = sum(sum(Xs .* mandible_mask, 1), 3) ./ bwProj;
% yc = sum(sum(Ys .* mandible_mask, 1), 3) ./ bwProj;
zc = sum(sum(Zs .* mandible_mask, 1), 3) ./ bwProj;

yc = 1:N;
xc = xc(:);
yc = yc(:);
zc = zc(:);

tmp = bwProj(:) > 256;

xc = xc(tmp);
yc = yc(tmp);
zc = zc(tmp);

validBW = (zc >= min(foramen_cords(3, :) - margin)) & (zc <= max(foramen_cords(3, :) + margin));

xc = xc(validBW);
yc = yc(validBW);
zc = zc(validBW);

for iter=1:20
    xc = linearSmooth7(xc);
%     yc = linearSmooth7(yc);
    zc = linearSmooth7(zc);
end

[arclen,~] = arclength(xc, yc, zc, 'spline');%曲线长度

[centerLine,~,~, ~]  = interparc(0:1/arclen:1,xc, yc, zc, 'spline');%插点

leftInd = 0.5*size(centerLine, 1) - margin;
rightInd = 0.5*size(centerLine, 1) + margin;

x = centerLine(:,1);
y = centerLine(:,2);
z = centerLine(:,3);
% [tt, r, s, k,t] = frenet(x',y',z');
[tt, r, s, tt_adf, r_adf, s_adf] = DR(x',y',z');
BB = s;
NN = r;
TT = tt;

sampleXCords = zeros(sampleWidth+1, sampleDepth+1, size(centerLine,1));
sampleYCords = zeros(sampleWidth+1, sampleDepth+1, size(centerLine,1));
sampleZCords = zeros(sampleWidth+1, sampleDepth+1, size(centerLine,1));

[XX, YY] = ndgrid(-sampleWidth/2:sampleWidth/2, -sampleDepth/2:sampleDepth/2);
[m1, n1] = size(XX);
for k=1:size(centerLine, 1)
    rotMatrix = [BB(k, :); NN(k, :); TT(k, :)];
    transCord = centerLine(k, :);
    cck = rotMatrix' * [XX(:)'; YY(:)'; zeros(1, m1*n1)];
    xk = cck(1, :);
    xk = xk(:)  + transCord(1);
    yk = cck(2, :);
    yk = yk(:)  + transCord(2);
    zk = cck(3, :);
    zk = zk(:)  + transCord(3);
    
    sampleXCords(:,:,k) = reshape(xk, [m1, n1]);
    sampleYCords(:,:,k) = reshape(yk, [m1, n1]);
    sampleZCords(:,:,k) = reshape(zk, [m1, n1]);
end

sampleXCords = sampleXCords(1:end-1, 1:end-1, :);
sampleYCords = sampleYCords(1:end-1, 1:end-1, :);
sampleZCords = sampleZCords(1:end-1, 1:end-1, :);

F = griddedInterpolant(Xs, Ys, Zs, double(vol));
sampledVol = F(sampleXCords,sampleYCords, sampleZCords);
sampledVol(isnan(sampledVol)) = 0;
leftSubPano = sampledVol(:, :, 1:leftInd);
rightSubPano = sampledVol(:, :, rightInd:end);

sampledVol1 = mat2gray(sampledVol, [0, 2500]);
figure, montage(sampledVol1(:,:,8:16:end))

F = griddedInterpolant(Xs, Ys, Zs, double(canal_mask));
sampledCanalMask = F(sampleXCords,sampleYCords, sampleZCords);
sampledCanalMask(isnan(sampledCanalMask)) = 0;
leftSubCanalMask = sampledCanalMask(:, :, 1:leftInd);
rightSubCanalMask = sampledCanalMask(:, :, rightInd:end);


