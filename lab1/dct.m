clear;
close all

load p1in.dat
load p1win.dat
load window.dat
load p1fft.dat
load fft.dat
load p1bin.dat
load p1all.dat
load melbin.dat

inFrameCnt = 74;
inDimCnt = 26;
outDimCnt = 12;
numCoeffs = 12;

inFeats = p1bin;
outFeats = zeros(74,12);

for i = 1:numCoeffs
    for idx = 1:inFrameCnt
        sum = 0;
        for j = 1:inDimCnt
            temp = inFeats(idx,j)*cos(pi*i*(j - 1 + 0.5)/26);
            sum = sum + temp;
        end
        outFeats(idx,i) = sqrt(2/26)*sum;
    end
end
