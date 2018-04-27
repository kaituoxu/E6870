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
load all.dat;
load p1submit.dat;
load alldata.dat;

a = [-1.2918, 1.41698,0.389452,0.664329,0.531871,0.171098, 0.0946618, 0.27412, 0.474134, 0.574907, 0.303789, 0.542182];
b = [-1.34883, 1.70244, 1.05621, 0.890976, 0.669948, 0.215092, 0.415701, 0.160527, 0.675014, 0.572366, 0.689619, 0.717665];
d = 0;
for i = 1:12
    c = (a(i) - b(i))^2;
    d = d + c;
end
d = sqrt(d);
%{
inFrameCnt = 74;
inDimCnt = 512;
outDimCnt = 26;
numBins = 26;
inFeats = p1fft;
outFeats = zeros(74,26);
samplePeriod = 5*(10^-5);
fmax = 1/ (2*samplePeriod);
melfmax = 1127 * log(1+fmax/700);
interval = melfmax / (numBins + 1);
fre = zeros(outDimCnt + 2,1);
for i = 1:28
    fre(i) = interval * (i - 1);
end

for i = 1:numBins
    for idx = 1:inFrameCnt
        sum = 0;
        for j = 1:2:inDimCnt
            f = (j - 1) / (2 * inDimCnt * samplePeriod);
            melf = 1127 * log(1 + f/700);
            x = sqrt(inFeats(idx,j)^2 + inFeats(idx,j + 1)^2);
            h = 0;
            if(melf >= fre(i) && melf <= fre(i + 1))
                h = (melf - fre(i))/(fre(i + 1) - fre(i));
            elseif(melf >= fre(i + 1) && melf <= fre(i + 2))
                h = (fre(i + 2) -  melf)/(fre(i + 2) - fre(i + 1));
            end
            sum = sum + x * h;
        end
        outFeats(idx,i) = sum;
    end
end
%}
