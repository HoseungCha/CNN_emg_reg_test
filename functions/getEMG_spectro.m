function [window_DB,spect_img] = getEMG_spectro(x,winsize,wininc,datawin,dispstatus)

% 파라미터 설정
if isempty(winsize)
    winsize = size(x,1);
end
if isempty(wininc)
    wininc = winsize;
end
if isempty(datawin)
    datawin = ones(winsize,1);
end
if isempty(datawin)
    dispstatus = 0;
end

datasize = size(x,1);
Nsignals = size(x,2);
numwin = floor((datasize - winsize)/wininc)+1;

% allocate spect_imgure memory
% spect_img = zeros(numwin,4);
window_DB = cell(numwin,1);
spect_img = cell(numwin,1);
if dispstatus
    h = waitbar(0,'Computing Waveform Length spect_imgures...');
end



st = 1;
en = winsize;
for i = 1:numwin
   if dispstatus
       waitbar(i/numwin);
   end
   curwin = x(st:en,:).*repmat(datawin,1,Nsignals);
   
   % 2차항으로 EMG feature agumentation
%     idx2agu = permn(1:4,2);
%     curwin_agu = zeros(length(curwin),size(idx2agu,1));
%     for i_aug = 1 : size(idx2agu,1)
%         curwin_agu(:,i_aug) = curwin(:,idx2agu(i_aug,1)).*curwin(:,idx2agu(i_aug,2));
%     end
%     
%     curwin = [curwin,curwin_agu];
        
   % spectrum to image (RGB)
   Nx = length(curwin);
   N_ch = size(curwin,2);
    nsc = floor(Nx/3);
    nov = floor((nsc/10)*7);
    nff = max(512,2^nextpow2(nsc));
%     nff = 2^nextpow2(nsc);
    ps_ = cell(1,N_ch);
%    tic
   window_DB{i} = curwin;
   for i_ch = 1 : N_ch
       [~,~,~,ps] = spectrogram(curwin(:,i_ch)',hamming(nsc),nov,nff,2048);
       ps_{i_ch} = 10*log10(abs(ps(1:185,:))+eps); % 736 Hz 까지 컷
       %We kept only the first 95 points in
% frequency of the spectrogram(0?736.54Hz) because themajority
% of the sEMG energy was observed within frequency range from 0
% to ∼700 Hz (Zhai et al., 2016 참고논문: Self-Recalibrating Surface EMG Pattern Recognition for Neuroprosthesis Control Based on Convolutional Neural Network
   end
   temp = mat2im(cell2mat(ps_),parula(numel(ps)));
   spect_img{i} = temp;
%    spect_img{i} = imresize(temp,[224 224]);
   
%    toc;
%    disp(i);
%    imshow(spect_img{i})
%    imshow(B)

   
   st = st + wininc;
   en = en + wininc;
end


if dispstatus
    close(h)
end


function c=a2c(a,p,cp)
%Function A2C: Computation of cepstral coeficients from AR coeficients.
%
%Usage: c=a2c(a,p,cp);
%   a   - vector of AR coefficients ( without a[0] = 1 )
%   p   - order of AR  model ( number of coefficients without a[0] )
%   c   - vector of cepstral coefficients (without c[0] )
%   cp  - order of cepstral model ( number of coefficients without c[0] )

%                              Made by PP
%                             CVUT FEL K331
%                           Last change 11-02-99

for n=1:cp,

  sum=0;

  if n<p+1,
    for k=1:n-1,
      sum=sum+(n-k)*c(n-k)*a(k);
    end;
    c(n)=-a(n)-sum/n;
  else
    for k=1:p,
      sum=sum+(n-k)*c(n-k)*a(k);
    end;
    c(n)=-sum/n;
  end;

end;