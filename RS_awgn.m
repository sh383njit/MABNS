
clear all;
clc;

M = 2;           % Modulation order
bps = log2(M);   % Bits per symbol
N = 7;          
K = 3;           
m=30; %msg length, cw length is N/K*m
%of AB comparison, m=56, K=3, N=7

pskModulator = comm.PSKModulator('ModulationOrder',M,'BitInput',true);
pskDemodulator = comm.PSKDemodulator('ModulationOrder',M,'BitOutput',true);
awgnChannel = comm.AWGNChannel('BitsPerSymbol',bps);
errorRate = comm.ErrorRate;

rsEncoder = comm.RSEncoder('CodewordLength',N,'MessageLength',K);
rsDecoder = comm.RSDecoder('CodewordLength',N,'MessageLength',K);

ebnoVec = [0 0.6 1.2 1.4 1.8]'; %(3:0.5:8)';
errorStats = zeros(length(ebnoVec),3);
ber=zeros(1,length(ebnoVec));

trans=10000;
for trans=1:trans
    for i = 1:length(ebnoVec)
        awgnChannel.EbNo = ebnoVec(i);
        reset(errorRate)
        while errorStats(i,2) < 100 && errorStats(i,3) < 1e7
            data = ones(m,1); %randi([0 1],m,1);                 % msg vector
%             size(data)
            
            cw = rsEncoder(data);                  % codeword of length 
%             size(cw)
            
            modData = pskModulator(cw);            % Modulate
            rxSig = awgnChannel(modData);               % Pass signal through AWGN
            rxData = pskDemodulator(rxSig);             % Demodulate
            decData = rsDecoder(rxData);                % RS decode
            errorStats(i,:) = errorRate(data,decData);  % Collect error statistics
       end
    end

%     berCurveFit = berfit(ebnoVec,errorStats(:,1));
%     berNoCoding = berawgn(ebnoVec,'psk',8,'nondiff');

    ber=ber+errorStats(:,1)';
end

ber=ber/trans

% semilogy(ebnoVec,errorStats(:,1),'b*', ...
% ebnoVec,berCurveFit,'c-',ebnoVec,berNoCoding,'r')
% ylabel('BER')
% xlabel('Eb/No (dB)')
% legend('Data','Curve Fit','No Coding')
% grid

