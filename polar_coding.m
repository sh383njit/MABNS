
clear all;
clc;

varn = [1 0.870964 0.758578 0.724436 0.660693 0.630957 0.60256]; 

K = 64;
N = 128;

trans=100;
L = 4; % list decoder length

for j=1:length(varn)
    chan = comm.AWGNChannel('NoiseMethod','Variance','Variance',varn(j));
%     bpskMod = comm.BPSKModulator;
%     bpskDemod = comm.BPSKDemodulator('DecisionMethod', 'Approximate log-likelihood ratio','Variance',varn(j));
    err=0;
    for i=1:trans
        msg = randi([0 1],K,1,'int8');
        enc = nrPolarEncode(msg,N);

%         mod = bpskMod(enc);
%         rSig = chan(mod);
%         rxLLR = bpskDemod(rSig);
        
        y = chan(double(enc));
        rxLLR=2*y/varn(j);

        rxBits = nrPolarDecode(rxLLR,K,N,L);

        numBitErrs = biterr(rxBits,msg);

        err=err+numBitErrs;
    end
    BER(j)=err/(trans*N)
end

BER


