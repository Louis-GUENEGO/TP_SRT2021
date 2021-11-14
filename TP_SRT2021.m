%% Initialisation
clear all;
close all;
clc

Nfft = 1024;

k = 8*188*10; % nombre de bits transmis
nb = 2; % nombre de bits par symbole
Ds = 250e3; % debit symbole
Ts = 1/Ds;
Fse = 4; % différence entre le débit symbole et la fréquence d'échantillonage
Fe = Fse * Ds; % fréquence d'échantillonage
Te = 1/Fe;
span = 8; % largeur du filtre
phi = 2*pi/360  * 0;
DeltaF = 100000;
tau = 0;%Fe*Ts .* [0.1, 0.2, 0.3,0.33,0.35,0.4];

qpsk_mod = comm.QPSKModulator('PhaseOffset',0,'BitInput',1);
qpsk_demod = comm.QPSKDemodulator('PhaseOffset',0,'BitOutput',1);

g = comm.RaisedCosineTransmitFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',Fse); % filtre de mise en forme
ga = comm.RaisedCosineReceiveFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'InputSamplesPerSymbol',Fse,'DecimationFactor',1); % filtre adapte

calculError = comm.ErrorRate('ReceiveDelay', 2*span, 'ComputationDelay', 2*span+2);

freqComp = comm.CoarseFrequencyCompensator('Modulation','QPSK', 'FrequencyResolution', 1000, 'SampleRate', Fe);
freqSync = comm.CarrierSynchronizer('Modulation', 'QPSK', 'ModulationPhaseOffset', 'Custom', 'CustomPhaseOffset', 0, 'SamplesPerSymbol', 1, 'DampingFactor', 0.707, 'NormalizedLoopBandwidth', 0.005);
delay = dsp.VariableFractionalDelay;

syncTemp = comm.SymbolSynchronizer('SamplesPerSymbol',Fse,'TimingErrorDetector','Gardner (non-data-aided)', 'DampingFactor', 0.707, 'NormalizedLoopBandwidth', 0.005, 'DetectorGain', 1);

eb_n0_dB = 100;%1:1:10; % rapport signal sur bruit en dB
eb_n0 = 10.^(eb_n0_dB/10);
sig = 1; % variance du signal
Eg = 1; % energie du filtre

Pb = qfunc ( sqrt (2*eb_n0 ) );
sigma = sig * Eg ./ (nb * eb_n0);

for toto = 1:length(tau)
    for bar = 1:length(DeltaF)
        for foo = 1:length(eb_n0)
            nbr_err = 0;
            calculError.reset;
            g.reset()
            ga.reset()
            syncTemp.reset()
            d0 = 1;
            r0 = 1;
            buffer_bk  = [];
            buffer_bkr = [];
            %p = 1;
            %while nbr_err < 100
                for p = 1:3
                
                
                %% Emetteur
                
                bk = double ((randn(1,k)>0).'); % generation des bits
                
                an = qpsk_mod.step(bk);
                
                dn = zeros(size(an));
                dn(1) = an(1) * d0;
                for n = 2:size(an)
                    dn(n) = an(n)*dn(n-1);
                end
                d0 = dn(end);
                dn_se = dn;
                
                sl = step(g,dn_se);
                
                %% Canal
                
                h = 1;
                x = (1:size(sl))';
                Fphi = exp(1i*phi + 1i*2*pi*DeltaF(bar)*x/Fe );
                sl_in = conv(h,delay.step(sl, tau(toto)).*Fphi);
                bruit = sqrt(sigma(foo)/2)* ( randn(size(sl_in)) + randn(size(sl_in)) *1i ) ; %bruit
                sl_in = sl_in + bruit;
                
                [yl,estFreqOffset] = freqComp.step(sl_in);
                
                
                
                %% Recepteur
                
                rn_se = step(ga,yl);
                rn_ = syncTemp.step(rn_se);
                rn = freqSync.step(rn_);
                dnr = zeros(size(rn));
                dnr(1) = rn(1) * conj(r0);
                for n = 2:size(rn)
                    dnr(n) = rn(n) * conj(rn(n-1));
                end
                r0 = rn(end);
                figure(1);
                plot(dnr, '*');
                xlim([-5 5]);
                ylim([-5 5]);
                title(['d_n (reception), iteration : ',num2str(p)]);
                xlabel('I');
                ylabel('Q');
                drawnow limitrate

                anr = qpsk_demod.step(dnr);
                bkr = anr;
                
                if p > 2 
                    error_cnt = calculError.step( bk, bkr );
                    nbr_err = error_cnt(2);
                end
                
                end
            %end
            
            TEB(foo, toto) = error_cnt(1)
            
        end
        
    end
    
    figure;
    xlim([-5 5])
    ylim([-5 5])
    plot(rn, '*');
    hold on
    plot(an, 'r*');
    title('Constellation de dnr et an')
    legend({'dnr', 'an'},'location','southwest');
    xlabel('I')
    ylabel('Q')
    
end




%% Figures

%constellation(qpsk_mod);
%constellation(qpsk_demod);


% figure
% plot(an, '*');
% title('Constellation de an')
% xlabel('I')
% ylabel('Q')
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

% figure;
% plot(rn(), '*');
% title('Constellation de rn')
% xlabel('I')
% ylabel('Q')
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

% figure;
% plot(dn, '*');
% title('Constellation de dn')
% xlabel('I')
% ylabel('Q')
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

figure;
pwelch(sl_in, hanning(Nfft), 0, Nfft, Fe, 'centered')
title('Densité de puissance de sl')

[pw, fpw] = pwelch(sl_in, hanning(Nfft), 0, Nfft, Fe, 'centered'); % diagramme de welch avec Nfft points et pas de recouvrement entre deux fenetres
[Ga, f] = ga.freqz(Nfft,'whole',Fe); %
figure('Name', 'Périodogramme de Welch et DSP théorique de sl(t)');
plot(fpw,10*log10(abs(pw)),f-1000000, 10*log10(fftshift(abs(Ga)))-67,'g');
xlabel('fréquence en Hz');
ylabel('dB / (rad/sample)');
legend({'Périodogramme de Welch','DSP théorique'},'Location','southwest');
title('Périodogramme de Welch et DSP théorique de sl(t)');


figure('Name', 'évolution du TEB en fonction de Eb/N0 en dB, avec plusieur valeurs de \Delta_f');
semilogy(eb_n0_dB, Pb, eb_n0_dB, TEB);
ylabel('TEB');
xlabel('Eb/N0 en dB');
title('évolution du TEB en fonction de Eb/N0 en dB');
legend({'théorique', '\tau = 0,1Ts', '\tau = 0,2Ts', '\tau = 0,3Ts', '\tau = 0,33Ts', '\tau = 0,35Ts', '\tau = 0,4Ts'},'location','southwest');
grid on

% figure
% pwelch(sl_in.^4, hanning(Nfft), 0, Nfft, Fe, 'centered')
















