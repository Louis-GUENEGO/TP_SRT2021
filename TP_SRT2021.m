%% Initialisation 
clear all;
close all;
clc

k = 8*188; % nombre de bits transmis
nb = 2; % nombre de bits par symbole
Ds = 250e3; % debit symbole
Ts = 1/Ds;
Fse = 8; % différence entre le débit symbole et la fréquence d'échantillonage
Fe = Fse * Ds; % fréquence d'échantillonage
Te = 1/Fe;
span = 8; % largeur du filtre
phi = 2*pi/360  * 0;
DeltaF = 100000;

g = comm.RaisedCosineTransmitFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',Fse); % filtre de mise en forme
ga = comm.RaisedCosineReceiveFilter('RolloffFactor',0.35,'FilterSpanInSymbols',span,'InputSamplesPerSymbol',Fse,'DecimationFactor',Fse); % filtre adapte

calculError = comm.ErrorRate('ReceiveDelay', 2*span, 'ComputationDelay', 2*span+2);

eb_n0_dB = 100; % rapport signal sur bruit en dB
eb_n0 = 10.^(eb_n0_dB/10);
sig = 1; % variance du signal
Eg = 1; % energie du filtre

Pb = qfunc ( sqrt (2*eb_n0 ) );
sigma = sig * Eg ./ (nb * eb_n0);

for foo = 1:length(eb_n0)
    nbr_err = 0;
    calculError.reset;
    g.reset()
    ga.reset()
    d0 = 1;
    r0 = 1;
    % while nbr_err < 100
        
        %% Emetteur

        bk = double ((randn(1,k)>0).'); % generation des bits
        qpsk_mod = comm.QPSKModulator('PhaseOffset',0,'BitInput',1);
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
        Fphi = exp(1i*phi + 1i*2*pi*DeltaF*x/Fe );
        sl_in = conv(h,sl.*Fphi);
        bruit = sqrt(sigma(foo)/2)* ( randn(size(sl_in)) + randn(size(sl_in)) *1i ) ; %bruit
        sl_in = sl_in + bruit;
        
        %% Recepteur

        rn_se = step(ga,sl_in);
        rn = rn_se; 
        dn = zeros(size(rn));
        dn(1) = rn(1) * conj(r0);
        for n = 2:size(an)
           dn(n) = rn(n) * conj(rn(n-1));
        end
        r0 = rn(end);
        qpsk_demod = comm.QPSKDemodulator('PhaseOffset',0,'BitOutput',1);
        anr = qpsk_demod.step(dn);
        bkr = anr;

        error_cnt = calculError.step( bk, bkr );
        nbr_err = error_cnt(2);
        
    % end
    
    TEB(foo) = error_cnt(1)
        
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

figure;
plot(rn, '*');
title('Constellation de rn')
xlabel('I') 
ylabel('Q') 
xlim([-1.5 1.5])
ylim([-1.5 1.5])

% figure;
% plot(dn, '*');
% title('Constellation de dn')
% xlabel('I') 
% ylabel('Q') 
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

figure;
Nfft = 512;
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


% figure('Name', 'évolution du TEB en fonction de Eb/N0 en dB');
% semilogy(eb_n0_dB, TEB, eb_n0_dB, Pb);
% ylabel('TEB');
% xlabel('Eb/N0 en dB');
% title('évolution du TEB en fonction de Eb/N0 en dB');
% legend({'Simulé', 'théorique'},'Location','southwest');
% grid on



















