function [Ps] = Ps(sig,f,nfig,splt)

if size(sig,1) > size(sig,2)
    l = size(sig,1);
else
    l = size(sig,1);
end

S = fft(sig);
Ps = S.*conj(S)/l;
Ps = Ps./max(Ps);

figure(nfig);
% fig = figure(nfig);
% set(fig,'Position',[0 0 720 360])
if splt == 1
    subplot(1,2,1)
else
    plot(f,Ps)
    hold on
    grid on
    ylabel('|s(f)|')
    xlabel('Frequenza [Hz]')
    axis([-0.2 max(f)/2 -0.1 1.1])
    if splt == 1
        subplot(1,2,2)
        plot(sig)
        grid on
        hold on
        ylabel('s(t)')
        xlabel('Tempo acquisizione [s]')
    end
end

