close all; clear all; clc;

%LAB2
% Parametry pomieszczenia
room_width = 20;   % szerokość pomieszczenia (m)
room_height = 12;  % wysokość pomieszczenia (m)

% Współrzędne nadajnika
tx_x = 5.05;
tx_y = 6.55;

% Parametry sygnału
Pt = 1e-3;        % Moc nadajnika (1 mW)
f = 2.4e9;        % Częstotliwość (2.4 GHz)
c = 3e8;          % Prędkość światła (m/s)
lambda = c / f;   % Długość fali (m)

% Definicja ścian jako odcinków(wektorów?) [x1 y1 x2 y2]
walls = [
    0, 0, 0, 12;          % Lewa ściana
    20, 0, 20, 12;        % Prawa ściana
    0, 0, 20, 0;          % Dolna ściana
    0, 12, 20, 12;        % Górna ściana
    13.05, 0, 13.05, 5;   % Ściana działowa dolna
    13.05, 7, 13.05, 12;  % Ściana działowa górna
];

% Siatka punktów odbiornika
[wspX, wspY] = meshgrid(0.1:0.1:20, 0.1:0.1:12);

% Inicjalizacja macierzy mocy odbieranej
Pr_matrix = zeros(size(wspX));

% Pętla po wszystkich pozycjach odbiornika
for idx = 1:numel(wspX)
    rx_x = wspX(idx);
    rx_y = wspY(idx);
    
    % Sprawdzenie warunku Line-of-Sight
    los = true;
    for w = 1:size(walls,1)
        % Współrzędne odcinka ściany
        x1 = walls(w,1); y1 = walls(w,2);
        x2 = walls(w,3); y2 = walls(w,4);
        
        % Sprawdzenie przecięcia z użyciem funkcji dwawektory
        result = dwawektory(tx_x, tx_y, rx_x, rx_y, x1, y1, x2, y2);
        if result == 1
            los = false;
            break;
        end
    end
    
    % Obsługa drzwi w ścianie działowej
    %{
    if ~los && tx_x < 13.05 && rx_x > 13.05
        y_at_wall = tx_y + (13.05 - tx_x)*(rx_y - tx_y)/(rx_x - tx_x);
        if y_at_wall >= 5 && y_at_wall <= 7
            los = true;
        end
    end
    %}
    
    if los
        % Obliczenie odległości
        d = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
        % Obliczenie mocy odbieranej (FSL)
        Pr = Pt * (lambda / (4 * pi * d))^2;
        Pr_matrix(idx) = Pr;
    else
        Pr_matrix(idx) = NaN; % Brak sygnału
    end
end

% Konwersja mocy do skali dBm
Pr_dBm = 10 * log10(Pr_matrix / 1e-3);

% Tworzenie mapy mocy sygnału
%figure;
pcolor(wspX, wspY, Pr_dBm);
shading interp;
colorbar;
title('Mapa mocy odbieranego sygnału (tylko LOS)');
xlabel('X [m]');
ylabel('Y [m]');
colormap jet;
