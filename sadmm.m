% Wyczyść przestrzeń roboczą i zamknij wszystkie figury
clear all; close all; clc;

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

% Współczynnik odbicia
reflection_coeff = 0.8;

% Definicja ścian jako odcinków [x1 y1 x2 y2]
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
    
    % Inicjalizacja pola całkowitego (amplituda zespolona)
    total_field = 0;
    
    % --- Bezpośrednia ścieżka (LOS) ---
    los = true;
    for w = 1:size(walls,1)
        x1 = walls(w,1); y1 = walls(w,2);
        x2 = walls(w,3); y2 = walls(w,4);
        result = dwawektory(tx_x, tx_y, rx_x, rx_y, x1, y1, x2, y2);
        if result == 1
            los = false;
            break;
        end
    end
    
    % Obsługa drzwi w ścianie działowej
    if ~los && tx_x < 13.05 && rx_x > 13.05
        y_at_wall = tx_y + (13.05 - tx_x)*(rx_y - tx_y)/(rx_x - tx_x);
        if y_at_wall >= 5 && y_at_wall <= 7
            los = true;
        end
    end
    
    if los
        d = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
        field = sqrt(Pt) * (lambda / (4 * pi * d)) * exp(-1j * 2 * pi * d / lambda);
        total_field = total_field + field;
    end
    
    % --- Odbicia pojedyncze ---
    for w = 1:size(walls,1)
        wall = walls(w,:);
        % Obliczenie obrazu nadajnika względem ściany
        [image_x, image_y] = mirrorPoint(tx_x, tx_y, wall(1), wall(2), wall(3), wall(4));
        
        % Punkt odbicia na ścianie
        [intersect_x, intersect_y] = lineSegmentIntersection(image_x, image_y, rx_x, rx_y, wall(1), wall(2), wall(3), wall(4));
        if isempty(intersect_x)
            continue;
        end
        
        % Sprawdzenie, czy droga odbita nie jest zasłonięta
        obstructed = false;
        % Od nadajnika do punktu odbicia
        for w2 = 1:size(walls,1)
            if w2 == w
                continue;
            end
            result = dwawektory(tx_x, tx_y, intersect_x, intersect_y, walls(w2,1), walls(w2,2), walls(w2,3), walls(w2,4));
            if result == 1
                obstructed = true;
                break;
            end
        end
        % Od punktu odbicia do odbiornika
        if ~obstructed
            for w2 = 1:size(walls,1)
                if w2 == w
                    continue;
                end
                result = dwawektory(rx_x, rx_y, intersect_x, intersect_y, walls(w2,1), walls(w2,2), walls(w2,3), walls(w2,4));
                if result == 1
                    obstructed = true;
                    break;
                end
            end
        end
        
        if obstructed
            continue;
        end
        
        % Obliczenie całkowitej drogi
        d1 = sqrt((tx_x - intersect_x)^2 + (tx_y - intersect_y)^2);
        d2 = sqrt((rx_x - intersect_x)^2 + (rx_y - intersect_y)^2);
        d_total = d1 + d2;
        
        % Obliczenie pola z uwzględnieniem odbicia
        field = reflection_coeff * sqrt(Pt) * (lambda / (4 * pi * d_total)) * exp(-1j * 2 * pi * d_total / lambda);
        total_field = total_field + field;
    end
    
    % Obliczenie mocy odbieranej
    Pr = abs(total_field)^2;
    Pr_matrix(idx) = Pr;
end

% Konwersja mocy do skali dBm
Pr_dBm = 10 * log10(Pr_matrix / 1e-3);

% Tworzenie mapy mocy sygnału
figure;
pcolor(wspX, wspY, Pr_dBm);
shading interp;
colorbar;
title('Mapa mocy odbieranego sygnału (uwzględniono odbicia)');
xlabel('X [m]');
ylabel('Y [m]');
colormap jet;

% --- Funkcje pomocnicze ---

% Funkcja obliczająca punkt lustrzany względem ściany
function [xm, ym] = mirrorPoint(x0, y0, x1, y1, x2, y2)
    % Równanie prostej ściany: Ax + By + C = 0
    A = y2 - y1;
    B = x1 - x2;
    C = x2 * y1 - x1 * y2;
    % Obliczenie punktu lustrzanego
    D = (A * x0 + B * y0 + C) / (A^2 + B^2);
    xm = x0 - 2 * A * D;
    ym = y0 - 2 * B * D;
end

% Funkcja znajdująca punkt przecięcia dwóch odcinków
function [xi, yi] = lineSegmentIntersection(x1, y1, x2, y2, x3, y3, x4, y4)
    % Sprawdzenie przecięcia prostych
    denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4);
    if denom == 0
        xi = []; yi = [];
        return;
    end
    xi = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / denom;
    yi = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / denom;
    % Sprawdzenie, czy punkt leży na obu odcinkach
    onSegment1 = min(x1,x2) <= xi && xi <= max(x1,x2) && min(y1,y2) <= yi && yi <= max(y1,y2);
    onSegment2 = min(x3,x4) <= xi && xi <= max(x3,x4) && min(y3,y4) <= yi && yi <= max(y3,y4);
    if ~(onSegment1 && onSegment2)
        xi = []; yi = [];
    end
end
