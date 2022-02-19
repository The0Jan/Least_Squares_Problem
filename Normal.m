x = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
y = [-32.9591 -20.7011 -12.6986 -5.1508 -1.6893 0.1266 0.0743 -0.8709 -1.7371 -3.9952 -4.8987];

% Podajemy szerokość macierzy A, czyli wartości większe o jeden niż
% wykorzystywane potęgi wielomianów.
power = [2 3 5 7 9];
power2 = [11 13 14 15];
x_probe = -5.5:0.1:5.5;


% Jak wygląda i jak się dopasowuje Normal
y_outcome = perform_test(x, y, power, x_probe, @normal);
draw_graph(x, y, x_probe, y_outcome, power, "normalny")

% Norma residuum dla Normal
[r_outcome, rw_outcome] = residuum_count(x, y, power, @normal);
draw_graph_norm(x, r_outcome, power, "normalny");


% Jak wygląda i jak się dopasowuje QR
y2_outcome = perform_test(x, y, power, x_probe, @qr);
draw_graph(x, y, x_probe, y2_outcome, power, "liniowych z R")

% Norma residuum dla QR
[r2_outcome, rw2_outcome] = residuum_count(x, y, power, @qr);
draw_graph_norm(x, r2_outcome, power, "liniowych z R");

% Rysowanie norm residuum od stopnia macierzy dla obu
draw_grap_2r(rw_outcome, rw2_outcome, power, "Normalny", "Liniowy z R")


% Funkcja do obliczania normy residuum (zarówno dla od próbki i stopnia
% wielomianu jak i tylko od stopnia wielomianu)
function [out, whole]= residuum_count(x, y, power, method)
    out = zeros(length(power), length(x));
    whole = zeros(length(power),1);
    for i = 1:length(power)
        a_vertices = method(x, y,power(i));
        for j = 1:length(x)
            out(i,j) = sqrt((y(j)-count_y(a_vertices, x(j)))^2);
            whole(i) = whole(i) + (y(j)-count_y(a_vertices, x(j)))^2;
        end
        whole(i) = sqrt(whole(i));      
    end
end

% Przeprowadzanie rozkładu QR
function out = qr(x, y, power)
    % Budowanie macierzy A
    A = zeros(length(x), power);
     for k = 1: length(x)
        for i = 1: power
            A(k,i) = x(k)^(i-1);
        end
     end
     % Przeprowadzanie algorytmu Grama-Schmidta
     Q = zeros(length(x), power);
     R = zeros(power, power);
     d = zeros(power);
     for i =1:power
        Q(:,i) = A(:,i);
        R(i,i) = 1;
        d(i) = transpose(Q(:,i))*Q(:,i);

        for j = i+1:power
            R(i,j) = (transpose(Q(:,i))*A(:,j))/d(i);
            A(:,j) = A(:,j) - R(i,j)*(Q(:,i));
        end
     end
     % Normowanie rozkładów
     for i=1:power
        dd = norm(Q(:,i));
        Q(:,i) = Q(:,i)/dd;
        R(i,i:power) = R(i,i:power)*dd;
     end
     Qy = transpose(Q)*transpose(y);
     out = my_solve(R,Qy);   
end


% Obliczanie za pomocą układu równań normalnych.
function out = normal(x, y, power)
    % Budowanie macierzy A
    A = zeros(length(x), power);
     for k = 1: length(x)
        for i = 1: power
            A(k,i) = x(k)^(i-1);
        end
     end
    
    AtA = transpose(A)*A;
    AtY = transpose(A)*transpose(y);
    out = my_solve(AtA, AtY);
end


% Obliczanie wektora Y
function out = count_y(a_vertix, x)
    out = 0;

    for i = 1:length(a_vertix)
        out = out + a_vertix(i)*(x^(i-1));
    end
end


% Wykonaj test dla odpowiednich rzędów wielomianów
function out = perform_test(x, y, power, tested_x, method)
        out = zeros(length(power), length(tested_x));
    for i = 1:length(power)
        a_vertices = method(x, y,power(i));

        for j = 1:length(tested_x)
            out(i,j) = count_y(a_vertices, tested_x(j));
            
        end
    end
end

% Porównywanie norm residuum od potęgi dla obu metod
function draw_grap_2r(rw1_outcome, rw2_outcome, power, method1, method2)
    figure
    legends = strings(2,1);
    power = power -1;
    semilogy(power,rw1_outcome,"--");
    legends(1) =  method1 ;
    hold on

    semilogy(power,rw2_outcome,'-.');
    legends(2) =  method2 ;
    hold off

    grid on
    lgd = legend( legends(:), "FontSize", 14);
    
    title("Zależność normy residuum od stopnia wielomianu. ", "FontSize", 18);
    xlabel('Stopień wielomianu', "FontSize", 16);
    ylabel('Norma residuum', "FontSize", 16);
end

% Rysowanie jak funkcja przechodzi przez punkty
function draw_graph(x_org, y_org, tested_x, outcome_y, power, name)
    figure
    legends = strings(length(power),1);

    for i= 1:length(power)
        plot(tested_x,outcome_y(i,:));
        legends(i) = "Stopień wielomianu: " + num2str(power(i)-1) ;
        hold on
    end

    legends(length(power)+1) = "Próbki pomiarowe";
    plot(x_org,y_org,'*');

    grid on
    lgd = legend( legends(:), "FontSize", 14);
    
    title("Zależność dopasowania funkcji f(x) od stopnia wielomianu dla: " + name, "FontSize", 18);
    xlabel('X', "FontSize", 16);
    ylabel('Y', "FontSize", 16);
end


% Rysowanie grafy dla błędu rozwiązania 
function draw_graph_norm(x_org, outcome_y, power, name)
    figure
    legends = strings(length(power),1);

    for i= 1:length(power)
        semilogy(x_org,outcome_y(i,:));
        legends(i) = "Stopień wielomianu: " + num2str(power(i)-1) ;
        hold on
    end

    grid on
    lgd = legend( legends(:), "FontSize", 14);
    
    title("Zależność normy residuum od stopnia wielomianu dla: " + name, "FontSize", 18);
    xlabel('Wartość X', "FontSize", 16);
    ylabel('Wielkość normy residuum', "FontSize", 16);
end


% Algorytm Gaussa
function out = my_solve(A, B)
    matrix = [A B];
    len = length(matrix);
    for j = 1:len-2
        [number, position] = max(abs(matrix(j:len -1,j)));
        matrix([j position+j-1], :) = matrix([position+j-1 j],:);
        for i = j+1:len-1
            m = matrix(i, j)/matrix(j, j);
            matrix(i, :) = matrix(i, :) - matrix(j, :)*m;
        end
    end    
    N = length(matrix) - 1;
    out = zeros(N,1);
    out(N) = matrix(N, N+1)/matrix(N,N);
    for i = N -1: -1: 1
        out(i) = (matrix(i, N+1) - matrix(i,i+1:N)*out(i+1:N))/matrix(i,i);
    end
end