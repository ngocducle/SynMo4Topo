%Here's the MATLAB conversion of the provided Python code for calculating and plotting edge states in a heterostructure:
% Minh's code converted to MATLAB
%```matlab
% Heterostructure Edge State Calculation

% Parameters
epsilon = 1e-3;
NE = 501;
epsilonE = 1e-4;

% Hamiltonian function
function H = Hamiltonian(k, delta, side)
    % Parameters of upper layer (InP)
    omega0 = 0.27815; U0 = 0.02232; v0 = 0.37602;
    
    % Parameters of lower layer (Sb3S2) & Interlayer coupling
    if strcmp(side, 'right')
        omegai = 0.27999; Ui = 0.02621; vi = 0.40590;
        Vi = 0.04511;
    elseif strcmp(side, 'left')
        omegai = 0.28010; Ui = 0.02142; vi = 0.40454;
        Vi = 0.04728;
    end
    
    V = Vi * exp(-1i * delta * pi);
    Vc = conj(V);
    
    H = [omega0+v0*k, U0, V, 0;
         U0, omega0-v0*k, 0, Vc;
         Vc, 0, omegai+vi*k, Ui;
         0, V, Ui, omegai-vi*k];
end

% FindEdgeStates function
function edgeStates = FindEdgeStates(delta, EE)
    function GH = PlaneWaves(Q, E, side)
        % Parameters of upper layer (InP)
        omega0 = 0.27815; U0 = 0.02232; v0 = 0.37602;
        
        % Parameters of lower layer (Sb3S2) & Interlayer coupling
        if strcmp(side, 'right')
            omegai = 0.27999; Ui = 0.02621; vi = 0.40590;
            Vi = 0.04511;
        elseif strcmp(side, 'left')
            omegai = 0.28010; Ui = 0.02142; vi = 0.40454;
            Vi = 0.04728;
        end
        
        V = Vi * exp(-1i * Q * pi);
        Vc = conj(V);
        
        QH = [omega0-E, U0, V, 0;
              U0, omega0-E, 0, Vc;
              Vc, 0, omegai-E, Ui;
              0, V, Ui, omegai-E];
        
        QS = diag([-v0, v0, -vi, vi]);
        
        [vecs, vals] = eig(QH, QS);
        vals = diag(vals);
        args = angle(vals * exp(1e-9 * 1i));
        args(args < 0) = args(args < 0) + 2*pi;
        [~, idx] = sort(args);
        GH = vecs(:, idx);
        GH = GH ./ vecnorm(GH);
    end

    Wa = arrayfun(@(E) PlaneWaves(delta, E, 'right'), EE, 'UniformOutput', false);
    Wa = cat(3, Wa{:});
    Wb = arrayfun(@(E) PlaneWaves(delta, E, 'left'), EE, 'UniformOutput', false);
    Wb = cat(3, Wb{:});
    
    WW = cat(2, Wb(:,1:2,:), -Wa(:,3:4,:));
    S = abs(arrayfun(@(i) det(WW(:,:,i)), 1:size(WW,3)));
    
    edgeStates = EE(find(diff(sign(diff(S))) == 2) + 1);
end

% Main code
Nq = 101;
qq = linspace(0.45, 0.55, Nq);
kk = linspace(-0.12, 0.12, 251);

bulk1 = zeros(1, Nq);
bulk2 = zeros(1, Nq);
EdgeStates = [];
Qedge = [];

for nq = 1:Nq
    q = qq(nq);
    EB1 = arrayfun(@(k) eig(Hamiltonian(k, q, 'left')), kk, 'UniformOutput', false);
    EB1 = cat(2, EB1{:});
    EB2 = arrayfun(@(k) eig(Hamiltonian(k, q, 'right')), kk, 'UniformOutput', false);
    EB2 = cat(2, EB2{:});
    
    Emax = min(min(EB1(2,:)), min(EB2(2,:)));
    bulk1(nq) = Emax;
    Emin = max(max(EB1(1,:)), max(EB2(1,:)));
    bulk2(nq) = Emin;
    
    EE = linspace(Emin-epsilonE, Emax+epsilonE, NE);
    edge = FindEdgeStates(q, EE);
    
    EdgeStates = [EdgeStates; edge'];
    Qedge = [Qedge; repmat(q, length(edge), 1)];
end

% Plotting
figure('Position', [100, 100, 900, 1200]);
plot(Qedge, EdgeStates, 'Color', [1 0.5 0], 'LineWidth', 4);
hold on;
plot(qq, bulk1, 'Color', [0.4 0.4 0.4], 'LineWidth', 4);
plot(qq, bulk2, 'Color', [0.4 0.4 0.4], 'LineWidth', 4);
fill([qq fliplr(qq)], [bulk1 ones(1,Nq)*0.234], [0.7 0.7 0.7]);
fill([qq fliplr(qq)], [bulk2 ones(1,Nq)*0.221], [0.7 0.7 0.7]);

xlabel('$vt/\Lambda$', 'Interpreter', 'latex', 'FontSize', 32);
ylabel('$\omega$ ($2\pi c/\Lambda$)', 'Interpreter', 'latex', 'FontSize', 32);
xlim([0.45 0.55]);
ylim([0.2225 0.2325]);
xticks([0.45 0.5 0.55]);
yticks(0.224 + (0:4)*0.002);
set(gca, 'FontSize', 32, 'TickLength', [0.02 0.02], 'LineWidth', 4);
legend('edge', 'bulk', 'FontSize', 32);
box on;

saveas(gcf, 'chiral_edge.png');
```

This MATLAB code replicates the functionality of the original Python script for calculating and plotting edge states in a heterostructure[1]. Key changes and considerations include:

1. Function definitions are adjusted to MATLAB syntax.
2. Array operations are modified to use MATLAB's built-in functions like `arrayfun` and `cat`.
3. The `eig` function in MATLAB returns eigenvalues and eigenvectors separately, so the output is adjusted accordingly.
4. Plotting commands are adapted to MATLAB's syntax, including figure size, line properties, and text formatting.
5. The `saveas` function is used to save the figure instead of `savefig`.

Note that this conversion focuses on the main functionality of calculating and plotting edge states. Some additional features from the original Python script, such as time tracking and file I/O operations, are not included in this MATLAB version[1].

Citations:
[1] https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/40621784/63e98310-0f4b-4c9f-a05e-4d14fa2f1514/heterostructure.py
