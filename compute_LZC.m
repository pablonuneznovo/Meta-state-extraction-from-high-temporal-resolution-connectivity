function [LZComplexity] = compute_LZC(sequence)
% Computes the Lempel Ziv complexity of a symbolic sequence
%
%       Input:
%               - sequence: array with a symbolic sequence
%
%       Output:
%               - LZComplexity: Lempel-Ziv complexity of the symbolic
%               sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Author: Daniel Abasolo
%  Modified by: Carlos Gomez Peña and Pablo Nuñez Novo
%  last update: 10/06/2020

numStates=numel(unique(sequence));

n = length(sequence);

b=n/(log(n)/log(numStates));

% Initialize values
c = 1;
S = sequence(1);
Q = sequence(2);

for i = 2:n
    SQ = [S,Q];
    SQ_pi = [SQ(1:(length(SQ)-1))];
    
    k = findstr(Q,SQ_pi);
    
    if length(k)==0

        c = c+1;
        if (i+1)>n
            break;
        else
            S = [S,Q];
            Q = sequence(i+1);
        end
    else

        if (i+1)>n
            break;
        else
            Q = [Q,sequence(i+1)];
        end
    end
end

LZComplexity = c/b;