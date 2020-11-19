%% ========================================================================
%  function to generate monomial moment matrix or gram matrix (GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%  output : generate three monomial moment matrices (N=7, N=5, N=5)
%% ========================================================================

function modGramMatrix
FilNam = sprintf('gram.mat'); %load ID, G1

% user-defined index to delete from the monomial moment matrix G1 
index0 = [6, 7];
index1 = [2, 5];
index2 = [1, 5];

% load file 
load(FilNam);
tmp1 = find(ID(:, index0(1)) ~= 0);
tmp2 = find(ID(:, index0(2)) ~= 0);
INDEX0 = unique([tmp1, tmp2]);
clear 'tmp1', 'tmp2'; 

tmp1 = find(ID(:, index1(1)) ~= 0);
tmp2 = find(ID(:, index1(2)) ~= 0);
INDEX1 = unique([tmp1, tmp2]);
clear 'tmp1', 'tmp2'; 

tmp1 = find(ID(:, index2(1)) ~= 0);
tmp2 = find(ID(:, index2(2)) ~= 0);
INDEX2 = unique([tmp1, tmp2]);
clear 'tmp1', 'tmp2';

% modify the monomial moment matrix
G0 = G1;
G0(INDEX0, :) = [];
G0(:, INDEX0) =[];

G2 = G1;
G2(INDEX1, :) = [];
G2(:, INDEX1) =[];

G3 = G1;
G3(INDEX2, :) = [];
G3(:, INDEX2) = [];

% generate whitenning matrix 
QQ1 = inv(chol(G0, 'lower'));
QQ2 = inv(chol(G2, 'lower'));
QQ3 = inv(chol(G3, 'lower'));

save(FilNam, 'ID', 'QQ1', 'QQ2', 'QQ3', 'INDEX0', 'INDEX1', 'INDEX2');