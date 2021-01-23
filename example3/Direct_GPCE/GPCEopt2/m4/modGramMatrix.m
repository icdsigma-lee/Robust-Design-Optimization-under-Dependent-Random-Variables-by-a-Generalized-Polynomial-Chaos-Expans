%% ========================================================================
%  function to generate monomial moment matrix or gram matrix (GPCE)
%  written by Dongjin Lee (dongjin-lee@uiowa.edu) 
%  output : generate three monomial moment matrices (N=7, N=5, N=5)
%% ========================================================================
function modGramMatrix
N = 7;
m = 4;
nA = nchoosek(N+m,m);
FilNam1 = sprintf('gram.mat'); 
FilNam2 = sprintf('frtst.mat'); %load ID, G1
% load file 
load(FilNam1);
clear ID;
% user-defined index to delete from the monomial moment matrix G1 
index0 = [6, 7];
index1 = [2, 5];
index2 = [1, 5];
cnt = 0;
for m0=1:m+1
    mm = m0-1; %total degree 
for i1=mm+1:-1:1
    for i2=mm+1:-1:1              
        for i3=mm+1:-1:1    
            for i4=mm+1:-1:1    
                for i5=mm+1:-1:1    
                	for i6=mm+1:-1:1   
                		for i7=mm+1:-1:1   
                            j1 = i1-1;
                            j2 = i2-1;
                            j3 = i3-1;
                            j4 = i4-1;
                            j5 = i5-1;
                            j6 = i6-1;
                            j7 = i7-1;
            if ((j1+j2+j3+j4+j5+j6+j7)==mm)
                cnt = cnt + 1;
                ID(cnt,:) = [j1 j2 j3 j4 j5 j6 j7];
            end 
                        end
                    end
                end 
         	end            
        end            
    end 
end 
end  


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
GRN = G1(1:nA,1:nA);
% modify the monomial moment matrix
GRN0 = GRN;
GRN0(INDEX0, :) = [];
GRN0(:, INDEX0) = [];

GRN1 = GRN;
GRN1(INDEX1, :) = [];
GRN1(:, INDEX1) = [];

GRN2 = GRN;
GRN2(INDEX2, :) = [];
GRN2(:, INDEX2) = [];

% generate whitenning matrix 
ORN0 = inv(chol(GRN0, 'lower'));
ORN1 = inv(chol(GRN1, 'lower'));
ORN2 = inv(chol(GRN2, 'lower'));
save(FilNam2, 'ID', 'ORN0', 'ORN1', 'ORN2', 'INDEX0', 'INDEX1', 'INDEX2');