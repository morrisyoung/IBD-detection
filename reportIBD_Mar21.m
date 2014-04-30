clear all

tic


% a substitute version of reading file
%fid = fopen('1000G_50I_100Million.trees', 'r');
% fid = fopen('test_10P_1Chrom_100Million.trees', 'r');
% fname = '10000N_100C_1M_1000m.trees';
fname = 'test_1000G_50I_100Million.trees';
fid = fopen(fname, 'r');
%fid = fopen('test_4P_1Chrom.trees', 'r');
i=1;
while 1
    nextline = fgetl(fid);
    if ~ischar(nextline)
        break
    end
    %disp(nextline);
    data{i} = nextline;
    i = i+1;
end


% data = readFile('1000G_25I_50Million.trees');
% data = readFile('1000G_50I_100Million.trees');
nt = length(data)-4;
N = 1000;
% N = 10000;

% n = 25;
% L = 0.5;

n = 50;
% n =250;
L = 1;

last_seg = zeros(n,n);
np = nchoosek(n,2);

%m = 0.0005;
%m = 0.00001;
m = 0.0005;
scale = 1e8;

% f = fopen(sprintf('report_cutoff%d.txt',round(m*scale)),'w');
f = fopen(sprintf(strcat(fname, '_result_matlab'),round(m*scale)),'w');

tot_len = 0;
tot_seg = 0;
mx = 0;
for i=1:(nt+1)
    if mod(i,100)==0
        fprintf('i=%d\n',i);
    end
    
    if i<=nt
        tree_str = data{i+3};
        x = strsplit(tree_str);
        y = strsplit(x{3},'_');
        pos = str2num(y{6});                
        
        tree = phytreeread(tree_str);
        
        name_strs = get(tree,'LeafNames');
        id = zeros(1,n);
        for j=1:n
            z = strsplit(name_strs{j},'.');
            id(j) = str2num(z{1});
        end
        
        Dtmp = pdist(tree,'SquareForm',true);
        D = zeros(n,n);
        D(id,id) = Dtmp;        
    else
        D = Inf*ones(n,n);
        pos = L*scale;
    end
    
    if i>1
        
        for n1=1:n
            for n2=(n1+1):n
                d = abs(D(n1,n2)-lastD(n1,n2));
                if d>0.008 % previous 0.002
                    len = pos-last_seg(n1,n2);
                    if len > (m*scale)
                        fprintf(f,'%d-%d [%d,%d]\n',n1,n2,last_seg(n1,n2),pos);
                        tot_len = tot_len + len;
                        tot_seg = tot_seg + 1;
                    end
                    last_seg(n1,n2) = pos;
                end
            end
        end
    end
    
    lastD = D;
end
fclose(f);
fprintf('Total len: %g, expected: %g\n',tot_len/(np*scale*L),(1+4*m*N)/(1+2*m*N)^2);
fprintf('Total seg: %g, expected: %g\n',tot_seg/np,2*N*L/(1+2*m*N)^2);

toc
