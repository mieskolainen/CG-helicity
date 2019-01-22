% Generate (semi-PDG style) SU(2) Clebsch-Gordan table in a latex format
% up to spin 4.
%
% Should match 1-to-1 tables given in PDG.
%
% mikael.mieskolainen@cern.ch, 2017

clear; close all;
addpath ./src

% Create different j1 x j2 combinations
Jmin = 1/2;
Jmax = 4;
step = 1/2;
JJ = [];
for k = Jmin:step:Jmax
    for l = 1/2:1/2:k
        JJ = [JJ; k l];
    end
end

fileID = fopen('CGtable.tex','w');
tic;

% Loop over (j1,j2) combinations
for k = 1:size(JJ,1)

    j1 = JJ(k,1);
    j2 = JJ(k,2);
    
    m1 = -j1:step:j1;
    m2 = -j1:step:j2;
    
    J  = 0:step:j1+j2;
    
    [j1a,j1b] = rat(j1);
    [j2a,j2b] = rat(j2);
    
    titlestring = sprintf('[j_1 \\otimes j_2 = %d/%d \\otimes %d/%d]', j1a,j1b, j2a,j2b);
    titlestring = regexprep(titlestring, '/1', ''); % Replace /1 with empty
    fprintf(fileID,'\n\\begin{table}\n');
    fprintf(fileID,'\\tiny\n');
    
    fprintf(fileID,'\\caption{');
    fprintf(fileID, '$%s$}\n', titlestring);
    
    fprintf(fileID,'\\begin{center}\n');
    fprintf(fileID,'\\begin{tabular}{|c|c|c|c|c|c|}\n');
    fprintf(fileID, '\\hline \n');
    fprintf(fileID, '   & $m_1$ & $m_2$ & $J$ & $M$ & $\\langle j_1 j_2 m_1 m_2 | j_1 j_2 J M \\rangle$ \\\\ \n');
    fprintf(fileID, '\\hline \n');
    
    ind = 1; % Table row index
    
    for m1_ = m1
        for m2_ = m2
            for J_ = J
                for M_ = -J_:step:J_
                    
                    % Get CG coefficient
                    [gg,a,b] = clebschgordan(j1,j2,m1_,m2_,J_,M_);
                    if abs(gg - 0) > 1e-9 % not zero
                        
                        if (sign(gg) > 0)
                            sgn = '';
                        else
                            sgn = '-';
                        end
                        
                        [m1a,m1b] = rat(m1_);
                        [m2a,m2b] = rat(m2_);
                        [Ja,Jb]   = rat(J_);
                        [Ma,Mb]   = rat(M_);
                        
                        if (abs(a-b) < 1e-9) % Check if equal
                        string = sprintf('$%d$ & $%d/%d$ & $%d/%d$ & $%d/%d$ & $%d/%d$ & $1$ \\\\ \n', ...
                               ind, m1a,m1b, m2a,m2b, Ja,Jb, Ma,Mb);
                        else
                        string = sprintf('$%d$ & $%d/%d$ & $%d/%d$ & $%d/%d$ & $%d/%d$ & $%s\\sqrt{%d/%d}$ \\\\ \n', ...
                               ind, m1a,m1b, m2a,m2b, Ja,Jb, Ma,Mb, sgn,a,b); 
                        end
                        ind = ind + 1;
                        string = regexprep(string, '/1\$', '$'); % Remove redundant /1$ signs
                        fprintf(fileID, '%s', string);
                    end
                end
            end
        end
    end
    fprintf(fileID, '\\hline \n');
    fprintf(fileID, '\\end{tabular}\n');
    fprintf(fileID, '\\end{center}\n');
    fprintf(fileID, '\\end{table}\n');
end

fclose(fileID);
toc;

%% Create the latex file and create pdf
mainID = fopen('SU2.tex','w');

fprintf(mainID,'\\documentclass[12pt]{article}\n');
fprintf(mainID,'\\twocolumn \n');
fprintf(mainID,'\\usepackage[left=10px,right=10px,top=10px,bottom=10px,paperwidth=8in,paperheight=85in]{geometry}\n');
%fprintf(mainID,'\\usepackage{multicol}\n');
fprintf(mainID,'\\begin{document}\n');
fprintf(mainID,'\\huge{$SU(2)$ Clebsch-Gordan} \\\\ \n');
fprintf(mainID,'\\large{mikael.mieskolainen@cern.ch} \\\\ \n');
fprintf(mainID,'\\vspace{5em} \\\\ \n');
fprintf(mainID,'\\small \n');
fprintf(mainID,'$|j_1 j_2 m_1 m_2\\rangle \\equiv | j_1 m_1 \\rangle \\otimes | j_2 m_2 \\rangle$ \n');
fprintf(mainID,'\\vspace{2em} \\\\ \n');
fprintf(mainID,'$|JM\\rangle = \\sum_{-j_1 \\leq m_1 \\leq j_1} \\sum_{-j_2 \\leq m_2 \\leq j_2} |j_1 j_2 m_1 m_2\\rangle \\langle j_1 j_2 m_1 m_2 | (j_1 j_2) J M\\rangle$ \n');
fprintf(mainID,'\\vspace{0.5em} \\\\ \n');
fprintf(mainID,'$|j_1 j_2 m_1 m_2 \\rangle = \\sum_{J}\\sum_{-J \\leq M \\leq J}  | J M\\rangle \\langle (j_1 j_2) J M |j_1 j_2 m_1 m_2\\rangle$  \n');
fprintf(mainID,'\\newpage\n');
%fprintf(mainID,'\\begin{multicols}{3}\n');
fprintf(mainID,'\\input{CGtable.tex}\n');
%fprintf(mainID,'\\end{multicols}\n');
fprintf(mainID,'\\end{document}\n');

fclose(mainID);

system('pdflatex SU2.tex');
system('rm SU2.aux SU2.log');
