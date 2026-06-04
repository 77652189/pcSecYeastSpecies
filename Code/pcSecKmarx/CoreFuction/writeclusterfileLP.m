function writeclusterfileLP(allname,subname)

subfileName = [subname,'.sh'];
fptr = fopen(subfileName,'w');
if fptr < 0
    error('Unable to open %s for writing.', subfileName);
end
soplexBin = getenv('SOPLEX_BIN');
if isempty(soplexBin)
    soplexBin = 'soplex';
end
if contains(soplexBin, ' ')
    soplexBin = ['"', soplexBin, '"'];
end
soplexReadMode = getenv('SOPLEX_READMODE');
if isempty(soplexReadMode)
    soplexReadMode = '0';
end
fprintf(fptr,'#!/bin/bash\n');
fprintf(fptr,'#SBATCH -n 20\n');
fprintf(fptr,'#SBATCH -o out.txt\n');
fprintf(fptr,'#SBATCH --time 0-2:00:00\n');
fprintf(fptr,'#SBATCH --mail-user=liulz23@mails.tsinghua.edu.cn\n');
fprintf(fptr,'#SBATCH --mail-type=end\n');
for i = 1:length(allname)
    fprintf(fptr,[soplexBin,' -s0 -g5 -t10000 -f1e-20 -o1e-20 -x -q -c --int:readmode=',soplexReadMode,' --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 ',allname{i},' > ', [allname{i},'.out'],' &\n']);
    if mod(i,50)==0
        fprintf(fptr,'wait;\n');
    end

%bash sub.sh
end
fprintf(fptr,'wait;\n');
fclose(fptr);
