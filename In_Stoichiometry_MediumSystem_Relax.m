clear
clc

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

LoadFolder='Parameter/Parameter.mat';
load(LoadFolder);

Damp=125;
BeadMass=Damp*BeadCsi;

InFolder='MediumSystem_Stoichiometry/In/';
mkdir(InFolder);
OutFolder='Out_Relax/';
mkdir([InFolder OutFolder]);

Replicates=1; %10
index=0;

for ns=1:8
    if ns==1
        L1=8; 
        L2=8;
        DR=log(L1/(L1-1))/4;
        Ratio=exp(0:DR:(6*DR));
    elseif ns==2
        L1=8; 
        L2=7; 
        DR=log(L1/(L1-1))/4;
        Ratio=exp((-6*DR):DR:(6*DR));
    elseif ns==3
        L1=10; 
        L2=10;  
        DR=log(L1/(L1-1))/3;
        Ratio=exp(0:DR:(6*DR));      
    elseif ns==4  
        L1=10; 
        L2=9;
        DR=log(L1/(L1-1))/3;
        Ratio=exp((-6*DR):DR:(6*DR));
    elseif ns==5  
        L1=12; 
        L2=12; 
        DR=log(L1/(L1-1))/2;
        Ratio=exp(0:DR:(6*DR)); 
    elseif ns==6  
        L1=12; 
        L2=11; 
        DR=log(L1/(L1-1))/2;
        Ratio=exp((-6*DR):DR:(6*DR));
    elseif ns==7  
        L1=14; 
        L2=14; 
        DR=log(L1/(L1-1))/1.5;
        Ratio=exp(0:DR:(6*DR));        
    elseif ns==8  
        L1=14; 
        L2=13; 
        DR=log(L1/(L1-1))/1.5;
        Ratio=exp((-6*DR):DR:(6*DR));
        Ratio=[Ratio,L1/L2];
    end
    NR=length(Ratio);
    
    for nr=1:NR
        
        ratio=Ratio(nr);
        NP=2500;
        np1=round(NP*ratio/(ratio+1)/L1);
        np2=round(NP/(ratio+1)/L2);
 
        for rep=(1:Replicates)
            index=index+1;
            Filename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep)];
            InitCondFilename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep) '.restart'];
            InFilename=['Relaxation_Index_' num2str(index) '.in'];
            OutFilename=Filename;
            InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,rep)
        end
    end
end

function []=InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,rep)

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp
Af=-7*kBT;
eps=0.15*kBT;
lambda=0.68;
Rc=2;
BondLength=4.5;
Rb=50;

TimeStep=Damp/50;
RunSteps=10^8;
TimeRecordInterval=RunSteps/100;
Thermo=RunSteps;

fid=fopen([InFolder InFilename],'w');

fprintf(fid, ['processors 1 * *\n\n']);

fprintf(fid, ['read_restart Out_Anneal/' InitCondFilename '\n\n']);

fprintf(fid, ['pair_style hybrid lj/cut/soft 1 1 ' num2str(5) ' soft ' num2str(Rc) '\n']);
fprintf(fid, ['pair_coeff 1 1 lj/cut/soft ' num2str(eps) ' ' num2str(3.5) ' ' num2str(lambda) '\n']);
fprintf(fid, ['pair_coeff 2 2 lj/cut/soft ' num2str(eps) ' ' num2str(3.5) ' ' num2str(lambda) '\n']);
fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Af) ' ' num2str(Rc) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 1.0 1.0 1.0\n\n']);

fprintf(fid, ['bond_style harmonic\n']);
fprintf(fid, ['bond_coeff 1 ' num2str(20*kBT/BondLength^2) ' ' num2str(BondLength) '\n\n']);

fprintf(fid, ['neighbor 5 bin \n']);
fprintf(fid, ['neigh_modify every 1 delay 0\n\n']);

fprintf(fid, ['fix 1 all nve\n']);
fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero yes\n']);

fprintf(fid, ['thermo ' num2str(Thermo) '\n']);
fprintf(fid, ['timestep ' num2str(TimeStep) '\n\n']);

if rep==1
fprintf(fid, ['dump 1 all movie ' num2str(TimeRecordInterval) ' Out_Relax/' OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
fprintf(fid, ['dump_modify 1 acolor 1 blue\n']); 
fprintf(fid, ['dump_modify 1 acolor 2 orange\n']); 
fprintf(fid, ['dump_modify 1 adiam 1 3\n']); 
fprintf(fid, ['dump_modify 1 adiam 2 3\n\n']);
end

fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);

fprintf(fid, ['write_restart Out_Relax/' OutFilename '.restart\n']);

fclose(fid); 

end

