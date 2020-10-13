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
OutFolder='Out_Anneal/';
mkdir([InFolder OutFolder]);

Replicates=1;
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
    
        LoadFolder='InitialState_MediumSystem_Stoichiometry/';
        load([LoadFolder 'L1_' num2str(L1) '_L2_' num2str(L2) ...
                        '_N1_' num2str(np1) '_N2_' num2str(np2) '.mat']);

        for rep=(1:Replicates)
            index=index+1;
            Filename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep)];
            InitCondFilename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '.initial'];
            if rep==1
                InitCondGenerate(InFolder,InitCondFilename,rep)
            end
            InFilename=['Anneal_Index_' num2str(index) '.in'];
            OutFilename=Filename;
            InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,rep)
        end
    end
end

function []=InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,rep)

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp
Ai=-0*kBT;
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
fprintf(fid, ['units nano \n']);
fprintf(fid, ['boundary p p p \n']);
fprintf(fid, ['atom_style bond\n\n']);

fprintf(fid, ['processors 1 * *\n\n']);

fprintf(fid, ['read_data ' InitCondFilename '\n\n']);

fprintf(fid, ['pair_style hybrid lj/cut/soft 1 1 ' num2str(5) ' soft ' num2str(Rc) '\n']);
fprintf(fid, ['pair_coeff 1 1 lj/cut/soft ' num2str(eps) ' ' num2str(3.5) ' ' num2str(lambda) '\n']);
fprintf(fid, ['pair_coeff 2 2 lj/cut/soft ' num2str(eps) ' ' num2str(3.5) ' ' num2str(lambda) '\n']);
fprintf(fid, ['pair_coeff 1 2 soft ' num2str(Ai) ' ' num2str(Rc) '\n']);
fprintf(fid, ['pair_modify shift yes\n']);
fprintf(fid, ['special_bonds lj/coul 1.0 1.0 1.0\n\n']);

fprintf(fid, ['bond_style harmonic\n']);
fprintf(fid, ['bond_coeff 1 ' num2str(20*kBT/BondLength^2) ' ' num2str(BondLength) '\n\n']);

fprintf(fid, ['variable A equal "ramp(' num2str(Ai) ',' num2str(Af) ')"\n\n']);

fprintf(fid, ['region wallx block -' num2str(Rb) ' ' num2str(Rb) ...
                                ' -' num2str(60) ' ' num2str(60) ...
                                ' -' num2str(60) ' ' num2str(60) ' open 3 open 4 open 5 open 6\n\n']);

fprintf(fid, ['neighbor 5 bin \n']);
fprintf(fid, ['neigh_modify every 1 delay 0\n\n']);

fprintf(fid, ['fix 1 all nve\n']);
fprintf(fid, ['fix 2 all langevin ' num2str(Temp) ' ' num2str(Temp) ' ' num2str(Damp) ' ' num2str(randi(10^7)) ' zero yes\n']);
fprintf(fid, ['fix 3 all adapt 1 pair soft a 1 2 v_A\n']);
fprintf(fid, ['fix 4 all wall/region wallx harmonic ' num2str(100*kBT/BeadSize^2) ' 0 ' num2str(2*BeadSize) '\n\n']);

fprintf(fid, ['thermo ' num2str(Thermo) '\n']);
fprintf(fid, ['timestep ' num2str(TimeStep) '\n\n']);

if rep==1
fprintf(fid, ['dump 1 all movie ' num2str(TimeRecordInterval) ' Out_Movie/' OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
fprintf(fid, ['dump_modify 1 acolor 1 blue\n']); 
fprintf(fid, ['dump_modify 1 acolor 2 orange\n']); 
fprintf(fid, ['dump_modify 1 adiam 1 3\n']); 
fprintf(fid, ['dump_modify 1 adiam 2 3\n\n']);
end

fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);
fprintf(fid, ['write_restart Out_Anneal/' OutFilename '.restart\n']);

fclose(fid); 

end

function []=InitCondGenerate(InFolder,InitCondFilename,rep)

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp

Natom_type=2; %number of atom types;
Natom=length(Monomer);
Nbond_type=1; %number of bond types;
Nbond=length(Bond); %number of DNA bonds;
Nangl_type=0;  %number of angle types;
Nangl=0; %number of angles;
Ndihe_type=0; %number of dihedral types;
Ndihe=0; %number of dihedrals;
Nimpr_type=0; %number of improper types;
Nimpr=0; %number of impropers;

xlo=-BoxSize(1)/2; xhi=BoxSize(1)/2; %x boundary
ylo=-BoxSize(2)/2; yhi=BoxSize(2)/2; %y boundary
zlo=-BoxSize(3)/2; zhi=BoxSize(3)/2; %z boundary

V=zeros(3,Natom);

fid=fopen([InFolder InitCondFilename],'w');
fprintf(fid,'LAMMPS chain data file\n\n');
fprintf(fid,'%d atoms\n', Natom);
fprintf(fid,'%d bonds\n', Nbond);
fprintf(fid,'%d angles\n', Nangl);
fprintf(fid,'%d dihedrals\n', Ndihe);
fprintf(fid,'%d impropers\n\n', Nimpr);
fprintf(fid,'%d atom types\n', Natom_type);
fprintf(fid,'%d bond types\n', Nbond_type);
fprintf(fid,'%d angle types\n', Nangl_type);
fprintf(fid,'%d dihedral types\n', Ndihe_type);
fprintf(fid,'%d improper types\n\n', Nimpr_type);
fprintf(fid,'%8.5f %8.5f xlo xhi\n', xlo, xhi);
fprintf(fid,'%8.5f %8.5f ylo yhi\n', ylo, yhi);
fprintf(fid,'%8.5f %8.5f zlo zhi\n\n', zlo, zhi);

fprintf(fid,'Masses\n\n');
fprintf(fid,'%d %8.5f\n',1,BeadMass);
fprintf(fid,'%d %8.5f\n\n',2,BeadMass);

fprintf(fid,'Atoms\n\n');
for i=1:Natom
    Mole_type=Atype(i);
    Atom_type=Atype(i);
    fprintf(fid,[num2str(i) ' ' num2str(Mole_type) ' ' num2str(Atom_type) ' ' ...
                 num2str(Monomer(1,i)) ' ' num2str(Monomer(2,i)) ' ' num2str(Monomer(3,i)) '\n']);
end
fprintf(fid,'\n');

fprintf(fid,'Velocities\n\n');
for i=1:Natom
    fprintf(fid,[num2str(i) ' ' num2str(V(1,i)) ' ' num2str(V(2,i)) ' ' num2str(V(3,i)) '\n']);
end
fprintf(fid,'\n');

fprintf(fid,'Bonds\n\n');
for i=1:Nbond
    Bond_type=1;
    fprintf(fid,'%d %d %d %d\n',...
            i,Bond_type,Bond(1,i),Bond(2,i));
end
fprintf(fid,'\n');
fclose(fid); 
end
