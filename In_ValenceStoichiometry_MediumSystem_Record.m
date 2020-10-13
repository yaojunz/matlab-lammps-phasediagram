clear
clc

global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

LoadFolder='Parameter/Parameter.mat';
load(LoadFolder);

Damp=10;
BeadMass=Damp*BeadCsi;

InFolder='MediumSystem_ValenceStoichiometry/In/';
mkdir(InFolder);
OutFolder='Out_Record0/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record1/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record2/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record3/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record4/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record5/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record6/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record7/';
mkdir([InFolder OutFolder]);
OutFolder='Out_Record8/';
mkdir([InFolder OutFolder]);

Replicates=1; %10
L1=14;
Ratio=14./(12:16);
for record=0:8
    index=0;
    for L2=12:16
        for ratio=Ratio
            NP=5000;
            np1=round(NP*ratio/(ratio+1)/L1);
            np2=round(NP/(ratio+1)/L2);    
            for rep=(1:Replicates)
                index=index+1;
                Filename=['L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep)];
                if record==0
                    InitCondFilename=['Out_Relax/L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep) '.restart'];
                else
                    InitCondFilename=['Out_Record' num2str(record-1) '/L1_' num2str(L1) '_L2_' num2str(L2) '_N1_' num2str(np1) '_N2_' num2str(np2) '_Rep' num2str(rep) '.restart'];
                end
                InFilename=['Record' num2str(record) '_Index_' num2str(index) '.in'];
                OutFilename=['Out_Record' num2str(record) '/' Filename];
                InFileGenerate(InFolder,InFilename,InitCondFilename,OutFilename,rep)
            end
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

TimeStep=Damp/20;
RunSteps=10^8;
TimeRecordInterval=RunSteps/100;
Thermo=RunSteps;

fid=fopen([InFolder InFilename],'w');

fprintf(fid, ['processors 1 * *\n\n']);

fprintf(fid, ['read_restart ' InitCondFilename '\n\n']);

fprintf(fid, ['mass 1 ' num2str(BeadMass) '\n']);
fprintf(fid, ['mass 2 ' num2str(BeadMass) '\n\n']);

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
fprintf(fid, ['dump 1 all movie ' num2str(TimeRecordInterval) ' ' OutFilename '.mpeg type type zoom 4.5 box yes 0.01 view 85 85 size 1000 400 shiny 0.5\n']);
fprintf(fid, ['dump_modify 1 acolor 1 blue\n']); 
fprintf(fid, ['dump_modify 1 acolor 2 orange\n']); 
fprintf(fid, ['dump_modify 1 adiam 1 3\n']); 
fprintf(fid, ['dump_modify 1 adiam 2 3\n\n']);
end

fprintf(fid, ['dump 2 all xyz ' num2str(TimeRecordInterval) ' ' OutFilename '.xyz\n\n']); 

fprintf(fid, ['run ' num2str(RunSteps) '\n\n']);

fprintf(fid, ['write_restart ' OutFilename '.restart\n\n']);

fclose(fid); 

end

% function []=InitCondGenerate(InFolder,InitCondFilename,rep)
% % ReadFolder=['SAWLCPeriodicBoundary\In\Out_Equilibrium\'];
% % ReadFilename=['SAWLC_Rep' num2str(rep) '.dump'];
% % fid=fopen([ReadFolder ReadFilename],'r');
% % NH=2; 
% % for nh=1:NH
% %     Header=fgetl(fid);
% % end
% % Data=textscan(fid,'%f %f %f %f',BeadNumberP);
% % Header=fgetl(fid);
% % X=Data{1,2}';
% % Y=Data{1,3}';
% % Z=Data{1,4}';
% % fclose(fid);
% 
% global Monomer Bond Atype BoxSize BeadSize BeadCsi BeadMass Temp kBT Damp
% 
% Natom_type=2; %number of atom types;
% Natom=length(Monomer);
% Nbond_type=1; %number of bond types;
% Nbond=length(Bond); %number of DNA bonds;
% Nangl_type=0;  %number of angle types;
% Nangl=0; %number of angles;
% Ndihe_type=0; %number of dihedral types;
% Ndihe=0; %number of dihedrals;
% Nimpr_type=0; %number of improper types;
% Nimpr=0; %number of impropers;
% 
% xlo=-BoxSize(1)/2; xhi=BoxSize(1)/2; %x boundary
% ylo=-BoxSize(2)/2; yhi=BoxSize(2)/2; %y boundary
% zlo=-BoxSize(3)/2; zhi=BoxSize(3)/2; %z boundary
% 
% % Ratio=max([max(abs(X)),max(abs(Y)),max(abs(Z))])/BoxSize
% % X=zeros(1,Natom);
% % Y=zeros(1,Natom);
% % Z=zeros(1,Natom);
% % V=(kBT/BeadMass)^0.5*normrnd(0,1,Natom,3);
% 
% V=zeros(3,Natom);
% 
% fid=fopen([InFolder InitCondFilename],'w');
% fprintf(fid,'LAMMPS chain data file\n\n');
% fprintf(fid,'%d atoms\n', Natom);
% fprintf(fid,'%d bonds\n', Nbond);
% fprintf(fid,'%d angles\n', Nangl);
% fprintf(fid,'%d dihedrals\n', Ndihe);
% fprintf(fid,'%d impropers\n\n', Nimpr);
% fprintf(fid,'%d atom types\n', Natom_type);
% fprintf(fid,'%d bond types\n', Nbond_type);
% fprintf(fid,'%d angle types\n', Nangl_type);
% fprintf(fid,'%d dihedral types\n', Ndihe_type);
% fprintf(fid,'%d improper types\n\n', Nimpr_type);
% fprintf(fid,'%8.5f %8.5f xlo xhi\n', xlo, xhi);
% fprintf(fid,'%8.5f %8.5f ylo yhi\n', ylo, yhi);
% fprintf(fid,'%8.5f %8.5f zlo zhi\n\n', zlo, zhi);
% 
% fprintf(fid,'Masses\n\n');
% fprintf(fid,'%d %8.5f\n',1,BeadMass);
% fprintf(fid,'%d %8.5f\n\n',2,BeadMass);
% 
% fprintf(fid,'Atoms\n\n');
% for i=1:Natom
%     Mole_type=Atype(i);
%     Atom_type=Atype(i);
%     fprintf(fid,[num2str(i) ' ' num2str(Mole_type) ' ' num2str(Atom_type) ' ' ...
%                  num2str(Monomer(1,i)) ' ' num2str(Monomer(2,i)) ' ' num2str(Monomer(3,i)) '\n']);
% end
% fprintf(fid,'\n');
% 
% fprintf(fid,'Velocities\n\n');
% for i=1:Natom
%     fprintf(fid,[num2str(i) ' ' num2str(V(1,i)) ' ' num2str(V(2,i)) ' ' num2str(V(3,i)) '\n']);
% end
% fprintf(fid,'\n');
% 
% fprintf(fid,'Bonds\n\n');
% for i=1:Nbond
%     Bond_type=1;
%     fprintf(fid,'%d %d %d %d\n',...
%             i,Bond_type,Bond(1,i),Bond(2,i));
% end
% fprintf(fid,'\n');
% % 
% % fprintf(fid,'Angles\n\n');
% % for i=1:Nangl
% %     Angl_type=1;
% %     fprintf(fid,'%d %d %d %d %d\n',...
% %             i,Angl_type,i,i+1,i+2);
% % end
% 
% end
