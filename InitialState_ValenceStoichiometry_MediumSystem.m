clear
clc

Folder='InitialState_MediumSystem_ValenceStoichiometry/';
mkdir(Folder)

Rescale=4.5;
BoxSize(1)=250;
BoxSize(2)=50;
BoxSize(3)=50;
BoxSize=BoxSize*2^(1/3);

L1=14;
Ratio=14./(12:16);
for L2=12:16
    for ratio=Ratio
        NP=5000;
        NP1=round(NP*ratio/(ratio+1)/L1);
        NP2=round(NP/(ratio+1)/L2);

        x=(0)*max(L1,L2)-(max(L1,L2)-1)/2;
        y=((0:15)-7.5)*2/3*1.12;
        z=((0:15)-7.5)*2/3*1.12;

        [X,Y,Z]=meshgrid(x,y,z);
        X=reshape(X,[],1);
        Y=reshape(Y,[],1);
        Z=reshape(Z,[],1);

        nx=length(x);
        ny=length(y);
        nz=length(z);

        Polymer=zeros(3,max(L1,L2),nx*ny*nz);

        for np=1:nx*ny*nz
            Polymer(1,:,np)=(0:max(L1,L2)-1)+X(np,1);
            Polymer(2,:,np)=0+Y(np,1);
            Polymer(3,:,np)=0+Z(np,1);
        end
        Polymer1=Polymer(:,1:L1,1:NP1);
        Polymer2=Polymer(:,1:L2,1:NP2);

        Polymer1=Polymer1*Rescale;
        Polymer2=Polymer2*Rescale;

        NM=L1*NP1+L2*NP2;
        Monomer=zeros(3,NM);
        NB=(L1-1)*NP1+(L2-1)*NP2;
        Bond=zeros(2,NB);
        Atype=zeros(1,NM);

        nm=0;
        nb=0;
        for np=1:NP1
            for l=1:L1
                nm=nm+1;
                Monomer(1,nm)=Polymer1(1,l,np);
                Monomer(2,nm)=Polymer1(2,l,np);
                Monomer(3,nm)=Polymer1(3,l,np);
                Atype(1,nm)=1;
                if l~=L1
                    nb=nb+1;
                    Bond(1,nb)=nm;
                    Bond(2,nb)=nm+1;
                end
            end
        end
        for np=1:NP2
            for l=1:L2
                nm=nm+1;
                Monomer(1,nm)=Polymer2(1,l,np);
                Monomer(2,nm)=Polymer2(2,l,np);
                Monomer(3,nm)=Polymer2(3,l,np);
                Atype(1,nm)=2;
                if l~=L2
                    nb=nb+1;
                    Bond(1,nb)=nm;
                    Bond(2,nb)=nm+1;
                end
            end
        end

        save([Folder 'L1_' num2str(L1) '_L2_' num2str(L2) ...
                     '_N1_' num2str(NP1) '_N2_' num2str(NP2) '.mat'],'Monomer','Bond','Atype','BoxSize');
    end
end

