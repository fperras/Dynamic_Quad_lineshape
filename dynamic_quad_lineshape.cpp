#include <stdlib.h>
#include <stdio.h>
#include "wigner_complex.h"
#include "unsupported/Eigen/FFT"
#include "unsupported/Eigen/MatrixFunctions"
using namespace std;
#include <vector>

int main()
{
    char buffer[128], keyword[64], crystal_file[128], input_filename[128];
    int i,j,k,t,counter=0;
    int nsites=0, np=1024, gamma_steps=1;
    std::vector< std::vector<double> > aa,bb,gg,CQ,etaQ,diso,RQ0,RQ2;
    std::vector<double> kex, null_vector;
    std::vector<int> njumps;
    double offset=0., vL=94.7e6, dwell = 0.00001, sw=100000., lb=250,I=3.5, wr = 12500. *2.*Pi, magic_angle = acos(sqrt(1./3.));;
    FILE *input, *fp;

    printf("This program is for the simulation of static and MAS\n");
    printf("NMR spectra of quadrupolar nuclei.\n");
    printf("It comes with no warranty.\n\n");
    printf("written by Frederic A Perras at Ames National Laboratory\n\n");
    printf("What is the name of your input file?\n");
    scanf("%s",input_filename);
    char output_filename[strlen(input_filename)-3];
    sprintf(output_filename,"%.*s",strlen(input_filename)-4,input_filename);

    input=fopen(input_filename,"r");
    if(input==NULL){
        FILE *error_file;
        error_file=fopen("error.txt","w");
        fprintf(error_file, "\nERROR: Input file  not found\n");
        fclose(error_file);
        exit(1);
    }

    //reading the input file
    while(fgets(buffer, sizeof(buffer), input) != NULL){
        sscanf(buffer,"%s",keyword);

        if(strcmp(keyword, "site")==0){
            int site=0, junk;
            sscanf(buffer,"%s %d",keyword, &site);
            site=site-1;
            if(site>nsites){
                FILE *error_file;
                error_file=fopen("errors.txt","a");
                fprintf(error_file, "\nERROR: Sites should be introduced in chronological order\n", crystal_file);
                fclose(error_file);
                exit(1);
            }
            if (site>(nsites-1)){
                diso.push_back(null_vector);
                CQ.push_back(null_vector);
                etaQ.push_back(null_vector);
                RQ0.push_back(null_vector);
                RQ2.push_back(null_vector);
                aa.push_back(null_vector);
                bb.push_back(null_vector);
                gg.push_back(null_vector);
                kex.push_back(0.);
                njumps.push_back(0);
                nsites++;
            }
            diso[site].push_back(0.);
            CQ[site].push_back(0.);
            etaQ[site].push_back(0.);
            RQ0[site].push_back(0.);
            RQ2[site].push_back(0.);
            aa[site].push_back(0.);
            bb[site].push_back(0.);
            gg[site].push_back(0.);

            sscanf(buffer,"%s %d %lf %lf %lf %lf %lf %lf",keyword, &junk, &diso[site][njumps[site]], &CQ[site][njumps[site]], &etaQ[site][njumps[site]], &aa[site][njumps[site]], &bb[site][njumps[site]], &gg[site][njumps[site]]);
            sprintf(keyword,"void");
            njumps[site]++;
        }
        if(strcmp(keyword, "offset")==0){
            sscanf(buffer,"%s %lf",keyword, &offset);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "spinrate")==0){
            sscanf(buffer,"%s %lf",keyword, &wr);
            wr=wr*2.*Pi;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "spin")==0){
            sscanf(buffer,"%s %lf",keyword, &I);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "Larmor")==0){
            sscanf(buffer,"%s %lf",keyword, &vL);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "kex")==0){
            int site,junk;
            sscanf(buffer,"%s %d",keyword, &site);
            site=site-1;
             if(site>nsites){
                FILE *error_file;
                error_file=fopen("errors.txt","a");
                fprintf(error_file, "\nERROR: Sites should be introduced in chronological order\n", crystal_file);
                fclose(error_file);
                exit(1);
            }
            if (site>(nsites-1)){
                diso.push_back(null_vector);
                CQ.push_back(null_vector);
                etaQ.push_back(null_vector);
                RQ0.push_back(null_vector);
                RQ2.push_back(null_vector);
                aa.push_back(null_vector);
                bb.push_back(null_vector);
                gg.push_back(null_vector);
                kex.push_back(0.);
                njumps.push_back(0);
                nsites++;
            }
            sscanf(buffer,"%s %d %lf",keyword, &junk, &kex[site]);
            kex[site]=kex[site]*2.*Pi;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "lb")==0){
            sscanf(buffer,"%s %lf",keyword, &lb);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "sw")==0){
            sscanf(buffer,"%s %lf",keyword, &sw);
            dwell=1./sw;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "np")==0){
            sscanf(buffer,"%s %d",keyword, &np);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "gamma_angles")==0){
            sscanf(buffer,"%s %d",keyword, &gamma_steps);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "crystal_file")==0){
            sscanf(buffer,"%s %s",keyword, crystal_file);
            sprintf(crystal_file,"%s.cry",crystal_file);
            sprintf(keyword,"void");
        }
    }//done reading input file

    double w0=vL*2*Pi;

    if(I==1){
    for(i=0;i<nsites;i++){
        njumps.push_back(njumps[i]);
        RQ0.push_back(null_vector);
        RQ2.push_back(null_vector);
        aa.push_back(null_vector);
        bb.push_back(null_vector);
        gg.push_back(null_vector);
        kex.push_back(kex[i]);
        diso.push_back(null_vector);
        for(j=0;j<njumps[i];j++){
            double eq=-2.*sqrt(6.)*Pi/(4.*I*(2*I-1.))*CQ[i][j];
            RQ0[i][j] = eq*etaQ[i][j]/2.;
            RQ2[i][j] = sqrt(1.5)*eq;
            diso[i][j]=2.*Pi*diso[i][j];
            aa[i][j]=aa[i][j]*Pi/180.;
            bb[i][j]=bb[i][j]*Pi/180.;
            gg[i][j]=gg[i][j]*Pi/180.;

            RQ0[i+nsites].push_back(-RQ0[i][j]);
            RQ2[i+nsites].push_back(-RQ2[i][j]);
            aa[i+nsites].push_back(aa[i][j]);
            bb[i+nsites].push_back(bb[i][j]);
            gg[i+nsites].push_back(gg[i][j]);
            diso[i+nsites].push_back(diso[i][j]);
    }}
    nsites=nsites*2;
    }

    else{
    for(i=0;i<nsites;i++){
        for(j=0;j<njumps[i];j++){
            double eq=-2.*sqrt(6./w0)*Pi/(4.*I*(2*I-1.))*CQ[i][j];
            RQ0[i][j] = eq/2.*etaQ[i][j]/sqrt(6.);
            RQ2[i][j] = eq/2.;
            diso[i][j]=2.*Pi*diso[i][j];
            aa[i][j]=aa[i][j]*Pi/180.;
            bb[i][j]=bb[i][j]*Pi/180.;
            gg[i][j]=gg[i][j]*Pi/180.;
    }}
    }

    //reading the crystal file
    FILE *cry;
    cry=fopen(crystal_file,"r");
    if(cry==NULL){
        FILE *error_file;
        error_file=fopen("errors.txt","a");
        fprintf(error_file, "\nERROR: crystal file '%s' not found\n", crystal_file);
        fclose(error_file);
        exit(1);
    }
    fgets(buffer, sizeof(buffer), cry);
    int n_orient;
    sscanf(buffer,"%d",&n_orient);
    double *al = (double *) malloc(n_orient*gamma_steps*sizeof(double));
    double *bt = (double *) malloc(n_orient*gamma_steps*sizeof(double));
    double *gm = (double *) malloc(n_orient*gamma_steps*sizeof(double));

    for(i=0;i<n_orient;i++){
        fgets(buffer, sizeof(buffer), cry);
        for(j=0;j<gamma_steps;j++){
            sscanf(buffer,"%lf %lf %lf",&al[i*gamma_steps+j],&bt[i*gamma_steps+j],&gm[i*gamma_steps+j]);
            gm[i*gamma_steps+j]=360.*j/gamma_steps;
    }}
    fclose(cry);

    //FID and spectrum initialization
    Eigen::VectorXcd total_FID(np);
    Eigen::VectorXcd spectrum(np);
    total_FID.setZero();
    total_FID(0)=n_orient*gamma_steps*nsites;

    #pragma omp parallel private(i,j,k,t)
    {
    for(int site=0;site<nsites;site++){
        //Exchange matrix
        Eigen::MatrixXcd K(njumps[site],njumps[site]);
        for(i=0;i<njumps[site];i++){
            K(i,i)=kex[site]*dwell*(1.-1.*njumps[site]);
            for(j=i+1;j<njumps[site];j++){
                K(i,j)=K(j,i)=kex[site]*dwell;
            }
        }

        //variable declarations
        double alpha, beta, gamma, vCS;
        Eigen::VectorXcd FID(np);
        Eigen::VectorXcd R2(5);
        Eigen::MatrixXcd D_MOL_MAS(5,5);
        Eigen::MatrixXcd D_MAS_LAB(5,5);
        Eigen::MatrixXcd D_PAS_MOL(5,5);
        Eigen::VectorXcd R_Q(5);
        R_Q(1) = R_Q(3) = 0.;
        Eigen::MatrixXcd Omega(njumps[site],njumps[site]);
        Eigen::MatrixXcd Prop(njumps[site],njumps[site]);
        Eigen::VectorXcd M(njumps[site]);
        Eigen::VectorXcd M_temp(njumps[site]);
        Omega.setIdentity();

        //powder averaging
        for(int orientation=omp_get_thread_num(); orientation<n_orient*gamma_steps; orientation=orientation+omp_get_num_threads()){
            beta=2.*Pi*bt[orientation]/360.;
            alpha=2.*Pi*al[orientation]/360.;
            gamma=2.*Pi*gm[orientation]/360.;
            counter++;
            printf("%d/%d\n",counter,n_orient*gamma_steps*nsites);

            FID.setZero();
            FID(0)=1.0/njumps[site];
            for(i=0;i<njumps[site];i++){
                M(i)=M_temp(i)=1./njumps[site];
            }
            D_MOL_MAS=D2_matrix(alpha, beta, gamma);

            for(t=1;t<np;t++){//time
                D_MAS_LAB=D2_matrix(wr*t*dwell, magic_angle,0.);

                for(i=0;i<njumps[site];i++){ //loop over the different jump orientations
                    //EFG tensor
                    R_Q(0) = R_Q(4) = RQ0[site][i];
                    R_Q(2) = RQ2[site][i];
                    D_PAS_MOL=D2_matrix(aa[site][i], bb[site][i], gg[site][i]);

                    //Rotating the PAS tensor into the lab frame
                    R2 = D_MAS_LAB*D_MOL_MAS*D_PAS_MOL*R_Q;

                    //Calculating the quadrupole shift
                    if(I==1)
                        vCS =diso[site][i]*vL/1000000.-offset*2.*Pi+  real(R2(2)); //I=1
                    else
                        vCS =diso[site][i]*vL/1000000.-offset*2.*Pi+  (4.*I*(I+1.)-3.)* (2*(real(R2(1))*real(R2(3))-imag(R2(1))*imag(R2(3)))+ (real(R2(0))*real(R2(4))- imag(R2(0))*imag(R2(4)))); //CT for I=N/2

                    Omega(i,i)=vCS*1i*dwell;
                }
                Prop=Omega+K;
                M_temp=Prop.exp()*M;
                FID(t)=0.;
                for(i=0;i<njumps[site];i++){
                    M(i)=M_temp(i);
                    FID(t)=FID(t)+M(i);
                }
                //The value for that data point is added to the total FID.
                total_FID(t) = total_FID(t)+FID(t)*exp(-dwell*t*lb);
                }//end of time loop
            } //end of powder averaging
        }//sites
    }//end of parallel bit

    free(al);
    free(bt);
    free(gm);

    //Fast Fourier Transform of the FID
    Eigen::FFT<double> fft;
    fft.fwd(spectrum,total_FID);

    //writing out the fid
    char out_file[128];
    sprintf(out_file,"%s.fid",output_filename);
    fp=fopen(out_file,"w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nTYPE=FID\nDATA\n",np,1./dwell);
    for(i=0;i<np;i++){
        fprintf(fp,"%f %f\n",real(total_FID(i)),imag(total_FID(i)));
    }
    fprintf(fp,"%s","END");
    fclose(fp);

    //writing the spectrum as a SPE file
    sprintf(out_file,"%s.spe",output_filename);
    fp=fopen(out_file,"w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nTYPE=SPE\nDATA\n",np,1./dwell,0.5/dwell+offset,vL/1000000.);
    //right edge of the spectrum
    k=0;
    for(i=np/2-1;i>=0;i--){
        fprintf(fp,"%f %f\n",real(spectrum(i)),imag(spectrum(i)));
        k++;
    }
    //left edge of the spectrum
    k=0;
    for(i=np-1;i>=np/2;i--){
        fprintf(fp,"%f %f\n",real(spectrum(i)),imag(spectrum(i)));
        k++;
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    return 0;
}
