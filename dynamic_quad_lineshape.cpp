#include <stdlib.h>
#include <stdio.h>
#include "wigner_complex.h"
#include "unsupported/Eigen/FFT"
#include "unsupported/Eigen/MatrixFunctions"
using namespace std;
#include <vector>

int main()
{
    int i,j,k,t,counter=0;
    int njumps=0;
    int np=1024, gamma_steps=1;
    double shift, diso=0.;
    char buffer[128], keyword[64], crystal_file[128];
    std::vector<double> aa,bb,gg;

    //Fit parameters
    double wr = 12500. *2.*Pi; //MAS rate
    double I=7./2.; //spin quantum number
    double w0=94.7e6 *2*Pi; //Larmor frequency
    double CQ=11.8e6 *2*Pi; //Quadrupole coupling constant
    double etaQ=0.0;  //Quadrupole asymmetry parameter
    double kex=1.0e10 *2*Pi; //exchange constant for dynamics
    double lb=250.; //Lorentzian line broadening in Hz


    //spectral width and time increment size
    double dwell = 0.0000025; //dwell time  This parameter is also the maxdt so it needs to be small such that the Hamiltonian is time independent.


    FILE *input, *fp;
    input=fopen("input.txt","r");
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
        if(strcmp(keyword, "CQ")==0){
            sscanf(buffer,"%s %lf",keyword, &CQ);
            CQ=CQ*2.*Pi;
            sprintf(keyword,"void");
        }
        if(strcmp(keyword, "diso")==0){
            sscanf(buffer,"%s %lf",keyword, &diso);
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
            sscanf(buffer,"%s %lf",keyword, &w0);
            w0=w0*2.*Pi;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "eta")==0){
            sscanf(buffer,"%s %lf",keyword, &etaQ);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "kex")==0){
            sscanf(buffer,"%s %lf",keyword, &kex);
            kex=kex*2.*Pi;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "lb")==0){
            sscanf(buffer,"%s %lf",keyword, &lb);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "dwell")==0){
            sscanf(buffer,"%s %lf",keyword, &dwell);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "np")==0){
            sscanf(buffer,"%s %d",keyword, &np);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "orientation")==0){
            aa.push_back(0.);
            bb.push_back(0.);
            gg.push_back(0.);
            sscanf(buffer,"%s %lf %lf %lf",keyword, &aa[njumps], &bb[njumps], &gg[njumps]);
            njumps++;
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

    double vL = w0/2./Pi;    //Larmor frequency in Hz
    double width = 1./dwell/vL*1.e6; //spectral width in ppm
    double leftmax = width/2.;
    double rightmax = -1.*leftmax;

    //EFG tensor
    Eigen::VectorXcd R_Q(5);
    CQ=CQ/(4.*I*(2*I-1.)); //Applying spin-specific scaling
    double eq=-sqrt(6./w0)*CQ;
    R_Q(0) = R_Q(4) = eq/2.*etaQ/sqrt(6.);
    R_Q(1) = R_Q(3) = 0.;
    R_Q(2) = eq/2.;

    //Exchange matrix
    Eigen::MatrixXcd K(njumps,njumps);
    for(i=0;i<njumps;i++){
        K(i,i)=kex*dwell*(1.-1.*njumps);
        for(j=i+1;j<njumps;j++){
            K(i,j)=K(j,i)=kex*dwell;
        }
    }

    for(i=0;i<njumps;i++){
        aa[i]=aa[i]*Pi/180.;
        bb[i]=bb[i]*Pi/180.;
        gg[i]=gg[i]*Pi/180.;
    }

    //FID and spectrum
    Eigen::VectorXcd total_FID(np);
    Eigen::VectorXcd spectrum(np);
    total_FID.setZero();

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
    total_FID(0)=n_orient*gamma_steps;

    for(i=0;i<n_orient;i++){
        fgets(buffer, sizeof(buffer), cry);
        for(j=0;j<gamma_steps;j++){
            sscanf(buffer,"%lf %lf %lf",&al[i*gamma_steps+j],&bt[i*gamma_steps+j],&gm[i*gamma_steps+j]);
            gm[i*gamma_steps+j]=360.*j/gamma_steps;
    }}
    fclose(cry);

    #pragma omp parallel private(i,j,k,t)
    {

    Eigen::VectorXcd FID(np);

    //variable declarations
    double alpha, beta, gamma, vCS;
    double magic_angle = acos(sqrt(1./3.));
    Eigen::VectorXcd R2(5);
    Eigen::MatrixXcd D_MOL_MAS(5,5);
    Eigen::MatrixXcd D_MAS_LAB(5,5);
    Eigen::MatrixXcd D_PAS_MOL(5,5);

    //Frequency matrix
    Eigen::MatrixXcd Omega(njumps,njumps);
    Eigen::MatrixXcd Prop(njumps,njumps);
    Eigen::VectorXcd M(njumps);
    Eigen::VectorXcd M_temp(njumps);
    Omega.setIdentity();

    //powder averaging
    int orientation;
    for(orientation=omp_get_thread_num(); orientation<n_orient*gamma_steps; orientation=orientation+omp_get_num_threads()){

    //for(cosbeta=-1.+omp_get_thread_num()*2./(beta_steps); cosbeta<0.999; cosbeta=cosbeta+omp_get_num_threads()*2./(beta_steps)){
        beta=2.*Pi*bt[orientation]/360.;
        alpha=2.*Pi*al[orientation]/360.;
        gamma=2.*Pi*gm[orientation]/360.;
      //  for(alpha=0.; alpha<1.999*Pi; alpha=alpha+Pi*2./(alpha_steps)){
          //  for(gamma=0.; gamma<1.999*Pi; gamma=gamma+2.*Pi/(gamma_steps)){

                counter++;
                printf("%d/%d\n",counter,n_orient*gamma_steps);

                FID.setZero();
                FID(0)=1.;
                for(i=0;i<njumps;i++){
                    M(i)=M_temp(i)=0.25;
                }
                D_MOL_MAS=D2_matrix(alpha, beta, gamma);

                for(t=1;t<np;t++){//time
                D_MAS_LAB=D2_matrix(wr*t*dwell, magic_angle,0.);

                for(i=0;i<njumps;i++){ //loop over the different jump orientations
                D_PAS_MOL=D2_matrix(aa[i], bb[i], gg[i]);

                //Rotating the PAS tensor into the lab frame
                R2 = D_MAS_LAB*D_MOL_MAS*D_PAS_MOL*R_Q;

                //Calculating the quadrupole shift
                vCS =(4.*I*(I+1.)-3.)* (2*(real(R2(1))*real(R2(3))-imag(R2(1))*imag(R2(3)))+ (real(R2(0))*real(R2(4))- imag(R2(0))*imag(R2(4))));

                Omega(i,i)=vCS*1i*dwell;
                }
                Prop=Omega+K;

                M_temp=Prop.exp()*M;

                FID(t)=0.;
                for(i=0;i<njumps;i++){
                    M(i)=M_temp(i);
                    FID(t)=FID(t)+M(i);
                }

                //The value for that data point is added to the total FID.
                total_FID(t) = total_FID(t)+FID(t)*exp(-dwell*t*lb);

            }//end of time loop
    }//}} //end of powder averaging

    }//end of parallel bit

    free(al);
    free(bt);
    free(gm);

    //Fast Fourier Transform of the FID
    Eigen::FFT<double> fft;
    fft.fwd(spectrum,total_FID);

    //writing out the fid
    fp=fopen("fid.fid","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nTYPE=FID\nDATA\n",np,1./dwell);
    for(i=0;i<np;i++){
        fprintf(fp,"%f %f\n",real(total_FID(i)),imag(total_FID(i)));
    }
    fclose(fp);

    //writing the spectrum as a SPE file
    fp=fopen("spec.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nTYPE=SPE\nDATA\n",np,1./dwell,0.5/dwell+diso*vL/1000000.,vL/1000000.);
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

    fclose(fp);
    return 0;
}
