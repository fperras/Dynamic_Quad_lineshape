#include <math.h>
#include <complex>
#include "Eigen/Dense"
static double Pi = 3.1415926535897932384626433;
using namespace std::complex_literals;

double d2mn(int m, int n, double beta){
/* reduced Wigner rotation matrix elements */
if(n==2){
    if(m==2){
        return pow(1.+cos(beta),2)/4.;
    }
    else if(m==1){
        return (1.+cos(beta))/2.*sin(beta);
    }
    else if(m==0){
        return sqrt(3./8.)*pow(sin(beta),2);
    }
    else if(m==-1){
        return (1.-cos(beta))/2.*sin(beta);
    }
    else if(m==-2){
        return pow(1.-cos(beta),2)/4.;
    }
}
else if(n==1){
    if(m==2){
        return -1.*(1.+cos(beta))/2.*sin(beta);
    }
    else if(m==1){
        return pow(cos(beta),2)-(1.-cos(beta))/2.;
    }
    else if(m==0){
        return sqrt(3./8.)*sin(2.*beta);
    }
    else if(m==-1){
        return (1.+cos(beta))/2.-pow(cos(beta),2);
    }
    else if(m==-2){
        return (1.-cos(beta))/2.*sin(beta);
    }
}
else if(n==0){
    if(m==2){
        return sqrt(3./8.)*pow(sin(beta),2);
    }
    else if(m==1){
        return -1.*sqrt(3./8.)*sin(2.*beta);
    }
    else if(m==0){
        return (3.*pow(cos(beta),2)-1)/2.;
    }
    else if(m==-1){
        return sqrt(3./8.)*sin(2.*beta);
    }
    else if(m==-2){
        return sqrt(3./8.)*pow(sin(beta),2);
    }
}
else if(n==-1){
    if(m==2){
        return -1.*(1.-cos(beta))/2.*sin(beta);
    }
    else if(m==1){
        return (1.+cos(beta))/2.-pow(cos(beta),2);
    }
    else if(m==0){
        return -1.*sqrt(3./8.)*sin(2.*beta);
    }
    else if(m==-1){
        return pow(cos(beta),2)-(1.-cos(beta))/2.;
    }
    else if(m==-2){
        return (1.+cos(beta))/2.*sin(beta);
    }
}
else if(n==-2){
    if(m==2){
        return pow(1.-cos(beta),2)/4.;
    }
    else if(m==1){
        return -1.*(1.-cos(beta))/2.*sin(beta);
    }
    else if(m==0){
        return sqrt(3./8.)*pow(sin(beta),2);
    }
    else if(m==-1){
        return -1.*(1.+cos(beta))/2.*sin(beta);
    }
    else if(m==-2){
        return pow(1.+cos(beta),2)/4.;
    }
}

else
    exit(0);

}

std::complex<double> D2mn(int m,int n, double alpha, double beta, double gamma)
/*   Wigner rotation matrix element */
{
    return d2mn(m,n, beta)*(cos(-alpha*m-gamma*n)+1i*sin(-alpha*m-gamma*n));
}

Eigen::MatrixXcd D2_matrix(double alpha, double beta, double gamma)
/*   Wigner rotation matrix */
{
    int i, j;
    Eigen::MatrixXcd D2(5,5);
    //Calculating the PAS to LAB Wigner rotation matrix
            for(i=-2;i<=2;i++){
                for(j=-2;j<=2;j++){
                    D2(i+2,j+2) = D2mn(i,j,gamma,beta,alpha);
            }}
    return D2;
}

