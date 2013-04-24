#include "mex.h"
#include "math.h"

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

/* summation function */
double sum(double a[], int size)
{
    double summed=0;
    int i;
    for (i=0;i<size;i++) 
    {
	summed=summed+a[i];
    }
    return summed;
}

void lfSource(double alpha[], double epsi[], double Tc, double fs, double Tp, double Te, double Ta, double EE)
{

   /*// Initialise */
   float TolFun=0.0000001;
   int MaxIter=100, count=1, i, count_alpha=1;
   double Tb=Tc-Te, omega_g=pi/Tp, eps0, change=1.0, change_alpha=1.0, f_eps, f_eps_prime, a0=0, A2, f_a, f_a_prime,e0,part1,part2,part3,partAtan,part4,A1,E0;

   /*// Solve epsilon using newton raphson method */
   eps0=1/Ta;
   while(count<=MaxIter&fabs(change)>TolFun)
   {
     f_eps = (eps0*Ta-1.0+exp(-eps0*Tb));
     f_eps_prime = (Ta-Tb*exp(-eps0*Tb));
     change=f_eps/f_eps_prime;
     eps0=eps0-change;
     eps0=eps0-(eps0*Ta-1+exp(-eps0*Tb))/(Ta-Tb*exp(-eps0*Tb));
     count++;
   }
   epsi[0]=eps0;

   /* Solve alpha - Do newton-raphson iterations */
   e0==-EE/(exp(a0*Te)*sin(omega_g*Te));

   A2 =(-EE/((epsi[0]*epsi[0])*Ta))*(1-exp(-epsi[0]*Tb)*(1+epsi[0]*Tb));

   while(count_alpha<=MaxIter&fabs(change_alpha)>TolFun)
   {
     part1=sqrt((a0*a0)+(omega_g*omega_g));
     partAtan=2*atan((sqrt((a0*a0)+(omega_g*omega_g))-a0)/omega_g);
     part2=sin(omega_g*Te-partAtan);
     part3=omega_g*exp(-a0*Te)-((A2/EE)*((a0*a0)+(omega_g*omega_g))*sin(omega_g*Te));
     part4=(sin(omega_g*Te)*(1-2*a0*A2/EE)-omega_g*Te*exp(-a0*Te));
     a0=a0-((part1*part2)+part3)/part4 ;

     part1=sqrt((a0*a0)+(omega_g*omega_g));
     partAtan=2*atan((sqrt((a0*a0)+(omega_g*omega_g))-a0)/omega_g);
     part2=sin(omega_g*Te-partAtan);
     part3=omega_g*exp(-a0*Te)-((A2/EE)*((a0*a0)+(omega_g*omega_g))*sin(omega_g*Te));
     part4=(sin(omega_g*Te)*(1-2*a0*A2/EE)-omega_g*Te*exp(-a0*Te));

     a0=a0-((part1*part2)+part3)/part4 ;

     count_alpha++;
   }
   
   alpha[0]=a0;
    /* print to see area balance  */
   /* E0=-EE/(exp(a0*Te)*sin(omega_g*Te)); */
   /* A1=(E0*exp(a0*Te)/sqrt((a0*a0)+(omega_g*omega_g)))*sin(omega_g*Te-(2*atan((sqrt((a0*a0)+(omega_g*omega_g))-a0)/omega_g)))+E0*omega_g/((a0*a0)+(omega_g*omega_g)); */
   /* printf("%d\n",count_alpha;)*/
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*// inputs*/
    double Tc, fs, Tp, Te, Ta, EE;

    /*// outputs*/
    double *alpha, *epsi;


    /*  check for proper number of arguments */
    if(nrhs<6)
    {
        mexErrMsgTxt("Six inputs needed.");
    }
    if(nlhs!=2)
    {
        mexErrMsgTxt("Two outputs required.");
    }

    /*  assign inputs  */
    Tc = mxGetScalar(prhs[0]);
    fs = mxGetScalar(prhs[1]);
    Tp = mxGetScalar(prhs[2]);
    Te = mxGetScalar(prhs[3]);
    Ta = mxGetScalar(prhs[4]);
    EE = mxGetScalar(prhs[5]);
      
    /* Assign pointers to each output. */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

    
    alpha   = mxGetPr(plhs[0]);
    epsi   = mxGetPr(plhs[1]);

  
    /*  call the C subroutine */
    lfSource(alpha, epsi, Tc, fs, Tp, Te, Ta, EE);
  
}
