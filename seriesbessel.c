#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

int nsum;
float *zeros;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

float bessj0(float x)
{
	float ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

float bessj1(float x)
{
	float ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

main()
{
     FILE *fp1,*fp2;
     float rf,rr,rc,*f,kappaT,deltax,a,*prefact,*topo,*J0rr,*J0rc,*J1rr,*J1rc;
     float *J0rf,*J1rf,*Intrf,*J0r0,*J1r0,*Intrr,*Intrc;
     int lattice_size_x,i,nsum;

     fp1=fopen("besselzerosintegral","r");
     fp2=fopen("volcanicconeoutput","w");
     kappaT=0.0001;  /* dimensionless – scaled to a^2 */
     lattice_size_x=300; 
     deltax=10.0;  
     a=lattice_size_x*deltax;
     rf=200.0/a;
     rr=500.0/a;
     rc=1500.0/a;
     topo=vector(1,lattice_size_x);
     f=vector(1,lattice_size_x);
     zeros=vector(1,300);
     J0rc=vector(1,300);
     J1rc=vector(1,300);
     Intrc=vector(1,300);
     J0rr=vector(1,300);
     J1rr=vector(1,300);
     Intrr=vector(1,300);
     J0rf=vector(1,300);
     J1rf=vector(1,300);
     Intrf=vector(1,300);
     prefact=vector(1,300);
     for (nsum=1;nsum<=300;nsum++)
      {fscanf(fp1,"%f %f %f %f",&zeros[nsum],&Intrf[nsum],&Intrr[nsum],&Intrc[nsum]);
       J0rc[nsum]=bessj0(rc*zeros[nsum]);
       J1rc[nsum]=bessj1(rc*zeros[nsum]);
       J0rr[nsum]=bessj0(rr*zeros[nsum]);
       J1rr[nsum]=bessj1(rr*zeros[nsum]);
       J0rf[nsum]=bessj0(rf*zeros[nsum]);
       J1rf[nsum]=bessj1(rf*zeros[nsum]);}
     for (i=1;i<=lattice_size_x;i++)
      {if (i*deltax/a<rf) 
        f[i]=1+(rf-rr)/(rc-rr);
       else if (i*deltax/a<rr)
        f[i]=1+(i*deltax/a-rr)/(rc-rr);
       else if (i*deltax/a<rc)
        f[i]=1+(rr-i*deltax/a)/(rc-rr);
       else f[i]=0;
       /* print initial condition first */
       fprintf(fp2,"%f %12.10f\n",i*deltax,f[i]);}
     for (i=1;i<=lattice_size_x;i++)
      topo[i]=0; 
     for (nsum=1;nsum<=300;nsum++)
      prefact[nsum]=2.0/zeros[nsum]*exp(-kappaT*zeros[nsum]*zeros[nsum])/
       (bessj1(zeros[nsum])*bessj1(zeros[nsum]));
     for (i=1;i<=lattice_size_x;i++)
      {for (nsum=1;nsum<=300;nsum++)
        topo[i]+=prefact[nsum]*bessj0(zeros[nsum]*i*deltax/a)*(rf*rf/(rc-rr)*
         J1rf[nsum]-2*rr*rr/(rc-rr)*J1rr[nsum]+(1+rr/(rc-rr))*
         rc*J1rc[nsum]+2*zeros[nsum]/(rc-rr)*Intrr[nsum]-zeros[nsum]/
         (rc-rr)*(Intrf[nsum]+Intrc[nsum]));
       fprintf(fp2,"%f %12.10f\n",i*deltax,topo[i]);} 
   fclose(fp1);
   fclose(fp2); 
}
