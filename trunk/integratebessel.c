#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

float bessj0(float x);

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

float func(float x)
{
        return x*x*bessj0(zeros[nsum]*x);
}

#define EPS 1.0e-6
#define JMAX 20

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

#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

float qsimp(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	return 0.0;
}
#undef EPS
#undef JMAX

main()
{    FILE *fp0,*fp1;
     double a,zero=0.0,s1,s2,s3,rf,rr,rc;

     fp0=fopen("besselzeros","r");
     fp1=fopen("besselzerosintegral","w");
     zeros=vector(1,300);
     a=3000.0;
     rf=200.0/a;
     rr=500.0/a;
     rc=1500.0/a;
     for (nsum=1;nsum<=300;nsum++)
      {fscanf(fp0,"%f",&zeros[nsum]);
       s1=qsimp(func,zero,rf);
       s2=qsimp(func,zero,rr);
       s3=qsimp(func,zero,rc);
       fprintf(fp1,"%f %12.10f %12.10f %12.10f\n",zeros[nsum],s1,s2,s3);}
     fclose(fp0);
     fclose(fp1);
}
