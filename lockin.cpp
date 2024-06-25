//path=D:\Roch\Bcc55\bin
//bcc32 -WDR lockin

////////////////
// LOCKIN.CPP //
////////////////

#include <math.h>

extern "C" __declspec(dllexport) int PSD_ref(

	double *ref,    //reference mesured       [V]
	double *nod,    //nodes table             [rational index]
	int sn,         //number of samples       [integer]

	double offset,  //user defined offset     [V]

	double *reft,   //reference pulsation     [rational index]
	double *refp    //reference phase         [rational index]	

	){
	
	int nodn,i;
	double x,y1,y2,sx,sy,sxx,sxy;

	//REMOVE USER DEFINED OFFSET
	for(i=0;i<sn;i++) ref[i]-=offset;
		
	//INTERPOLATE RISING EDGE                                            
	nodn=0;
	for(i=1;i<sn;i++) 
		if(ref[i-1]<=0 && ref[i]>=0)
			nod[nodn++]=(double)(i++-1)+ref[i-1]/(ref[i-1]-ref[i]);   

	//PERIOD AND PHASE
	sx=0.0;
	sy=0.0;
	sxx=0.0;
	sxy=0.0;
	for(i=0;i<nodn;i++){
		sx+=i;
		sy+=nod[i];
		sxx+=i*i;
		sxy+=i*nod[i];
	}
	*reft=(sxy-sx*sy/nodn)/(sxx-sx*sx/nodn);
	*refp=(sx*sxy-sxx*sy)/(sx*sx-nodn*sxx);

	return 0;
}

extern "C" __declspec(dllexport) int PSD_sig(

	double *sig,    //signal mesured                    [V]
	double dt,      //sampling intervals                [s]
	int sn,         //number of samples                 [integer]

	double reft,    //reference pulsation               [rational index]
	double refp,    //reference phase                   [rational index]	

	double tc,      //Low pass filter time constant     [s]
	double *vx,     //Low-pass filters outputs (PHASE)  [Vpp]
	double *vy,     //Low-pass filters outputs (QUADR)  [Vpp]
	int fn          //number of LP filters              [integer]
		
	){     

	int i,j;
	double t=dt/(dt+tc);
	double refw=6.2831853071795864/reft;

	for(i=0;i<sn;i++){	
		
		vx[0]+=2.0*(sin(refw*(i-refp))*sig[i]-vx[0])*t;
		for(j=1;j<fn;j++) vx[j]+=(vx[j-1]-vx[j])*t;
		
		vy[0]+=2.0*(cos(refw*(i-refp))*sig[i]-vy[0])*t;
		for(j=1;j<fn;j++) vy[j]+=(vy[j-1]-vy[j])*t;
		
		sig[i]=vx[fn-1];
	}

	return 0;
}
