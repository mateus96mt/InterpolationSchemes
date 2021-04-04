#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double calculateError(const double *analyticSolution, const double *numericSolution, int nx, int errorType) {

    double *difAnalyticToNum = (double*) malloc((nx)*sizeof(double));;
    double sumNum = 0.0, sumDem = 0.0, maxDifAnalyticToNum = 0.0, maxAnalytic = 0.0;

    for (int i = 0; i < nx; i++) {
        difAnalyticToNum[i] = abs(analyticSolution[i] - numericSolution[i]);

        if (difAnalyticToNum[i] > maxDifAnalyticToNum) {
            maxDifAnalyticToNum = difAnalyticToNum[i];
        }

        if (analyticSolution[i] > maxAnalytic) {
            maxAnalytic = analyticSolution[i];
        }

        switch (errorType) {

            case 1:
                sumNum += abs(analyticSolution[i] - numericSolution[i]);
                sumDem += abs(analyticSolution[i]);
                break;

            case 2:
                sumNum += pow(analyticSolution[i] - numericSolution[i], 2);
                sumDem += pow(analyticSolution[i], 2);
                break;

            case 3:
                break;
        }

    }

    switch (errorType) {

        case 1:
            return sumNum / sumDem;

        case 2:
            return sqrt(sumNum / sumDem);

        case 3:
            return maxDifAnalyticToNum / maxAnalytic;

        default:
            return 0.0;
    }
}

int main (int c, char *argv[]) {
   double *x,*phi_in,*phi_exata,*phi_next,*phi,*f, param;
   double dx,dt,CFL,a,b,d,tf,k,f0, fiu,numa,numb,dena,denb,x0,xN;
   int j,op,fluxo,N,n;
   FILE *fp, *fp1; 

   /*
   printf("Informe o esquema upwind:\n");
   printf("(1) ADBQUICKEST \n");
   printf("(2) TOPUS \n");
   printf("(3) FSFL \n");
   printf("(4) SDPUS \n");
   printf("(5) EPUS \n");
   
   scanf("%d",&fluxo);
   printf("Parametro livre: \n");
   scanf("%lf",&param);
   

   printf("\nInforme o valor do CFL: \n");
   scanf("%lf",&CFL);
   */
   
   fluxo = atoi(argv[2]);
   param = atof(argv[3]);
   CFL = atof(argv[4]);
   int errorType = atoi(argv[5]);

//    printf("saida: %s\n", argv[1]);
//    printf("fluxo: %d\n", fluxo);
//    printf("param: %lf\n", param);
//    printf("CFL: %.3lf\n", CFL);
   
   
   if (fluxo == 1) {
       /*parametros do esquema ADBQUICKEST*/
     	numa = (2.0-(3.0*fabs(CFL))+(CFL*CFL));
     	dena = (7.0-(6.0*CFL)-(3.0*fabs(CFL))+(2.0*CFL*CFL));
     	a      = numa/dena;
     	numb = ((CFL*CFL)-4.0-(3.0*fabs(CFL))+(6*CFL));
     	denb = (-5.0-(3.0*fabs(CFL))+(2.0*CFL*CFL)+(6.0*CFL));
     	b    = numb/denb;
   }
   
   /*parametros da malha*/
   dx    = 0.005;
   N = 400;
   
   op = 1;
   if (op==1){
        tf = 1.0;
        x0    = 0.0;
        xN    = 2.0;
   }
   
   dt = CFL*dx; 
   k = tf/dt;

//    printf("dt = %f\n", dt);
//    printf("iterations = %d\n", (int)k);
   
   x         = (double*) malloc((N+2)*sizeof(double));
   phi_in    = (double*) malloc((N+2)*sizeof(double));
   phi_exata = (double*) malloc((N+2)*sizeof(double));
   phi_next  = (double*) malloc((N+2)*sizeof(double));
   phi       = (double*) malloc((N+2)*sizeof(double));
   f         = (double*) malloc((N+2)*sizeof(double));

   for(j=0;j<N;j++){
      x[j]        = 0.0;
      phi_in[j]   = 0.0;
      phi_exata[j]= 0.0;
      phi[j]      = 0.0;
      phi_next[j] = 0.0;
      f[j]        = 0.0;

   }
   
   x[0] = x0;   

   for(j=0;j<N;j++){
      x[j+1]=x[j]+dx;
   }

   char *fileName = malloc(sizeof(char) * 1024);
   sprintf(fileName,"%s/FINAL.data", argv[1]);
   fp1=fopen(fileName, "w");
//    if (fluxo==1){
//        sprintf(fileName,"%s/Eqdv_ADBQUICKEST.dat", argv[1]);
//        fp1=fopen(fileName,"w");
//    }
//    if (fluxo==2){ /*saida TOPUS*/
//        sprintf(fileName,"%s/Eqdv_TOPUS_%.2f.dat", argv[1], param);
//        fp1=fopen(fileName,"w");
//    }
//    if (fluxo==3){ /*saida FSFL*/
//        sprintf(fileName,"%s/Eqdv_FSFL_%.2f.dat", argv[1], param);
//        fp1=fopen(fileName,"w");
//    }
//    if (fluxo==4){ /*saida SDPUS*/
//        sprintf(fileName,"%s/Eqdv_SDPUS_%.2f.dat", argv[1], param);
//        fp1=fopen(fileName,"w");
//    }
//    if (fluxo==5){ /*saida EPUS*/
//        sprintf(fileName,"%s/Eqdv_EPUS_%.2f.dat", argv[1], param);
//        fp1=fopen(fileName,"w");
//    }
   
   
   /*CONDICAO INICIAL*/
   if (op==1)
   {  
      for (j=0;j<N;j++)
      { 
         if ((x[j] >= 0.0) && (x[j] < 0.2))
         {
            phi_in[j] = exp(-log(50.0)*((x[j]-0.15)/(0.05))*((x[j]-0.15)/(0.05)));
         }
         else if ((x[j] > 0.3) && (x[j] < 0.4))
         {
            phi_in[j] = 1.0;
         }
         else if ((x[j] > 0.5) && (x[j] < 0.55))
         {
            phi_in[j] = 20.0*x[j]-10.0;
         }
         else if ((x[j] >= 0.55) && (x[j] <0.6))
         {
            phi_in[j] = 12.0-20.0*x[j];
         }
         else if ((x[j] > 0.7) && (x[j] < 0.8))
         {
            phi_in[j] = sqrt(1.0-((x[j]-0.75)/(0.05))*((x[j]-0.75)/(0.05)));
        }
         else if ((x[j] > 0.8) && (x[j] <= 2.0))
         {
            phi_in[j] = 0.0;
         }
      }
   }
   
   
   /* Soluca Exata */
    if (op==1)
   {  
      for (j=0;j<N;j++)
      {
          if ((x[j]-tf < 0.2) && (x[j]-tf >=0.0))
          {
             phi_exata[j] = exp(-log(50.0)*((x[j]-tf-0.15)/(0.05))*((x[j]-tf-0.15)/(0.05)));
          }

          if ((x[j]-tf>0.3) && (x[j]-tf<0.4))
          {
             phi_exata[j] = 1.0;
          }
          if ((x[j]-tf > 0.5) && (x[j]-tf <0.55))
          {
             phi_exata[j] = 20.0*(x[j]-tf)-10.0;
          }
          if ((x[j]-tf >= 0.55) && (x[j]-tf <0.6))
          {
             phi_exata[j] = 12.0-20.0*(x[j]-tf);
          }
          if ((x[j]-tf > 0.7) && (x[j]-tf < 0.8))
          {
             phi_exata[j] = sqrt(1.0-((x[j]-tf-0.75)/(0.05))*((x[j]-tf-0.75)/(0.05)));
          }
          if ((x[j]-tf >0.8) && (x[j]-tf <=2.0))
          {
             phi_exata[j] = 0.0;
          }
      }
   } 
   
   
   /* Soluca Numerica*/
   for (j=0;j<=N;j++)
   {
      phi[j] = phi_in[j];
   }

   
   for (n=1;n<=(int)k;n++) {
       /*ADBQUICKEST*/
        if (fluxo==1)
        {   
            phi[N] = 0.0;
            for (j=1;j<N;j++)
            {
                if ((phi[j+1]-phi[j-1]) == 0.0)
                {
                   f[j] = phi[j]; 
                }
                else
                {
                   fiu = (phi[j]-phi[j-1])/(phi[j+1]-phi[j-1]);
                   /*phiU = phi[j]; phiR = phi[j-1]; phiD = phi[j+1]*/

                   if ((fiu >= a) && (fiu <= b)) {
                   
                   	f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * (fiu + 0.5 * (1.0 -fabs(CFL)) *(1.0-fiu) -(1.0/6.0) *(1.0-pow(CFL,2.0)) * (1.0-2.0*fiu));

                   } 
                   else if ((fiu >= 0.0) && (fiu < a)) {
                         f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * ((2.0 - CFL) * fiu);
                   
                   }
                   else if ((fiu > b) && (fiu <= 1.0)) {
                         f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * (1.0 - CFL + CFL * fiu);
                   }
                   else {
                      f[j] = phi[j];
                   }

                }
                f[0] = 0.0;
                phi_next[0] = 0.0;
            }
        }
      
         //Esquema TOPUS
         else if (fluxo==2)
         {
             phi[N] = 0.0;
            for (j=1;j<N;j++)
            {
                if ((phi[j+1]-phi[j-1]) == 0.0)
                {
                   f[j] = phi[j]; 
                }
                else
                {
                    fiu = (phi[j]-phi[j-1])/(phi[j+1]-phi[j-1]);
                   /*phiU = phi[j]; phiR = phi[j-1]; phiD = phi[j+1]*/

                   if ((fiu >= 0.0) && (fiu <= 1.0)) {
                   
                    
                       
                   	f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * (
                        (param * pow(fiu, 4))
                   + (((-2 * param) + 1) * pow(fiu, 3))
                   + ((((5 * param) - 10) / 4) * pow(fiu, 2))
                   + (((-param + 10) / 4) * fiu)
                        
                    );

                   } 
                   
                   else {
                      f[j] = phi[j];
                   }
                }
                f[0] = 0.0;
                phi_next[0] = 0.0;
            }
        }
        
        //Esquema FSFL
        else if (fluxo==3)
         {
             phi[N] = 0.0;
            for (j=1;j<N;j++)
            {
                if ((phi[j+1]-phi[j-1]) == 0.0)
                {
                   f[j] = phi[j]; 
                }
                else
                {
                 fiu = (phi[j]-phi[j-1])/(phi[j+1]-phi[j-1]);
                   /*phiU = phi[j]; phiR = phi[j-1]; phiD = phi[j+1]*/

                   if ((fiu >= 0.0) && (fiu <= 1.0)) {
                   
                    
                       
                   	f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * ((((-2 * param) + 4) * pow(fiu, 4))
                   + (((4 * param) - 8) * pow(fiu, 3))
                   + ((((-5 * param) + 8) / 2) * pow(fiu, 2))
                   + (((param + 2) / 2) * fiu));

                   } 
                   
                   else {
                      f[j] = phi[j];
                   }
                }
                f[0] = 0.0;
                phi_next[0] = 0.0;
            }
        }
        
        //Esquema SDPUS
        else if (fluxo==4)
         {
             phi[N] = 0.0;
            for (j=1;j<N;j++)
            {
                if ((phi[j+1]-phi[j-1]) == 0.0)
                {
                   f[j] = phi[j]; 
                }
                else
                {
                 fiu = (phi[j]-phi[j-1])/(phi[j+1]-phi[j-1]);
                   /*phiU = phi[j]; phiR = phi[j-1]; phiD = phi[j+1]*/

                   if ((fiu >= 0.0) && (fiu <= 1.0)) {
                   
                    
                       
                   	f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * (((-24.0 + (4 * param)) * pow(fiu, 6))
                   + ((68.0 - (12.0 * param)) * pow(fiu, 5))
                   + ((-64.0 + (13 * param)) * pow(fiu, 4))
                   + ((20.0 - (6 * param)) * pow(fiu, 3))
                   + (param * pow(fiu, 2))
                   + fiu);

                   } 
                   
                   else {
                      f[j] = phi[j];
                   }
                }
                f[0] = 0.0;
                phi_next[0] = 0.0;
            }
        }
        
        //Esquema EPUS
        else if (fluxo==5)
         {
             phi[N] = 0.0;
            for (j=1;j<N;j++)
            {
                if ((phi[j+1]-phi[j-1]) == 0.0)
                {
                   f[j] = phi[j]; 
                }
                else
                {
                 fiu = (phi[j]-phi[j-1])/(phi[j+1]-phi[j-1]);
                   /*phiU = phi[j]; phiR = phi[j-1]; phiD = phi[j+1]*/

                   if ((fiu >= 0.0) && (fiu <= 1.0)) {
                   
                    
                       
                   	f[j] = phi[j-1] + (phi[j+1] - phi[j-1]) * ((-4 * (param - 24.0) * pow(fiu, 8))
                   + (16.0 * (param - 23.0) * pow(fiu, 7))
                   + ((528.0 - (25 * param)) * pow(fiu, 6))
                   + (((19.0 * param) - 336.0) * pow(fiu, 5))
                   + ((80.0 - (7.0 * param)) * pow(fiu, 4))
                   + (param * pow(fiu, 3))
                   + fiu);

                   } 
                   
                   else {
                      f[j] = phi[j];
                   }
                }
                f[0] = 0.0;
                phi_next[0] = 0.0;
            }
        }

        for (j=1;j<N;j++)
        
            phi_next[j] = phi[j]+CFL*(f[j-1]-f[j]);    
       
        /* atualiza valores */
	for (j=0;j<=N;j++)
            phi[j] = phi_next[j];
        
   } /* fim do loop de tempo */
     
     		 
   fprintf(fp1,"jump line %d\n", 1);
   fprintf(fp1,"jump line %d\n", 2);
   for(j=0;j<=N;j++)
         fprintf(fp1,"%1.5f %1.5f %1.5f\n", x[j], phi_exata[j], phi_next[j]);
     
   
   printf("%.10f\n", calculateError(phi_exata, phi_next, N+2, errorType));
   
   free(x);
   free(phi_in);  
   free(phi);
   free(phi_exata);   
   free(phi_next);  
   free(f);  
   
   fclose(fp1);
   
   return (1);
}
