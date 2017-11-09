/*------------------------------------------------------------------------------*/
/***		Part 1 (of 3)	S-function definition of Block Parameter		  ***/
/*------------------------------------------------------------------------------*/
/*		‡@ Name of the S-function Block                                         */
/*------------------------------------------------------------------------------*/
#define S_FUNCTION_NAME 		IM_AG_Fit_Sfunc

/*------------------------------------------------------------------------------*/
/*		‡A Number of Inputs														*/
/*------------------------------------------------------------------------------*/
#define INPUT_PORT_WIDTH		8

/*------------------------------------------------------------------------------*/
/*		‡B Number of Outputs													*/
/*------------------------------------------------------------------------------*/
#define OUTPUT_PORT_WIDTH		30

/*------------------------------------------------------------------------------*/
/*		‡C Sample time of S-function block execution [sec]						*/
/*		?iIf it's continuous, set CONTINUOUS_SAMPLE_TIME?j                      */
/*------------------------------------------------------------------------------*/
#define SAMPLE_TIME				0.001

/*------------------------------------------------------------------------------*/
/***		Part 1 End of Definition of S-function block parameters			  ***/
/*------------------------------------------------------------------------------*/

#define S_FUNCTION_LEVEL 2

/*
 * Input Parameters
 */
#include "simstruc.h"
#include <stdio.h>

#define 	NUM_PARAMS		5
#define 	Rs1             (ssGetSFcnParam(S,0))
#define 	Rr1             (ssGetSFcnParam(S,1))
#define 	Ls1             (ssGetSFcnParam(S,2))
#define 	Lr1             (ssGetSFcnParam(S,3))
#define 	Msr1            (ssGetSFcnParam(S,4))

#define 	Rs              ((float)mxGetPr(Rs1)[0])
#define 	Rr              ((float)mxGetPr(Rr1)[0])
#define 	Ls              ((float)mxGetPr(Ls1)[0])
#define 	Lr              ((float)mxGetPr(Lr1)[0])
#define 	Msr             ((float)mxGetPr(Msr1)[0])




/*------------------------------------------------------------------------------*/
/***		Part 2 (of 3)	Description of User's Code						  ***/
/*------------------------------------------------------------------------------*/
/*		‡@ Redifinition of Data Type and Macro                                  */
/*------------------------------------------------------------------------------*/
typedef unsigned long		uint32;
typedef unsigned short		uint16;
typedef unsigned char		uchar8;
typedef unsigned char		BOOL;
typedef signed   long		sint32;
typedef signed   short		sint16;
typedef signed   char		schar8;

/*------------------------------------------------------------------------------*/
/*		‡B Global Variable														*/
/*------------------------------------------------------------------------------*/
#define CHROM_NB            5          //Number of chromozomes i.e. of parameters to be identified
#define POPULATION_SIZE     300
#define AVG_KING_FERT       1.13       //Average childs number of the best individual (1 to 2)
#define P_CROSS             0.2        //Chromozome crossing probability
#define CROS_MIX_RATE       0.1        //Not used Stenght of the cossing of two chromo (0 strongest to 1 no chross)
#define P_MUT               0.0005
#define RANDOM_MAX          20000      //max value of int random number used for float rand --> resolution     
#define MAX_GEN             500000     //calculation maximum time: 30s max
#define AVG_MIN_END         0.07       //stop algo when fitness average < AVG_MIN_END
#define DISP_MIN_END        0.3        //stop algo when dispertion < DISP_MIN_END
#define PMUT_END            (80/SAMPLE_TIME)        //stop mutations progressivley at x sec/SAMPLE_TIME 
#define FIT_ARR_SIZE        3          //number n of the fitness sigma 0 to n sum

/********************************************************************************/
/*************************	Structures declaration 	*****************************/
/********************************************************************************/
typedef struct 
{
    float           Psel;                   //selection probability
    float           Fit;                    //fitness value 
    float           Chrom[CHROM_NB];        //target of the identification
} Individual;

typedef struct 
{
    float           Real;                   //selection probability
    float           Imag;                   //fitness value 
} Complex;

struct
{
    float AllelThdMin[CHROM_NB];
    float AllelThdMax[CHROM_NB];
} ChromThd;

/********************************************************************************/
/*************************	Global variables		*****************************/
/********************************************************************************/
Individual Population[POPULATION_SIZE];
float Average=3000, Dispersion=300;
static float Fit[FIT_ARR_SIZE];
static float Fit2[FIT_ARR_SIZE];
static float Fit3[FIT_ARR_SIZE];
long int Generation=0;
long int mut=0, cros=0, Term=0;

/*------------------------------------------------------------------------------*/
/*		‡C Function																*/
/*------------------------------------------------------------------------------*/
//Fonction de trie: meilleur individu (fitness min) en premiere positiom
void SortPopulation_Best_First(Individual* Pop)
{
    int i,j, indMin;
    Individual tempo;

   for (i = 0; i < POPULATION_SIZE - 1; i++)
   {
      indMin = i;
      
      for (j = i + 1; j < POPULATION_SIZE; j++)
      {
          if (Pop[j].Fit < Pop[indMin].Fit)
             indMin = j;
      }
       if (indMin != i)
       {
          tempo      = Pop[i];
          Pop[i] = Pop[indMin];
          Pop[indMin] = tempo;
       }
   }
}

//Calculate the selection probability of each indiviudal function of its rank and of the fertility of the king
void SelectionProbaCalc(Individual* Pop)
{
    int i;
    
    for(i=0; i<POPULATION_SIZE; i++)
    {
        Pop[i].Psel = (AVG_KING_FERT - (i*(2*AVG_KING_FERT - 2)/(POPULATION_SIZE-1))) / POPULATION_SIZE;
    }
}

//generate a random chromosome within initialization Min and Max chromosome limits
void NewChromosomeRandom(float* ChromTmp)
{   
    int i,tmp;
    
    for(i=0; i<CHROM_NB; i++)
    {
        tmp = rand()%(RANDOM_MAX);
        ChromTmp[i] = (((float)tmp/(float)(RANDOM_MAX)) * (ChromThd.AllelThdMax[i]-ChromThd.AllelThdMin[i])) + ChromThd.AllelThdMin[i];
    }
}

//Cross Chomozome A and B and return two new chromozomes
void ChromoCross(float* ChromA, float* ChromB, float* ChromARes, float* ChromBRes)
{
    int i,tmp;
    float prob;
   
    for(i=0; i<CHROM_NB; i++)
    {
        //Random number between 0 and 1
        tmp = rand()%(RANDOM_MAX);
        prob = ((float)tmp/(float)(RANDOM_MAX)) * 1;
        //prob = CROS_MIX_RATE;
        ChromARes[i] = prob*ChromA[i] + (1-prob)*ChromB[i];
        ChromBRes[i] = prob*ChromB[i] + (1-prob)*ChromA[i];
    }
}
//Cross Chomozome A and B and return two new chromozomes
void ChromoCross2(float* ChromA, float* ChromB, float* ChromARes, float* ChromBRes)
{
    int i,tmp;
    float prob=1;
    
    //Random number between 0 and 1
    tmp = rand()%(RANDOM_MAX);
    //prob = ((float)tmp/(float)(RANDOM_MAX)) * 1.2;
    for(i=0; i<CHROM_NB; i++)
    {

        if((i%2)==0)
        {
            ChromARes[i] = prob*ChromB[i];
            ChromBRes[i] = ChromA[i];
        }
        else
        {
            ChromARes[i] = ChromA[i];
            ChromBRes[i] = prob*ChromB[i];
        }
    }
}
//random +1 or -1
int Coin_Toss()
{
    int res;
    res = rand()% 2;  
    if(res == 0)
        res = -1;
    return(res);
}

//Calculates the absolute value of a float
float Abs(float A)
{
	if(A>=0) return A;
	else return(-A);	
}

//Cross Chomozome A and B and return two new chromozomes
void ChromoMutate(float* Chrom, float* ChromRes)
{
    int i,tmp, Allel_Sel;
    float prob;
    
    for(i=0; i<CHROM_NB; i++)
    {
        //Select a random Allel to be mutate
        Allel_Sel = i;//rand()% (CHROM_NB); 
        //Select a number between 0 and 1
        tmp = rand()%(RANDOM_MAX);
        prob = ((float)tmp/(float)(RANDOM_MAX)) * 1;
        ChromRes[Allel_Sel] = Chrom[Allel_Sel] + prob*Coin_Toss()*Abs(prob*Chrom[Allel_Sel]-ChromThd.AllelThdMax[Allel_Sel]); 
        ChromRes[Allel_Sel] = Abs(ChromRes[Allel_Sel]);
    }
}

//Prability selection of individual
int Proba_Sel(float Proba)
{
    float Randf;
    int tmp, res;
        
    //random float beween 0 and 1 
    tmp = rand()%(RANDOM_MAX);
    Randf = ((float)tmp/(float)(RANDOM_MAX)) * 1;
    
    if(Randf < Proba)
        res = 1;
    else
        res = 0;
    return (res);
}

// Initialization function
void Initialize_Pop(Individual* Pop)
{
    int i;
    
    mut = 0;
    cros = 0;
    Average=3000; 
    Dispersion=300;
    Generation=0;
    Term=0;
    
    ChromThd.AllelThdMax[0] = 4;        //R1
    ChromThd.AllelThdMin[0] = 0.1;      //R1
    ChromThd.AllelThdMax[1] = 40;       //X1
    ChromThd.AllelThdMin[1] = 0.1;      //X1
    ChromThd.AllelThdMax[2] = 4;        //R2
    ChromThd.AllelThdMin[2] = 0.1;      //R2
    ChromThd.AllelThdMax[3] = 40;       //X2
    ChromThd.AllelThdMin[3] = 0.1;      //X2
    ChromThd.AllelThdMax[4] = 400;      //Xm
    ChromThd.AllelThdMin[4] = 1;        //Xm
    
    for(i=0; i<POPULATION_SIZE; i++)
    {
        Pop[i].Psel = 1;
        Pop[i].Fit = 2;
        NewChromosomeRandom(Pop[i].Chrom);   
    }
    for(i=0; i<FIT_ARR_SIZE; i++)
    {
        Fit[i] = 0;
        Fit2[i]= 0;
        Fit3[i]= 0;
    }
    
}
float Maxx(float a, float b)
{
    if(a>=b)
        return a;
    else
        return b;
}

Complex ComplexDiv(Complex a, Complex b)
{
    Complex Res;
    float Den;
    
    Den = b.Real*b.Real + b.Imag*b.Imag;
    Res.Real = (a.Real*b.Real + a.Imag*b.Imag)/Den;
    Res.Imag = (a.Imag*b.Real - a.Real*b.Imag)/Den;
    
    return (Res);
}

float Sigma(float* Input, float new, int n)
{
    float Sum=0;
    int i;
    
    for(i=0; i<(n-1); i++)
    {
        Input[i] = Input[i+1]; 
        Sum += Input[i];
    }
    Input[n-1] = new;
    Sum += Input[i];
    return (Sum);
}

int Check_Chromo(float* Chrom)
{
    int i,test;
    test=0;
    for(i=0; i<CHROM_NB; i++)
    {
        if((Chrom[i]>ChromThd.AllelThdMax[i])||(Chrom[i]<=ChromThd.AllelThdMin[i]))
            test = 1;
    }
    return (test);
}

float StandardDev(Individual* Pop, float Avg)
{
    int i;
    float SD=0;
    
    for(i=0; i<POPULATION_SIZE; i++)
    {
        SD += (Pop[i].Fit - Avg) * (Pop[i].Fit - Avg);
    }
    SD = sqrt(SD/POPULATION_SIZE);
    
    return (SD);
}
float LowPassDo(float input, float Cut, int Id)
{
    static float out_old[10]={5,0,0,0,0,0,0,0,0,0};
    float out=0, div;
    
    div = (float)SAMPLE_TIME / Cut;
    out = div*input + (1- div)*out_old[Id];
    out_old[Id] = out;
    
    return (out);
}
float LowPassDo2(float input, float Cut)
{
    static float out_old=5;
    float out, div;
    
    div = (float)SAMPLE_TIME / Cut;
    out = div* input + (1- div)*out_old;
    out_old = out;
    
    return (out);
}
float Avoid_zero(float input, float min)
{
    float lim;
    
    if(Abs(input)<=min)
        lim=min;
    else
        lim=input;
    return(lim);
}


/*------------------------------------------------------------------------------*/
/***		Part 2 End of User Code description     						  ***/
/*------------------------------------------------------------------------------*/



/*------------------------------------------------------------------------------*/
/***		Part 3 (of 3) Description of the actual work of S-function block  ***/
/*------------------------------------------------------------------------------*/
/*		‡@ Initiation part 1 of S-function block's inner state					*/
/*		?iExecuted only once when Simulation starts.?j							*/
/*------------------------------------------------------------------------------*/
void My_mdlStart(SimStruct *S)
{
    Initialize_Pop(Population);
}
/*------------------------------------------------------------------------------*/
/*		‡A Initiation part 2 of S-function block's inner state					*/
/*		?iExecuted at the start and the reset of simulation?j					*/
/*		  How to reset S-function?: ($Revison 1.04)                             */
/*					Place S-function block in Enabled Subsystem, and 			*/
/*					put some signal to reset the block.                         */
/*------------------------------------------------------------------------------*/
void My_Initialize(SimStruct *S)
{
	Initialize_Pop(Population);
}
/*------------------------------------------------------------------------------*/
/*		‡B Description of the output of S-function block						*/
/*		The signals in input port can be referred as u[0] u[1] ...				*/
/*		and output as y[0] y[1] ....                                            */
/*------------------------------------------------------------------------------*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
	const real_T *u = (const real_T*) ssGetInputPortSignal(S,0);
	real_T       *y = ssGetOutputPortSignal(S,0);
    
	//Inner state variables used for local calculations
	static Individual Population_new[POPULATION_SIZE];
    static float Avg_fil=5, Disp_fil=5, Fit_abs_best=10;
    static Individual AlainJupe;
    Individual Father, Mother, Son, Daughter, Xman;
    float Vs, Is, Slip, Cosphi,TmpAv,mutend=1,Req, Xeq,R22S2,Den,Iest,Fit_new,IestBest,XeqdReq, XrdRr_mes, Fit_new2, XrdRr_best;
    float Fitnewbest, Fitnew2best;
    int i, popnew_ctr, j, brake,check;
    Complex Iest_complex, Vest_complex;
    static int discard=0;
    //Theoretical valus
    float Xs, Xr, Xm, ws, R22S2_c, Den_c, Xeq_c, Req_c, Iest_c;
    float XeqdReq_c;
    Complex Iest_complex_c;
    //wr estimation
    float Id, Iq, W,wr_c, wrest_best, wr_mes,Msrest, Trest,Phird_lim,Phird_lim_c,wrest;
    static float Phird_est=0;
    static float Phird_est_c=0;
    float Fit_new3=0, Fitnew3best;
    
	//Measurements treatment
	Vs= (float)u[0];	
	Is= (float)u[1];
    Slip= (float)u[2];
    
    //Theoretical calculations
    ws= (float)u[3];
    XrdRr_mes= (float)u[4];
    //wr estimation
    Id= (float)u[5];
    Iq= (float)u[6];
    wr_mes = (float)u[7];
            
    if (1)//if((Generation < 1000) ||( ((Avg_fil>AVG_MIN_END) || (Disp_fil>DISP_MIN_END)) && (Generation < MAX_GEN) && (Term==0) ))
    {
        Average = 0; 
        Dispersion = 0;
        for(i=0; i<POPULATION_SIZE; i++)
        {
            //Permanent regime IM model based calculations
            R22S2 = (Population[i].Chrom[2]*Population[i].Chrom[2]/(Slip*Slip));
            Den = R22S2 + ((Population[i].Chrom[3] + Population[i].Chrom[4])*(Population[i].Chrom[3] + Population[i].Chrom[4]));
            Req = Population[i].Chrom[0] + ( (Population[i].Chrom[4]*Population[i].Chrom[4]*Population[i].Chrom[2]/Slip) / Den);
            Xeq = Population[i].Chrom[1] + ( (R22S2*Population[i].Chrom[4] + Population[i].Chrom[3]*Population[i].Chrom[4]*(Population[i].Chrom[3] + Population[i].Chrom[4]))
                                            / Den);
            Vest_complex.Real = Vs;
            Vest_complex.Imag = 0;
            Iest_complex.Real = Req;
            Iest_complex.Imag = Xeq;
            Iest_complex = ComplexDiv(Vest_complex, Iest_complex);
            Iest = sqrt(Iest_complex.Real*Iest_complex.Real + Iest_complex.Imag*Iest_complex.Imag);
            //Iest = Iest_complex.Real;
            
            //Cos phi based calculations
            XeqdReq = Xeq / Req;
            
            //wr etimation based alculations
            ws = Avoid_zero(ws, 0.01);
            Msrest = Population[i].Chrom[4] / ws;
            Trest = (Msrest + Population[i].Chrom[3]/ws) / Population[i].Chrom[2];
            Phird_est=((Msrest*SAMPLE_TIME*Id)+(Trest-SAMPLE_TIME)*Phird_est)/Trest;
            Phird_lim= Avoid_zero(Phird_est, 0.001);
            wrest = (Iq*Msrest/Trest)/Phird_lim;
            
            Fit_new = ((Iest/Is) - 1)*((Iest/Is) - 1);
            //Fit_new2 is not using cos phis as in paper but tan phi = Xeq/Req 
            Fit_new2 = ((XeqdReq/XrdRr_mes) - 1)*((XeqdReq/XrdRr_mes) - 1);
            Fit_new3 = ((wrest/wr_mes) - 1)*((wrest/wr_mes) - 1);
            
            if(i==0)
            {
                IestBest = Iest;
                XrdRr_best = XeqdReq;
                Fitnewbest = Fit_new;
                Fitnew2best = Fit_new2;
                wrest_best = wrest;
                Fitnew3best = Fit_new3;
            }
            
            //if any chromozome implosible discard by setting big fitness
            check = Check_Chromo(Population[i].Chrom);
            if(check == 1)
            {
                Population[i].Fit = 20;
                discard++;
            }
            else
            {
                Population[i].Fit = Sigma(Fit, Fit_new, FIT_ARR_SIZE)
                                    + 0.1*Sigma(Fit2, Fit_new2, FIT_ARR_SIZE)
                                    + 0*Sigma(Fit3, Fit_new3, FIT_ARR_SIZE);
            }
            //Population[i].Fit = Population[i].Chrom[0]*dW + Population[i].Chrom[1] - Telec;
            //Population[i].Fit = Abs(Population[i].Fit);
            Average += Population[i].Fit;
        }
        Average = Average / (POPULATION_SIZE);
        Dispersion = StandardDev(Population, Average);
        
        Generation ++;
        SortPopulation_Best_First(Population);        

        //Calculate the probability of selection of each individual based on its fitness
        SelectionProbaCalc(Population);
        popnew_ctr=0;
        brake=0;
        j=0;
        //Create new population     
        while((popnew_ctr<POPULATION_SIZE))
        {
            for(i=0; i<POPULATION_SIZE; i++)
            {
                if (Population[i].Fit < Fit_abs_best)
                {
                    Fit_abs_best = Population[i].Fit;
                    AlainJupe = Population[i];
                }
                if(Proba_Sel(Population[i].Psel)==1)
                {
                    Population_new[popnew_ctr]=Population[i];
                    popnew_ctr++;
                    j++;
                }
                else    
                {   //Mutate
                    //mutate only weaker individuals
                    mutend = Maxx(((PMUT_END - Generation)/PMUT_END),0.0000001);
                    if((Proba_Sel(P_MUT*mutend)==1)&&(popnew_ctr<POPULATION_SIZE))
                    {
                        ChromoMutate(Population_new[i].Chrom, Xman.Chrom);
                        Population_new[popnew_ctr] = Xman;
                        popnew_ctr++;
                        mut++;
                        j=0;    //Don't cross mutations 
                    }
                }
                if(j>=2)
                {
                    j=0;
                    //Coss two parents with a probability
                    if((Proba_Sel(P_CROSS*(2-mutend))==1)&&(popnew_ctr<POPULATION_SIZE))
                    {
                        ChromoCross(Population_new[popnew_ctr-2].Chrom, Population_new[popnew_ctr-1].Chrom, Son.Chrom, Daughter.Chrom);
                        Population_new[popnew_ctr] = Son;
                        popnew_ctr++;
                        Population_new[popnew_ctr] = Daughter;
                        popnew_ctr++;
                        cros++;
                    }
                }  
            }
            brake++;
        } 
        //God save the king
        for(i=1; i<POPULATION_SIZE; i++)
        {
            Population[i] = Population_new[i-1];
        }
    }
    //else 
    {
        Term = 1;
    }
    if(Average>1000) 
        TmpAv = 1000;
    else
        TmpAv = Average;
    
    if(Generation >= 1000)
    {
        Avg_fil = LowPassDo(Average, (float)0.2, 0);
        Disp_fil = LowPassDo2(Dispersion, (float)0.2);
    }
    
    
    //Theoretical calculations
    Xm = ws*Msr;    
    Xr = ws*(Lr-Msr);
    Xs = ws*(Ls-Msr);
    R22S2_c = (Rr*Rr/(Slip*Slip));
    Den_c = R22S2_c + ((Xr + Xm)*(Xr + Xm));
    Req_c = Rs + ( (Xm*Xm*Rr/Slip)/ Den_c );
    Xeq_c = Xs + ( (R22S2_c*Xm + Xr*Xm*(Xr + Xm))/ Den_c);
    Iest_complex_c.Real = Req_c;
    Iest_complex_c.Imag = Xeq_c;
    Iest_complex_c = ComplexDiv(Vest_complex, Iest_complex_c);
    Iest_c = sqrt(Iest_complex_c.Real*Iest_complex_c.Real + Iest_complex_c.Imag*Iest_complex_c.Imag);
    //Iest_c = Iest_complex_c.Real;
    XeqdReq_c = Xeq_c / Req_c;
    
    //wr theoretical
    Phird_est_c=((Msr*SAMPLE_TIME*Id)+((Lr/Rr)-SAMPLE_TIME)*Phird_est_c)/(Lr/Rr);
    Phird_lim_c= Avoid_zero(Phird_est_c, 0.001);
    wr_c = (Iq*Msr*Rr/Lr)/Phird_lim_c;
    
	/* Outouts */
	y[0]=(float)Population[0].Chrom[0]; 								
	y[1]=(float)LowPassDo(y[0],(float)0.2, 1);      //Rs;
    y[2]=(float)Population[0].Chrom[2]; 								
	y[3]=(float)LowPassDo(y[2],(float)0.2, 2);      //Rr;
    y[8]=(float)Population[0].Chrom[4]/ws;
	y[9]=(float)LowPassDo(y[8],(float)0.2, 5);      //Xm;
    y[4]=(float)y[8] + Population[0].Chrom[1]/ws;	
	y[5]=(float)LowPassDo(y[4],(float)0.2, 3);      //Xs;
    y[6]=(float)y[8] + Population[0].Chrom[3]/ws;
    y[7]=(float)LowPassDo(y[6],(float)0.2, 4);      //Xr;
    y[10]=(float)cros;                              //Population[5].Psel;//Chrom[3];
	y[11]=(float)mut;                               //Chrom[4];
	y[12]=(float)Term*10000;
    y[13]=(float)TmpAv;
    y[14]=(float)IestBest;
    y[15]=(float)Iest_c;
    y[16]=(float)Dispersion;
    y[17]=(float)XeqdReq_c;
    y[18]=(float)XrdRr_best;
    y[19]=(float)Avg_fil;
    y[20]=(float)Disp_fil;
    y[21]=(float)wrest_best;
    y[22]=(float)wr_c;
    y[23]=(float)Fitnew3best;
    y[24]=(float)AlainJupe.Chrom[0];
    y[25]=(float)AlainJupe.Chrom[2];
    y[28]=(float)AlainJupe.Chrom[4]/ws;
    y[26]=(float)y[28]+AlainJupe.Chrom[1]/ws;
    y[27]=(float)y[28]+AlainJupe.Chrom[3]/ws;
    y[29]=(float)AlainJupe.Fit;
    
	/****** End ******/
	
}
/*------------------------------------------------------------------------------*/
/***		Part 3 That was all for your coding for S-function!				  ***/
/***			 ?iYou don't have to rewrite following codes?j						  ***/
/*------------------------------------------------------------------------------*/
/*		Appendix : What you are doing in following scripts	($Revision 1.03)				  		*/
/*		‡@ Model Update at the process of Integration in S-Function Block
                                        : mdlUpdate(SimStruct *S, int_T tid)	*/
/*		‡A Calculation of Derivatives in S-Function Block	: mdlDerivatives(SimStruct *S)			*/
/*		‡B Termination of the Simulation	: mdlTerminate(SimStruct *S)			*/
/*------------------------------------------------------------------------------*/



/* Error handling
 * --------------
 *
 * You should use the following technique to report errors encountered within
 * an S-function:
 *
 *       ssSetErrorStatus(S,"Error encountered due to ...");
 *       return;
 *
 * Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
 * It cannot be a local variable. For example the following will cause
 * unpredictable errors:
 *
 *      mdlOutputs()
 *      {
 *         char msg[256];         {ILLEGAL: to fix use "static char msg[256];"}
 *         sprintf(msg,"Error due to %s", string);
 *         ssSetErrorStatus(S,msg);
 *         return;
 *      }
 *
 * See matlabroot/simulink/src/sfuntmpl_doc.c for more details.
 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    /* See sfuntmpl_doc.c for more details on the macros below */

    ssSetNumSFcnParams(S, NUM_PARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 1)) return;
	ssSetInputPortWidth(S, 0, INPUT_PORT_WIDTH);
    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    /*
     * Set direct feedthrough flag (1=yes, 0=no).
     * A port has direct feedthrough if the input is used in either
     * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
     * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
     */
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, OUTPUT_PORT_WIDTH);

	/* ssSetNumSampleTimes ‚ð2ˆÈ?ã‚É‚·‚é‚±‚Æ‚Æ?C?”?s‰º‚ÌssSetSampletime‚ð	*/
	/* •¡?”Žw’è‚·‚é‚Å?CS-Functionƒuƒ?ƒbƒN‚Ì‰‰ŽZŽüŠú‚ð‰Â•Ï‚É‚Å‚«‚é?D			*/
	/* Commented : Hisafumi Asai $Revision 1.01 							*/
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetOptions(S, 0);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}



#define MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */
#if defined(MDL_INITIALIZE_CONDITIONS)
  /* Function: mdlInitializeConditions ========================================
   * Abstract:
   *    In this function, you should initialize the continuous and discrete
   *    states for your S-function block.  The initial states are placed
   *    in the state vector, ssGetContStates(S) or ssGetRealDiscStates(S).
   *    You can also perform any other initialization activities that your
   *    S-function may require. Note, this routine will be called at the
   *    start of simulation and if it is present in an enabled subsystem
   *    configured to reset states, it will be call when the enabled subsystem
   *    restarts execution to reset the states.
   */
  static void mdlInitializeConditions(SimStruct *S)
  {
		My_Initialize(S);
  }
#endif /* MDL_INITIALIZE_CONDITIONS */



#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {
		/* Added : Hisafumi Asai $Revision 1.03 */
		My_mdlStart(S);
  }
#endif /*  MDL_START */



/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector, ssGetY(S).
 */
#if 0
/* ?ã‚ÉˆÚ“®‚µ‚½	*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T *u = (const real_T*) ssGetInputPortSignal(S,0);
    real_T       *y = ssGetOutputPortSignal(S,0);
}
#endif


#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
  /* Function: mdlUpdate ======================================================
   * Abstract:
   *    This function is called once for every major integration time step.
   *    Discrete states are typically updated here, but this function is useful
   *    for performing any tasks that should only take place once per
   *    integration step.
   */
  static void mdlUpdate(SimStruct *S, int_T tid)
  {
  }
#endif /* MDL_UPDATE */



#define MDL_DERIVATIVES  /* Change to #undef to remove function */
#if defined(MDL_DERIVATIVES)
  /* Function: mdlDerivatives =================================================
   * Abstract:
   *    In this function, you compute the S-function block's derivatives.
   *    The derivatives are placed in the derivative vector, ssGetdX(S).
   */
  static void mdlDerivatives(SimStruct *S)
  {
  }
#endif /* MDL_DERIVATIVES */



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
}


/*======================================================*
 * See sfuntmpl_doc.c for the optional S-function methods *
 *======================================================*/

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
