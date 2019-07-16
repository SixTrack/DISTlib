#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinterface.h"
#include "outputdist.h"



void gensixcanonical(){

	int counter = 0;
    int type = dist->incoordtype;
    double tc[dim];
    double normalized[dim], cancoord[dim];

    for(int i =0; i< dist->totincoord; i++){
    	for(int k=0; k<6; k++){
    		tc[k]=dist->incoord[i]->coord[k];
    	}
    	

        if(type==0 || type==3){
        	action2normalized(tc, normalized);
        	normalized2canonical(normalized, cancoord);
        }
        else if(type==1){
	        for(int k=0; k<6; k++){
	    		normalized[k] = tc[k];
	        	normalized2canonical(tc, cancoord);
	    	}
	    }
	    else if(type==2){
			for(int k=0; k<6; k++){
    			cancoord[k]=tc[k];
    		}

	    }


        if(particle_within_limits_physical(cancoord)==1 && particle_within_limits_normalized(normalized)){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->coord[p]   = cancoord[p];
            }
            dist->outcoord[i]->mass  = dist->incoord[i]->mass;
            dist->outcoord[i]->a     = dist->incoord[i]->a;
            dist->outcoord[i]->z     = dist->incoord[i]->z;
            counter++;
        }
    }
    dist->totoutcoord=counter;
    dist->isDistrcalculated=1;

}

/*If emittance is defined it converts to canonical coordinates */
void action2normalized(double acangl[6], double normalized[6]){
    
    normalized[0]= sqrt(acangl[0]/2)*cos(acangl[1]);
    normalized[1]=-sqrt(acangl[0]/2)*sin(acangl[1]);
    normalized[2]= sqrt(acangl[2]/2)*cos(acangl[3]);
    normalized[3]=-sqrt(acangl[2]/2)*sin(acangl[3]);
    normalized[4]= sqrt(acangl[4]/2)*cos(acangl[5]);
    normalized[5]=-sqrt(acangl[4]/2)*sin(acangl[5]); // used to devide with 1000 here before.. 
}

void normalized2canonical(double normalized[6], double cancoord[6]){
    normalized[0] = sqrt(dist->emitt->e1)*normalized[0];
    normalized[1] = sqrt(dist->emitt->e1)*normalized[1];
    normalized[2] = sqrt(dist->emitt->e2)*normalized[2];
    normalized[3] = sqrt(dist->emitt->e2)*normalized[3];
    normalized[4] = sqrt(dist->emitt->e3)*normalized[4];
    normalized[5] = sqrt(dist->emitt->e3)*normalized[5];

   
    if(dist->incoordtype==3) {
        double lindp = 0;
        double lindeltas=0;
        double deltap = normalized[4];
        double deltas = normalized[5];
        double *xap;
        double det = (dist->tas[4][4]*dist->tas[5][5] - dist->tas[4][5]*dist->tas[5][4]);
        for(int i=0; i<4;i++){
            lindeltas = lindeltas+dist->tas[4][i]*normalized[i];
            lindp=lindp+dist->tas[5][i]*normalized[i];
        } 

        lindp = deltap - lindp;
        lindeltas = deltas - lindeltas;

        xap = (double*)malloc(2*sizeof(double));
        solve2by2eq(dist->tas[4][4], dist->tas[4][5], lindeltas, dist->tas[5][4], dist->tas[5][5], lindp, xap );
        normalized[4] = xap[0];
        normalized[5] = xap[1];

    }
    mtrx_vector_mult_pointer(dim,dim, dist->tas, normalized,cancoord);

}


/*Checks if the particle is within the physical limit set by the user*/
int particle_within_limits_physical(double *physical){
    
    if(dist->cuts2apply->isset_p==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->physical[i]->isset==1){
            if(physical[i] > dist->cuts2apply->physical[i]->min && physical[i] < dist->cuts2apply->physical[i]->max) return 0;
        }   
    }
    
    return 1;

}

/*Checks if the particle is within the normalized limit set by the user*/
int particle_within_limits_normalized(double *normalized){
    
    if(dist->cuts2apply->isset_n==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->normalized[i]->isset==1){
            if(normalized[i] > dist->cuts2apply->normalized[i]->min && normalized[i] < dist->cuts2apply->normalized[i]->max) return 0;
        }
    }
    return 1;
}
/*
void createcoordinates(int index,  double start, double stop, int length, int type){

    dist->coord[index-1]->start = start;
    dist->coord[index-1]->stop = stop;
    dist->coord[index-1]->length = length;
    dist->coord[index-1]->type = type;
//  dist->coord[index-1]->coordtype = coordtype;
    

    if(type ==0){ //Constant value 
        dist->coord[index-1]->values = (double)malloc((length)sizeof(double));
        dist->coord[index-1]->values[0] = start;
        dist->coord[index-1]->length = 1; //if it is a constant the length should always be 1.
    }

    if(type > 0){ //Allocate space for the array
        if(dist->disttype==1){
            if(dist->totallength>0){
                int tmp;
                tmp=dist->totallength;
                length=&tmp;
                }
            else
                printf("You need to set a totallength for disttype 1!");
        }
        dist->coord[index-1]->values = (double)malloc((length)*sizeof(double));
     //   memcpy(dist->coord[index-1]->values , start, sizeof(double));   //not sure 
    }
    if(type==1){ //Linearly spaced intervalls
    
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
    }
    if(type==2){ //Exponentially spaced
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
        for(int i=0;i <length; i++){
               
            dist->coord[index-1]->values[i] = exp(dist->coord[index-1]->values[i]);
        }
    }
    if(type==3){ //Spaced with  ^2
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
        for(int i=0;i <length; i++){
            dist->coord[index-1]->values[i] = pow(dist->coord[index-1]->values[i],2);
        
        }

    }
    if(type==4){ // uniform random 
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
        for(int i=0;i <length; i++){
            dist->coord[index-1]->values[i] = rand_uni(start, stop);
            //printf("%f \n", dist->coord[index-1]->values[i] );
        }
    }

    if(type==5){ // Gaussian random (Here start is mean and stop is the standard deviation)
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
        for(int i=0;i <length; i++){
            dist->coord[index-1]->values[i] = randn(start, stop);
        }
    }

    if(type==6){ // Rayleigh distribution
        createLinearSpaced(length, start, stop,dist->coord[index-1]->values);
        for(int i=0;i <length; i++){
            dist->coord[index-1]->values[i] = randray(start, stop);

        }
    }
}*/