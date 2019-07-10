#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinterface.h"
#include "outputdist.h"



void gensixcanonical(){
                    printf("hjeegensixcanonicalgensixcanonicaleeeree  %d \n", dist->incoordtype);
    if(dist->incoordtype==0){ // action angle
        generatefromaction();
    }
    else if(dist->incoordtype==1){ //normalized coordinates

        generatefromnormalized();
    }
    else if(dist->incoordtype==2){ //physical coord
        generatefromphysical();
    }
    else if(dist->incoordtype==3){ // mixed coord
        generatefrommixed();
    }
}

void generatefromnormalized(){
    int counter = 0;
    double tc[dim];
    double normalized[dim], cancoord[dim];

    for(int i =0; i< dist->totincoord; i++){
        normalized[0]=dist->incoord[i]->normalized[0];
        normalized[1]=dist->incoord[i]->normalized[1];
        normalized[2]=dist->incoord[i]->normalized[2];
        normalized[3]=dist->incoord[i]->normalized[3];
        normalized[4]=dist->incoord[i]->normalized[4];
        normalized[5]=dist->incoord[i]->normalized[5];
        normalized2canonical(normalized, cancoord);

        if(particle_within_limits_physical(cancoord)==1){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->physical[p] = cancoord[p];
            }

            // Not nescessary at the moment but might be in the future.
            dist->outcoord[i]->mass  = dist->incoord[i]->mass;
            dist->outcoord[i]->a     = dist->incoord[i]->a;
            dist->outcoord[i]->z     = dist->incoord[i]->z;
            counter++;
        }
    }
    dist->totoutcoord=counter;
    dist->isDistrcalculated=1;
}

void generatefromaction(){
    int counter = 0;
    double tc[dim];
    double normalized[dim], cancoord[dim];

    for(int i =0; i< dist->totincoord; i++){

        tc[0]=dist->incoord[i]->action[0];
        tc[1]=dist->incoord[i]->action[1];
        tc[2]=dist->incoord[i]->action[2];
        tc[3]=dist->incoord[i]->action[3];
        tc[4]=dist->incoord[i]->action[4];
        tc[5]=dist->incoord[i]->action[5];
        action2normalized(tc, normalized);
        normalized2canonical(normalized, cancoord);
        if(particle_within_limits_physical(tc)==1){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->physical[p]   = cancoord[p];
                dist->outcoord[counter]->normalized[p] = normalized[p];
            }
            // Not nescessary at the moment but might be in the future.
            dist->outcoord[i]->mass  = dist->incoord[i]->mass;
            dist->outcoord[i]->a     = dist->incoord[i]->a;
            dist->outcoord[i]->z     = dist->incoord[i]->z;
            counter++;
        }
    }
    dist->totoutcoord=counter;
    dist->isDistrcalculated=1;
}

void generatefrommixed(){

}


void generatefromphysical(){
    int counter = 0;
    double tc[dim];
    double tmp_n[dim];

    for(int i =0; i< dist->totincoord; i++){
        tc[0]=dist->incoord[i]->physical[0];
        tc[1]=dist->incoord[i]->physical[1];
        tc[2]=dist->incoord[i]->physical[2];
        tc[3]=dist->incoord[i]->physical[3];
        tc[4]=dist->incoord[i]->physical[4];
        tc[5]=dist->incoord[i]->physical[5];
        
        if(particle_within_limits_physical(tc)==1){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->physical[p] = tc[p];
            }
            // Not nescessary at the moment but might be in the future.
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
    normalized[5]=-sqrt(acangl[4]/2/1000)*sin(acangl[5]);
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