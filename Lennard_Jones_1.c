#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double L; //Length of the box
int    N; // Number of particles
double rho; // Number Density
double dt; // Timestep
double m; // Mass of each particke
double K0; // Initial Kinetic Energy
double K; // Kinetic Energy
double U; //Potential Energy
double H; // Total Energy
double T; //Temperature
double rc = 3.0; //Potential cutoff
double rch = 0.3; // Hard Sphere Cut-off
double epsilon = 1;
int    numSteps; 
double latticex; // X-Position of a particular particle in Lattice
double latticey; // Y-Position of a particular particle in Lattice
double latticez; // Z-Position of a particular particle in Lattice
int static count = 0;

struct particle{
	double rx, ry, rz;
	double vx, vy, vz;
	double ax, ay, az;
};

double NormalRandomGen(void){
	double U1, U2, norVar;
	
	U1 = rand();
	U2 = rand();
	U1 /= RAND_MAX;
	U2 /= RAND_MAX;
	norVar = sqrt(-2*log(U1))*cos(2*M_PI*U2); // Box-Muller Transform
	// https://www.baeldung.com/cs/uniform-to-normal-distribution
	
	return norVar;
}

void systemInit(struct particle *atoms){
	int i=0, j=0, k=0, numpart;
	double lsp = L/((int)(cbrt(N))-1); // Lattice Spacing
	double U1, U2, norVar;
	double K0 = 0;
	
		
	for(numpart=0; numpart<N; numpart++){
		/* Position Initialization */
		latticex = i*lsp - 0.5*L;
		latticey = j*lsp - 0.5*L;
		latticez = k*lsp - 0.5*L;
		
		//printf("%d %.3f %.3f %.3f \n", numpart, latticex, latticey, latticez);
		
		atoms[numpart].rx = latticex;
		atoms[numpart].ry = latticey;
		atoms[numpart].rz = latticez;
		
		i += 1;
		if(i*lsp > L){
			j += 1;
			
			i = 0;
			
			if(j*lsp > L){
				k+= 1;
				
				i = 0;
				j = 0;
			}
		}
	
	/* Velocity Initialization */
	atoms[numpart].vx = sqrt(T)*NormalRandomGen()/m;
	atoms[numpart].vy = sqrt(T)*NormalRandomGen()/m;
	atoms[numpart].vz = sqrt(T)*NormalRandomGen()/m;
	
	K0 += 0.5*m*((atoms[numpart].vx)*(atoms[numpart].vx)+
		 atoms[numpart].vy*atoms[numpart].vy+
		 atoms[numpart].vz*atoms[numpart].vz);
	}	
	// printf("Initial Energy is %f \n", K0);
	/* for(i=0; i<N; i++){
		printf("%d %.2f %.2f %.2f %.2f %.2f %.2f \n", i+1, atoms[i].rx,atoms[i].ry,atoms[i].rz,atoms[i].vx,atoms[i].vy,atoms[i].vz);
	} */
	
	for(numpart=0; numpart<N; numpart++){
		atoms[numpart].ax = 0;
		atoms[numpart].ay = 0;
		atoms[numpart].az = 0;
	}
}


double Periodic(double r){
	//printf("%f \n", r);
	while(r > 0.5*L){
		r -= L;
	}
	while(r<-0.5*L){
		r += L;
	} 
	
	//printf("%f \n", r);
	return r;
}


double calcuForces(struct particle *atoms){
	double r, r2, r2i, r6i;
	double eij, fij, temp;
	int i, j, k = 0;
	double dx, dy, dz;
	
	U = 0;
	
	for(i=0; i<N-1; i++){
		for(j=i+1; j<N; j++){
			dx = Periodic(atoms[i].rx-atoms[j].rx);
			dy = Periodic(atoms[i].ry-atoms[j].ry);
			dz = Periodic(atoms[i].rz-atoms[j].rz);
			
			//printf("%d %.3f %.3f %.3f \n", k+1, dx, dy, dz);
			count += 1;
			
			r2 = dx*dx + dy*dy + dz*dz;
			//printf("%.3f ", r2);
			r = sqrt(r2);
			
			if(r> L/1.732){
				if((atoms[i].rx-atoms[j].rx)>L/2 || (atoms[j].rx-atoms[i].rx)>L/2){
					temp = atoms[j].rx - L;
					dx = atoms[i].rx- temp;
				}
				
				if((atoms[i].ry-atoms[j].ry)>L/2 || (atoms[j].ry-atoms[i].ry)>L/2){
					temp = atoms[j].ry - L;
					dx = atoms[i].ry - temp;
				}
				
				if((atoms[i].rz-atoms[j].rz)>L/2 || (atoms[j].rz-atoms[i].rz)>L/2){
					temp = atoms[j].rz - L;
					dx = atoms[i].rz- temp;
				}
				
				r2 = dx*dx + dy*dy + dz*dz;
				r = sqrt(r2);
			}
			
			
			
			if(r2<rc){
				if(r2>rch){
					r2i = 1/(r*r);
					r6i = r2i*r2i*r2i;
					eij = 4*epsilon*r6i*(r6i-1);
					fij = 48*epsilon*r2i*r6i*(r6i-0.5);
					//printf("%d %.3f %.3f %.3f \n",count, r2, eij, fij);
				
					U += eij;
				
					//Position Change
					atoms[i].ax = atoms[i].ax + fij*dx/m;
	            	atoms[i].ay = atoms[i].ay + fij*dy/m;
	            	atoms[i].az = atoms[i].az + fij*dz/m;
	            	atoms[j].ax = atoms[j].ax - fij*dx/m;
	            	atoms[j].ay = atoms[j].ay - fij*dy/m;
	            	atoms[j].az = atoms[j].az - fij*dz/m; 
	        	}
			}
			else{
				//printf("\n");
			}
		}
	}
	
	return U;
}


double motionInteg(struct particle *atoms){
	int i;
	double alpha;
	
	K = 0;
	H = 0;
	
	//velocity change: 1st half step
	for(i=0; i<N; i++){
		
		// U = calcuForces(atoms);
		
		/* for(i=0; i<N; i++){
		printf("%d %.2f %.2f %.2f %.2f %.2f %.2f \n", i+1, atoms[i].rx,atoms[i].ry,atoms[i].rz,atoms[i].vx,atoms[i].vy,atoms[i].vz);
	} */
	
		atoms[i].vx += 0.5*atoms[i].ax*dt;
		atoms[i].vy += 0.5*atoms[i].ay*dt;
		atoms[i].vz += 0.5*atoms[i].az*dt;
		
	// printf("%d %.2f %.2f %.2f %.2f %.2f %.2f \n", i+1, atoms[i].rx,atoms[i].ry,atoms[i].rz,atoms[i].vx,atoms[i].vy,atoms[i].vz);
		
		
	K += 0.5*m*(atoms[i].vx*atoms[i].vx 
	   + atoms[i].vy*atoms[i].vy
	   + atoms[i].vz*atoms[i].vz);
	   
	   //printf("%.3f \n", K);
	}
	
	//printf("%f \n", K);
	
	//Temperature Correction by Feedback Method
	alpha = sqrt(K0/K);
	/* for(i=0; i<N; i++){
		atoms[i].vx *= alpha;
		atoms[i].vy *= alpha;
		atoms[i].vz *= alpha;
	} */
	
	// Position Change
	for(i=0; i<N; i++){
		atoms[i].rx += atoms[i].vx*dt;
		atoms[i].ry += atoms[i].vy*dt;
		atoms[i].rz += atoms[i].vz*dt;
	}
	
	U = calcuForces(atoms);
	
	for(i=0; i<N; i++){
		atoms[i].vx += 0.5*atoms[i].ax*dt;
		atoms[i].vy += 0.5*atoms[i].ay*dt;
		atoms[i].vz += 0.5*atoms[i].az*dt;
	}
	
	H = U + K;
	
	printf("%f %f %f \n", U, K, H);
	
	return H;
}


void simulation(void){
	struct particle atoms[N];
	double sumH =0 , sumH2 = 0, fluctH = 0, Hbar;
	int i;
	
	systemInit(atoms);
	
	for(i=0; i<numSteps; i++){
		
		H = motionInteg(atoms); 
		
		sumH  = sumH + H;
		Hbar = sumH/(i+1);
		sumH2 = sumH2 + H*H;
      	fluctH = sqrt((sumH2/(i+1)) - Hbar*Hbar)/Hbar; 
      	
      	//printf("%d %f %f %f %f \n",i, H, sumH, sumH2, fluctH);
	}
	printf("The Total Energy is %f \n",H);
}


void main(void){
	N = 1000;
	rho = 1;
	dt = 0.0001;
	T = 300;
	m = 1;
	numSteps = 10;
	
	L = cbrt(N/rho);
	//printf("L = %f \n", L);
	
	simulation();
	
}
