#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

using namespace std;
const int N = 16;         // Number of particles
double L = 4;            // Declare & assign the linear size of the system
double runMax = 8000;      // Declare & assign the maximum number of runs
double sigma = 1.0;       // sigma (in reduced units)
//double a = sigma;       // a = lattice constant
double delta_r = sigma / 50;     // => moves particles randomly within +/- sigma/2
double r[N][2];   //r[N particles][x and y-coordinates]: Declare array of positions;  3d: [N][3]
double v[N][2];   //v[N particles][x and y-coordinates]: Declare array of velocities; 3d: [N][3]
double acc[N][2];         // Declare arrays of accelerations   3d: [N][3]
double vcm[2];
double Vmax = 1.0;        // Declare & assign the maximum initial velocity
double dt = 0.005;
double force =0;
double Potential=0;
double Energy=0;
double Temperature=0;
double vSquared = 0;
double Center_Mass_Velocity(); 

void Initialize()
{


    // Initialize positions
    int n = int(ceil(pow(N, 1.0/2)));  // n=square root(N)=number of atoms in each row/column
    double a = L / int(ceil(pow(N, 1.0/2)));        // a = lattice spacing
                                       // See photo of calculations in Telegram
                                       // CompPhys953_Gen or _Edu    [3d: 1.0/3]

   ofstream write("positions2.dat"); // Use ofstream (output file stream) to create a file
                                    // named "positions.dat" to write data to.
                                    // "cout" writes to screen. "write" writes to file.

    int p = 0;                 // Counter: number of particles Initialized
    for (int x = 0; x < n; x++)
      for (int y = 0; y < n; y++)   // 3d: for (int z = 0; z < n; z++)
        {
           if (p < N) // Continue placing particles
           {
 //          r[particle p][x-coordinate]
             r[p][0] = (x+0.5)*a + 2*(rand()/double(RAND_MAX)-0.5)*delta_r; // Displace x of
                                             // particle p randomly within +/- sigma/2
			 write << r[p][0] << '\t';
 //          r[particle p][y-coordinate]
             r[p][1] = (y+0.5)*a + 2*(rand()/double(RAND_MAX)-0.5)*delta_r; // Displace y of
                                             // particle p randomly within +/- sigma/2
			 write << r[p][1] << '\t';
			 write << '\n';
           }
           ++p;       // Increment particle counter by 1
        }
    write.close();

    // Initialize velocities
    for (int i = 0; i < N; i++)        // Initialize all velocities
        for (int j = 0; j < 2; j++)    // Vx and Vy        3d: [< 3]
            v[i][j] = Vmax * (2*rand()/double(RAND_MAX)-1);  // particle velocity is given values
   Center_Mass_Velocity();                                                            // between -Vmax and +Vmax
}
template <typename T>
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

void Accelerate()
{
    Potential = 0;
    force = 0;
   for (int i = 0; i < N; i++)         // set accelerations to zero
      for (int j = 0; j < 2; j++)      // 3d: < 3
         acc[i][j] = 0;

   for (int i = 0; i < N-1; i++)
      for (int j = i+1; j < N; j++)   // Sum over all j>i pairs
      {
          double rij[2];              // 3d: [3]    // position of i relative to j
          double rSquared = 0;
         for (int k = 0; k < 2; k++) {		// 3d: [3]
				rij[k] = r[i][k] - r[j][k];
				if (fabs(rij[k]) > L / 2)
					rij[k] -= sgn(rij[k]) * L;
    			rSquared += rij[k] * rij[k];
			}
           //Differentiate LJ -V(x) with respect to x to get the force
            Potential += 4 * (pow(rSquared, -6) - pow(rSquared, -3));
            force = 24 * (2 * pow(rSquared, -7) - pow(rSquared, -4));
         for (int k = 0; k < 2; k++) // 3d: [3]
	 {
           acc[i][k] += rij[k]*force;
           acc[j][k] -= rij[k]*force;
         }
       }
}
double Center_Mass_Velocity() {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 2; j++)
			vcm[j] += v[i][j] / N;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < 2; j++)
			v[i][j] -= vcm[j];
}

void pbc()
{
     for(int i=0; i<N; i++)
     {
         if (r[i][0] < 0) r[i][0] += L;
         if (r[i][0] > L) r[i][0] -= L;
         if (r[i][1] < 0) r[i][1] += L;
         if (r[i][1] > L) r[i][1] -= L;
     }
}
void velocityVerlet(double dt) {
    Accelerate();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 2; j++) {    //[2d [2]]
            r[i][j] += v[i][j] * dt + 0.5 * acc[i][j] * dt * dt;
            v[i][j] +=0.5 * acc[i][j] * dt;
        }
    Accelerate();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 2; j++)      //[2d [2]]
            v[i][j] += 0.5 * acc[i][j] * dt;
            pbc();
	
}



double instantaneousTemperature() {
	vSquared = 0;
	
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 2; j++)                     //[2d [2]]
			vSquared += v[i][j] * v[i][j];
	double T = vSquared / (2 * (N - 1));
	return T;
}

////////////////////     MAIN    ////////////////////////
int main(int argc, char *argv[]) {

   srand (time(NULL));

   Initialize();
    ofstream position1("position1.txt");
    ofstream position2("position2.txt");
    ofstream position3("position3.txt");
    ofstream temperature("temperature.txt");
    ofstream energy("energy.txt");
	//std::ofstream xyz_file ("test.xyz");
    
    for (int i = 0; i < runMax; i++){
        Energy=0;
        Potential=0;
        velocityVerlet(dt);
        /////////////////////////////////////////////////////////////////////////////////
     
        
       /* xyz_file << N << '\n' << i << '\n';
        
        for (int p = 0; p < N; p++) {
            xyz_file << 1 << '\t';
            for (int q = 0; q < 2; q++) {
                
                xyz_file << r[p][q] << '\t';
            }
            
            xyz_file << 0 << '\t';
            for (int q = 0; q < 2; q++) {
                xyz_file << v[p][q] << '\t';
            }
            xyz_file << 0 << '\n';
            
        } */
        
        
        /*if ((i >= 0) && (i % 10 == 0) && (i <= 20)) {
			for (int p = 0; p < N; p++) {
				for (int q = 0; q < 2; q++) {
					position1 << r[p][q] << '\t';
				}
				position1 << '\n';
			}
		}*/


		if ((i >= 40) && (i % 10 == 0) && (i <= 800)) {
			for (int p = 0; p < N; p++) {
				for (int q = 0; q < 2; q++) {
					position2 << r[p][q] << '\t';
				}
				position2 << '\n';
			}
		}


		if ((i >= 6000) && (i % 10 == 0) && (i <= 6600)) {
			for (int p = 0; p < N; p++) {
				for (int q = 0; q < 2; q++) {
					position3 << r[p][q] << '\t';
				}
				position3 << '\n';
			}
		}
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        Temperature = instantaneousTemperature();
        //temperature << i*dt  << '\t' << Temperature << '\n';
        
        // Energy
        Energy=Potential + vSquared / 2;
        //energy << i*dt << '\t' << Energy << '\n';
        //////////////////////////////////////////////////////////////////////
        
        
        
    }
    return 0;
}
