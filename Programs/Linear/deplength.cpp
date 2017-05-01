#include "pol_montecarlo.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

// Remember all the length and time units are in mfp=lsa.

const double lambda=0.5145;	// micro
const double n_water=1.336;
const double dia=4.66*2;	// diameter in micro
const double x=pi*dia/lambda*n_water;
const double mr=1.19, mi=0.0;	// m=mr-i*mi is the refractive index
const double L=200;		// the thickness of the slab
const double lsa=27.63;		// micro, to produce lt=314 micro in Fig 4 in Lenke (2002) in J. Opt. A: Pure Appl. Opt 4, 293-298
const double verdet=0.01;       // angle of Faraday rotation per one lsa, = delta_n/2*(2*pi*n_water/lambda)*lsa.
const int photons=100*4*10;
const double theta_step = 0.002; // degree

// find the new normal direction after scattering
// Input: Incident direction s0, incident normal n0, outgoing direction s1
// Output: Outgoing normal n1
void normalDirection(const double s0[3], const double n0[3], const double s1[3], double n1[3])
{
	double n_len;

	// n1 = s0 x s1
	n1[0] = s0[1]*s1[2]-s0[2]*s1[1];
	n1[1] = -(s0[0]*s1[2]-s0[2]*s1[0]);
	n1[2] = s0[0]*s1[1]-s0[1]*s1[0];

	n_len = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);

	if (n_len < 1e-8) {	// NEED TO CHECK THE OPTIMAL VALUE
		if (s0[0]*s1[0]+s0[1]*s1[1]+s0[2]*s1[2] > 0) {
			n1[0] = n0[0];
			n1[1] = n0[1];
			n1[2] = n0[2];
		}
		else {
			n1[0] = -n0[0];
			n1[1] = -n0[1];
			n1[2] = -n0[2];
		}
	}
	else {
		n1[0] /= n_len;
		n1[1] /= n_len;
		n1[2] /= n_len;
	}
}


// evaluate the rotation matrix R(phi) when the coordinate system rotates from 
// (m0, n0, s0) to (m1, n1, s1)
void rotationMatrix(const double m0[3], const double n0[3], const double n1[3], double* cosphi, double* sinphi)
{
	*cosphi = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
	*sinphi = -(m0[0]*n1[0] + m0[1]*n1[1] + m0[2]*n1[2]);
}


// Faraday rotation of (cosphi, sinphi) over a distance of d measured in lsa.
// (m0, n0, s0) to (m1, n1, s1)
void FaradayRotation(const double d, double* cosphi, double* sinphi)
{
	double delta = verdet*d;
	*cosphi = (*cosphi)*cos(delta) - (*sinphi)*sin(delta);
	*sinphi = (*cosphi)*sin(delta) + (*sinphi)*cos(delta);
}


// c = a x b
void axb(const double a[3], const double b[3], double c[3])
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}


int main ()
{

	//Print to file
	ofstream fp_out;
	fp_out.open("Linearoutput.txt",ios::out); 

	dcmplx m(mr, mi);
	scatterer sct(x, m, 50001, 50000);
	photonPacket *ph = new photonPacket(&sct);
	double csca = sct.csca();

	// incident E field
	const dcmplx Ein1=dcmplx(1, 0), Ein2=dcmplx(0, 0); //linearly polarized
	//const dcmplx Ein1=dcmplx(1/sqrt(2), 0), Ein2=dcmplx(0, 1/sqrt(2)); //circularly polarized

	const double m0[3]={1,0,0}, n0[3]={0,1,0}, s0[3]={0,0,1};

	// double sn[3]={0.14106735979665894,0,-0.99};
	// TotalNumOfPaths does not include the paths which scatters only once 
	// It counts both forward and backward paths
	// I0: single scattering contribution, I_NoCBS and I_CBS with scattering# >= 2
	
	int N = 1; 76; //The number of thetaB angles measured
	int N2 = 1; 73; //The number of phiB angles measured
	double TotalNumOfPaths=0, I0=0, I_NoCBS[N][N2], I_NoCBS1[N][N2], I_NoCBS2[N][N2], I_CBS[N][N2], I_CBS1[N][N2], I_CBS2[N][N2];
	
	for (int j=0; j<N; j++) {
	
		for(int i=0; i<N2; i++){
	
			I_NoCBS[j][i]  = 0;
			I_CBS[j][i]    = 0;
			I_NoCBS1[j][i] = 0;
			I_CBS1[j][i]   = 0;
			I_NoCBS2[j][i] = 0;
			I_CBS2[j][i]   = 0;
		}
	}
	// the number of photons escape from the front or back surface
	int Reflection=0, Transmission=0;

	for (int i = 1; i <= photons; i++) {
		//launch a linear polarized light
		ph->launch ();
		//launch a circularly polarized light
		//ph->launch (Ein1, Ein2, m0[0], m0[1], m0[2], s0[0], s0[1], s0[2]);

		dcmplx T[2][2]; dcmplx TRev[2][2];

		double m1[3], n1[3], s1[3], s1rev[3], n1rev[3];
		double x1, y1, z1, sinphi0, cosphi0;
		double l, k = 1;

		double nplus = 1.01;
		double nminus = 1.82;
		double ndelta = nplus - nminus;

		double d, dRev;


		dcmplx S1_F1, S2_F1, S1_Fn, S2_Fn;
		dcmplx S1_R1, S2_R1, S1_Rn, S2_Rn;



		// T matrix before scattering
		T[0][0]=1;
		T[0][1]=0;
		T[1][0]=0;
		T[1][1]=1;

		TRev[0][0] = 1;
		TRev[0][1] = 0;
		TRev[1][0] = 0;
		TRev[1][1] = 1;

		if (i % 50 == 0) fprintf(stderr, "%d\n", i);

		while (1) {

			l = ph->move ();

			d = (ndelta / 2.0) * k * l;
			dRev = (ndelta/ 2.0) * k * -l;
			

			if (ph->z < 0) { Reflection++; break; }
			if (ph->z > L) { Transmission++; break; }

			// current m, n, s vectors before scattering
			// nrev = -n, srev = -s in the reversed path
			double m[3], n[3], s[3], nrev[3], srev[3];

			m[0] = ph->P[0][0]; m[1] = ph->P[0][1]; m[2] = ph->P[0][2];
			n[0] = ph->P[1][0]; n[1] = ph->P[1][1]; n[2] = ph->P[1][2];
			s[0] = ph->P[2][0]; s[1] = ph->P[2][1]; s[2] = ph->P[2][2];
			nrev[0] = -n[0]; nrev[1] = -n[1]; nrev[2] = -n[2];
			srev[0] = -s[0]; srev[1] = -s[1]; srev[2] = -s[2];

			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//
			// scattering event occurs
			//
			// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			ph->scatter ();

			// current scattering angles and rotation matrix
			double cosphi, sinphi;
			double cosphin, sinphin;
			
			FaradayRotation(d, &cosphi, &sinphi);	
			rotationMatrix(m, n, ph->P[1], &cosphi, &sinphi);
			

			if (ph->nsct==1) {
				// propagation direction after first scattering event
				m1[0] = ph->P[0][0]; m1[1] = ph->P[0][1]; m1[2] = ph->P[0][2];
				n1[0] = ph->P[1][0]; n1[1] = ph->P[1][1]; n1[2] = ph->P[1][2];
				s1[0] = ph->P[2][0]; s1[1] = ph->P[2][1]; s1[2] = ph->P[2][2];

				// used in the reversed path
				n1rev[0] = -n1[0]; n1rev[1] = -n1[1]; n1rev[2] = -n1[2];
				s1rev[0] = -s1[0]; s1rev[1] = -s1[1]; s1rev[2] = -s1[2];

				// rotation matrix R(phi0)
				sinphi0 = sinphi;
				cosphi0 = cosphi;

				// Scattering amplitude at the 1st scattering site in the foward path
				S1_F1 = ph->S1;
				S2_F1 = ph->S2;

				// Location of the first scattering site
				x1 = ph->x; y1 = ph->y; z1 = ph->z;
			}



			for(int j=0; j<N; j++) {
				for(int i=0; i<N2; i++) {
					double phiB=(i*5)*pi/180;
					double thetaB=(j*theta_step)*(pi/180);
					//double thetaB=(.1)*(pi/180);
					double sn[3]={-cos(phiB)*sin(thetaB),sin(phiB)*sin(thetaB),-cos(thetaB)};


					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//
					// Calculate E in the forward direction
					//
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//
					// Eout = S^(n)*R(phi_n-1)*T*S^(1)*R(phi0)*Ein, if nsct >= 2
					// Eout0 = S^(n)*R(phi_n-1)*Ein,                if nsct = 1  
					// Aout is the amplitude of light emerging from the medium
					dcmplx Eout0[2], Eout[2];
					double Aout;

					// the normal direction nn for scattering from s --> sn
					// and the rotation angle phi_(n-1)
					// mn = nn x sn
					// the outgoing electric field is E1*mn + E2*nn
					double mn[3], nn[3], cosphiNminusOne, sinphiNminusOne;
					
					FaradayRotation(d, &cosphiNminusOne, &sinphiNminusOne);
					normalDirection(s, n, sn, nn);
					rotationMatrix(m, n, nn, &cosphiNminusOne, &sinphiNminusOne);
					axb(nn, sn, mn);

					// the scattering stuff in the forward path (S^(n)) at the nth scatterer
					sct.phasef_lu(sn[0]*s[0]+sn[1]*s[1]+sn[2]*s[2], &S1_Fn, &S2_Fn);

					// the escaping without further scattering probability
					// the factor 2 accounts for the electric field
					double e = exp(ph->z / (2*sn[2]) );

					if (ph->nsct == 1) {
						/* the extra factor e/sqrt(csca): e for escaping probability,
						   and 1/sqrt(csca) scales Eout to be the electric field 
						   over the unit solid angle in sn direction */

						Eout0[0] = S2_Fn * (cosphiNminusOne * Ein1 + sinphiNminusOne * Ein2);
						Eout0[1] = S1_Fn * (-sinphiNminusOne * Ein1 + cosphiNminusOne * Ein2);

						Aout = sqrt(norm(Eout0[0]) + norm(Eout0[1]));
						Eout0[0] /= Aout;
						Eout0[1] /= Aout;

						Aout *= e/sqrt(csca);

						// factor 2 counts the contribution from the backward path
						I0 += 2*Aout*Aout;
					}

					if (ph->nsct >= 2) {
						dcmplx tmp1, tmp2;
						tmp1 = T[0][0]*S2_F1*(cosphi0*Ein1 + sinphi0*Ein2) + T[0][1]*S1_F1*(-sinphi0*Ein1 + cosphi0*Ein2);
						tmp2 = T[1][0]*S2_F1*(cosphi0*Ein1 + sinphi0*Ein2) + T[1][1]*S1_F1*(-sinphi0*Ein1 + cosphi0*Ein2);

						/* the extra factor e/sqrt(csca): e for escaping probability, 
						   and 1/sqrt(csca) scales Eout to be the electric field 
						   over the unit solid angle in sn direction */

						Eout[0] = S2_Fn*(cosphiNminusOne*tmp1 + sinphiNminusOne*tmp2);
						Eout[1] = S1_Fn*(-sinphiNminusOne*tmp1 + cosphiNminusOne*tmp2);

						Aout = sqrt(norm(Eout[0]) + norm(Eout[1]));
						Eout[0] /= Aout;
						Eout[1] /= Aout;

						Aout *= e/sqrt(csca);

#if 0
						// comparison with the point estimator
						double td, deposit, Q[3][3];
						dcmplx Ed1, Ed2;
						ph->pointestimator(&td, &deposit, Q, &Ed1, &Ed2, 0, sn[0], sn[1], sn[2]);
						printf("\n**** %d\n", ph->nsct);
						printf("z, e:  %f, %f\n", ph->z, e);
						printf("Eout1: %g %g\n", Eout[0].real(), Eout[0].imag());
						printf("Eout2: %g %g\n", Eout[1].real(), Eout[1].imag());
						printf("Amplitude: %g\n", Aout);
						printf("Ed1: %g %g\n", Ed1.real(), Ed1.imag());
						printf("Ed2: %g %g\n", Ed2.real(), Ed2.imag());
						printf("Amplitude: %g\n", sqrt(deposit));
						printf("Q:  %f %f %f\n", Q[0][0], Q[0][1], Q[0][2]);
						printf("    %f %f %f\n", Q[1][0], Q[1][1], Q[1][2]);
						printf("    %f %f %f\n", Q[2][0], Q[2][1], Q[2][2]);
						printf("mn: %f %f %f\n", mn[0], mn[1], mn[2]);
						printf("nn: %f %f %f\n", nn[0], nn[1], nn[2]);
#endif
					}


					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//
					// Calculate E in the reverse direction
					//
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//
					dcmplx EoutRev[2];
					double AoutRev;
					// first and second rotation angles phi_n and phi_n_prime in the reverse path: light moves from the incident direction and scattered from s0 into -s_(n-1)=-s[3]=srev direction and then the coordinate system is aligned to {m, -n, -s} which is the reversal of that in the forward path.
					// nPrime the normal for the scattering from s0 -> -s=srev
					// with rotation angle phin
					double nPrime[3] /*sinphin, cosphin*/;
		
					FaradayRotation(dRev, &sinphin, &cosphin),
					normalDirection(s0, n0, srev, nPrime);
					rotationMatrix(m0, n0, nPrime, &cosphin, &sinphin);

					// the scattering stuff for the reverse path (S^(n)) at the nth scatterer
					sct.phasef_lu(s0[0]*srev[0]+s0[1]*srev[1]+s0[2]*srev[2],&S1_Rn,&S2_Rn);

					// mPrime = nPrime x -s
					// rotation angle phinprime from (m'', n'', -s) -> (m, -n, -s)
					double mPrime[3], sinphinprime, cosphinprime;

					axb(nPrime, srev, mPrime);
					rotationMatrix(mPrime, nPrime, nrev, &cosphinprime, &sinphinprime);

					// n'' is the normal for the scattering from -s1 --> sn
					// with rotation angle phi1'
					double nPrimePrime[3], cosphi1prime, sinphi1prime;
					
					FaradayRotation(dRev, &cosphi1prime, &sinphi1prime);
					normalDirection(s1rev, n1rev, sn, nPrimePrime);
					rotationMatrix(m1, n1rev, nPrimePrime, &cosphi1prime, &sinphi1prime);

					// the scattering stuff for the reverse path (S^(1)) at the 1st site 
					sct.phasef_lu(s1rev[0]*sn[0]+s1rev[1]*sn[1]+s1rev[2]*sn[2],&S1_R1,&S2_R1);

					// m'' = n'' x sn
					// the outgoing field in the reverse path is E1*m'' + E2*n''
					double mPrimePrime[3];
					axb(nPrimePrime, sn, mPrimePrime);

					// the escaping probability without further scattering
					// the fcator 2 accounts for the electric field
					double e1=exp(z1/(2*sn[2]));

					if (ph->nsct == 1) {
						/* the extra factor e/sqrt(csca): e for escaping probability,
						   and 1/sqrt(csca) scales Eout to be the electric field 
						   over the unit solid angle in sn direction */

						EoutRev[0] = S2_Fn * (cosphi1prime * Ein1 
									+ sinphi1prime * Ein2);

						EoutRev[1] = S1_Fn * (-sinphi1prime * Ein1 
									+ cosphi1prime * Ein2);

						AoutRev = sqrt(norm(Eout0[0]) + norm(Eout0[1]));
						EoutRev [0] /= AoutRev;
						EoutRev [1] /= AoutRev;

						AoutRev *= e/sqrt(csca);

						// factor 2 counts the contribution from the backward path
						I0 += 2*Aout*Aout;
					}


					if (ph->nsct >= 2) {

						dcmplx tmp1, tmp2;

						tmp1 = TRev[0][0]*S2_F1*(cosphi1prime * Ein1 + sinphi1prime * Ein2) 
							+ TRev[0][1]*S1_F1*(-sinphi1prime * Ein1 + cosphi1prime * Ein2);

						tmp2 = TRev[1][0]*S2_F1*(cosphi1prime * Ein1 + sinphi1prime * Ein2) 
							+ TRev[1][1]*S1_F1*(-sinphi1prime * Ein1 + cosphi1prime * Ein2);

						/* the extra factor e/sqrt(csca): e for escaping probability, 
						   and 1/sqrt(csca) scales Eout to be the electric field 
						   over the unit solid angle in sn direction */

						EoutRev [0] = S2_Fn*(cosphiNminusOne*tmp1 + sinphiNminusOne*tmp2);
						EoutRev [1] = S1_Fn*(-sinphiNminusOne*tmp1 + cosphiNminusOne*tmp2);

						AoutRev = sqrt(norm(Eout[0]) + norm(Eout[1]));
						EoutRev [0] /= AoutRev;
						EoutRev [1] /= AoutRev;

						AoutRev *= e/sqrt(csca);

					}


					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					//
					// Compute I_NoCBS and I_CBS
					//
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					if (ph->nsct >= 2) {
						double Dphi = ((2*pi*n_water)/lambda)*lsa*((s0[0]+sn[0])*(ph->x-x1)+(s0[1]+sn[1])*(ph->y-y1)+(s0[2]+sn[2])*(ph->z-z1));
						dcmplx f(cos(Dphi), sin(Dphi));

						//printf("Dphi: %f\n", Dphi);
						//printf("Aout: %g\n", Aout);
						//printf("AoutRev: %g\n", AoutRev);

						// I_NoCBS (discount the reversal path probability)
						I_NoCBS[j][i] += Aout*Aout + Aout*Aout;
						I_NoCBS1[j][i] += Aout*Aout*norm(Eout[0]*mn[0]+Eout[1]*nn[0]) + Aout*Aout*norm(EoutRev[0]*mPrimePrime[0]+EoutRev[1]*nPrimePrime[0]);
						I_NoCBS2[j][i] += Aout*Aout*norm(Eout[0]*mn[1]+Eout[1]*nn[1]) + Aout*Aout*norm(EoutRev[0]*mPrimePrime[1]+EoutRev[1]*nPrimePrime[1]);

						// I_CBS
						dcmplx Ex, Ey, Ez;
						Ex=(Eout[0]*mn[0]+Eout[1]*nn[0])*Aout + (EoutRev[0]*mPrimePrime[0]+EoutRev[1]*nPrimePrime[0])*Aout*f;
						Ey=(Eout[0]*mn[1]+Eout[1]*nn[1])*Aout + (EoutRev[0]*mPrimePrime[1]+EoutRev[1]*nPrimePrime[1])*Aout*f;
						Ez=(Eout[0]*mn[2]+Eout[1]*nn[2])*Aout + (EoutRev[0]*mPrimePrime[2]+EoutRev[1]*nPrimePrime[2])*Aout*f;

						I_CBS[j][i] += norm(Ex)+norm(Ey)+norm(Ez);
						// light is leaving in -z direction, so (Ex, Ey, Ez) -> (Ex, -Ey, -Ez) in the coordinate system (xhat, -yhat, -zhat).
						I_CBS1[j][i] += norm(Ex);
						I_CBS2[j][i] += norm(Ey);
					}

				}// end for i loop
			}// end for j loop


			// update T matrix if needed
			if ( ph->nsct >=2 ) {
				// T matrix
				// T = SRT,  R = | cosphi, sinphi|
				//               |-sinphi, cosphi|;
				dcmplx Tnew[2][2]; dcmplx TRevnew[2][2];

				Tnew[0][0] = ph->S2*(cosphi*T[0][0] + sinphi*T[1][0]);
				Tnew[0][1] = ph->S2*(cosphi*T[0][1] + sinphi*T[1][1]);
				Tnew[1][0] = ph->S1*(-sinphi*T[0][0] + cosphi*T[1][0]);
				Tnew[1][1] = ph->S1*(-sinphi*T[0][1] + cosphi*T[1][1]);

				T[0][0] = Tnew[0][0]; T[0][1] = Tnew[0][1];
				T[1][0] = Tnew[1][0]; T[1][1] = Tnew[1][1];


				TRevnew[0][0] = ph->S2*(cosphin*TRev[0][0] + sinphin*TRev[1][0]);
				TRevnew[0][1] = ph->S2*(cosphin*TRev[0][1] + sinphin*TRev[1][1]);
				TRevnew[1][0] = ph->S1*(-sinphin*TRev[0][0] + cosphin*TRev[1][0]);
				TRevnew[1][1] = ph->S1*(-sinphin*TRev[0][1] + cosphin*TRev[1][1]);

				TRev[0][0] = TRevnew[0][0]; TRev[0][1] = TRevnew[0][1];
				TRev[1][0] = TRevnew[1][0]; TRev[1][1] = TRevnew[1][1];



				TotalNumOfPaths += 2;
			}

			// printf("Sct#, I_NoCBS, I_CBS: %g %g %g\n", TotalNumOfPaths, I_NoCBS[j][i], I_CBS[j][i]);

			// if (ph->nsct == 2) break;
		} // end while

	}       
	/*
	   printf("\nTotal number of photons used: %d\n", photons);
	   printf("Total number of paths (including forward and backward): %g\n", TotalNumOfPaths); 
	   printf("I Single backscattering:     %g\n", I0);
	   printf("I_NoCBS, I_NoCBS1, I_NoCBS2: %g %g %g\n", I_NoCBS[j][i], I_NoCBS1[j][i], I_NoCBS2[j][i]);
	   printf("I_CBS,   I_CBS1,   I_CBS2:   %g %g %g\n", I_CBS[j][i], I_CBS1[j][i], I_CBS2[j][i]);
	   printf("I_CBS/I_NoCBS:   %g\n", I_CBS[j][i]/I_NoCBS[j][i]);
	   printf("I_CBS1/I_NoCBS1: %g\n", I_CBS1[j][i]/I_NoCBS1[j][i]);
	   printf("I_CBS2/I_NoCBS2: %g\n", I_CBS2[j][i]/I_NoCBS2[j][i]);
	 */
	fp_out << "\nTotal number of photons used: " << photons << endl;
	fp_out << "Total number of paths (including forward and backward): " << TotalNumOfPaths << endl; 
	fp_out << "I Single backscattering: " << I0 << endl;
	fp_out << "Transmission percentage: " << 100*((Transmission+0.0)/photons) << endl;


	for(int j=0; j<N; j++){
		for(int i=0; i<N2; i++){
			//printf("ThetaB: %g", j*0.002);
			//printf(" phiB: %i", i*5);
			//printf(" I_CBS2/I_NoCBS2: %g\n", I_CBS2[j][i]/I_NoCBS2[j][i]);

			fp_out /*<< "ThetaB: "   */ << j*theta_step   << " ";
			fp_out /*<< " phiB: "    */ << i*5            << " ";
			fp_out /*<< " I_CBS: "   */ << I_CBS[j][i]    << " "; 
			fp_out /*<< " I_NoCBS: " */ << I_NoCBS[j][i]  << " ";
			fp_out /*<< " I_CBS1: "  */ << I_CBS1[j][i]   << " ";
			fp_out /*<< " I_NoCBS1: "*/ << I_NoCBS1[j][i] << " ";
			fp_out /*<< " I_CBS2: "  */ << I_CBS2[j][i]   << " ";
			fp_out /*<< " I_NoCBS2: "*/ << I_NoCBS2[j][i] << endl;

		}
	}

	delete ph;

	fp_out.close();
}







