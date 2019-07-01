

			// LambdaTwistAndEllipticP3P.c 

			// Author: M. Q. Rieck 
			// Date: 6/30/2019  

// This code implements and tests two P3P solvers, the "Lambda 
// Twist" method of Mikael Persson and Klaus Nordberg, and my own 
// elliptic-curve-based method. The Lambda Twist code here is a
// combination of tested formulas in their paper, "Lambda Twist: 
// an accurate fast robust perspective three point (P3P) solver"
// and their C++ source code, translated to C. 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Comment out if you want main() to use the Elliptic Curve method 
// instead of the Lambda Twist method:
//#define TEST_LAMBDA_TWIST

// How many trials should main() execute? 
#define NUMBER_TRIALS 1000000

// A couple levels of debugging capability:
//#define DEBUG
//#define DEBUG2
#define DEBUG2_CUTOFF 1

// Can restrict methods to reject trials outside a "attack angle range"
#define ATTACK_ANGLE_MIN 0
#define ATTACK_ANGLE_MAX 45

// Use the Newton-Raphson method to find root of cubic polynomial,
// provided by Persson and Nordberg, instead of my algebraic method? 
#define USE_CUBICK_FOR_CUBICS
// Parameters for Newton-Raphson method
#define KLAS_P3P_CUBIC_SOLVER_ITER 70 
#define NUMERIC_LIMIT 1e-40

#define FALSE 0 
#define TRUE 1 
#define DOT(u,v)	(u[0]*v[0]+u[1]*v[1]+u[2]*v[2])
#define LEN(v)		sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
#define ROTF(i)		((i+1)%3)
#define ROTB(i)		((i+2)%3)
#define SQR(x)		((x)*(x))
#define CUBE(x)		((x)*(x)*(x))
#define QUART(x)	(((x)*(x))*((x)*(x)))
#define SIGN(x)		((x) < 0 ? -1 : 1)

#define NO_ACCEPTABLE_ESTIMATES		-1
#define NO_SYSTEM_SOLUTIONS			-2
#define SINGULARITY_TOO_SMALL		-3
#define NO_ELIGIBLE_ROTATION		-4
#define REJECT_TESTING_PARAMETERS	-999

// Globals used for timimg tests
clock_t start_time, end_time, exec_time;
double getTime() { return (double)exec_time/CLOCKS_PER_SEC; }

// Lower case versions are now used instead, which can be set dynamically
int attack_angle_min = ATTACK_ANGLE_MIN;
int attack_angle_max = ATTACK_ANGLE_MAX;
void set_attack_angle_min(int degrees) { attack_angle_min = degrees; }
void set_attack_angle_max(int degrees) { attack_angle_max = degrees; }

// Normalize a vector 
void normalize(double v[3], double w[3]) {
	double len = LEN(v); 
	w[0] = v[0]/len; w[1] = v[1]/len; w[2] = v[2]/len; 
}

// Cross product of two vectors 
void crossProduct(double u[3], double v[3], double w[3]) {
	for(int i=0; i<3; i++) 
		w[i] = u[ROTF(i)]*v[ROTB(i)]-u[ROTB(i)]*v[ROTF(i)];
}

// Display a matrix 
void showMatrix(double m[3][3]) {
	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) 
			printf("%lf ", m[i][j]); 
		printf("\n");
	}
}

// Display a vector 
void showVector(double v[3]) {
	for(int i=0; i<3; i++) 
		printf("%lf ", v[i]); 
	printf("\n");
}

// Copy a matrix 
void copyMatrix(double m[3][3], double mp[3][3]) {
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) 
			mp[i][j] = m[i][j];
}

// Add two vectors 
void addVectors(double u[3], double v[3], double w[3]) {
	w[0] = u[0] + v[0]; 
	w[1] = u[1] + v[1]; 
	w[2] = u[2] + v[2]; 
}

// Subtract two vectors 
void subVectors(double u[3], double v[3], double w[3]) {
	w[0] = u[0] - v[0]; 
	w[1] = u[1] - v[1]; 
	w[2] = u[2] - v[2]; 
}

// Subtract two matrices 
void subMatrices(double m1[3][3], double m2[3][3], double m3[3][3]) {
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			m3[i][j] = m1[i][j] - m2[i][j];
}

// Multiply two matrices 
void multMatrices(double m1[3][3], double m2[3][3], double m3[3][3]) {
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++) {
			m3[i][j] = 0;
			for(int k=0; k<3; k++)
				m3[i][j] += m1[i][k]*m2[k][j]; 
		}
}

// Multiply a matrix and a vector 
void multMatrixVector(double m[3][3], double u[3], double v[3]) {
	for(int i=0; i<3; i++) {
		v[i] = 0; 
		for(int j=0; j<3; j++) 
			v[i] += m[i][j]*u[j]; 
	}
}

// Transpose a matrix 
void transposeMatrix(double m[3][3], double mt[3][3]) {
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			mt[i][j] = m[j][i];  
}

// Invert a matrix 
void invertMatrix(double m[3][3], double mi[3][3]) {
	double det; 
	det = 
	    m[0][0]*m[1][1]*m[2][2] 
	  + m[0][1]*m[1][2]*m[2][0] 
	  + m[0][2]*m[1][0]*m[2][1] 
	  - m[0][0]*m[1][2]*m[2][1] 
	  - m[0][1]*m[1][0]*m[2][2] 
	  - m[0][2]*m[1][1]*m[2][0]; 
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			mi[i][j] = (m[ROTF(j)][ROTF(i)]*m[ROTB(j)][ROTB(i)]
			  - m[ROTF(j)][ROTB(i)]*m[ROTB(j)][ROTF(i)]) / det;
}

// Return a random real number in a specified range 
double randomReal(double a, double b) {
   return a + (float)rand()/(float)(RAND_MAX/(b-a));
}

// Create the matrix for a random rotation
void randomRotation(double m[3][3]) {
	double a, b, c, d, theta, phi, psi, pi, w[3];
  	pi = 4*atan(1); // goofy way to get pi in C 
	theta = randomReal(-pi, pi);  // longitude 
	phi = acos(randomReal(-1,1)); // lattitude wrt north pole 
	psi = randomReal(-pi, pi);    // rotation amount
	w[0] = sin(phi)*cos(theta); 
	w[1] = sin(phi)*sin(theta); 
	w[2] = cos(phi);
	a = cos(psi/2); 
	b = w[0]*sin(psi/2); 
	c = w[1]*sin(psi/2); 
	d = w[2]*sin(psi/2); 
	m[0][0] = a*a+b*b-c*c-d*d; 
	m[0][1] = 2*(b*c-a*d); 
	m[0][2] = 2*(b*d+a*c); 
	m[1][0] = 2*(b*c+a*d); 
	m[1][1] = a*a-b*b+c*c-d*d; 
	m[1][2] = 2*(c*d-a*b); 
	m[2][0] = 2*(b*d-a*c); 
	m[2][1] = 2*(c*d+a*b); 
	m[2][2] = a*a-b*b-c*c+d*d;	
}

// Rotate two given unit vectors to two other given unit vectors, 
// assuming the angle between the vector pairs is the same
int rotateIntoPosition(double u[3], double up[3], double v[3], 
  double vp[3], double m[3][3]) {
  	int extra; 
  	double a, b, c, d, theta, phi, psi, pi, uu[3], vv[3], w[3], 
  	  tempVec1[3], tempVec2[3], tempMat1[3][3], tempMat2[3][3], 
  	  m1[3][3], m2[3][3];
  	extra = FALSE;
	pi = 4*atan(1);
  	uu[0] = u[0]; uu[1] = u[1]; uu[2] = u[2];
  	vv[0] = v[0]; vv[1] = v[1]; vv[2] = v[2];
  	subVectors(uu, up, tempVec1);
	subVectors(vv, vp, tempVec2);
	crossProduct(tempVec1, tempVec2, w); 
	while (LEN(w) < .001) {
		theta = randomReal(-pi, pi);
		phi = acos(randomReal(-1,1)); 
		psi = randomReal(-pi, pi); 
		w[0] = sin(phi)*cos(theta); 
		w[1] = sin(phi)*sin(theta); 
		w[2] = cos(phi);
		a = cos(psi/2); 
		b = w[0]*sin(psi/2); 
		c = w[1]*sin(psi/2); 
		d = w[2]*sin(psi/2); 
		m1[0][0] = a*a+b*b-c*c-d*d; 
		m1[0][1] = 2*(b*c-a*d); 
		m1[0][2] = 2*(b*d+a*c); 
		m1[1][0] = 2*(b*c+a*d); 
		m1[1][1] = a*a-b*b+c*c-d*d; 
		m1[1][2] = 2*(c*d-a*b); 
		m1[2][0] = 2*(b*d-a*c); 
		m1[2][1] = 2*(c*d+a*b); 
		m1[2][2] = a*a-b*b-c*c+d*d;	
#ifdef DEBUG
		printf("random rotate\n");	
		showMatrix(m1);
		transposeMatrix(m1, tempMat1);
		multMatrices(m1, tempMat1, tempMat2);
		showMatrix(tempMat2);
#endif
		multMatrixVector(m1, u, uu); 
		multMatrixVector(m1, v, vv); 
		extra = TRUE; 
		subVectors(uu, up, tempVec1);
		subVectors(vv, vp, tempVec2);
		crossProduct(tempVec1, tempVec2, w); 
	}
	normalize(w, w);
#ifdef DEBUG 
	printf("w = %lf %lf %lf\n", w[0], w[1], w[2]);
	showVector(uu);
	showVector(up);
	showVector(w);
#endif 
	theta = acos((DOT(uu,up) - SQR(DOT(uu,w)))
	  / (1-SQR(DOT(uu,w)))); 
	crossProduct(uu, up, tempVec1); 
	if (DOT(tempVec1,w) < 0) theta = -theta; 
	a = cos(theta/2); 
	b = w[0]*sin(theta/2); 
	c = w[1]*sin(theta/2); 
	d = w[2]*sin(theta/2); 
	m[0][0] = a*a+b*b-c*c-d*d; 
	m[0][1] = 2*(b*c-a*d); 
	m[0][2] = 2*(b*d+a*c); 
	m[1][0] = 2*(b*c+a*d); 
	m[1][1] = a*a-b*b+c*c-d*d; 
	m[1][2] = 2*(c*d-a*b); 
	m[2][0] = 2*(b*d-a*c); 
	m[2][1] = 2*(c*d+a*b); 
	m[2][2] = a*a-b*b-c*c+d*d;
	if (extra) {
		multMatrices(m, m1, m2);
#ifdef DEBUG
		showMatrix(m2);
		printf("copying\n");
#endif
		copyMatrix(m2, m);
	}
	multMatrixVector(m, u, tempVec1);
	multMatrixVector(m, v, tempVec2);
#ifdef DEBUG
	printf("theta = %lf\n", theta);	
	printf("a b c d = %lf %lf %lf %lf\n", a, b, c, d);	
	showMatrix(m); 
	showVector(tempVec1);
	showVector(tempVec2);
#endif 
	return extra; 	// was a random rotation needed? 
}

// Rotate so first unit vector becomes vertical, and projection 
// of other two unit vectors onto the xy-axis become symmetrically
// placed about the x-axis 
int rotateTriple(double u[3], double v[3], double w[3], 
  double m[3][3]) {
  	double uu[3], vv[3], ww[3], nv[3], nw[3], n[3], e2[3], e3[3];
  	normalize(u, uu); 
  	normalize(v, vv); 
  	normalize(w, ww); 
  	crossProduct(vv, uu, nv); 
  	normalize(nv, nv); 
  	crossProduct(uu, ww, nw);
  	normalize(nw, nw); 
  	addVectors(nv, nw, n); 
  	normalize(n, n); 
#ifdef DEBUG
  	showVector(uu); 
  	showVector(vv); 
  	showVector(ww); 
  	showVector(nv); 
  	showVector(nw); 
  	showVector(n); 
#endif 
  	e3[2] = e2[1] = 1; 
  	e3[0] = e3[1] = e2[0] = e2[2] = 0; 
  	return rotateIntoPosition(uu, e3, n, e2, m);
}

// Newton-Raphson method for obtaining a real root of a cubic polynomial
// (Mostly copy-and-pasted code from the original Lambda Twist code 
// of Persson and Nordberg)
double cubick(double b, double c, double d) {
#if 1
    double r0;
    if (b*b  >= 3.0*c) {
        double v = sqrt(b*b-3.0*c);
        double t1 = (-b - v)/(3.0);
        double k = ((t1+b)*t1+c)*t1+d;
        if (k > 0.0) {
            r0 = t1 - sqrt(-k/(3.0*t1 + b));
        } else {
            double t2 = (-b + v)/(3.0);
            k = ((t2+b)*t2+c)*t2+d;
            r0 = t2 + sqrt(-k/(3.0*t2 + b));
        }
    }
    else{
        r0 = -b/3.0;
        if (fabs(((3.0*r0+2.0*b)*r0+c))<1e-4) r0 += 1;
    }
#endif
    double fx,fpx;    
    for (unsigned int cnt = 0; cnt < KLAS_P3P_CUBIC_SOLVER_ITER; ++cnt)
    {
        fx=(((r0+b)*r0+c)*r0+d);
        if (cnt < 7 || fabs(fx) > NUMERIC_LIMIT) {
            fpx=((3.0*r0+2.0*b)*r0+c);
            r0-= fx/fpx;
        }
        else
            break;
    }
    return r0;
}

// Return a real root of a cubic polynomial (mostly algebraic method)
double realRootOfCubic(double a, double b, double c, double d) {
	double p, q, del, sqt, root, rt1, rt2, rt3, temp;
	int flip = FALSE; 

	if (fabs(a) < fabs(d)) {
		temp = a; a = d; d = temp; 
		temp = b; b = c; c = temp; 
		flip = TRUE; 
	}
// Depressed cubic is t^3 + p t + q where t = x + b / (3*a)
	p = (3*a*c-SQR(b))/(3*SQR(a));
	q = (2*CUBE(b)-9*a*b*c+27*SQR(a)*d)/(27*CUBE(a));
	del = 4*CUBE(p) + 27*SQR(q); 
	if (del >= 0) {
	   sqt = sqrt(del/108);
	   root = cbrt(-q/2+sqt)+cbrt(-q/2-sqt);
#ifdef DEBUG
       printf("only real root of depressed cubic: %lf\n", root);
#endif
	} else {
       rt1 =  2*sqrt(-p/3)*cos(acos( 3*q*sqrt(-3/p)/(2*p))/3); 
       rt2 = -2*sqrt(-p/3)*cos(acos(-3*q*sqrt(-3/p)/(2*p))/3); 
       rt3 = -rt1-rt2; 
#ifdef DEBUG
       printf("roots of depressed cubic: %lf %lf %lf\n", 
       	 rt1, rt2, rt3);
#endif
       root = rt1;
       if (rt2 > root) root = rt2;
       if (rt3 > root) root = rt3;
	}
	root = root - b/(3*a); 
#ifdef DEBUG
	printf("real root of cubic: %lf\n", root);
#endif
    if (flip) return 1/root; else return root; 
}

// Get roots of a quadratic polynomial 
// (Mostly copy-and-pasted code from the original Lambda Twist code 
// of Persson and Nordberg)
int root2real(double b, double c, 
  double* r1, double* r2){
    double v = b*b - 4.0*c;
    if (v < 0){
        *r1 = *r2 = 0.5*b;
        return FALSE;
    }
    double y = sqrt(v);
    if (b < 0){
        *r1= 0.5*(-b+y);
        *r2= 0.5*(-b-y);
    } else {
        *r1= 2.0*c/(-b+y);
        *r2= 2.0*c/(-b-y);
    }
    return TRUE;
}

// Obtain eigenvalues and eigenvectors for special symmetic matrix 
// (Mostly copy-and-pasted code from the original Lambda Twist code 
// of Persson and Nordberg)
void eigwithknown(double M[3][3], double E[3][3], double L[3]) {
	double b, c, e1, e2, mx0011, prec_0, prec_1, e, tmp, a1, a2, 
	rnorm, tmp2, a21, a22, rnorm2, MT[3][3], v1[3], v2[3], v3[3];

	transposeMatrix(M, MT);
	crossProduct(MT[0], MT[1], v3);
	normalize(v3, v3);
	b = -M[0][0] - M[1][1] - M[2][2];
	c = -SQR(M[0][1]) - SQR(M[0][2]) - SQR(M[1][2])
	  + M[0][0]*(M[1][1] + M[2][2]) + M[1][1]*M[2][2];
	root2real(b, c, &e1, &e2); 
	if (fabs(e1) < fabs(e2)) {
		tmp = e1;
		e1 = e2;
		e2 = tmp;
	}
	L[0] = e1; L[1] = e2, L[2] = 0;
	mx0011 = -M[0][0]*M[1][1];
	prec_0 = M[0][1]*M[1][2] - M[0][2]*M[1][1];
	prec_1 = M[0][1]*M[0][2] - M[0][0]*M[1][2];
	tmp = 1.0 / (e1*(M[0][0] + M[1][1]) + mx0011 
		- e1*e1 + SQR(M[0][1]));
	a1 = -(e1*M[0][2] + prec_0)*tmp;
	a2 = -(e1*M[1][2] + prec_1)*tmp;
	rnorm = 1.0 / sqrt(SQR(a1) + SQR(a2) + 1.0); 
	a1 *= rnorm;
	a2 *= rnorm;
	v1[0] = a1; v1[1] = a2; v1[2] = rnorm;
	tmp2 = 1.0 / (e2*(M[0][0] + M[1][1]) + mx0011
	 	- SQR(e2) + SQR(M[0][1]));
	a21 = -(e2*M[0][2] + prec_0)*tmp2;
	a22 = -(e2*M[1][2] + prec_1)*tmp2;
	rnorm2 = 1.0 / sqrt(SQR(a21) + SQR(a22) + 1.0);
	a21 *= rnorm2;
	a22 *= rnorm2;
	v2[0] = a21; v2[1] = a22; v2[2] = rnorm2;
	E[0][0] = v1[0]; E[0][1] = v2[0]; E[0][2] = v3[0]; 
	E[1][0] = v1[1]; E[1][1] = v2[1]; E[1][2] = v3[1]; 
	E[2][0] = v1[2]; E[2][1] = v2[2]; E[2][2] = v3[2];
}

// Some code to test Persson and Nordberg's code 
void eigwithknownTest() {
	double M[3][3], E[3][3], L[3], ET[3][3], P[3][3], 
	V[3], S[3][3], T[3][3];
	/*  
		It only works when the first eigenvalue is
		negative, and second eigenvalue is positive 
		with a smaller absolute value 
	*/
	while(TRUE) {
		printf("First eigenvalue (negative): ");
		scanf("%lf", &(L[0]));
		printf("Second eigenvalue (positive, smaller): ");
		scanf("%lf", &(L[1]));
		randomRotation(E); 
		printf("Random rotation:\n");
		showMatrix(E); 
		L[2] = 0; 
		S[0][0] = L[0]; S[0][1] = 0; S[0][2] = 0;
		S[1][0] = 0; S[1][1] = L[1]; S[1][2] = 0;
		S[2][0] = 0; S[2][1] = 0; S[2][2] = L[2];
		transposeMatrix(E, ET);
		multMatrices(E, S, T);
		multMatrices(T, ET, M); 
		printf("Symmetric matrix for testing:\n");
		showMatrix(M);
		eigwithknown(M, E, L);
		printf("Resulting orthogonal matrix:\n");  
		showMatrix(E);
		printf("Resulting eigenvalues:\n");  
		showVector(L);
		transposeMatrix(E, ET);
		multMatrices(E, S, T);
		multMatrices(T, ET, M); 
		printf("Should get same symmetric matrix:\n");
		showMatrix(M);
		printf("\n");  
	}
}

// Solve the P3P problem using the Lambda Twist method  
// (Written by blending the contents of Persson and Nordberg's
// research paper and their C++ code)
double lambdaTwistSolver(double v[3][3]) {
	int i, j, k, numRoots, valid; 
	double pi, a12, a23, a31, b12, b23, b31, blob, s12sq, 
	 s23sq, s31sq, p0, p1, p2, p3, gamma, ssq, tmp, ss,  
	 w0, w1, den, a, b, c, tau1, tau2, tau, d, l1, l2, l3, 
	 error, minerror, 
	 vn[3][3], vl[3], sl[3], s[3][3], sn[3][3], vndp[3],
	 sndp[3], D0[3][3], D00[3][3], 
	 M[3][3], E[3][3], L[3], Ls[8][3], T[3][3], ET[3][3], 
	 Mp[3][3], S[3][3], guess[3][3], diff[3][3],
	 tempMat1[3][3], tempMat2[3][3], direction[3], 
	 triangleNormal[3], tempVec[3], 

	errorLimit = 1000,
	tol = 0;

	for(i=0; i<3; i++) {
		vl[i] = LEN(v[i]);
		normalize(v[i], vn[i]);
		for(j=0; j<3; j++) {
			s[i][j] = v[ROTF(i)][j]-v[ROTB(i)][j];
		}
		sl[i] = LEN(s[i]);
		normalize(s[i], sn[i]); 
	}
// reject (not a miss) if attack angle is outside acceptable range  
	addVectors(vn[0], vn[1], tempVec);
	addVectors(tempVec, vn[2], tempVec); 
	normalize(tempVec, direction);
	crossProduct(sn[0], sn[1], triangleNormal);
	normalize(triangleNormal, triangleNormal);
	if (fabs(DOT(direction, triangleNormal)) < 
	  cos(attack_angle_max*atan(1)/45)) return REJECT_TESTING_PARAMETERS;
	if (fabs(DOT(direction, triangleNormal)) > 
	  cos(attack_angle_min*atan(1)/45)) return REJECT_TESTING_PARAMETERS;
// continue for acceptable setups
	start_time = clock();  
	for(i=0; i<3; i++) {
		vndp[i] = DOT(vn[ROTF(i)],vn[ROTB(i)]);
		sndp[i] = DOT(sn[ROTF(i)],sn[ROTB(i)]);
	}
	a23 = DOT(s[0],s[0]); 
	a31 = DOT(s[1],s[1]); 
	a12 = DOT(s[2],s[2]); 
	b23 = vndp[0];
	b31 = vndp[1];
	b12 = vndp[2];
	// They changed notation in their code, using c's 
	// instead of paper's b's; I'm going back to b's 
	blob = b12*b23*b31 - 1.0; 
	s23sq = 1.0 - b23*b23; //SQR(b23); 
	s31sq = 1.0 - b31*b31; //SQR(b31); 
	s12sq = 1.0 - b12*b12; //SQR(b12); 
	p0 = a12*(a12*s23sq-a23*s12sq); 
	p1 =  - ( a23*(a31-a23)*s12sq - a12*a12*s23sq 
	          - 2.0*a12*(blob*a23 + a31*s23sq) );  
	p2 =  2.0*blob*a23*a31 + a31*(2.0*a12+a31)*s23sq 
	         + a23*(a23-a12)*s31sq; 
	p3 =  -a31*(a23*s31sq-a31*s23sq);
#ifdef DEBUG
	printf("v[0] = %lf %lf %lf\n", v[0][0], v[0][1], v[0][2]); 
	printf("v[1] = %lf %lf %lf\n", v[1][0], v[1][1], v[1][2]); 
	printf("v[2] = %lf %lf %lf\n", v[2][0], v[2][1], v[2][2]);
	printf("distances to control points = %lf %lf %lf\n", 
	  vl[0], vl[1], vl[2]);
	printf("vn[0] = %lf %lf %lf\n", vn[0][0], vn[0][1], vn[0][2]); 
	printf("vn[1] = %lf %lf %lf\n", vn[1][0], vn[1][1], vn[1][2]); 
	printf("vn[2] = %lf %lf %lf\n", vn[2][0], vn[2][1], vn[2][2]);
	printf("vndp = %lf %lf %lf\n", vndp[0], vndp[1], vndp[2]); 
	printf("sidelengths = %lf %lf %lf\n", sl[0], sl[1], sl[2]);
	printf("squared sidelengths = %lf %lf %lf\n", 
	  SQR(sl[0]), SQR(sl[1]), SQR(sl[2]));
	printf("sn[0] = %lf %lf %lf\n", sn[0][0], sn[0][1], sn[0][2]); 
	printf("sn[1] = %lf %lf %lf\n", sn[1][0], sn[1][1], sn[1][2]); 
	printf("sn[2] = %lf %lf %lf\n", sn[2][0], sn[2][1], sn[2][2]); 
	printf("sndp = %lf %lf %lf\n", sndp[0], sndp[1], sndp[2]); 
	printf("a23 a31 a12 = %lf %lf %lf\n", a23, a31, a12); 
	printf("b23 b31 b12 = %lf %lf %lf\n", b23, b31, b12); 
	printf("s23sq s31sq s12sq blob = %lg %lg %lg %lg\n", 
	  s23sq, s31sq, s12sq, blob);
	printf("p0 p1 p2 p3 = %lg %lg %lg %lg\n", p0, p1, p2, p3);
#endif
#ifdef USE_CUBICK_FOR_CUBICS
	gamma = cubick(p2/p3, p1/p3, p0/p3); 
#else 
	gamma = realRootOfCubic(1, p2/p3, p1/p3, p0/p3); 
#endif
#ifdef DEBUG
	printf("p0/p3 p1/p3 p2/p3 gamma = %lf %lf %lf %lf\n", 
	  p0/p3, p1/p3, p2/p3, gamma);
	printf("Check root %lf: %lf\n", gamma, 
	  p0 + p1*gamma + p2*SQR(gamma) + p3*CUBE(gamma)); 
#endif
	D0[0][0] = a23*(1.0 + gamma); 
	D0[1][0] = D0[0][1] = -a23*b12;
	D0[2][0] = D0[0][2] = -a23*b31*gamma;
	D0[1][1] = a23 - a12 - a31*gamma;
	D0[2][1] = D0[1][2] = b23*(a12 + a31*gamma);
	D0[2][2] = (a23 - a31)*gamma - a12;
#ifdef DEBUG
	printf("D0:\n");
	showMatrix(D0); 
	printf("determinant = %lf\n", determinant(D0));
#endif
	eigwithknown(D0, E, L); // they use A and V instead of M and E 
	/* ssq is square of s in the paper; */ 
	ssq = (tmp = -L[1]/L[0]) >= 0 ? tmp : 0.0; 
#ifdef DEBUG
	printf("inintial ssq = %lf\n", ssq);
#endif
	if (ssq < 1) {
		tmp = L[0];
		L[0] = L[1]; 
		L[1] = tmp; 
		tmp = E[0][0]; 
		E[0][0]	= E[0][1]; 
		E[0][1] = tmp;
		tmp = E[1][0]; 
		E[1][0]	= E[1][1]; 
		E[1][1] = tmp;
		tmp = E[2][0]; 
		E[2][0]	= E[2][1]; 
		E[2][1] = tmp;
		ssq = (tmp = -L[1]/L[0]) >= 0 ? tmp : 0.0; 
#ifdef DEBUG
		printf("Swapping.\n");
		printf("updated ssq = %lf\n", ssq);
#endif
	}
#ifdef DEBUG
	printf("eigenvalues: ");
	showVector(L); 
	printf("E:\n");
	showMatrix(E); 
	transposeMatrix(E, ET);
	multMatrices(E, ET, T);
	printf("E ET:\n");
	showMatrix(T); 
	S[0][0] = L[0]; S[0][1] = 0; S[0][2] = 0;
	S[1][0] = 0; S[1][1] = L[1]; S[1][2] = 0;
	S[2][0] = 0; S[2][1] = 0; S[2][2] = L[2];
	printf("S:\n");
	showMatrix(S); 
	transposeMatrix(E, ET);
	multMatrices(E, S, T);
	multMatrices(T, ET, D00); 
	printf("Should get D0 again:\n");
	showMatrix(D00);
#endif
	valid = 0;
	minerror = NO_SYSTEM_SOLUTIONS; 
	{
		ss = sqrt(ssq);
		den = ss*E[0][1] - E[0][0]; 
		w0 = (E[1][0] - ss*E[1][1]) / den;
		w1 = (E[2][0] - ss*E[2][1]) / den; 
		a = (a31-a12)*SQR(w1) + 2.0*a12*b31*w1 - a12; 
		b = 2.0*( (a31-a12)*w0*w1 - a31*b12*w1 + a12*b31*w0 ); 
		c = (a31-a12)*SQR(w0) - 2*a31*b12*w0 + a31; 
#ifdef DEBUG
		printf("ss den w0 w1 = %lf %lf %lf %lf\n", 
		  ss, den, w0, w1);
		printf("a b c = %lf %lf %lf\n", a, b, c);
		printf("b*b - 4.0*a*c = %lf\n", b*b - 4.0*a*c); 
#endif
		if (b*b - 4.0*a*c >= 0) {
			minerror = NO_ACCEPTABLE_ESTIMATES; 
			root2real(b/a, c/a, &tau1, &tau2); 
#ifdef DEBUG
			printf("tau1 tau2 = %lf %lf\n", tau1, tau2); 
#endif			
			if (tau1 > 0) {
			  d = a23 / (tau1*(tau1-2*b23)+1.0); 
			  l2 = sqrt(d);
			  l3 = tau1*l2;
			  l1 = w0*l2 + w1*l3; 
#ifdef DEBUG
			  printf("(tau1>0) tau1 l1 l2 l3 = %lf %lf %lf %lf\n", 
			    tau1, l1, l2, l3); 
#endif
			  if (l1 >= 0) {
				Ls[valid][0] = l1; 
				Ls[valid][1] = l2; 
				Ls[valid][2] = l3; 
				valid++; 
				guess[0][0] = l1 * vn[0][0];
				guess[0][1] = l1 * vn[0][1];
				guess[0][2] = l1 * vn[0][2];
				guess[1][0] = l2 * vn[1][0];
				guess[1][1] = l2 * vn[1][1];
				guess[1][2] = l2 * vn[1][2];
				guess[2][0] = l3 * vn[2][0];
				guess[2][1] = l3 * vn[2][1];
				guess[2][2] = l3 * vn[2][2];
				/* coordinate differences */ 
				subMatrices(guess, v, diff);
				/* want relative error so divide by actual distances */ 
				diff[0][0] /= vl[0]; diff[0][1] /= vl[0]; diff[0][2] /= vl[0]; 
				diff[1][0] /= vl[1]; diff[1][1] /= vl[1]; diff[1][2] /= vl[1]; 
				diff[2][0] /= vl[2]; diff[2][1] /= vl[2]; diff[2][2] /= vl[2]; 
				transposeMatrix(diff, tempMat1);
				multMatrices(diff, tempMat1, tempMat2);
				error = sqrt(tempMat2[0][0]+tempMat2[1][1]+tempMat2[2][2]);
				if (error < errorLimit &&
				  (minerror < 0 || error < minerror))
				  	minerror = error; 
			  }
			}
			if (tau2 > 0) {
			  d = a23 / (tau2*(tau2-2*b23)+1.0); 
			  l2 = sqrt(d); 
			  l3 = tau2*l2; 
			  l1 = w0*l2 + w1*l3; 
#ifdef DEBUG
			  printf("(tau2>0) tau2 l1 l2 l3 = %lf %lf %lf %lf\n", 
			    tau2, l1, l2, l3); 
#endif
			  if (l1 >= 0) {
				Ls[valid][0] = l1; 
				Ls[valid][1] = l2; 
				Ls[valid][2] = l3; 
				valid++; 
				guess[0][0] = l1 * vn[0][0];
				guess[0][1] = l1 * vn[0][1];
				guess[0][2] = l1 * vn[0][2];
				guess[1][0] = l2 * vn[1][0];
				guess[1][1] = l2 * vn[1][1];
				guess[1][2] = l2 * vn[1][2];
				guess[2][0] = l3 * vn[2][0];
				guess[2][1] = l3 * vn[2][1];
				guess[2][2] = l3 * vn[2][2];
				/* coordinate differences */ 
				subMatrices(guess, v, diff);
				/* want relative error so divide by actual distances */ 
				diff[0][0] /= vl[0]; diff[0][1] /= vl[0]; diff[0][2] /= vl[0]; 
				diff[1][0] /= vl[1]; diff[1][1] /= vl[1]; diff[1][2] /= vl[1]; 
				diff[2][0] /= vl[2]; diff[2][1] /= vl[2]; diff[2][2] /= vl[2]; 
				transposeMatrix(diff, tempMat1);
				multMatrices(diff, tempMat1, tempMat2);
				error = sqrt(tempMat2[0][0]+tempMat2[1][1]+tempMat2[2][2]);
				if (error < errorLimit &&
				  (minerror < 0 || error < minerror))
				  	minerror = error; 
			  }
			}
		}
	}
	{
		ss = -ss;
		den = ss*E[0][1] - E[0][0]; 
		w0 = (E[1][0] - ss*E[1][1]) / den;
		w1 = (E[2][0] - ss*E[2][1]) / den; 
		a = (a31-a12)*SQR(w1) + 2.0*a12*b31*w1 - a12; 
		b = 2.0*( (a31-a12)*w0*w1 - a31*b12*w1 + a12*b31*w0 ); 
		c = (a31-a12)*SQR(w0) - 2*a31*b12*w0 + a31; 
#ifdef DEBUG
		printf("ss den w0 w1 = %lf %lf %lf %lf\n", 
		  ss, den, w0, w1);
		printf("a b c = %lf %lf %lf\n", a, b, c);
		printf("b*b - 4.0*a*c = %lf\n", b*b - 4.0*a*c); 
#endif
		if (b*b - 4.0*a*c >= 0) {
			minerror = NO_ACCEPTABLE_ESTIMATES; 
			root2real(b/a, c/a, &tau1, &tau2); 
#ifdef DEBUG
			printf("tau1 tau2 = %lf %lf\n", tau1, tau2); 
#endif			
			if (tau1 > 0) {
			  d = a23 / (tau1*(tau1-2*b23)+1.0); 
			  l2 = sqrt(d);
			  l3 = tau1*l2;
			  l1 = w0*l2 + w1*l3; 
#ifdef DEBUG
			  printf("(tau1>0) tau1 l1 l2 l3 = %lf %lf %lf %lf\n", 
			    tau1, l1, l2, l3); 
#endif
			  if (l1 >= 0) {
				Ls[valid][0] = l1; 
				Ls[valid][1] = l2; 
				Ls[valid][2] = l3; 
				valid++; 
				guess[0][0] = l1 * vn[0][0];
				guess[0][1] = l1 * vn[0][1];
				guess[0][2] = l1 * vn[0][2];
				guess[1][0] = l2 * vn[1][0];
				guess[1][1] = l2 * vn[1][1];
				guess[1][2] = l2 * vn[1][2];
				guess[2][0] = l3 * vn[2][0];
				guess[2][1] = l3 * vn[2][1];
				guess[2][2] = l3 * vn[2][2];
				/* coordinate differences */ 
				subMatrices(guess, v, diff);
				/* want relative error so divide by actual distances */ 
				diff[0][0] /= vl[0]; diff[0][1] /= vl[0]; diff[0][2] /= vl[0]; 
				diff[1][0] /= vl[1]; diff[1][1] /= vl[1]; diff[1][2] /= vl[1]; 
				diff[2][0] /= vl[2]; diff[2][1] /= vl[2]; diff[2][2] /= vl[2]; 
				transposeMatrix(diff, tempMat1);
				multMatrices(diff, tempMat1, tempMat2);
				error = sqrt(tempMat2[0][0]+tempMat2[1][1]+tempMat2[2][2]);
				if (error < errorLimit &&
				  (minerror < 0 || error < minerror))
				  	minerror = error; 
			  }
			}
			if (tau2 > 0) {
			  d = a23 / (tau2*(tau2-2*b23)+1.0); 
			  l2 = sqrt(d); 
			  l3 = tau2*l2; 
			  l1 = w0*l2 + w1*l3; 
#ifdef DEBUG
			  printf("(tau2>0) tau2 l1 l2 l3 = %lf %lf %lf %lf\n", 
			    tau2, l1, l2, l3); 
#endif
			  if (l1 >= 0) {
				Ls[valid][0] = l1; 
				Ls[valid][1] = l2; 
				Ls[valid][2] = l3; 
				valid++; 
				guess[0][0] = l1 * vn[0][0];
				guess[0][1] = l1 * vn[0][1];
				guess[0][2] = l1 * vn[0][2];
				guess[1][0] = l2 * vn[1][0];
				guess[1][1] = l2 * vn[1][1];
				guess[1][2] = l2 * vn[1][2];
				guess[2][0] = l3 * vn[2][0];
				guess[2][1] = l3 * vn[2][1];
				guess[2][2] = l3 * vn[2][2];
				/* coordinate differences */ 
				subMatrices(guess, v, diff);
				/* want relative error so divide by actual distances */ 
				diff[0][0] /= vl[0]; diff[0][1] /= vl[0]; diff[0][2] /= vl[0]; 
				diff[1][0] /= vl[1]; diff[1][1] /= vl[1]; diff[1][2] /= vl[1]; 
				diff[2][0] /= vl[2]; diff[2][1] /= vl[2]; diff[2][2] /= vl[2]; 
				transposeMatrix(diff, tempMat1);
				multMatrices(diff, tempMat1, tempMat2);
				error = sqrt(tempMat2[0][0]+tempMat2[1][1]+tempMat2[2][2]);
				if (error < errorLimit &&
				  (minerror < 0 || error < minerror))
				  	minerror = error; 
			  }
			}
		}
	}
#ifdef DEBUG2
	if (minerror < 0 || minerror > DEBUG2_CUTOFF) {
    	printf("minimum error = %lf\n", minerror);
    	printf("distance = %lf\n", 
    	  (LEN(v[0])+LEN(v[1])+LEN(v[2]))/3); 
    	printf("size = %lf\n\n", 
    	  (LEN(s[0])+LEN(s[1])+LEN(s[2]))/3); 
	}
#endif	
	end_time = clock(); 
	exec_time = end_time - start_time;  
	return minerror;
}

// Obtain all real roots of a depressed quartic polynomial 
int allRealRootsOfDepressedQuartic(double p, double q, double r, 
  double shift, double roots[]) {

  	int count = 0; 
    double rootOfCubic, part1, part2, part3, part4, part5;
#ifdef USE_CUBICK_FOR_CUBICS
    rootOfCubic = cubick(p, (SQR(p)-4*r)/4, -SQR(q)/8);
#else
    rootOfCubic = 
      realRootOfCubic(8, 8*p, 2*SQR(p)-8*r, -SQR(q));
#endif
    part1 = sqrt(2*rootOfCubic);
    part2 = sqrt(2/rootOfCubic); 
    part3 = -2*(p+rootOfCubic);
#ifdef DEBUG
    printf("selected root of cubic = %lf\n", rootOfCubic); 
    printf("shift = %lf\n", shift); 
    printf("part1 part2 part3 = %lf %lf %lf\n", part1, part2, part3);
#endif
    part4 = part3-q*part2;
    part5 = part3+q*part2;
    if (part4 >= 0) {
      part4 = sqrt(part4);
      roots[count++] = shift + (part1+part4)/2; 
      roots[count++] = shift + (part1-part4)/2; 
    } 
    if (part5 >= 0) {
      part5 = sqrt(part5);
      roots[count++] = shift + (-part1+part5)/2; 
      roots[count++] = shift + (-part1-part5)/2; 
    } 
    return count; 
}

// Solve a system of two equations, one involving a quartic 
// "quinomial" and the other linear; 
int solveSystem(double A, double B, double C, double D, double E, 
  double F, double G, double H, double U[4], double V[4] ) {
	int numSolns, swap; 
	double p, q, r, shift, temp, tempVec[3]; 
	swap = (fabs(F) < fabs(G));
	if (swap) {
		temp = B; B = C; C = temp;
		temp = F; F = G; G = temp;  
	}
	shift = -H/(2*F);
	p = (C*SQR(F) - D*F*G + B*SQR(G) - SQR(A*H))
	  / (2*SQR(A*F));
	q = (C*SQR(F) - B*SQR(G))*H / (2*SQR(A*F)*F);
	r = (8*E*SQR(F*G) + 2*C*SQR(F*H) + 2*D*F*G*SQR(H)
	  + 2*B*SQR(G*H) + SQR(A*H)*SQR(H)) 
	  / (16*SQR(A*F)*SQR(F));
	numSolns = allRealRootsOfDepressedQuartic(p, q, r, shift, U);
#ifdef DEBUG
	printf("number of system solutions = %d\n", numSolns); 
	printf("shift p q r = %lf %lf %lf %lf\n", shift, p, q, r);
#endif
	for (int i=0; i<numSolns; i++) {
		V[i] = -(F*U[i]+H)/G; 
		if (swap) { 
			temp = U[i]; U[i] = V[i]; V[i] = temp; 
#ifdef DEBUG
		printf("U V = %lf %lf\n", U[i], V[i]);
#endif
		}
	}
	return numSolns;
}

// Solve the P3P problem using the elliptic-curve method
double ellipticP3PSolver(double v[3][3]) {
	int i, j, k, rotChoice, tripleProdSign, tempInt, 
	  randomRotate, numSolns; 
	double a, b, c, d, A, B, C, D, E, F, G, H, rootOfCubic, 
	  tripleProd, mu0, nu0, alpha1, alpha2, alpha1sq, 
	  alpha2sq, alphaSqDiff, alphaSqSum, eta, beta, singDot, 
	  singLen0, singLen1, mu1, mu2, nu1, nu2, temp,   
	  lambda0, lambda1, lambda2, lambda, error, minerror, 
	  tol1, tol2, tol3, tol4, tol5, tol6, tol7, 
	  vn[3][3], vl[3], sl[3], s[3][3], sn[3][3], cn[3][3], 
	  list[3], pi[3], singularity[2][3], vnr[3][3], snr[3][3], 
	  cnr[3][3], vnd[3][3], dcn[3], a1[3], a2[3], nr[3], n[3], 
	  guess[3][3], tempMat1[3][3], tempMat2[3][3], tempVec[3], 
	  U[4], V[4], X[4], Y[4], Z[4], dp[3], rotMat[3][3], 
	  invRotMat[3][3], defMat[3][3], invDefMat[3][3],
	  diff[3][3], vndp[3], sndp[3], cndp[3], mu00[3], nu00[3], 
	  alpha10[3], alpha20[3], eta0[3], eligible[3], eta0prime[3],
	  direction[3], triangleNormal[3], mu10[4], mu20[4], nu10[4], 
	  nu20[4], errors[4], 

	errorLimit = 1000;
	tol1 = tol2 = tol3 = tol4 = tol5 = tol6 = tol7 = 0;

	// Phase One - Find some basic knowable geometric 
	// quantities (also compute some unknowable ones to 
	// test the accuracy of the various computations 
	// during debugging)

	for(i=0; i<3; i++) {
		vl[i] = LEN(v[i]);
		normalize(v[i], vn[i]);
		for(j=0; j<3; j++) {
			s[i][j] = v[ROTF(i)][j]-v[ROTB(i)][j];
		}
		sl[i] = LEN(s[i]);
		normalize(s[i], sn[i]); 
	}
// reject (not a miss) if attack angle is outside acceptable range  
	addVectors(vn[0], vn[1], tempVec);
	addVectors(tempVec, vn[2], tempVec); 
	normalize(tempVec, direction);
	crossProduct(sn[0], sn[1], triangleNormal);
//	crossProduct(s[0], s[1], triangleNormal);
	normalize(triangleNormal, triangleNormal);
	if (fabs(DOT(direction, triangleNormal)) < 
	  cos(attack_angle_max*atan(1)/45)) return REJECT_TESTING_PARAMETERS;
	if (fabs(DOT(direction, triangleNormal)) > 
	  cos(attack_angle_min*atan(1)/45)) return REJECT_TESTING_PARAMETERS;
// continue for acceptable setups 
	start_time = clock(); 
	for(i=0; i<3; i++) {
		vndp[i] = DOT(vn[ROTF(i)],vn[ROTB(i)]);
		sndp[i] = DOT(sn[ROTF(i)],sn[ROTB(i)]);
		crossProduct(vn[ROTF(i)], vn[ROTB(i)], cn[i]);
//		crossProduct(v[ROTF(i)], v[ROTB(i)], cn[i]);
		normalize(cn[i], cn[i]);
	}
	for(i=0; i<3; i++) {
		cndp[i] = DOT(cn[ROTF(i)],cn[ROTB(i)]);
	}
	crossProduct(vn[0],vn[1],tempVec);
	tripleProd = DOT(tempVec, vn[2]); 
	tripleProdSign = SIGN(tripleProd); 	
#ifdef DEBUG
	printf("v[0] = %lf %lf %lf\n", v[0][0], v[0][1], v[0][2]); 
	printf("v[1] = %lf %lf %lf\n", v[1][0], v[1][1], v[1][2]); 
	printf("v[2] = %lf %lf %lf\n", v[2][0], v[2][1], v[2][2]); 
	printf("distances to control points = %lf %lf %lf\n", 
	  vl[0], vl[1], vl[2]);
	printf("vn[0] = %lf %lf %lf\n", vn[0][0], vn[0][1], vn[0][2]); 
	printf("vn[1] = %lf %lf %lf\n", vn[1][0], vn[1][1], vn[1][2]); 
	printf("vn[2] = %lf %lf %lf\n", vn[2][0], vn[2][1], vn[2][2]);
	printf("vndp = %lf %lf %lf\n", vndp[0], vndp[1], vndp[2]); 
	printf("sidelengths = %lf %lf %lf\n", sl[0], sl[1], sl[2]);
	printf("sn[0] = %lf %lf %lf\n", sn[0][0], sn[0][1], sn[0][2]); 
	printf("sn[1] = %lf %lf %lf\n", sn[1][0], sn[1][1], sn[1][2]); 
	printf("sn[2] = %lf %lf %lf\n", sn[2][0], sn[2][1], sn[2][2]); 
	printf("sndp = %lf %lf %lf\n", sndp[0], sndp[1], sndp[2]); 
	printf("cn[0] = %lf %lf %lf\n", cn[0][0], cn[0][1], cn[0][2]); 
	printf("cn[1] = %lf %lf %lf\n", cn[1][0], cn[1][1], cn[1][2]); 
	printf("cn[2] = %lf %lf %lf\n", cn[2][0], cn[2][1], cn[2][2]); 
	printf("cndp = %lf %lf %lf\n", cndp[0], cndp[1], cndp[2]); 
	printf("triple product sign = %d\n", tripleProdSign); 
#endif

	// Phase Two - Compute some basic parameters for the 
	// algorithm, and decide which view line to rotate to 
	// a vertical position to simplify the equations

	tempInt = -1;
	for(int i=0; i<3; i++) {
		mu00[i] = sqrt((1+cndp[i])/2); 
		nu00[i] = tripleProdSign*sqrt((1-cndp[i])/2);
		if (fabs(1-SQR(sndp[i])) > tol1) {
			alpha10[i] = (sndp[(i+2)%3] - 
			  sndp[i]*sndp[(i+1)%3]) / (1-SQR(sndp[i])); 
			alpha20[i] = (sndp[(i+1)%3] - 
			  sndp[i]*sndp[(i+2)%3]) / (1-SQR(sndp[i])); 
			eta0prime[i] = SQR(2*mu00[i]*nu00[i]*
			  (SQR(alpha10[i])-SQR(alpha20[i]))); 
			eta0[i] = 1 - eta0prime[i];
#ifdef DEBUG
			printf("i mu00 nu00 alpha10 alpha20 eta0 = %d %lf %lf %lf %lf %lf\n", 
				i, mu00[i], nu00[i], alpha10[i], alpha20[i], eta0[i]);
#endif
			if (eta0[i] < tol2 || fabs(mu00[i]) < tol3 || 
			  fabs(nu00[i]) < tol3 || fabs(alpha10[i]) < tol4 ||
			  fabs(alpha20[i]) < tol4
/* ? */
			)
				eligible[i] = FALSE; 
			else {
				eligible[i] = TRUE; 
				tempInt = i; 			
			}
		} else { 
			eligible[i] = FALSE; 
		}
	}
	if (tempInt < 0) {
#ifdef DEBUG
      printf("Error: no eligible rotations"); 
#endif
	  return NO_ELIGIBLE_ROTATION; 
	}	
	rotChoice = tempInt;
/*
	if (eligible[(tempInt+1)%3] && 
	 eta0[(tempInt+1)%3] > eta0[rotChoice])
	  rotChoice = (tempInt+1)%3; 
	if (eligible[(tempInt+2)%3] &&
	 eta0[(tempInt+2)%3] > eta0[rotChoice])
	  rotChoice = (tempInt+2)%3;
*/
	if (eligible[(tempInt+1)%3] && 
	 fabs(eta0prime[(tempInt+1)%3]) < fabs(eta0prime[rotChoice]))
	  rotChoice = (tempInt+1)%3; 
	if (eligible[(tempInt+2)%3] &&
	 fabs(eta0prime[(tempInt+2)%3]) < fabs(eta0prime[rotChoice]))
	  rotChoice = (tempInt+2)%3;
	mu0 = mu00[rotChoice]; 
	nu0 = nu00[rotChoice]; 
	alpha1 = alpha10[rotChoice]; 
	alpha2 = alpha20[rotChoice]; 
	eta = eta0[rotChoice]; 
	alpha1sq = SQR(alpha1); 
	alpha2sq = SQR(alpha2);
	alphaSqDiff = alpha1sq - alpha2sq; 
	alphaSqSum = alpha1sq + alpha2sq; 
	beta = alphaSqSum - SQR(alphaSqDiff); 
#ifdef DEBUG
	printf("rotChoice = %d\n", rotChoice);
	printf("tripleProd = %lf\n", tripleProd);
	printf("tripleProdSign = %d\n", tripleProdSign);
	printf("mu0 nu0 = %lf %lf\n", mu0, nu0);
	printf("alpha1 alpha2 = %lf %lf\n", alpha1, alpha2);
	printf("eta beta = %lf %lf\n", eta, beta);
#endif

    // Phase Three - Compute rotation matrix, its inverse 
    // and the images of some vectors; Also compute the 
    // a few vectors under the deforming transformation 
	// This phase alse computes the two singularities for 
	// the original quartic in X, Y and Z, thereby 
	// obtaining the deformation matrix and its inverse   

	singularity[0][0] = 1 - sqrt(eta); 
	singularity[0][1] = -2*alphaSqDiff*mu0*nu0;
	singularity[0][2] = 0; 
	singularity[1][0] = 1 + sqrt(eta); 
	singularity[1][1] = -2*alphaSqDiff*mu0*nu0;
	singularity[1][2] = 0; 
    singLen0 = LEN(singularity[0]); 
    singLen1 = LEN(singularity[1]); 
    if (fabs(singLen0) < tol5 || fabs(singLen1) < tol5) {
#ifdef DEBUG
      printf("Error: singularity too small"); 
#endif
      return SINGULARITY_TOO_SMALL; 
    }
	normalize(singularity[0], singularity[0]); 
	normalize(singularity[1], singularity[1]); 
    singDot = DOT(singularity[0],singularity[1]);
#ifdef DEBUG
    printf("lengths = %lf %lf\n", singLen0, singLen1); 
    printf("first singularity = %lf %lf %lf\n", singularity[0][0], 
    	singularity[0][1], singularity[0][2]);
    printf("second singularity = %lf %lf %lf\n", singularity[1][0], 
    	singularity[1][1], singularity[1][2]);
    printf("dot prod of these = %lf\n", singDot); 
#endif 
    defMat[0][0] = singularity[0][0]; 
    defMat[1][0] = singularity[0][1]; 
    defMat[2][0] = singularity[0][2]; 
    defMat[0][1] = singularity[1][0]; 
    defMat[1][1] = singularity[1][1]; 
    defMat[2][1] = singularity[1][2]; 
    defMat[0][2] = defMat[1][2] = 0; 
    defMat[2][2] = 1; 
    invertMatrix(defMat, invDefMat); 
    randomRotate = rotateTriple(vn[rotChoice], vn[(rotChoice+1)%3], 
      vn[(rotChoice+2)%3], rotMat);
    invertMatrix(rotMat, invRotMat); 
    multMatrixVector(rotMat, vn[0], vnr[0]); 
    multMatrixVector(rotMat, vn[1], vnr[1]); 
    multMatrixVector(rotMat, vn[2], vnr[2]); 
    multMatrixVector(rotMat, sn[0], snr[0]); 
    multMatrixVector(rotMat, sn[1], snr[1]); 
    multMatrixVector(rotMat, sn[2], snr[2]); 
    multMatrixVector(rotMat, cn[0], cnr[0]); 
    multMatrixVector(rotMat, cn[1], cnr[1]); 
    multMatrixVector(rotMat, cn[2], cnr[2]); 
    multMatrixVector(invDefMat, vnr[0], vnd[0]); 
    multMatrixVector(invDefMat, vnr[1], vnd[1]); 
    multMatrixVector(invDefMat, vnr[2], vnd[2]); 
    crossProduct(vnd[(rotChoice+1)%3], vnd[(rotChoice+2)%3], dcn);
    normalize(dcn, dcn);
#ifdef DEBUG 
    printf("Deformation matrix:\n");
    showMatrix(defMat); 
    printf("Its inverse:\n");
    showMatrix(invDefMat); 
    printf("Rotation matrix:\n");
    showMatrix(rotMat);
    printf("Its matrix:\n");
    showMatrix(invRotMat);    
    printf("Rotated (unit) view vectors:\n");
    showVector(vnr[0]); 
    showVector(vnr[1]); 
    showVector(vnr[2]); 
    printf("Rotated side direction vectors:\n");
    showVector(snr[0]); 
    showVector(snr[1]); 
    showVector(snr[2]); 
    printf("Rotated containment normal vectors:\n");
    showVector(cnr[0]); 
    showVector(cnr[1]); 
    showVector(cnr[2]); 
    printf("Deformed (unit) view vectors:\n");
    showVector(vnd[0]); 
    showVector(vnd[1]); 
    showVector(vnd[2]); 
    printf("Deformed special containment normal vector:\n");
    showVector(dcn);
#endif    

    // Phase Four - Solve system consisting of a quartic 
    // quinomial equation and a linear equation 

    A = 2*sqrt(2)*alphaSqDiff*mu0*eta; 
    B = 2*(eta-2*beta*mu0*mu0 - (1-2*alphaSqSum*mu0*mu0)*sqrt(eta)); 
    C = 2*(eta-2*beta*mu0*mu0 + (1-2*alphaSqSum*mu0*mu0)*sqrt(eta));
    D = -32*SQR(alphaSqDiff)*beta*SQR(mu0*nu0)*SQR(mu0);
    E = (1-SQR(alpha1+alpha2))*(1-SQR(alpha1-alpha2))*mu0*mu0; 
	A /= (singLen0*singLen1); 
	B /= SQR(singLen0); 
	C /= SQR(singLen1); 
	D /= (singLen0*singLen1); 
	F = dcn[0]; 
	G = dcn[1]; 
	H = dcn[2]; 
#ifdef DEBUG
	printf("A B C D E after = %lf %lf %lf %lf %lf\n", A, B, C, D, E);
	printf("F G H = %lf %lf %lf\n", F, G, H);
#endif
    numSolns = solveSystem(A, B, C, D, E, F, G, H, U, V); 
    for(i=0; i<numSolns; i++) {
		tempVec[0] = defMat[0][0]*U[i]+defMat[0][1]*V[i]+defMat[0][2];
		tempVec[1] = defMat[1][0]*U[i]+defMat[1][1]*V[i]+defMat[1][2];
		tempVec[2] = defMat[2][0]*U[i]+defMat[2][1]*V[i]+defMat[2][2];
		normalize(tempVec, tempVec); 
		X[i] = tempVec[0]; Y[i] = tempVec[1]; Z[i] = tempVec[2];
#ifdef DEBUG
		printf("X Y Z = %lf %lf %lf\n", X[i], Y[i], Z[i]);
#endif
	}

    // Phase Five - Process each potential solution, rotating 
    // back to the original coordinates to discover a possible 
    // plane through the origin that is parallel to the control- 
    // points plane; Attempt to translate this plane to see if 
    // it really could be moved to the control-points plane 

	errors[0] = errors[1] = errors[2] = errors[3] = 
		minerror = NO_ACCEPTABLE_ESTIMATES;
	if (numSolns == 0) minerror = NO_SYSTEM_SOLUTIONS;
	else 
	  for(i=0; i<numSolns; i++) {
        if (fabs(Z[i]) > tol6) {
		  mu10[i] = mu1 = ( (SQR(nu0)-SQR(mu0)-1)*SQR(nu0*X[i])
		   + (SQR(mu0)-SQR(nu0)-1)*SQR(mu0*Y[i])
		   + 2*mu0*nu0*X[i]*Y[i]
		   + 2*(1+SQR(alpha1)-SQR(alpha2))*SQR(mu0*nu0) )
		   / (alpha1*SQR(2*mu0*nu0)*Z[i]); 
		  mu20[i] = mu2 = ( (SQR(nu0)-SQR(mu0)-1)*SQR(nu0*X[i])
		   + (SQR(mu0)-SQR(nu0)-1)*SQR(mu0*Y[i])
		   - 2*mu0*nu0*X[i]*Y[i]
		   + 2*(1-SQR(alpha1)+SQR(alpha2))*SQR(mu0*nu0) )
		   / (alpha2*SQR(2*mu0*nu0)*Z[i]); 
    	  nu10[i] = nu1 = (nu0*X[i]-mu0*Y[i])/(2*alpha1*mu0*nu0); 
    	  nu20[i] = nu2 = (nu0*X[i]+mu0*Y[i])/(2*alpha2*mu0*nu0); 
		  a1[0] = mu0*nu1; a1[1] = -nu0*nu1; a1[2] = mu1; 
		  a2[0] = mu0*nu2; a2[1] =  nu0*nu2; a2[2] = mu2; 
		  crossProduct(a1, a2, nr);
		  normalize(nr, nr); 
		  multMatrixVector(invRotMat, nr, n); 
		  /*normalize(n, n);*/ 
		  dp[0] = DOT(vn[0],n); 
		  dp[1] = DOT(vn[1],n); 
		  dp[2] = DOT(vn[2],n);
#ifdef DEBUG
		  printf("mu1 mu2 nu1 nu2 = %lf %lf %lf %lf\n", 
			mu1, mu2, nu1, nu2); 
		  printf("Vectors a1, a2, nr, n:\n");
		  showVector(a1); 
		  showVector(a2); 
		  showVector(nr); 
		  showVector(n); 
		  printf("dp's = %lf %lf %lf\n", dp[0], dp[1], dp[2]);
#endif
		  if (
		    SQR(dp[1])+SQR(dp[2])-2*dp[1]*dp[2]*vndp[0] > tol7 && 
		    SQR(dp[2])+SQR(dp[0])-2*dp[2]*dp[0]*vndp[1] > tol7 && 
		    SQR(dp[0])+SQR(dp[1])-2*dp[0]*dp[1]*vndp[2] > tol7 ) {
			  lambda0 = sl[0]*dp[1]*dp[2] / 
			    sqrt(SQR(dp[1])+SQR(dp[2])-2*dp[1]*dp[2]*vndp[0]); 
			  lambda1 = sl[1]*dp[2]*dp[0] / 
			    sqrt(SQR(dp[2])+SQR(dp[0])-2*dp[2]*dp[0]*vndp[1]); 
			  lambda2 = sl[2]*dp[0]*dp[1] / 
			    sqrt(SQR(dp[0])+SQR(dp[1])-2*dp[0]*dp[1]*vndp[2]); 
			  lambda = (lambda0 + lambda1 + lambda2) / 3;
			  for(j=0; j<3; j++)
			    for(k=0; k<3; k++)
				  guess[j][k] = tripleProdSign*lambda*vn[j][k]/dp[j]; 
			  // coordinate differences 
			  subMatrices(guess, v, diff);
			  // want relative error so divide by actual distances  
			  diff[0][0] /= vl[0]; diff[0][1] /= vl[0]; diff[0][2] /= vl[0]; 
			  diff[1][0] /= vl[1]; diff[1][1] /= vl[1]; diff[1][2] /= vl[1]; 
			  diff[2][0] /= vl[2]; diff[2][1] /= vl[2]; diff[2][2] /= vl[2]; 
			  transposeMatrix(diff, tempMat1);
			  multMatrices(diff, tempMat1, tempMat2);
			  errors[i] = error = 
			    sqrt(tempMat2[0][0]+tempMat2[1][1]+tempMat2[2][2]);
			  if (error < errorLimit && 
			   (minerror < 0 || error < minerror)) 
			  	 minerror = error; 
#ifdef DEBUG
		  printf("lambdas = %lf %lf %lf\n", lambda0, lambda1, lambda2);
		  printf("lambda = %lf\n", lambda);
		  printf("Candidate solution:\n");
		  showMatrix(guess); 
		  printf("Difference with actual control points:\n");
		  showMatrix(diff); 
		  printf("error = %lf\n", error); 
#endif
		  }
	  } 
    }
#ifdef DEBUG
    if (minerror < -1)
      printf("Error: no solutions to algebraic system"); 
    else if (minerror < 0)
      printf("Error: no acceptable estimates"); 
#endif
#ifdef DEBUG2
	if (minerror < 0 || minerror > DEBUG2_CUTOFF) {
    	printf("minimum error = %lf\n", minerror);
    	printf("numSolns = %d\n", numSolns);
    	printf("errors = %lg %lg %lg %lg\n", 
    	  errors[0], errors[1], errors[2], errors[3]);
    	printf("distance = %lf\n", 
    	  (LEN(v[0])+LEN(v[1])+LEN(v[2]))/3); 
    	printf("size = %lf\n", 
    	  (LEN(s[0])+LEN(s[1])+LEN(s[2]))/3); 
		printf("rotation choice = %d\n", rotChoice);
		printf("eta = %lf\n", eta);
    	printf("dot prod of singularities = %lf\n", singDot); 
    	printf("square of that = %lf\n", SQR(singDot)); 
    	printf("eta' list: %lf %lf %lf\n", 
    	  eta0prime[0], eta0prime[1], eta0prime[2]);
    	printf("mu0's: %lf %lf %lf\n", mu00[0], mu00[1], mu00[2]);
    	printf("nu0's: %lf %lf %lf\n", nu00[0], nu00[1], nu00[2]);
    	printf("alpha1's: %lf %lf %lf\n", alpha10[0], alpha10[1], alpha10[2]);
    	printf("alpha2's: %lf %lf %lf\n", alpha20[0], alpha20[1], alpha20[2]);
    	printf("randomRotate = %d\n", randomRotate); 
    	printf("Z's = %lf %lf %lf %lf\n", Z[0], Z[1], Z[2], Z[3]);
    	printf("mu1's = %lf %lf %lf %lf\n", mu10[0],mu10[1],mu10[2],mu10[3]);
    	printf("mu2's = %lf %lf %lf %lf\n", mu20[0],mu20[1],mu20[2],mu20[3]);
    	printf("nu1's = %lf %lf %lf %lf\n", nu10[0],nu10[1],nu10[2],nu10[3]);
    	printf("nu2's = %lf %lf %lf %lf\n\n", nu20[0],nu20[1],nu20[2],nu20[3]);
	}
#endif
	end_time = clock();
	exec_time = end_time - start_time; 
	return minerror; 
}

// Helpful interface to lambdaTwistSolver from Python

double solveLT(double v00, double v01, double v02, double v10, 
  double v11, double v12, double v20, double v21, double v22) {
	double v[3][3]; 
	v[0][0] = v00; v[0][1] = v01; v[0][2] = v02; 
	v[1][0] = v10; v[1][1] = v11; v[1][2] = v12; 
	v[2][0] = v20; v[2][1] = v21; v[2][2] = v22; 
	return lambdaTwistSolver(v);
}

// Helpful interface to ellipticP3Psolver from Python

double solveEP(double v00, double v01, double v02, double v10, 
  double v11, double v12, double v20, double v21, double v22) {
	double v[3][3]; 
	v[0][0] = v00; v[0][1] = v01; v[0][2] = v02; 
	v[1][0] = v10; v[1][1] = v11; v[1][2] = v12; 
	v[2][0] = v20; v[2][1] = v21; v[2][2] = v22; 
	return ellipticP3PSolver(v);
}

// Test one of the two methods for solving the P3P problem 
int main() {

	int i, count, hits = 0, misses = 0, index, 
	  failure[5], errorCounts[17], numErrorRanges;
	double controlPts[3][3], error, sum = 0, sum2 = 0, 
	  timeTotal = 0, avgerror, minerror = -1, maxerror = -1, 
	  stddev, errorRanges[17];

	srand(time(0));
	for(i=0; i<5; i++) failure[i] = 0; 
	for(i=0; i<17; i++) errorCounts[i] = 0; 
	errorRanges[1]  = 1; 
	errorRanges[2]  = .1; 
	errorRanges[3]  = .01; 
	errorRanges[4]  = .001; 
	errorRanges[5]  = .0001; 
	errorRanges[6]  = .00001; 
	errorRanges[7]  = .000001; 
	errorRanges[8]  = .0000001; 
	errorRanges[9]  = .00000001; 
	errorRanges[10] = .000000001; 
	errorRanges[11] = .0000000001; 
	errorRanges[12] = .00000000001; 
	errorRanges[13] = .000000000001; 
	errorRanges[14] = .0000000000001;
	errorRanges[15] = .00000000000001;
	errorRanges[16] = .000000000000001;
	numErrorRanges = 17;

	for(count = 0; count < NUMBER_TRIALS; count++) {

		controlPts[0][0] = 1500. + randomReal(0,1); 
		controlPts[0][1] = 700.  + randomReal(0,1); 
		controlPts[0][2] = 800.  + randomReal(0,1);
		controlPts[1][0] = 1500. + randomReal(0,1);
		controlPts[1][1] = 700.  + randomReal(0,1);
		controlPts[1][2] = 800.  + randomReal(0,1);
		controlPts[2][0] = 1500. + randomReal(0,1);
		controlPts[2][1] = 700.  + randomReal(0,1);
		controlPts[2][2] = 800.  + randomReal(0,1);

		controlPts[0][0] = 20. + randomReal(0,1); 
		controlPts[0][1] = 30.  + randomReal(0,1); 
		controlPts[0][2] = 00.  + randomReal(0,1);
		controlPts[1][0] = 20. + randomReal(0,1);
		controlPts[1][1] = 30.  + randomReal(0,1);
		controlPts[1][2] = 00.  + randomReal(0,1);
		controlPts[2][0] = 20. + randomReal(0,1);
		controlPts[2][1] = 30.  + randomReal(0,1);
		controlPts[2][2] = 00.  + randomReal(0,1);

#ifdef TEST_LAMBDA_TWIST
		error = lambdaTwistSolver(controlPts);
#else 
		error = ellipticP3PSolver(controlPts);
#endif 
		if (error >= 0) {
			sum += error; 
			sum2 += error*error; 
			if (minerror < 0 || error < minerror) 
			 minerror = error; 
			if (error > maxerror) maxerror = error; 
			for (i=numErrorRanges-1; i>0; i--) {
				if (error < errorRanges[i]) {
					errorCounts[i]++; 
					break; 
				}
			}
			if (i==0) errorCounts[0]++; 
			hits++; 
			/*printf("error = %lf\n", error);*/
		} else {
			index = (int)(0.5 - error); 
			if (index > 3) index = 0; 
			failure[index]++; 
			if (index > 0) misses++;
			/*printf("failure code = %d\n", index);*/
		}
		timeTotal += getTime();
	}
	avgerror = sum/hits;
	stddev = sqrt(sum2/hits - avgerror*avgerror);
	printf("average error = %lg\n", avgerror);
	printf("standard deviation = %lg\n", stddev);
	printf("minimum error = %lg\n", minerror);
	printf("maximum error = %lg\n", maxerror);
	printf("number of hits = %d\n", hits);
	printf("number of misses %d\n", misses);
	printf("success rate = %f%%\n", 100.*hits/(hits+misses));
	printf("number rejected test cases (outside testing parameters) = %d\n", 
		failure[0]); 
	printf("number failures due to no acceptable estimates = %d\n", 
		failure[1]); 
	printf("number failures due to no algebraic system solutions = %d\n", 
		failure[2]); 
#ifndef TEST_LAMBDA_TWIST	
	printf("number failures due to small singularity = %d\n", 
		failure[3]); 
	printf("number failures due to no eligible rotations = %d\n", 
		failure[4]); 
#endif 
	for (i=0; i<numErrorRanges; i++) 
		printf("Count for range #%2d: %d\n", i, errorCounts[i]); 
	printf("Execution time: %lf\n", timeTotal); 
}
