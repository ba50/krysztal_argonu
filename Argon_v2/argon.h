#define K_B 8.13e-3f
#define PI 3.1416f

struct Vector{
	double x, y, z, r;
};

struct parameters{
	int N;
	int n;
	double m;
	double e;
	double R;
	double f;
	double L;
	double a;
	double tau;
	int So;
	int Sd;
	int Sout;
	int Sxyz;
};

struct state{
	double V_p;
	double V_s;
	double P;
	double P_mean;
	double H;
	double H_mean;
	double T;
	double T_mean;
	double E_kin;
	Vector *r;
	Vector *p;
	Vector *F_p;
	Vector *F_s;
};

double znak(){
	double znak;
	double tmp = (double) rand() / RAND_MAX;

	if (tmp < 0.5){
		znak = 1;
	}
	else{
		znak = -1;
	}

	return znak;
}

double r0_1(){
	return (double)rand() / RAND_MAX;
}