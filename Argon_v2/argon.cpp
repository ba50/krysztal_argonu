#include <fstream>
#include <string>

#include "argon.h"

void main(int argc, char *argv[]){

	int i, j, k, licznik_xyz, licznik_out, licznik_p;
	double time, r_ij;

	FILE *XYZ_file, *p_file, *OUT_file;

	std::ifstream IN_file;

	std::string buffer;
	std::string IN_filename = argv[1];
	std::string OUT_filename = argv[2];
	std::string XYZ_filename = argv[3];

	parameters param;
	state stat;
	Vector sumaP;

	if (OUT_filename == "-"){
		OUT_filename = "out.dat";
	}

	/* Wczytywanie danych */

	IN_file.open(IN_filename);

	IN_file >> param.n;
	std::getline(IN_file, buffer);

	IN_file >> param.m;
	std::getline(IN_file, buffer);

	IN_file >> param.e;
	std::getline(IN_file, buffer);

	IN_file >> param.R;
	std::getline(IN_file, buffer);

	IN_file >> param.f;
	std::getline(IN_file, buffer);

	IN_file >> param.L;
	std::getline(IN_file, buffer);

	IN_file >> param.a;
	std::getline(IN_file, buffer);

	IN_file >> stat.T;
	std::getline(IN_file, buffer);

	IN_file >> param.tau;
	std::getline(IN_file, buffer);

	IN_file >> param.So;
	std::getline(IN_file, buffer);

	IN_file >> param.Sd;
	std::getline(IN_file, buffer);

	IN_file >> param.Sout;
	std::getline(IN_file, buffer);

	IN_file >> param.Sxyz;
	std::getline(IN_file, buffer);

	IN_file.close();

	param.L = param.L * param.a*(param.n - 1);

	/* Okreslanie stanu poczatkowego ukladu */

	srand(time_t(NULL));

	param.N = (int) pow(param.n, 3);

	stat.r = new Vector[param.N];
	stat.p = new Vector[param.N];
	stat.F_p = new Vector[param.N];
	stat.F_s = new Vector[param.N];

	for (i = 0; i < param.n; i++){
		for (j = 0; j < param.n; j++){
			for (k = 0; k < param.n; k++){
				stat.r[i + j*param.n + k*param.n*param.n].x = (i - (param.n - 1) / 2.f)*param.a + (j - (param.n - 1) / 2) * param.a / 2.f + (k - (param.n - 1) / 2)*param.a / 2.f;
				stat.r[i + j*param.n + k*param.n*param.n].y = (j - (param.n - 1) / 2)*param.a*sqrt(3.f) / 2 + (k - (param.n - 1) / 2)*param.a*sqrt(3.f) / 6;
				stat.r[i + j*param.n + k*param.n*param.n].z = (k - (param.n - 1) / 2)*param.a * sqrt(2.f / 3.f);
			}
		}
	}

	for (i = 0; i < param.N; i++){
		stat.r[i].r = sqrt(pow(stat.r[i].x, 2) + pow(stat.r[i].y, 2) + pow(stat.r[i].z, 2));
	}

	fopen_s(&XYZ_file, XYZ_filename.c_str(), "w");

	fprintf(XYZ_file, "%d\n", param.N);
	fprintf(XYZ_file, "STEP:\t%d\n", param.N);
	for (i = 0; i < param.N; i++){
		fprintf(XYZ_file, "Ar\t%f\t%f\t%f\n", stat.r[i].x, stat.r[i].y, stat.r[i].z);
	}
	fprintf(XYZ_file, "\n", param.N);

	for (i = 0; i < param.N; i++){
		stat.p[i].x = znak()*sqrt(2*param.m*(-K_B*stat.T*log(r0_1()) / 2));
		stat.p[i].y = znak()*sqrt(2*param.m*(-K_B*stat.T*log(r0_1()) / 2));
		stat.p[i].z = znak()*sqrt(2*param.m*(-K_B*stat.T*log(r0_1()) / 2));

		stat.p[i].r = sqrt(pow(stat.p[i].x, 2) + pow(stat.p[i].y, 2) + pow(stat.p[i].z, 2));
	}

	sumaP.x = 0;
	sumaP.y = 0;
	sumaP.z = 0;
	

	for (i = 0; i < param.N; i++){
		sumaP.x += stat.p[i].x;
		sumaP.y += stat.p[i].y;
		sumaP.z += stat.p[i].z;
	}

	fopen_s(&p_file, "p.dat", "w");

	for (i = 0; i < param.N; i++){
		stat.p[i].x -= sumaP.x / (double) param.N;
		stat.p[i].y -= sumaP.y / (double) param.N;
		stat.p[i].z -= sumaP.z / (double) param.N;
		fprintf(p_file, "%f\t%f\t%f\n", stat.p[i].x, stat.p[i].y, stat.p[i].z);
	}

	fclose(p_file);

	/* Obliczanie potencjalu, sil oraz cisnienia w stanie poczatkowym */

	for (i = 0; i < param.N; i++){
		stat.F_s[i].x = 0;
		stat.F_s[i].y = 0;
		stat.F_s[i].z = 0;

		stat.F_p[i].x = 0;
		stat.F_p[i].y = 0;
		stat.F_p[i].z = 0;
	}

	stat.V_p = 0;
	stat.V_s = 0;

	for (i = 0; i < param.N; i++){

		if (stat.r[i].r >= param.L){

			stat.V_s += param.f*pow(stat.r[i].r - param.L, 2) / 2;

			stat.F_s[i].x += param.f*(param.L - stat.r[i].r)* stat.r[i].x / stat.r[i].r;
			stat.F_s[i].y += param.f*(param.L - stat.r[i].r)* stat.r[i].y / stat.r[i].r;
			stat.F_s[i].z += param.f*(param.L - stat.r[i].r)* stat.r[i].z / stat.r[i].r;
		}

		for (j = i; j < param.N; j++){
			if (i != j){
				r_ij = sqrt(pow(stat.r[i].x - stat.r[j].x, 2) + pow(stat.r[i].y - stat.r[j].y, 2) + pow(stat.r[i].z - stat.r[j].z, 2));

				stat.V_p += param.e*(pow(param.R / r_ij, 12) - 2 * pow(param.R / r_ij, 6));

				stat.F_p[i].x += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].x - stat.r[j].x) / pow(r_ij, 2);
				stat.F_p[j].x -= stat.F_p[i].x;

				stat.F_p[i].y += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].y - stat.r[j].y) / pow(r_ij, 2);
				stat.F_p[j].y -= stat.F_p[i].y;

				stat.F_p[i].z += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].z - stat.r[j].z) / pow(r_ij, 2);
				stat.F_p[j].z -= stat.F_p[i].z;
			}

		}
	}


	fopen_s(&OUT_file, OUT_filename.c_str(), "w");

	fprintf(OUT_file, "TIME\tH\tV\tT\tP\n");

	licznik_out = 0;
	licznik_xyz = 0;
	licznik_p = 0;

	stat.T = 0;
	stat.P = 0;
	stat.H = 0;

	stat.T_mean = 0;
	stat.P_mean = 0;
	stat.H_mean = 0;

	time = 0;
	for (int s = 0; s < (param.So + param.Sd); s++){

		for (i = 0; i < param.N; i++){

			stat.p[i].x += (stat.F_p[i].x + stat.F_s[i].x)*param.tau / 2.0;
			stat.p[i].y += (stat.F_p[i].y + stat.F_s[i].y)*param.tau / 2.0;
			stat.p[i].z += (stat.F_p[i].z + stat.F_s[i].z)*param.tau / 2.0;

			stat.r[i].x += stat.p[i].x*param.tau / param.m;
			stat.r[i].y += stat.p[i].y*param.tau / param.m;
			stat.r[i].z += stat.p[i].z*param.tau / param.m;
			stat.r[i].r = sqrt(pow(stat.r[i].x, 2) + pow(stat.r[i].y, 2) + pow(stat.r[i].z, 2));

			stat.F_s[i].x = 0;
			stat.F_s[i].y = 0;
			stat.F_s[i].z = 0;

			stat.F_p[i].x = 0;
			stat.F_p[i].y = 0;
			stat.F_p[i].z = 0;
		}

		stat.P = 0;
		stat.V_s = 0;
		stat.V_p = 0;
		stat.E_kin = 0;
		
		for (i = 0; i < param.N; i++){

			if (stat.r[i].r >= param.L){

				stat.V_s += param.f*pow(stat.r[i].r - param.L, 2) / 2;

				stat.F_s[i].x += param.f*(param.L - stat.r[i].r)* stat.r[i].x / stat.r[i].r;
				stat.F_s[i].y += param.f*(param.L - stat.r[i].r)* stat.r[i].y / stat.r[i].r;
				stat.F_s[i].z += param.f*(param.L - stat.r[i].r)* stat.r[i].z / stat.r[i].r;

				stat.F_s[i].r = sqrt(pow(stat.F_s[i].x, 2) + pow(stat.F_s[i].y, 2) + pow(stat.F_s[i].z, 2));
			}

			stat.P += stat.F_s[i].r / (4 * PI*pow(param.L, 2));

			for (j = i; j < param.N; j++){

				if (i != j){

					r_ij = sqrt(pow(stat.r[i].x - stat.r[j].x, 2) + pow(stat.r[i].y - stat.r[j].y, 2) + pow(stat.r[i].z - stat.r[j].z, 2));

					stat.V_p += param.e*(pow(param.R / r_ij, 12) - 2 * pow(param.R / r_ij, 6));

					stat.F_p[i].x += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].x - stat.r[j].x) / pow(r_ij, 2);
					stat.F_p[j].x -= stat.F_p[i].x;

					stat.F_p[i].y += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].y - stat.r[j].y) / pow(r_ij, 2);
					stat.F_p[j].y -= stat.F_p[i].y;

					stat.F_p[i].z += 12 * param.e*(pow(param.R / r_ij, 12) - pow(param.R / r_ij, 6))*(stat.r[i].z - stat.r[j].z) / pow(r_ij, 2);
					stat.F_p[j].z -= stat.F_p[i].z;
				}
			}

		}

		for (i = 0; i < param.N; i++){
			stat.p[i].x += (stat.F_p[i].x + stat.F_s[i].x)*param.tau / 2.0;
			stat.p[i].y += (stat.F_p[i].y + stat.F_s[i].y)*param.tau / 2.0;
			stat.p[i].z += (stat.F_p[i].z + stat.F_s[i].z)*param.tau / 2.0;

			stat.p[i].r = sqrt(pow(stat.p[i].x, 2) + pow(stat.p[i].y, 2) + pow(stat.p[i].z, 2));

			stat.E_kin += pow(stat.p[i].r, 2) / (2.0 * param.m);
		}

		stat.T = 2 * stat.E_kin / (3 * param.N*K_B);

		stat.H = stat.E_kin + stat.V_p + stat.V_s;

		if (s / param.Sout == licznik_out){

			fprintf(OUT_file, "%e\t%e\t%e\t%e\t%e\n", time, stat.H, stat.V_s+stat.V_p, stat.T, stat.P);

			licznik_out++;
		}

		if (s / param.Sxyz == licznik_xyz){
			fprintf(XYZ_file, "%d\n", param.N);
			fprintf(XYZ_file, "STEP:\t%d\n", s + 2);
			for (i = 0; i < param.N; i++){
				fprintf(XYZ_file, "Ar\t%f\t%f\t%f\n", stat.r[i].x, stat.r[i].y, stat.r[i].z);
			}
			fprintf(XYZ_file, "\n");
			licznik_xyz++;
		}

		if ((param.Sd + param.So)*licznik_p / 100 == s){
			printf("%d%%\r", licznik_p);
			licznik_p++;
		}

		if (s >= param.So){
			stat.T_mean += stat.T / param.Sd;
			stat.P_mean += stat.P / param.Sd;
			stat.H_mean += stat.H / param.Sd;
		}

		time += param.tau;
	}
	printf("100%%\n\n");

	printf("T_mean = %e\n", stat.T_mean);
	printf("P_mean = %e\n", stat.P_mean);
	printf("H_mean = %e\n", stat.H_mean);

	fclose(XYZ_file);
	fclose(OUT_file);

	system("pause");
}