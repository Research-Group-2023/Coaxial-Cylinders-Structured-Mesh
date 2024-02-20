#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <windows.h>

using namespace std;

#include "Variables.h"

void Arrays_Creation(int N, int M)
{

	r_cell = new double[N + 2];
	eps_cell = new double[M + 2];

	r_m = new double[150 * N];
	eps_m = new double[150 * N];

	dr = new double[N + 2];
	deps = new double[M + 2];

	V = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) V[i] = new double[M + 2];
	U = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) U[i] = new double[M + 2];

	v = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) v[i] = new double[M + 2];
	u = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) u[i] = new double[M + 2];

	ap_r = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) ap_r[i] = new double[M + 2];
	ap_eps = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) ap_eps[i] = new double[M + 2];

	PP = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) PP[i] = new double[M + 2];

	P = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) P[i] = new double[M + 2];

	psi = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) psi[i] = new double[M + 2];

	visc = new double* [N + 2]; for (int i = 0; i <= N + 1; i++) v[i] = new double[M + 2];

	sign = new int* [N + 2]; for (int i = 0; i <= N + 1; i++) sign[i] = new int[M + 2];

}

void Arrays_Remove()
{

	delete[] r_cell;
	delete[] eps_cell;

	delete[] r_m;
	delete[] eps_m;

	delete[] dr;
	delete[] deps;

	for (int i = 0; i <= N + 1; i++) delete[] V[i]; delete[] V;
	for (int i = 0; i <= N + 1; i++) delete[] U[i]; delete[] U;

	for (int i = 0; i <= N + 1; i++) delete[] v[i]; delete[] v;
	for (int i = 0; i <= N + 1; i++) delete[] u[i]; delete[] u;

	for (int i = 0; i <= N + 1; i++) delete[] ap_r[i]; delete[] ap_r;
	for (int i = 0; i <= N + 1; i++) delete[] ap_eps[i]; delete[] ap_eps;

	for (int i = 0; i <= N + 1; i++) delete[] PP[i]; delete[] PP;

	for (int i = 0; i <= N + 1; i++) delete[] P[i]; delete[] P;

	for (int i = 0; i <= N + 1; i++) delete[] visc[i]; delete[] visc;

	for (int i = 0; i <= N + 1; i++) delete[] psi[i]; delete[] psi;

	for (int i = 0; i <= N + 1; i++) delete[] sign[i]; delete[] sign;

}

void Redistricting()
{

	for (int i = 0; i <= N + 1; i++)
	{
		for (int j = 0; j <= M + 1; j++)
		{

			v[i][j] = V[i][j];
			u[i][j] = U[i][j];

		}
	}

}

void Read_From_Save()
{

	ifstream file_Save("Save/Save_80/Save_Konf_2/Re = 1/Mesh.DAT");
	/* Считывание начальных условий и сетки*/
	{

		file_Save >> _time;

		file_Save >> N >> M;

		file_Save >> R0 >> R1 >> alfa >> omega_0 >> omega_1 >> Re >> Pi;

		Arrays_Creation(N, M);

		for (int i = 0; i <= N + 1; i++)
		{

			file_Save >> i >> dr[i] >> r_cell[i];

		}

		for (int j = 0; j <= M + 1; j++)
		{

			file_Save >> j >> deps[j] >> eps_cell[j];

		}

		for (int i = 0; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{

				file_Save >> sign[i][j];

			}

		}

		file_Save.close();

	}

	ifstream file_Save1("Save/Save_80/Save_Konf_2/Re = 1/V_U_P.DAT");
	/* Считывание V, U, P */
	{

		for (int i = 0; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{

				file_Save1 >> V[i][j] >> v[i][j] >> U[i][j] >> u[i][j] >> P[i][j];

				i = i;

			}

		}

		file_Save1.close();

	}

}

void Topology(int N, int M)
{

	for (int i = 0; i <= N + 1; i++)
	{
		for (int j = 0; j <= M + 1; j++)
		{

			sign[i][j] = Computational_Cell;

			if (i == 0) sign[i][j] = Inner_Wall;

			if (i == N + 1) sign[i][j] = Outer_Wall;

			if (j == 0 || j == M + 1) sign[i][j] = Gluing_Region;



		}


	}

	if (Obstacle == true)
	{

		int i_1 = round(N * r_1);
		int i_2 = round(N * r_2);
		int j_1 = round(M * eps_1);
		int j_2 = round(M * eps_2);

		for (int i = 0; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{


				if (i >= i_1 && i <= i_2 && j >= j_1 && j <= j_2 /*|| i >= i_1 && i <= i_2 && j >= 58 && j <= 60 || i >= i_1 && i <= i_2 && j >= 37 && j <= 39 || i >= i_1 && i <= i_2 && j >= 77 && j <= 79*/ )
					sign[i][j] = Inner_Wall;

			}
		}

	}

}

void Initial_Conditions()
{
	if (Start_From_Save == false)
	{
		Arrays_Creation(N, M);

		Topology(N, M);

		/* Построение регулярной сетки */
		if (Regular_Mesh == true)
		{

			for (int i = 0; i <= N + 1; i++)
			{

				dr[i] = dr_1;

			}

			for (int j = 0; j <= M + 1; j++)
			{

				deps[j] = deps_1;

			}

			r_cell[0] = R0 - 0.5 * dr[0];

			for (int i = 1; i <= N + 1; i++)
			{

				r_cell[i] = r_cell[i - 1] + 0.5 * (dr[i - 1] + dr[i]);

			}

			eps_cell[0] = 0.5 * deps[0];

			for (int j = 1; j <= M + 1; j++)
			{

				eps_cell[j] = eps_cell[j - 1] + 0.5 * (deps[j - 1] + deps[j]);

			}

		}

		/* Построение нерегулярной сетки */
		if (Regular_Mesh == false)
		{

			for (int j = 0; j <= M + 1; j++)
			{

				deps[j] = deps_1;

			}

			dr[0] = R0 * deps_1 / (1. + deps_1);

			r_cell[0] = R0 - 0.5 * dr[0];

			dr[1] = deps_1 * R0;

			r_cell[1] = R0 + 0.5 * dr[1];

			for (int i = 2; i <= N + 1; i++)
			{
				dr[i] = deps_1 * (r_cell[i - 1] + 0.5 * dr[i - 1]);

				r_cell[i] = r_cell[i - 1] + 0.5 * (dr[i] + dr[i - 1]);

			}

			eps_cell[0] = 0.5 * deps[0];

			for (int j = 1; j <= M + 1; j++)
			{

				eps_cell[j] = eps_cell[j - 1] + 0.5 * (deps[j - 1] + deps[j]);

			}

		}

		/* Установившиеся начальные условия */
		if (Steady_Cond == true) for (int i = 0; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{

				v[i][j] = 0.0;
				V[i][j] = 0.0;

				u[i][j] = ((omega_1 * R1 * R1 - omega_0 * R0 * R0) * r_cell[i] + R0 * R0 * R1 * R1 * (omega_0 - omega_1) / r_cell[i]) / (R1 * R1 - R0 * R0);
				U[i][j] = ((omega_1 * R1 * R1 - omega_0 * R0 * R0) * r_cell[i] + R0 * R0 * R1 * R1 * (omega_0 - omega_1) / r_cell[i]) / (R1 * R1 - R0 * R0);

				P[i][j] = Re / pow((R1 * R1 - R0 * R0), 2) * (pow((omega_1 * R1 * R1 - omega_0 * R0 * R0), 2) * (r_cell[i] * r_cell[i] - R0 * R0) * 0.5 + 2 * R0 * R0 * R1 * R1 * (omega_0 - omega_1) * (omega_1 * R1 * R1 - omega_0 * R0 * R0) * log(r_cell[i] / R0) - 0.5 * (omega_0 - omega_1) * (omega_0 - omega_1) * pow(R0 * R1, 4) * (1.0 / (r_cell[i] * r_cell[i]) - 1.0 / (R0 * R0)));
				visc[i][j] = 1;
			}

		}

		/* Произвольные начальные условия */
		if (Steady_Cond == false) for (int i = 0; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{

				v[i][j] = 0.0;
				V[i][j] = 0.0;

				u[i][j] = 0.0;
				U[i][j] = 0.0;

				P[i][j] = 0.0;
				visc[i][j] = 1;
			}

		}

	}
	else Read_From_Save();

	for (int i = 0; i <= N + 1; i++)
	{
		for (int j = 0; j <= M + 1; j++)
		{

			PP[i][j] = 0.0;

			ap_r[i][j] = 0.0;
			ap_eps[i][j] = 0.0;

		}

	}

}

void Boundary_conditions()
{

	for (int i = 0; i <= N + 1; i++)
	{

		U[i][0] = U[i][M];
		U[i][M + 1] = U[i][1];

		ap_eps[i][0] = ap_eps[i][M];
		ap_eps[i][M + 1] = ap_eps[i][1];

		P[i][0] = P[i][M];
		P[i][M + 1] = P[i][1];

		V[i][0] = V[i][M];
		V[i][M + 1] = V[i][1];

		ap_r[i][0] = ap_r[i][M];
		ap_r[i][M + 1] = ap_r[i][1];

	}

	for (int j = 0; j <= M + 1; j++)
	{

		U[0][j] = omega_0 * R0 - dr[0] / dr[1] * (U[1][j] - omega_0 * R0);
		U[N + 1][j] = omega_1 * R1 - dr[N + 1] / dr[N] * (U[N][j] - omega_1 * R1);

		P[0][j] = P[1][j];
		P[N + 1][j] = P[N][j];

	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M + 1; j++)
		{

			if (sign[i][j] == Computational_Cell)
			{

				if (sign[i + 1][j] == Outer_Wall)
				{

					U[i + 1][j] = omega_1 * R1 - dr[i + 1] / dr[i] * (U[i][j] - omega_1 * R1);
					P[i + 1][j] = P[i][j];

				}

				if (sign[i - 1][j] == Inner_Wall && sign[i - 1][j + 1] == Inner_Wall)
				{

					U[i - 1][j] = omega_0 * R0 - dr[i - 1] / dr[i] * (U[i][j] - omega_0 * R0);
					P[i - 1][j] = P[i][j];

				}

				if (sign[i][j + 1] == Inner_Wall && sign[i + 1][j + 1] == Inner_Wall)
				{

					V[i][j + 1] = -deps[j + 1] / deps[j] * V[i][j];
					P[i][j + 1] = P[i][j];

				}

				if (sign[i + 1][j - 1] == Inner_Wall && sign[i][j - 1] == Inner_Wall)
				{

					V[i][j - 1] = -deps[j - 1] / deps[j] * V[i][j];
					P[i][j - 1] = P[i][j];

				}

			}

		}

	}

}

void Calculation_Velocity_Vr()
{

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			double rp, rn, rs, re, rw;
			double drp, drn, drs, dre, drw;
			double depsp, depsn, depss, depse, depsw;

			double visc_s, visc_n, visc_w, visc_e, visc_p;

			double Vp, Vn, Vs;
			double Up, Ue, Uw, Un, Us;

			double VN, VS, VE, VW;

			if (sign[i][j] == Computational_Cell && sign[i + 1][j] == Computational_Cell)
			{

				/* Присвоение значений узлам КО */
				{

					rp = r_cell[i] + 0.5 * dr[i];
					rn = r_cell[i + 1];
					rs = r_cell[i];
					re = r_cell[i] + 0.5 * dr[i];
					rw = r_cell[i] + 0.5 * dr[i];

					visc_s = visc[i][j];
					visc_n = visc[i + 1][j];
					visc_e = 0.5 * (0.5 * (visc[i][j] + visc[i + 1][j]) + 0.5 * (visc[i][j + 1] + visc[i + 1][j + 1]));
					visc_w = 0.5 * (0.5 * (visc[i][j] + visc[i + 1][j]) + 0.5 * (visc[i][j - 1] + visc[i + 1][j - 1]));
					visc_p = 0.5 * (visc_s + visc_n);

					drp = 0.5 * (dr[i] + dr[i + 1]);
					drn = dr[i + 1];
					drs = dr[i];
					dre = drp;
					drw = drp;

					depsp = deps[j];
					depsn = deps[j];
					depss = deps[j];
					depse = 0.5 * (deps[j] + deps[j + 1]);
					depsw = 0.5 * (deps[j] + deps[j - 1]);

					VN = V[i + 1][j];
					VS = V[i - 1][j];
					VE = V[i][j + 1];
					VW = V[i][j - 1];

					Vp = V[i][j];
					Vn = 0.5 * (V[i][j] + V[i + 1][j]);
					Vs = 0.5 * (V[i][j] + V[i - 1][j]);

					Ue = U[i][j] + 0.5 * dr[i] * (U[i + 1][j] - U[i][j]) / dre;
					Uw = U[i][j - 1] + 0.5 * dr[i] * (U[i + 1][j - 1] - U[i][j - 1]) / drw;
					Un = 0.5 * (U[i + 1][j] + U[i + 1][j - 1]);
					Us = 0.5 * (U[i][j] + U[i][j - 1]);
					Up = 0.5 * (Ue + Uw);
				}

				/* ************* */
				double R = Re * rn * depsn * Vn;
				double an_r = fmax(-R, 0.0) + visc_n*depsn * rn / drn;

				/* ************* */
				R = Re * rs * depss * Vs;
				double as_r = fmax(+R, 0.0) + visc_s*depss * rs / drs;

				/* ************* */
				R = Re * dre * Ue;
				double ae_r = fmax(-R, 0.0) + visc_e*dre / (re * depse);

				/* ************* */
				R = Re * drw * Uw;
				double aw_r = fmax(+R, 0.0) + visc_w*drw / (rw * depsw);

				/* ************* */
				double S1 = rp * depsp * (visc_n - visc_n) * (Vn - Vs) / drp + (visc_e - visc_w) * (Un - Us) - Up * drp * (visc_e / re - visc_w / rw) / rp;

				/* ************ */
				double S = Re * pow(Up, 2) * drp * depsp - 2.0 * visc_p * drp / rp * (Ue - Uw) - visc_p * Vp * drp * depsp / rp + v[i][j] * Re * rp * drp * depsp / dt + S1;

				/* ************ */
				ap_r[i][j] = an_r + as_r + ae_r + aw_r + Re * drp * rp * depsp / dt;

				/* ************ */
				V[i][j] = (an_r * VN + as_r * VS + ae_r * VE + aw_r * VW + S - rp * depsp * (P[i + 1][j] - P[i][j])) / ap_r[i][j];

				//if (V[i][j] != V[i][j])
				i = i;

			}

		}

	}

}

void Calculation_Velocity_Veps()
{
	double rp, rn, rs, re, rw;
	double drp, drn, drs, dre, drw;
	double depsp, depsn, depss, depse, depsw;

	double visc_s, visc_n, visc_w, visc_e, visc_p;

	double Vp, Vn, Vs, Ve, Vw;
	double Up, Ue, Uw;

	double UN, US, UE, UW;

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{

			if (sign[i][j] == Computational_Cell && (sign[i][j + 1] == Computational_Cell || sign[i][j + 1] == Gluing_Region))
			{

				/* Присвоение значений узлам КО */
				{
					rp = r_cell[i];
					rn = r_cell[i] + 0.5 * dr[i];
					rs = r_cell[i] - 0.5 * dr[i];
					re = r_cell[i];
					rw = r_cell[i];

					drp = dr[i];
					drn = 0.5 * (dr[i] + dr[i + 1]);
					drs = 0.5 * (dr[i] + dr[i - 1]);
					dre = dr[i];
					drw = dr[i];

					depsp = 0.5 * (deps[j] + deps[j + 1]);
					depsn = 0.5 * (deps[j] + deps[j + 1]);
					depss = 0.5 * (deps[j] + deps[j + 1]);
					depse = deps[j + 1];
					depsw = deps[j];

					visc_s = 0.5 * (0.5 * (visc[i][j] + visc[i][j + 1]) + 0.5 * (visc[i - 1][j] + visc[i - 1][j + 1]));
					visc_n = 0.5 * (0.5 * (visc[i][j] + visc[i][j + 1]) + 0.5 * (visc[i + 1][j] + visc[i + 1][j + 1]));
					visc_e = visc[i][j+1];
					visc_w = visc[i][j];
					visc_p = 0.5 * (visc_w + visc_e);

					UN = U[i + 1][j];
					US = U[i - 1][j];
					UE = U[i][j + 1];
					UW = U[i][j - 1];

					Up = U[i][j];
					Ue = 0.5 * (U[i][j] + U[i][j + 1]);
					Uw = 0.5 * (U[i][j] + U[i][j - 1]);

					Vp = 0.25 * (V[i][j] + V[i][j + 1] + V[i - 1][j] + V[i - 1][j + 1]);   //!!!!!!!!!!!!!!
					Vn = 0.5 * (V[i][j] + V[i][j + 1]);
					Vs = 0.5 * (V[i - 1][j] + V[i - 1][j + 1]);

					Ve = 0.5 * (V[i][j + 1] + V[i - 1][j + 1]);
					Vw = 0.5 * (V[i][j] + V[i - 1][j]);
				}

				/* ************* */
				double R = Re * rn * depsn * Vn;
				double an_eps = fmax(-R, 0.0) + visc_n*depsn * rn / drn;

				/* ************* */
				R = Re * rs * depss * Vs;
				double as_eps = fmax(+R, 0.0) + visc_s*depss * rs / drs;

				/* ************* */
				R = Re * dre * Ue;
				double ae_eps = fmax(-R, 0.0) + visc_e*dre / (re * depse);

				/* ************* */
				R = Re * drw * Uw;
				double aw_eps = fmax(+R, 0.0) + visc_w*drw / (rw * depsw);

				/* ************ */
				double S1 = (visc_n - visc_s) * (Ve - Vw) - Up * depsp * (visc_n - visc_s) + drp * (visc_e - visc_w) * (Ue / re - Uw / rw) / depsp + 2 * Vp * drp * (visc_e / re - visc_w / rw);

				/* ************ */
				double S = Re * Up * Vp * drp * depsp - visc_p * Up * drp * depsp / rp + 2.0 * visc_p * drp / rp * (Ve - Vw) + u[i][j] * Re * rp * drp * depsp / dt + S1;

				/* ************ */
				ap_eps[i][j] = an_eps + as_eps + ae_eps + aw_eps + Re * drp * rp * depsp / dt;

				/* ************ */
				U[i][j] = (an_eps * UN + as_eps * US + ae_eps * UE + aw_eps * UW + S - rp * depsp * (P[i][j + 1] - P[i][j])) / ap_eps[i][j];

				if (U[i][j] != U[i][j])
					i = i;

			}

			if (j == 20 || j == 21)
				i = i;

		}

	}

}

void Calculation_Pressure_P()
{

	double Epress = 1.0;

	Boundary_conditions();

	/* Обнуление */
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{

				PP[i][j] = 0.0;

			}

		}

	}


	double rp, rn, rs;
	double drp, dre, drw;
	double depsp, depsn, depss;

	double PPN, PPS, PPE, PPW;

	int iterP = 0;
	int ii, jj, ii1, jj1;

	while (Epress > 0.00001)
	{
		Epress = 0.0;

		iterP++;

		maxPP = 0.0;

		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= M; j++)
			{

				if (sign[i][j] == Computational_Cell)
				{

					rp = r_cell[i];
					rn = r_cell[i] + 0.5 * dr[i];
					rs = r_cell[i] - 0.5 * dr[i];

					drp = dr[i];
					dre = dr[i];
					drw = dr[i];

					depsp = deps[j];
					depsn = deps[j];
					depss = deps[j];

					PPN = PP[i + 1][j];
					PPS = PP[i - 1][j];
					PPE = PP[i][j + 1];
					PPW = PP[i][j - 1];

					double CN, CS, CW, CE;

					if (sign[i + 1][j] == Outer_Wall) CN = 0.0; else CN = rn * rn * depsn * depsp / ap_r[i][j];

					if (sign[i - 1][j] == Inner_Wall)  CS = 0; else CS = rs * rs * depss * depsp / ap_r[i - 1][j];

					if (sign[i][j + 1] == Inner_Wall)  CE = 0; else CE = dre * drp / ap_eps[i][j];

					if (sign[i][j - 1] == Inner_Wall)  CW = 0; else CW = drw * drp / ap_eps[i][j - 1];

					double CP = CN + CS + CE + CW;

					double pp1 = (CN * PPN + CS * PPS + CE * PPE + CW * PPW - (rn * V[i][j] * depsp - rs * V[i - 1][j] * depsp + U[i][j] * drp - U[i][j - 1] * drp)) / CP;

					if (Epress < abs(PP[i][j] - pp1))
					{

						Epress = abs(PP[i][j] - pp1);
						ii1 = i;
						jj1 = j;

						i = i;

					}

					PP[i][j] = pp1;

					if (abs(PP[i][j]) > maxPP)
					{

						maxPP = abs(PP[i][j]);

						ii = i;

						jj = j;

					}

					if (PP[i][j] != PP[i][j])
						i = i;

				}

			}

		}

		/* Граничные условия для ппроавки */
		{

			for (int i = 0; i <= N + 1; i++)
			{

				PP[i][0] = PP[i][M];
				PP[i][M + 1] = PP[i][1];

			}

		}

	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			drp = dr[i];
			depsp = deps[j];

			if (sign[i][j] == Computational_Cell)
				P[i][j] = P[i][j] + relax * PP[i][j];

			if (sign[i][j] == Computational_Cell && (sign[i][j + 1] == Computational_Cell || sign[i][j + 1] == Gluing_Region))
				U[i][j] = U[i][j] - drp * (PP[i][j + 1] - PP[i][j]) / ap_eps[i][j];

			if (sign[i][j] == Computational_Cell && sign[i + 1][j] == Computational_Cell)
				V[i][j] = V[i][j] - drp * depsp * (PP[i + 1][j] - PP[i][j]) / ap_r[i][j];

		}

	}

	Boundary_conditions();

}

void Write()
{

	int j_1 = round((eps_1 + eps_2) * 0.5 * M);


	if (iterGlob == 1)
	{
		ofstream file_E_vel("E_vel.DAT", ios_base::trunc);
		file_E_vel << fixed << setprecision(6) << "Re = " << Re << "\t" << N << "\t" << M << endl;

	}

	ofstream file_E_vel("E_vel.DAT", ios_base::app);
	file_E_vel << fixed << setprecision(6) << _time << "\t" << E_vel << endl;

	if (iterGlob == 1 || iterGlob % 100 == 0 || abs(_time - _time_final) < dt)
	{
		/* Файлы для записи */
		ofstream file_1("U.DAT", ios_base::trunc);
		ofstream file_2("V.DAT", ios_base::trunc);
		ofstream file_3("P.DAT", ios_base::trunc);

		file_1 << fixed << setprecision(6) << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;
		file_2 << fixed << setprecision(6) << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;
		file_3 << fixed << setprecision(6) << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

		for (int i = 0; i <= N; i++)
		{

			file_1 << fixed << setprecision(2) << r_cell[i] << fixed << setprecision(6) << "\t" << U[i][j_1] << endl;
			file_2 << fixed << setprecision(2) << r_cell[i] + 0.5 * dr[i] << fixed << setprecision(6) << "\t" << V[i][j_1] << endl;
			file_3 << fixed << setprecision(2) << r_cell[i] << fixed << setprecision(6) << "\t" << P[i][j_1] << endl;

		}

		file_1.close();
		file_2.close();
		file_3.close();
		file_E_vel.close();

	}

	cout << fixed << setprecision(4) << "time = " << _time << "\t" << " Re = " << Re << endl;
	cout << fixed << setprecision(7) << "Mesh:  " << N << "*" << M << "\t E_vel = " << E_vel << endl;
	if (iterGlob == 1) cin.get();
	cout << " r " << "\t" << "    U  " << "\t" << "\t" << "    V  " << "\t" << "\t" << "    P  " << endl;

	for (int i = 1; i <= N; i++)
	{

		cout << fixed << setprecision(2) << r_cell[i] << "\t";
		cout << fixed << setprecision(6) << U[i][j_1] << "\t" << V[i][j_1] << "\t" << P[i][j_1] << "\t" << endl;

	}

	cout << endl;

}

void WriteEnd()
{
	ofstream file_4("P_field.DAT", ios_base::trunc);

	file_4 << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

	for (int j = 0; j < M; j++)
	{

		double angle = (eps_cell[j]);

		for (int i = 0; i <= N + 1; i++)
		{
			double x = r_cell[i] * cos(angle);

			double y = r_cell[i] * sin(angle);

			file_4 << fixed << setprecision(6) << x << "\t" << y << "\t" << P[i][j] << "\t" << i << endl;

		}

		file_4 << j << endl;

	}

	ofstream file_5("U_field.DAT", ios_base::trunc);

	file_5 << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

	for (int j = 0; j < M; j++)
	{

		double angle = (eps_cell[j] + 0.5 * deps[j]);

		for (int i = 0; i <= N + 1; i++)
		{

			double x = r_cell[i] * cos(angle);

			double y = r_cell[i] * sin(angle);

			file_5 << fixed << setprecision(6) << x << "\t" << y << "\t" << U[i][j] << "\t" << i << endl;

		}

		file_5 << j << endl;

	}

	ofstream file_6("V_field.DAT", ios_base::trunc);

	file_6 << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

	for (int j = 0; j < M; j++)
	{
		double angle = eps_cell[j];

		for (int i = 0; i <= N + 1; i++)
		{

			double x = (r_cell[i] + 0.5 * dr[i]) * cos(angle);

			double y = (r_cell[i] + 0.5 * dr[i]) * sin(angle);

			file_6 << fixed << setprecision(6) << x << "\t" << y << "\t" << V[i][j] << "\t" << i << endl;

		}

		file_6 << j << endl;

	}

	file_4.close();
	file_5.close();
	file_6.close();

}

void BLN()
{

	const int NN = 720;

	double R1 = 1.0;

	double R0 = 0.2;

	ofstream file_7("Documents/BLN/BLN_R1.bln", ios_base::trunc);
	ofstream file_9("Documents/BLN/BLN_obstacle.bln", ios_base::trunc);

	file_7 << 4 * NN + 1 << " 0 " << endl;
	file_9 << M + 2 << " 1 " << endl;

	/* Отсечение внешней области */
	{

		for (int i = 0; i <= 2 * NN; i++)
		{

			double y = pow(R1 * R1 - (-1.0 + i * 1.0 / NN) * (-1.0 + i * 1.0 / NN), 0.5);

			file_7 << -1 + i * 1.0 / NN << " \t " << y << endl;

		}

		for (int i = 2 * NN - 1; i >= 0; i--)
		{

			double y = -pow(R1 * R1 - (-1.0 + i * 1.0 / NN) * (-1.0 + i * 1.0 / NN), 0.5);

			file_7 << -1 + i * 1.0 / NN << " \t " << y << endl;

		}

	}

	/* Отсечение внутренней области */
	{
		int i_1 = round(N * r_1);
		int i_2 = round(N * r_2);
		int j_1 = round(M * eps_1);
		int j_2 = round(M * eps_2);

		for (int j = 0; j <= M + 1; j++)
		{
			for (int i = 0; i <= N + 1; i++)
			{

				if (sign[i][j] == Inner_Wall && sign[i + 1][j] == Computational_Cell || j == 0 && i == 0 || sign[i][j + 1] == Inner_Wall && i == i_2)
				{

					double angle = (eps_cell[j] + 0.5 * deps[j]);

					file_9 << fixed << setprecision(6) << (r_cell[i]) * cos(angle) << " \t " << r_cell[i] * sin(angle) << endl;

					i = i;

				}

			}


		}

	}

	file_7.close();
	file_9.close();

}

void Mesh()
{
	BLN();

	ofstream file_7("Mesh_grf.DAT", ios_base::trunc);

	double h = 2.0 * Pi / 200;

	for (int i = N; i >= 0; i--)
	{

		for (int j = 0; j <= 200; j++)
		{

			double angle = (j * h + h) * 180 / Pi;

			file_7 << fixed << setprecision(6) << angle << "\t" << r_cell[i] + 0.5 * dr[i] << endl;

		}

		file_7 << endl;

	}

	for (int j = 0; j <= M; j++)
	{

		double angle = (eps_cell[j]) * 180 / Pi;

		file_7 << fixed << setprecision(6) << angle << "\t" << r_cell[0] + 0.5 * dr[0] << endl;

		file_7 << fixed << setprecision(6) << angle << "\t" << r_cell[N] + 0.5 * dr[N] << endl << endl;

	}

	file_7.close();

}

void Stream_Function()
{

	double psi_bound = 0.0;

	double Estream = 1.0;

	/* Начальное условие */
	{

		for (int i = 1; i <= N + 1; i++)
		{
			for (int j = 0; j <= M + 1; j++)
			{
				psi[i][j] = 0.0;
			}
		}

	}

	/* Интегрирование методом трапеций */
	double SS = 0.5 * 0.5 * dr[1] * (omega_0 * R0 + U[1][4]);
	SS += 0.5 * 0.5 * dr[N] * (U[N][4] + omega_1 * R1);

	for (int i = 2; i <= N; i++)
	{

		SS += 0.5 * 0.5 * (dr[i] + dr[i - 1]) * (U[i][4] + U[i - 1][4]);

	}

	SS = -SS;

	int iterPsi = 0;

	int ii, jj;

	while (Estream >= 0.0000001)
	{

		iterPsi++;

		Estream = 0.0;

		/* Граничные условия */
		{

			for (int i = 0; i <= N + 1; i++)
			{

				psi[i][0] = psi[i][M];
				psi[i][M + 1] = psi[i][1];

			}

			for (int j = 0; j <= M + 1; j++)
			{

				psi[0][j] = 0.0;
				psi[N][j] = SS;

			}

		}

		/* Вычисление функции тока */
		{

			for (int i = 1; i < N; i++)
			{
				for (int j = 1; j <= M; j++)
				{
					if (sign[i][j] == Computational_Cell && sign[i + 1][j] == Computational_Cell)
					{

						double rp = r_cell[i] + 0.5 * dr[i];
						double rn = r_cell[i + 1];
						double rs = r_cell[i];
						double re = r_cell[i] + 0.5 * dr[i];
						double rw = r_cell[i] + 0.5 * dr[i];

						double drp = 0.5 * (dr[i] + dr[i + 1]);
						double drn = dr[i + 1];
						double drs = dr[i];

						double depsp = 0.5 * (deps[j] + deps[j + 1]);
						double depse = deps[j + 1];
						double depsw = deps[j];

						double PsiP = psi[i][j];
						double PsiN = psi[i + 1][j];
						double PsiS = psi[i - 1][j];
						double PsiE = psi[i][j + 1];
						double PsiW = psi[i][j - 1];

						double Un = U[i + 1][j];
						double Us = U[i][j];
						double Ve = V[i][j + 1];
						double Vw = V[i][j];

						double S = -((rn * Un - rs * Us) / drp - (Ve - Vw) / depsp);
						double an = rn / (drp * drn);
						double as = rs / (drp * drs);
						double ae = 1.0 / (rp * depsp * depse);
						double aw = 1.0 / (rp * depsp * depsw);
						double ap = an + as + ae + aw;

						double tmp_psi = (an * PsiN + as * PsiS + ae * PsiE + aw * PsiW - S) / ap;

						if (Estream < abs(tmp_psi - psi[i][j]))
						{

							Estream = abs(tmp_psi - psi[i][j]);
							ii = i;
							jj = j;

						}

						psi[i][j] = tmp_psi;

						if (j == 21)
							i = i;

					}

					if (j == 20 || j == 21)
						i = i;

				}

			}

		}

	}

	/* Вывод в файл */
	{

		ofstream file_8("Psi_field.DAT", ios_base::trunc);

		file_8 << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

		for (int j = 1; j <= M + 1; j++)
		{

			double angle = (eps_cell[j] + 0.5 * deps[j]);

			for (int i = 0; i <= N; i++)
			{
				double x = (r_cell[i] + 0.5 * dr[i]) * cos(angle);

				double y = (r_cell[i] + 0.5 * dr[i]) * sin(angle);

				file_8 << fixed << setprecision(6) << x << "\t" << y << "\t" << psi[i][j] << "\t" << i << endl;

				if (j == 21)
					i = i;

			}

			file_8 << j << endl;

		}

		file_8.close();

		int j_1 = round((eps_1 + eps_2) * 0.5 * M);

		ofstream file_9("Psi.DAT", ios_base::trunc);

		for (int i = 0; i <= N; i++)
		{

			file_9 << fixed << setprecision(6) << r_cell[i] + 0.5 * dr[i] << "\t" << psi[i][21] << endl;

		}

		file_9.close();

	}

}

void Dissipative_Function()
{

	double dVr_dr, dVr_deps;
	double dVeps_dr, dVeps_deps;
	double F = 0.0, Vp;
	int iterP = 0;

	ofstream file_10("Q_field.DAT", ios_base::trunc);
	ofstream file_11("F.DAT", ios_base::trunc);

	file_10 << "time = " << _time << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

	for (int j = 0; j <= M; j++)
	{

		for (int i = 0; i <= N; i++)
		{

			if (sign[i][j] == Computational_Cell && j != 0)
			{

				dVr_dr = (V[i][j] - V[i - 1][j]) / (dr[i]);
				dVr_deps = ((V[i][j + 1] + V[i - 1][j + 1]) / 2.0 - (V[i][j - 1] + V[i - 1][j - 1]) / 2.0) / (2 * deps[j]);
				dVeps_dr = ((U[i + 1][j] + U[i + 1][j - 1]) / 2.0 - (U[i - 1][j] + U[i - 1][j - 1]) / 2.0) / (2 * dr[i]);
				dVeps_deps = (V[i][j] - V[i][j - 1]) / (deps[j]);
				Vp = (V[i][j] + V[i - 1][j]) / 2.0;

				double Q = 0.5 * (pow(2 * dVr_dr, 2) + 2 * r_cell[i] * r_cell[i] * pow(dVr_deps / pow(r_cell[i], 2) + dVeps_dr, 2) + pow(2 * (dVeps_deps + Vp / r_cell[i]), 2));

				F = F + Q * r_cell[i] * dr[i] * deps[j];

				i = i;

				/* Запись в файл */
				{

					double angle = eps_cell[j];

					double x = r_cell[i] * cos(angle);

					double y = r_cell[i] * sin(angle);

					file_10 << fixed << setprecision(6) << x << "\t" << y << "\t" << Q << "\t" << i << endl;

				}


			}

		}

		file_10 << j << endl;

	}

	file_11 << fixed << setprecision(6) << 32.0 * F / Re << "\t" << Re << endl;

	file_10.close();

}

void Write_Save()
{
	if (iterGlob == 1 || iterGlob % 200 == 0 || abs(_time - _time_final) < dt)
	{

		/* Запись начальных условий и сетки */
		{

			ofstream file_Save("Save/Mesh.DAT", ios_base::trunc);

			file_Save << _time << endl;

			file_Save << N << "\t" << M << endl;

			file_Save << R0 << "\t" << R1 << "\t" << alfa << "\t" << omega_0 << "\t" << omega_1 << "\t" << Re << "\t" << Pi << endl;

			for (int i = 0; i <= N + 1; i++)
			{

				file_Save << i << "\t" << dr[i] << "\t" << r_cell[i] << endl;

			}

			for (int j = 0; j <= M + 1; j++)
			{

				file_Save << j << "\t" << deps[j] << "\t" << eps_cell[j] << endl;

			}

			for (int i = 0; i <= N + 1; i++)
			{
				for (int j = 0; j <= M + 1; j++)
				{

					file_Save << sign[i][j] << endl;

				}

			}

			file_Save.close();

		}

		/* Запись V, U, P */
		{

			ofstream file_Save("Save/V_U_P.DAT", ios_base::trunc);

			for (int i = 0; i <= N + 1; i++)
			{
				for (int j = 0; j <= M + 1; j++)
				{

					file_Save << V[i][j] << "\t" << v[i][j] << "\t" << U[i][j] << "\t" << u[i][j] << "\t" << P[i][j] << endl;

				}

			}

			file_Save.close();

		}

	}

}

void Development()
{

	double VEL, vel;

	/* *************** */

	E_vel = 0.0;

	/* *************** */

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{

			VEL = 0.5 * sqrt(pow(U[i][j] + U[i - 1][j], 2.0) + pow(V[i][j] + V[i][j - 1], 2.0));

			vel = 0.5 * sqrt(pow(u[i][j] + u[i - 1][j], 2.0) + pow(v[i][j] + v[i][j - 1], 2.0));

			if (E_vel < fabs((VEL - vel) / VEL)) E_vel = fabs((VEL - vel) / VEL);

		}

	}

}

double Vapp(double eps, double r)
{

	for (int i = 0; i < N + 1; i++)
	{

		if ((r_cell[i] + 0.5 * dr[i]) <= r && (r_cell[i + 1] + 0.5 * dr[i + 1]) >= r)
		{

			I = i;

			break;

		}
	}

	for (int j = 0; j <= M + 1; j++)
	{

		if (eps_cell[j] <= eps && eps_cell[j + 1] >= eps)
		{

			J = j;

			break;

		}

	}

	double V_1 = V[I][J] + (eps - eps_cell[J]) / (eps_cell[J + 1] - eps_cell[J]) * (V[I][J + 1] - V[I][J]);

	double V_2 = V[I + 1][J] + (eps - eps_cell[J]) / (eps_cell[J + 1] - eps_cell[J]) * (V[I + 1][J + 1] - V[I + 1][J]);

	return V_1 + (r - (r_cell[I] + 0.5 * dr[I])) / dr[I + 1] * (V_2 - V_1);

}

double Uapp(double eps, double r)
{

	for (int i = 0; i <= N + 1; i++)
	{

		if (r_cell[i] <= r && r_cell[i + 1] >= r)
		{

			I = i;

			break;

		}
	}

	for (int j = 0; j <= M + 1; j++)
	{

		if ((eps_cell[j] + 0.5 * deps[j]) <= eps && (eps_cell[j + 1] + 0.5 * deps[j + 1]) >= eps)
		{

			J = j;

			break;

		}

	}

	double U_1 = U[I][J] + (eps - (eps_cell[J] + 0.5 * deps[J])) / ((eps_cell[J + 1] + 0.5 * deps[J]) - (eps_cell[J] + 0.5 * deps[J])) * (U[I][J + 1] - U[I][J]);

	double U_2 = U[I + 1][J] + (eps - (eps_cell[J] + 0.5 * deps[J])) / ((eps_cell[J + 1] + 0.5 * deps[J]) - (eps_cell[J] + 0.5 * deps[J])) * (U[I + 1][J + 1] - U[I + 1][J]);

	return U_1 + (r - r_cell[I]) / dr[I + 1] * (U_2 - U_1);

}

void Flow_Evolution()
{

	/* Начально условие */
	if (iterGlob == 1)
	{

		count_marker = N + 1;

		int init_J = 10;

		r_m[0] = R0;
		eps_m[0] = eps_cell[init_J];

		r_m[N + 1] = R1;
		eps_m[N + 1] = eps_cell[init_J];

		for (int i = 1; i <= N; i++)
		{

			r_m[i] = r_cell[i];
			eps_m[i] = eps_cell[init_J];

			i = i;

		}

	}

	if (iterGlob % 500 == 0 || iterGlob == 1)
	{

		ofstream f1("Flow Evolution/Flow Evolution " + to_string(iterGlob / 100) + ".DAT", ios_base::trunc);
		f1 << "time = " << tt << "\t" << "Re = " << Re << "\t" << N << "\t" << M << endl;

	}

	/* Основной цикл */
	{

		for (int i = 1; i < count_marker; i++)
		{

			double r_tmp = r_m[i] + Vapp(eps_m[i], r_m[i]) * dt;

			double eps_tmp = eps_m[i] + Uapp(eps_m[i], r_m[i]) / r_m[i] * dt;

			r_m[i] = r_tmp;

			eps_m[i] = eps_tmp;

			if (iterGlob % 500 == 0 || iterGlob == 1)
			{

				ofstream f1("Flow Evolution/Flow Evolution " + to_string(iterGlob / 100) + ".DAT", ios_base::app);

				double angle = (eps_m[i]/* - omega_1 * _time*/);

				double x = (r_m[i] * cos(angle));

				double y = (r_m[i] * sin(angle));

				f1 << fixed << setprecision(6) << r_m[i] * cos(angle) << "\t" << r_m[i] * sin(angle) << "\t";
				f1 << fixed << setprecision(6) << r_m[i] << "\t" << eps_m[i] << endl;

				i = i;

			}

		}

		if (iterGlob % 45000 == 0 || iterGlob == 1)
			int i = 0;

		for (int i = 1; i < count_marker; i++)
		{

			double dx = (r_m[i] * cos(eps_m[i]) - r_m[i - 1] * cos(eps_m[i - 1])) * (r_m[i] * cos(eps_m[i]) - r_m[i - 1] * cos(eps_m[i - 1]));
			double dy = (r_m[i] * sin(eps_m[i]) - r_m[i - 1] * sin(eps_m[i - 1])) * (r_m[i] * sin(eps_m[i]) - r_m[i - 1] * sin(eps_m[i - 1]));

			double ds = sqrt(dx + dy);

			if (ds > 0.03)
			{

				/* смещение от i + 1 до count_marker вправо на 1 */
				for (int j = count_marker; j >= i + 1; j--)
				{

					r_m[j] = r_m[j - 1];
					eps_m[j] = eps_m[j - 1];

				}

				/* вставка нового маркера */
				r_m[i] = 0.5 * (r_m[i - 1] + r_m[i]);
				eps_m[i] = 0.5 * (eps_m[i - 1] + eps_m[i]);

				count_marker++;

			}

		}

	}

}

int main()
{
	Initial_Conditions();
	Mesh();

	iterGlob = 0;
	tt = 0.0;

	do
	{

		iterGlob++;

		do
		{

			Boundary_conditions();
			Calculation_Velocity_Vr();
			Calculation_Velocity_Veps();
			Calculation_Pressure_P();

		} while (maxPP >= 0.0001);

		Development();
		Redistricting();
		Write();
		Write_Save();
		//Flow_Evolution();

		_time += dt;
		tt += dt;

		if (_time > 25.0) break;

	} while (E_vel > 0.000001); /*(_time < 85.0)*/;

	Stream_Function();
	Dissipative_Function();
	WriteEnd();
	Arrays_Remove();

}
