//#define _CRT_SECURE_NO_WARNINGS

#pragma once

/* ��������� ������� */
double R0 = 0.2;
double R1 = 1.0;
double alfa = R0 / R1;

double omega_0 = 0.0;
double omega_1 = 1.0;

double Re = 1.0;
double Pi = 3.14159;

int count_marker = 0;

/* ���������� ��������� �� ���������������� ���������� r */
int N = 40;

/* ���������� ��������� �� ���������������� ���������� eps */
int M = 80; //36//154

/* ��� �� r � eps*/
double dr_1 = (1.0 - alfa) / N;
double deps_1 = 2.0 * Pi / M;

/* ������� ������� ���������*/
double* r_cell, * eps_cell;
double* dr, * deps;

/* ��������� ��������� ������������ */
double* r_m, * eps_m;
int I, J;

/* ������� �������� �� ������� ����*/
double** V;
double** U;

/* ������� �������� �� ���������� ����*/
double** v;
double** u;

/* ������� �������� �� ������� ���� */
double** P;

/* ������� �������� �� ���������� ���� */
double** p;

/* ������� �������� �� ������� ����*/
double** VISC;

/* ������� �������� �� ���������� ����*/
double** visc;

/* ������� �������� �������� �� ������� ���� */
double** PP;

double maxPP = 1.0;

/* ��������� ������������ �� �������� */
double** ap_r;
double** ap_eps;

/* ������� ���� */
double** psi;

/* ����� */
double _time = 0.0;

/* ��� �� ������� */
double dt = 0.0001;

double _time_final = 25.01;
double _time_write = _time_final;

/* ��������� ������� */
bool Steady_Cond = false;

/* ����� */
bool Regular_Mesh = true;

/* ����� � ����� */
bool Start_From_Save = false;

/* �������� ������������ */
bool CVel;
double E_vel;

/* ������� ����������� */
bool Obstacle = true;
double r_1 = 0.0;
double r_2 = 0.0;
double eps_1 = 0.23; //0.73
double eps_2 = 0.25; //0.75

/* ����������� ���������� */
double relax = 0.3;

int iterGlob;
double tt;

/* �������������� ������� */
int** sign;
enum Topological_Attribute { Inner_Wall, Outer_Wall, Gluing_Region, Computational_Cell };

