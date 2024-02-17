//#define _CRT_SECURE_NO_WARNINGS

#pragma once

/* Параметры расчета */
double R0 = 0.2;
double R1 = 1.0;
double alfa = R0 / R1;

double omega_0 = 0.0;
double omega_1 = 1.0;

double Re = 1.0;
double Pi = 3.14159;

int count_marker = 0;

/* Количество разбиений по пространственной координате r */
int N = 40;

/* Количество разбиений по пространственной координате eps */
int M = 80; //36//154

/* Шаг по r и eps*/
double dr_1 = (1.0 - alfa) / N;
double deps_1 = 2.0 * Pi / M;

/* Массивы системы координат*/
double* r_cell, * eps_cell;
double* dr, * deps;

/* Параметры маркерной визуализации */
double* r_m, * eps_m;
int I, J;

/* Массивы скорости на текущем шаге*/
double** V;
double** U;

/* Массивы скорости на предыдущем шаге*/
double** v;
double** u;

/* Массивы давления на текущем шаге */
double** P;

/* Массивы давления на предыдущем шаге */
double** p;

/* Массивы вязкости на текущем шаге*/
double** VISC;

/* Массивы вязкости на предыдущем шаге*/
double** visc;

/* Массивы поправки давления на текущем шаге */
double** PP;

double maxPP = 1.0;

/* Расчетные коэффициенты по скорости */
double** ap_r;
double** ap_eps;

/* Функция тока */
double** psi;

/* Время */
double _time = 0.0;

/* Шаг по времени */
double dt = 0.0001;

double _time_final = 25.01;
double _time_write = _time_final;

/* Начальные условия */
bool Steady_Cond = false;

/* Сетка */
bool Regular_Mesh = true;

/* Старт с файла */
bool Start_From_Save = false;

/* Параметр установления */
bool CVel;
double E_vel;

/* Наличие препятствия */
bool Obstacle = true;
double r_1 = 0.0;
double r_2 = 0.0;
double eps_1 = 0.23; //0.73
double eps_2 = 0.25; //0.75

/* Коэффициент релаксации */
double relax = 0.3;

int iterGlob;
double tt;

/* Топологический признак */
int** sign;
enum Topological_Attribute { Inner_Wall, Outer_Wall, Gluing_Region, Computational_Cell };

