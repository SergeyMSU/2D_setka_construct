#pragma once
#include <vector>
#include "Yzel.h"
#include "Header.h"

using namespace std;


class Setka
{
public:
	struct Geo_param geo;           // Геометрические параметры сетки



	vector<Yzel*> All_Yzel;
	vector<Luch*> All_Luch;
	vector<Luch*> A_Luch;
	vector<Luch*> B_Luch;
	vector<Luch*> C_Luch;
	vector<Luch*> D_Luch;
	vector<Luch*> E_Luch;
	vector<Luch*> F_Luch;  // Он такой один
	vector<Yzel*> Yzel_zav_1;   // Узлы с зависимостью 1 (т.е. опредляются только координатами основных узлов)
	vector<Luch*> G_Luch;  // Он такой один

	vector<Cell*> All_Cell;

	

	Setka();
	void Set_geo();
	void Construct_initial();

	void Print_yzel();    // Печать всех узлов
	void Print_cell();    // Печать всех ячеек

};

