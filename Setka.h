#pragma once
#include <vector>
#include "Yzel.h"
#include "Header.h"

using namespace std;


class Setka
{
public:
	struct Geo_param geo;           // √еометрические параметры сетки



	vector<Yzel*> All_Yzel;
	vector<Luch*> All_Luch;
	vector<Luch*> A_Luch;
	vector<Luch*> B_Luch;
	vector<Luch*> C_Luch;
	vector<Luch*> D_Luch;
	vector<Luch*> E_Luch;
	vector<Luch*> H_Luch;
	vector<Luch*> G_Luch;

	vector<Cell*> All_Cell;

	

	Setka();
	void Set_geo();
	void Construct_initial();

	void Print_yzel();    // ѕечать всех узлов
	void Print_cell();    // ѕечать всех €чеек

};

