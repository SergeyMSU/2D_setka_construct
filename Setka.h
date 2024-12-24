#pragma once
#include <vector>
#include "Yzel.h"
#include "Header.h"

using namespace std;


class Setka
{
public:
	struct Geo_param geo;           // �������������� ��������� �����



	vector<Yzel*> All_Yzel;
	vector<Luch*> All_Luch;
	vector<Luch*> A_Luch;
	vector<Luch*> B_Luch;
	vector<Luch*> C_Luch;
	vector<Luch*> D_Luch;
	vector<Luch*> E_Luch;
	vector<Luch*> F_Luch;  // �� ����� ����
	vector<Yzel*> Yzel_zav_1;   // ���� � ������������ 1 (�.�. ����������� ������ ������������ �������� �����)
	vector<Luch*> G_Luch;  // �� ����� ����

	vector<Cell*> All_Cell;

	

	Setka();
	void Set_geo();
	void Construct_initial();

	void Print_yzel();    // ������ ���� �����
	void Print_cell();    // ������ ���� �����

};

