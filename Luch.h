#pragma once
using namespace std;
#include <string>
#include <vector>
#include "Header.h"
class Yzel;


class Luch
{
public:
	string type;            // ��� ����
	double the, phi; 

	vector<Yzel*> Yzels;          // ��� ����� ����
	vector<Yzel*> Yzels_opor;     // ������� ����� ����  (������� ����� ������ �������!)
	vector<Luch*> Luch_soseds;    // ����-������ (�� ��� ����� ������ �������� �����)
};

