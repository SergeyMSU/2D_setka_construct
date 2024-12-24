#pragma once
using namespace std;
#include <string>
#include <vector>
#include "Header.h"
class Yzel;


class Luch
{
public:
	string type;            // Тип луча
	double the, phi; 

	vector<Yzel*> Yzels;          // Все точки луча
	vector<Yzel*> Yzels_opor;     // Опорные точки луча  (крайние точки всегда опорные!)
	vector<Luch*> Luch_soseds;    // Лучи-соседи (по ним легко искать соседние точки)
};

