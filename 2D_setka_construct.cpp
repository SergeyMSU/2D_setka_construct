// 2D_setka_construct.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Header.h"

int main()
{
    Setka S = Setka();

    S.Set_geo();
    S.Construct_initial();
    S.Print_yzel();
    S.Print_cell();


}
