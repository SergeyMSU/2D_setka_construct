#include "Setka.h"
#include "Header.h"
#include <iostream>
#include <algorithm>
#include <fstream>


using namespace std;

Setka::Setka()
{

}

void Setka::Set_geo()
{
    this->geo.phi = 0.0;
    this->geo.tetta0 = const_pi / 18;
    this->geo.tetta1 = const_pi - const_pi / 12;
    this->geo.tetta2 = const_pi * 120.0 / 180.0;
    this->geo.N1 = 20;
    this->geo.N2 = 10;
    this->geo.N3 = 15; //10
    this->geo.N4 = 10;
    this->geo.N5 = 8;

    this->geo.M0 = 14;//14;
    this->geo.M1 = 15;//15;
    this->geo.M11 = 3;
    this->geo.M2 = 10;//10;
    this->geo.M3 = 20;
    this->geo.M4 = 30;
    this->geo.MF = 5;//4;

    this->geo.R0 = 1.0;
    this->geo.R1 = 3.0;
    this->geo.R2 = 6.0;
    this->geo.R3 = 9.0;
    this->geo.R4 = 12.0;
    this->geo.R5 = 20.0;
    this->geo.L6 = -10.0;
    this->geo.L7 = -20.0;

    this->geo.dd1 = 4.0;    // � ���� ������� ���������� �������������� ������
    this->geo.dd2 = 2.0;    // � ���� ������� ���������� �������������� ������

    this->geo.dd3 = 2.2;
    this->geo.dd4 = 2.2;

    this->geo.dd5 = 2.3;
    this->geo.dd6 = 3.0;

    this->geo.dd7 = 3.8;
    this->geo.dd8 = 5.3;
}

void Setka::Construct_initial()
{
    double x, y, z, the, r;

    // ��������� ���� �� �-�����
    for (int i = 0; i < this->geo.N1; i++)  
    {
        the = this->geo.tetta0 + i * (const_pi/2.0 - this->geo.tetta0)/(this->geo.N1 - 1);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "A";
        this->All_Luch.push_back(A);
        this->A_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;
        
        // ����� �� R1 �� R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;


        // ����� �� R2 (TS) �� R3 (HP)
        for (int j = 0; j < this->geo.M2; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R3 (HP) 
        r = this->geo.R3;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // ����� �� R3 (HP) �� R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R4;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // ����� �� R4 (BS) �� R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R5;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;

    }

    // ��������� ���� �� B-�����
    for (int i = 0; i < this->geo.N2; i++)
    {
        the = const_pi / 2.0 + (i + 1) * (this->geo.tetta2 - const_pi / 2.0) / (this->geo.N2);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "B";
        this->All_Luch.push_back(A);
        this->B_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;

        // ����� �� R1 �� R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;

        // ����� �� R2 (TS) �� R2 ++
        for (int j = 0; j < this->geo.M11; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������ �������� ����
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = cos(the) * this->geo.dd3;
        tt1 = sin(the) * this->geo.dd3;
        t2 = 0.0;
        tt2 = 1.0 * this->geo.dd4;
        x0 = x;        // ��� ���������� � ��������� �������� �����
        y0 = y;
        x1 = x;
        y1 = this->geo.R3;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // ����� �� R2 (TS) �� R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R3 (HP) 
        r = this->geo.R3;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // ����� �� R3 (HP) �� R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R4;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // ����� �� R4 (BS) �� R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R5;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;

    }

    // ��������� ���� �� C-�����
    for (int i = 0; i < this->geo.N3; i++)
    {
        the = this->geo.tetta2 + (i + 1) * (this->geo.tetta1 - this->geo.tetta2) / (this->geo.N3);
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "C";
        this->All_Luch.push_back(A);
        this->C_Luch.push_back(A);

        for (int j = 0; j < this->geo.M0; j++)
        {
            r = this->geo.R0 + j * (this->geo.R1 - this->geo.R0) / this->geo.M0;
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R1 
        r = this->geo.R1;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz1 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz1);
        A->Yzels_opor.push_back(Yz1);
        this->All_Yzel.push_back(Yz1);
        Yz1->luch = A;

        // ����� �� R1 �� R2 (TS)
        for (int j = 0; j < this->geo.M1; j++)
        {
            r = this->geo.R1 + (j + 1) * (this->geo.R2 - this->geo.R1) / (this->geo.M1 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R2 (TS) 
        r = this->geo.R2;
        x = r * cos(the);
        y = r * sin(the) * cos(A->phi);
        z = r * sin(the) * sin(A->phi);
        auto Yz2 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz2);
        A->Yzels_opor.push_back(Yz2);
        this->All_Yzel.push_back(Yz2);
        Yz2->luch = A;

        // ����� �� R2 (TS) �� R2 ++
        for (int j = 0; j < this->geo.M11; j++)
        {
            r = this->geo.R2 + (j + 1) * (this->geo.R3 - this->geo.R2) / (this->geo.M2 + 1);
            x = r * cos(the);
            y = r * sin(the) * cos(A->phi);
            z = r * sin(the) * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������ �������� ����
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        if (i == 0)
        {
            t1 = (4 * cos(the) - 6 * 1.0) / 10.0 * this->geo.dd5;
            tt1 = sin(the) * this->geo.dd5;
        }
        else
        {
            t1 = (cos(the) - 1.0) / 2.0 * this->geo.dd5;
            tt1 = sin(the) * this->geo.dd5;
        }
        t2 = -1.0 * this->geo.dd6;
        tt2 = 0.0;
        x0 = x;        // ��� ���������� � ��������� �������� �����
        y0 = y;
        x1 = this->geo.L6;
        y1 = y;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // ����� �� R2 (TS) �� L6
        for (int j = 0; j < this->geo.N4; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.N4 + 1);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ����� �� L6 �� L7
        for (int j = 0; j < this->geo.N5; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.N4 + 1);
            x = this->geo.L6 + (j) * (this->geo.L7 - this->geo.L6)/(this->geo.N5 - 1);
            y = y;
            z = z;
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;

            if(j == 0) A->Yzels_opor.push_back(Yz);
            if(j == this->geo.N5 - 1) A->Yzels_opor.push_back(Yz);
        }
    }

    // ��������� ���� �� D-�����
    for (int i = 0; i < this->geo.N4 + this->geo.N5 - 3; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "D";
        this->All_Luch.push_back(A);
        this->D_Luch.push_back(A);

        double x, y, z;

        auto yz = this->C_Luch[0]->Yzels[this->geo.M0 + 1 + this->geo.M1 + 1 + this->geo.M11 + 3 + i];

        x = yz->coord[0][0];
        y = yz->coord[0][1];
        z = yz->coord[0][2];
        y = sqrt(y * y + z * z);
        A->Yzels.push_back(yz);
        A->Yzels_opor.push_back(yz);

        // ������ �������� ����
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = 0.0;
        tt1 = 1.0 * 2;
        t2 = 0.0;
        tt2 = 1.0 * 2;
        x0 = x;        // ��� ���������� � ��������� �������� �����
        y0 = y;
        x1 = x;
        y1 = this->geo.R3;

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // ����� �� R2 (TS) �� R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11; j++)
        {
            double s = 1.0 * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            //x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            //y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            x = x0 + (x1 - x0) * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            y = y0 + (y1 - y0) * (j + 1) / (this->geo.M2 - this->geo.M11 + 1);
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R3 (HP) 
        r = this->geo.R3;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // ����� �� R3 (HP) �� R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R4;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // ����� �� R4 (BS) �� R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            x = x0;
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R5;
        x = x0;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;


    }

    // ��������� ���� �� E-�����
    for (int i = 0; i < 4; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "E";
        this->All_Luch.push_back(A);
        this->E_Luch.push_back(A);

        auto yz1 = this->B_Luch[size(this->B_Luch) - 1]->Yzels_opor[2];
        auto yz2 = this->D_Luch[0]->Yzels_opor[1];

        double x, y, r;

        x = yz1->coord[0][0] + (i + 1) * (yz2->coord[0][0] - yz1->coord[0][0]) / 5;
        y = yz1->coord[0][1] + (i + 1) * (yz2->coord[0][1] - yz1->coord[0][1]) / 5;
        z = yz1->coord[0][2] + (i + 1) * (yz2->coord[0][2] - yz1->coord[0][2]) / 5;

        r = sqrt(y * y + z * z);

        // ������� ����� �� ������� R3 (HP) 
        r = this->geo.R3;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz3 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz3);
        A->Yzels_opor.push_back(Yz3);
        this->All_Yzel.push_back(Yz3);
        Yz3->luch = A;

        // ����� �� R3 (HP) �� R4 (BS)
        for (int j = 0; j < this->geo.M3; j++)
        {
            r = this->geo.R3 + (j + 1) * (this->geo.R4 - this->geo.R3) / (this->geo.M3 + 1);
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R4;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz4 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz4);
        A->Yzels_opor.push_back(Yz4);
        this->All_Yzel.push_back(Yz4);
        Yz4->luch = A;


        // ����� �� R4 (BS) �� R5
        for (int j = 0; j < this->geo.M4; j++)
        {
            r = this->geo.R4 + (j + 1) * (this->geo.R5 - this->geo.R4) / (this->geo.M4 + 1);
            y = r * cos(A->phi);
            z = r * sin(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        // ������� ����� �� ������� R4  
        r = this->geo.R5;
        y = r * cos(A->phi);
        z = r * sin(A->phi);
        auto Yz5 = new Yzel(x, y, z);
        A->Yzels.push_back(Yz5);
        A->Yzels_opor.push_back(Yz5);
        this->All_Yzel.push_back(Yz5);
        Yz5->luch = A;


    }

    // ��������� H-����
    for (int i = 0; i < this->geo.MF - 1; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "H";
        this->All_Luch.push_back(A);
        this->H_Luch.push_back(A);

        int nn = this->B_Luch.size();
        auto yz1 = this->B_Luch[nn - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];
        auto yz3 = this->B_Luch[nn - 2]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];
        auto yz2 = this->D_Luch[0]->Yzels[1 + i];
        auto yz4 = this->D_Luch[1]->Yzels[1 + i];

        A->Yzels.push_back(yz1);
        A->Yzels_opor.push_back(yz1);

        // ������ �������� ����
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        t1 = (yz1->coord[0][0] - yz3->coord[0][0]) * this->geo.dd7;
        tt1 = (yz1->coord[0][1] - yz3->coord[0][1]) * this->geo.dd7;
        t2 = (yz4->coord[0][0] - yz2->coord[0][0]) * this->geo.dd8;
        tt2 = (yz4->coord[0][1] - yz2->coord[0][1]) * this->geo.dd8;
        x0 = yz1->coord[0][0];        // ��� ���������� � ��������� �������� �����
        y0 = yz1->coord[0][1];
        x1 = yz2->coord[0][0];        // ��� ���������� � ��������� �������� �����
        y1 = yz2->coord[0][1];

        double a1, b1, c1, d1, a2, b2, c2, d2;

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;


        for (int j = 0; j < 4; j++)
        {
            double s = 1.0 * (j + 1) / (5);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;
            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        A->Yzels.push_back(yz2);
        A->Yzels_opor.push_back(yz2);
    }

    // ��������� G-����
    for (int i = 0; i < 4; i++)
    {
        the = 0.0;
        auto A = new Luch();
        A->phi = this->geo.phi;
        A->the = the;
        A->type = "G";
        this->All_Luch.push_back(A);
        this->G_Luch.push_back(A);

        int k = this->H_Luch.size();

        auto yz1 = this->H_Luch[k - 1]->Yzels[i + 1];
        auto yz2 = this->E_Luch[i]->Yzels_opor[0];
        auto yz3 = this->H_Luch[k - 2]->Yzels[i + 1];

        auto yz4 = this->B_Luch[this->B_Luch.size() - 1]->Yzels[this->geo.M0 + 1 + geo.M1 + 1 + geo.M2];
        auto yz5 = this->B_Luch[this->B_Luch.size() - 1]->Yzels[this->geo.M0 + 1 + geo.M1 + 1 + geo.M2 - 1];

        double x, y, r;

        A->Yzels.push_back(yz1);
        A->Yzels_opor.push_back(yz1);

        // ������ �������� ����
        double x0, y0, x1, y1, t1, t2, tt1, tt2;

        // ��� ������� �������������
        double l1 = sqrt(kv((yz1->coord[0][0] - yz3->coord[0][0])) +
            kv((yz1->coord[0][1] - yz3->coord[0][1])));
        double l3 = sqrt(kv((yz4->coord[0][0] - yz5->coord[0][0])) +
            kv((yz4->coord[0][1] - yz5->coord[0][1])));

        double a1, b1, c1, d1, a2, b2, c2, d2;
        bool goto_aa1;

    aa1:
        goto_aa1 = false;
        t1 = (yz1->coord[0][0] - yz3->coord[0][0]) * this->geo.dd1;
        tt1 = (yz1->coord[0][1] - yz3->coord[0][1]) * this->geo.dd1;
        t2 = 0.0;
        tt2 = 1.0 * this->geo.dd2;
        x0 = yz1->coord[0][0];        // ��� ���������� � ��������� �������� �����
        y0 = yz1->coord[0][1];
        x1 = yz2->coord[0][0];
        y1 = yz2->coord[0][1];

        

        a1 = x0;
        b1 = t1;
        c1 = -2 * t1 - t2 - 3.0 * x0 + 3.0 * x1;
        d1 = t1 + t2 + 2.0 * x0 - 2.0 * x1;

        a2 = y0;
        b2 = tt1;
        c2 = -2 * tt1 - tt2 - 3.0 * y0 + 3.0 * y1;
        d2 = tt1 + tt2 + 2.0 * y0 - 2.0 * y1;

        // ����� �� R2 (TS) �� R3 (HP)
        for (int j = 0; j < this->geo.M2 - this->geo.M11 - this->geo.MF + 1; j++)
        {
            double s = 1.0 * (j + 1.0) / (this->geo.M2 - this->geo.M11 - this->geo.MF + 1 + 1.0);
            x = a1 + b1 * s + c1 * s * s + d1 * s * s * s;
            y = a2 + b2 * s + c2 * s * s + d2 * s * s * s;

            // ������ �������������� �������������
            if (i == 0 && j == 0) 
            {
                double l2 = sqrt(kv(x - x0) + kv(y - y0));
                if (fabs(l2 - l1) / l1 * 100.0 > 2.0)
                {
                    if (l2 < l1)
                    {
                        this->geo.dd1 *= 1.01;
                    }
                    else
                    {
                        this->geo.dd1 *= 0.99;
                    }
                    //cout << "1 " << l1 << " " << l2 << " " << this->geo.dd1 << endl;
                    //system("pause");
                    goto_aa1 = true;
                    break;
                }

                double xx, yy, ss;
                ss = 1.0 * (this->geo.M2 - this->geo.M11 - this->geo.MF + 1) / (this->geo.M2 - this->geo.M11 - this->geo.MF + 1 + 1.0);
                xx = a1 + b1 * ss + c1 * ss * ss + d1 * ss * ss * ss;
                yy = a2 + b2 * ss + c2 * ss * ss + d2 * ss * ss * ss;
                double l4 = sqrt(kv(xx - x1) + kv(yy - y1));
                if (fabs(l4 - l3) / l1 * 100.0 > 2.0)
                {
                    if (l4 < l3)
                    {
                        this->geo.dd2 *= 1.01;
                    }
                    else
                    {
                        this->geo.dd2 *= 0.99;
                    }
                    //cout << "2 " << l3 << " " << l4 << " " << this->geo.dd2 << endl;
                    //system("pause");
                    goto_aa1 = true;
                    break;
                }

            }


            z = y * sin(A->phi);
            y = y * cos(A->phi);
            auto Yz = new Yzel(x, y, z);
            A->Yzels.push_back(Yz);
            this->All_Yzel.push_back(Yz);
            Yz->luch = A;
        }

        if (goto_aa1 == true)
        {
            goto aa1;
        }


        A->Yzels.push_back(yz2);
        A->Yzels_opor.push_back(yz2);
    }




    // ��������� ������ 
    // /////////////////////////////////////////////////////////////////////////////////////////////
    // � - ����
    int n1 = this->A_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->A_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->A_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->A_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->A_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->A_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }

    // B - ����

    n1 = this->B_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->B_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->B_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->B_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->B_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->B_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }

    // C - ����

    n1 = this->C_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->C_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->C_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->C_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->C_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->C_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }

    // D - ����

    n1 = this->D_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->D_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->D_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->D_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->D_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->D_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }

    // E - ����

    n1 = this->E_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->E_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->E_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->E_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->E_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->E_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }


    // G - ����

    n1 = this->G_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->G_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->G_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->G_Luch[i]->Yzels[j + 1]);
            A->Yzels.push_back(this->G_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->G_Luch[i + 1]->Yzels[j]);
            this->All_Cell.push_back(A);
        }
    }


    // H - ����

    n1 = this->H_Luch.size();

    for (int i = 0; i < n1 - 1; i++)
    {
        int n2 = this->H_Luch[i]->Yzels.size();
        for (int j = 0; j < n2 - 1; j++)
        {
            auto A = new Cell();
            A->Yzels.push_back(this->H_Luch[i]->Yzels[j]);
            A->Yzels.push_back(this->H_Luch[i + 1]->Yzels[j]);
            A->Yzels.push_back(this->H_Luch[i + 1]->Yzels[j + 1]);
            A->Yzels.push_back(this->H_Luch[i]->Yzels[j + 1]);
            this->All_Cell.push_back(A);
        }
    }


    // BC - �������
    for (int j = 0; j < this->geo.M0 + this->geo.M1 + this->geo.M11 + 1; j++)
    {
        int k = this->B_Luch.size();
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->C_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->C_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // AB - �������
    for (int j = 0; j < this->B_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->A_Luch.size();
        auto A = new Cell();
        A->Yzels.push_back(this->A_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->A_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->B_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->B_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // BE - �������
    for (int j = 0; j < this->E_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->B_Luch.size();
        int nn = this->geo.M0 + this->geo.M1 + this->geo.M2 + 2;
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->E_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->E_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // ED - �������
    for (int j = 0; j < this->E_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->E_Luch.size();
        int nn = this->geo.M2 - this->geo.M11 + 1;
        auto A = new Cell();
        A->Yzels.push_back(this->E_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->E_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j]);
        this->All_Cell.push_back(A);
    }


    // BG - �������
    for (int j = 0; j < this->G_Luch[0]->Yzels.size() - 1; j++)
    {
        int k = this->B_Luch.size();
        int nn = this->geo.M0 + this->geo.M1 + this->geo.M11 + this->geo.MF;
        auto A = new Cell();
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j]);
        A->Yzels.push_back(this->B_Luch[k - 1]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->G_Luch[0]->Yzels[j + 1]);
        A->Yzels.push_back(this->G_Luch[0]->Yzels[j]);
        this->All_Cell.push_back(A);
    }

    // GD - �������
    for (int j = 0; j < this->G_Luch[0]->Yzels.size() - 2; j++)
    {
        int k = this->E_Luch.size();
        int nn = this->geo.MF - 1;
        auto A = new Cell();
        A->Yzels.push_back(this->G_Luch[k - 1]->Yzels[j]);
        A->Yzels.push_back(this->G_Luch[k - 1]->Yzels[j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j + 1]);
        A->Yzels.push_back(this->D_Luch[0]->Yzels[nn + j]);
        this->All_Cell.push_back(A);
    }

    // ������ ������� ������ �����

    if (true)
    {
        //1
        int k = this->B_Luch.size();
        auto yz1 = this->B_Luch[k - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1];
        auto yz2 = this->B_Luch[k - 1]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2];
        auto yz3 = this->H_Luch[0]->Yzels[1];
        auto yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1];

        auto A = new Cell();
        A->Yzels.push_back(yz1);
        A->Yzels.push_back(yz2);
        A->Yzels.push_back(yz3);
        A->Yzels.push_back(yz4);
        this->All_Cell.push_back(A);



        //3
        for (int i = 0; i < 3; i++)
        {
            yz1 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1 + i];
            yz2 = this->H_Luch[0]->Yzels[1 + i];
            yz3 = this->H_Luch[0]->Yzels[2 + i];
            yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + i];

            A = new Cell();
            A->Yzels.push_back(yz1);
            A->Yzels.push_back(yz2);
            A->Yzels.push_back(yz3);
            A->Yzels.push_back(yz4);
            this->All_Cell.push_back(A);
        }

        //4 
        yz1 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 1 + 3];
        yz2 = this->H_Luch[0]->Yzels[4];
        yz3 = this->D_Luch[0]->Yzels[1];
        yz4 = this->C_Luch[0]->Yzels[this->geo.M0 + this->geo.M1 + this->geo.M11 + 2 + 3];

        A = new Cell();
        A->Yzels.push_back(yz1);
        A->Yzels.push_back(yz2);
        A->Yzels.push_back(yz3);
        A->Yzels.push_back(yz4);
        this->All_Cell.push_back(A);
    }


}

void Setka::Print_yzel()
{
    ofstream fout;
    fout.open("All_yzel.txt");

    for (auto& i : this->All_Yzel)
    {
        fout << i->coord[0][0] << " " << i->coord[0][1] << " " << i->coord[0][2] << endl;
    }

}

void Setka::Print_cell()
{
    int ll = this->All_Cell.size();
    ofstream fout;
    fout.open("Setka_all_cell.txt");
    fout << "TITLE = \"HP\" ";
    fout << " VARIABLES = \"X\", \"Y\"  ZONE T= \"HP\", N=  " << 4 * ll;
    fout << " , E= " << 4 * ll;
    fout << " , F=FEPOINT, ET=LINESEG  " << endl;
    for (auto& i : this->All_Cell)
    {
        for (auto& j : i->Yzels)
        {
            fout << j->coord[0][0] << " " << j->coord[0][1] << endl;
        }
    }
    for (int i = 0; i < ll; i++)
    {
        fout << 4 * i + 1 << " " << 4 * i + 2 << endl;
        fout << 4 * i + 2 << " " << 4 * i + 3 << endl;
        fout << 4 * i + 3 << " " << 4 * i + 4 << endl;
        fout << 4 * i + 4 << " " << 4 * i + 1 << endl;
    }

    fout.close();
}
