#include <iostream>
#include <Windows.h>
#include "КЛАССЫ.h"
using namespace std;
//
unsigned short chislo_iteraciy = 0;
pair <double, double> Base_1_2_Mer_Algoritm(double a, double b, double epsilon, double r, unsigned short N, unsigned short m, double epsilon_granichnoe)
{
    pair <double, double> min_xy;
    if (N == 1)
    {
        HINSTANCE load_function = LoadLibrary(L"ФУНКЦИИ.dll");
        typedef double (*sh) (double);
        sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc");
        //
        pair<double, double>* start = new pair<double, double>;
        start->first = a;
        pair<double, double>* end = new pair<double, double>;
        end->first = b;
        //
        start->second = ShekelFunc(a);
        end->second = ShekelFunc(b);
        //
        Otrezok* otrezok = new Otrezok(start, end, N);
        //
        Otrezki* interval = new Otrezki(otrezok, r, N);
        //
        pair<double, double> pred_i_sled_shag;
        pred_i_sled_shag.first = a;
        pred_i_sled_shag.second = b;
        //
        while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
        {
            min_xy.first = pred_i_sled_shag.second;
            min_xy.second = ShekelFunc(min_xy.first);
            //
            pred_i_sled_shag.first = pred_i_sled_shag.second;
            //
            typedef double (*shag) (double, double, double, double, double, unsigned short);
            shag Shag = (shag)GetProcAddress(load_function, "Shag");
            pred_i_sled_shag.second = Shag(interval->Getm(), interval->GetX_Rmax().first, interval->GetX_Rmax().second, interval->GetY_Rmax().first, interval->GetY_Rmax().second, N);
            FreeLibrary(load_function);
            if (pred_i_sled_shag.second >= b)
            {
                srand(time(0));
                pred_i_sled_shag.second = (b - epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe / 2.0) - (b - epsilon_granichnoe));
            }
            if (pred_i_sled_shag.second <= a)
            {
                srand(time(0));
                pred_i_sled_shag.second = (a + epsilon_granichnoe / 2.0) + (double)rand() / RAND_MAX * ((a + epsilon_granichnoe) - (a + epsilon_granichnoe / 2.0));
            }
            //
            pair<double, double>* promejutochnaya_tochka = new pair<double, double>;
            promejutochnaya_tochka->first = pred_i_sled_shag.second;
            promejutochnaya_tochka->second = ShekelFunc(pred_i_sled_shag.second);
            interval->Add(promejutochnaya_tochka, N, -1.0, 0, -1, a, b, epsilon_granichnoe);
            chislo_iteraciy++;
        }
        return min_xy;
    }
    else
    {
        HINSTANCE load_function = LoadLibrary(L"ФУНКЦИИ.dll");
        typedef double (*grsh) (double, double);
        grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc");
        //
        Hilbert_Curve_2D Curve_2D = Hilbert_Curve_2D(a, b, a, b, m, Top);
        Hilbert_Curve_2D Curve_2D_PI_Na_Dva = Hilbert_Curve_2D(a, b, a, b, m, Left);
        Hilbert_Curve_2D Curve_2D_Minus_PI_Na_Dva = Hilbert_Curve_2D(a, b, a, b, m, Right);
        //
        pair<double, double>* start = new pair<double, double>;
        start->first = a;
        start->second = GrishaginFunc(Curve_2D.GetHead()->Get_x1(), Curve_2D.GetHead()->Get_x2());
        pair<double, double>* end = new pair<double, double>;
        end->first = b;
        end->second = GrishaginFunc(Curve_2D.GetEnd()->Get_x1(), Curve_2D.GetEnd()->Get_x2());
        pair<double, double>* start_PI_Na_Dva = new pair<double, double>;
        start_PI_Na_Dva->first = a;
        start_PI_Na_Dva->second = GrishaginFunc(Curve_2D_PI_Na_Dva.GetHead()->Get_x1(), Curve_2D_PI_Na_Dva.GetHead()->Get_x2());
        pair<double, double>* end_PI_Na_Dva = new pair<double, double>;
        end_PI_Na_Dva->first = b;
        end_PI_Na_Dva->second = GrishaginFunc(Curve_2D_PI_Na_Dva.GetEnd()->Get_x1(), Curve_2D_PI_Na_Dva.GetEnd()->Get_x2());
        pair<double, double>* start_Minus_PI_Na_Dva = new pair<double, double>;
        start_Minus_PI_Na_Dva->first = a;
        start_Minus_PI_Na_Dva->second = GrishaginFunc(Curve_2D_Minus_PI_Na_Dva.GetHead()->Get_x1(), Curve_2D_Minus_PI_Na_Dva.GetHead()->Get_x2());
        pair<double, double>* end_Minus_PI_Na_Dva = new pair<double, double>;
        end_Minus_PI_Na_Dva->first = b;
        end_Minus_PI_Na_Dva->second = GrishaginFunc(Curve_2D_Minus_PI_Na_Dva.GetEnd()->Get_x1(), Curve_2D_Minus_PI_Na_Dva.GetEnd()->Get_x2());
        //
        Otrezok* otrezok = new Otrezok(start, end, N);
        Otrezok* otrezok_PI_Na_Dva = new Otrezok(start_PI_Na_Dva, end_PI_Na_Dva, N);
        Otrezok* otrezok_Minus_PI_Na_Dva = new Otrezok(start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, N);
        //
        Otrezki* interval = new Otrezki(otrezok, r, N);
        Otrezki* interval_PI_Na_Dva = new Otrezki(otrezok_PI_Na_Dva, r, N);
        Otrezki* interval_Minus_PI_Na_Dva = new Otrezki(otrezok_Minus_PI_Na_Dva, r, N);
        //
        pair<double, double> pred_i_sled_shag;
        pred_i_sled_shag.first = a;
        pred_i_sled_shag.second = b;
        pair<double, double> pred_i_sled_shag_PI_Na_Dva;
        pred_i_sled_shag_PI_Na_Dva.first = a;
        pred_i_sled_shag_PI_Na_Dva.second = b;
        pair<double, double> pred_i_sled_shag_Minus_PI_Na_Dva;
        pred_i_sled_shag_Minus_PI_Na_Dva.first = a;
        pred_i_sled_shag_Minus_PI_Na_Dva.second = b;
        //
        unsigned short schetchick = 0;
        while (max(abs(pred_i_sled_shag.second - pred_i_sled_shag.first), max(abs(pred_i_sled_shag_PI_Na_Dva.second - pred_i_sled_shag_PI_Na_Dva.first), abs(pred_i_sled_shag_Minus_PI_Na_Dva.second - pred_i_sled_shag_Minus_PI_Na_Dva.first))) > epsilon)
        {
            if (pred_i_sled_shag.second == b)
            {
                min_xy.first = pred_i_sled_shag.second;
                min_xy.second = min(end->second, min(end_PI_Na_Dva->second, end_Minus_PI_Na_Dva->second));
            }
            else
            {
                unsigned int number = pred_i_sled_shag.second / ((b - a) / pow(2.0, m * N));
                if (number <= pow(2.0, m * N) / 2.0)
                {
                    CenterSquare* Curr = Curve_2D.GetHead();
                    while (number != 0)
                    {
                        number--;
                        Curr = Curr->GetNext();
                    }
                    min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                    min_xy.first = pred_i_sled_shag.second;
                }
                else
                {
                    CenterSquare* Curr = Curve_2D.GetEnd();
                    while (pow(2.0, m * N) - number - 1 != 0)
                    {
                        number++;
                        Curr = Curr->GetPrevious();
                    }
                    min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                    min_xy.first = pred_i_sled_shag.second;
                }
                number = pred_i_sled_shag_PI_Na_Dva.second / ((b - a) / pow(2.0, m * N));
                if (number <= pow(2.0, m * N) / 2.0)
                {
                    CenterSquare* Curr = Curve_2D_PI_Na_Dva.GetHead();
                    while (number != 0)
                    {
                        number--;
                        Curr = Curr->GetNext();
                    }
                    if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.second)
                    {
                        min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        min_xy.first = pred_i_sled_shag_PI_Na_Dva.second;
                    }
                }
                else
                {
                    CenterSquare* Curr = Curve_2D_PI_Na_Dva.GetEnd();
                    while (pow(2.0, m * N) - number - 1 != 0)
                    {
                        number++;
                        Curr = Curr->GetPrevious();
                    }
                    if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.second)
                    {
                        min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        min_xy.first = pred_i_sled_shag_PI_Na_Dva.second;
                    }
                }
                number = pred_i_sled_shag_Minus_PI_Na_Dva.second / ((b - a) / pow(2.0, m * N));
                if (number <= pow(2.0, m * N) / 2.0)
                {
                    CenterSquare* Curr = Curve_2D_Minus_PI_Na_Dva.GetHead();
                    while (number != 0)
                    {
                        number--;
                        Curr = Curr->GetNext();
                    }
                    if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.second)
                    {
                        min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        min_xy.first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
                    }
                }
                else
                {
                    CenterSquare* Curr = Curve_2D_Minus_PI_Na_Dva.GetEnd();
                    while (pow(2.0, m * N) - number - 1 != 0)
                    {
                        number++;
                        Curr = Curr->GetPrevious();
                    }
                    if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.second)
                    {
                        min_xy.second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        min_xy.first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
                    }
                }
            }
            //
            pred_i_sled_shag.first = pred_i_sled_shag.second;
            pred_i_sled_shag_PI_Na_Dva.first = pred_i_sled_shag_PI_Na_Dva.second;
            pred_i_sled_shag_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
            //
            typedef double (*shag) (double, double, double, double, double, unsigned short);
            shag Shag = (shag)GetProcAddress(load_function, "Shag");
            pred_i_sled_shag.second = Shag(interval->Getm(), interval->GetX_Rmax().first, interval->GetX_Rmax().second, interval->GetY_Rmax().first, interval->GetY_Rmax().second, N);
            pred_i_sled_shag_PI_Na_Dva.second = Shag(interval_PI_Na_Dva->Getm(), interval_PI_Na_Dva->GetX_Rmax().first, interval_PI_Na_Dva->GetX_Rmax().second, interval_PI_Na_Dva->GetY_Rmax().first, interval_PI_Na_Dva->GetY_Rmax().second, N);
            pred_i_sled_shag_Minus_PI_Na_Dva.second = Shag(interval_Minus_PI_Na_Dva->Getm(), interval_Minus_PI_Na_Dva->GetX_Rmax().first, interval_Minus_PI_Na_Dva->GetX_Rmax().second, interval_Minus_PI_Na_Dva->GetY_Rmax().first, interval_Minus_PI_Na_Dva->GetY_Rmax().second, N);
            FreeLibrary(load_function);
            if (pred_i_sled_shag.second >= b)
            {
                srand(time(0));
                pred_i_sled_shag.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
            }
            if (pred_i_sled_shag.second <= a)
            {
                srand(time(0));
                pred_i_sled_shag.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
            }
            if (pred_i_sled_shag_PI_Na_Dva.second >= b)
            {
                srand(time(0));
                pred_i_sled_shag_PI_Na_Dva.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
            }
            if (pred_i_sled_shag_PI_Na_Dva.second <= a)
            {
                srand(time(0));
                pred_i_sled_shag_PI_Na_Dva.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
            }
            if (pred_i_sled_shag_Minus_PI_Na_Dva.second >= b)
            {
                srand(time(0));
                pred_i_sled_shag_Minus_PI_Na_Dva.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
            }
            if (pred_i_sled_shag_Minus_PI_Na_Dva.second <= a)
            {
                srand(time(0));
                pred_i_sled_shag_Minus_PI_Na_Dva.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
            }
            //
            pair<double, double>* promejutochnaya_tochka = new pair<double, double>;
            promejutochnaya_tochka->first = pred_i_sled_shag.second;
            pair<double, double>* promejutochnaya_tochka_PI_Na_Dva = new pair<double, double>;
            promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second;
            pair<double, double>* promejutochnaya_tochka_Minus_PI_Na_Dva = new pair<double, double>;
            promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
            unsigned short flag = 0;
            CenterSquare* Curr;
            CenterSquare* Curr1;
            CenterSquare* Curr2;
            unsigned int number = (pred_i_sled_shag.second) / ((b - a) / pow(2.0, m * N));
            pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
            if (number <= pow(2.0, m * N) / 2.0)
            {
                Curr = Curve_2D.GetHead();
                while (number != 0)
                {
                    number--;
                    Curr = Curr->GetNext();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                }
                else
                {
                    promejutochnaya_tochka->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                }
            }
            else
            {
                Curr = Curve_2D.GetEnd();
                while (pow(2.0, m * N) - number - 1 != 0)
                {
                    number++;
                    Curr = Curr->GetPrevious();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                }
                else
                {
                    promejutochnaya_tochka->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                }
            }
            number = (pred_i_sled_shag_PI_Na_Dva.second) / ((b - a) / pow(2.0, m * N));
            pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
            if (number <= pow(2.0, m * N) / 2.0)
            {
                Curr1 = Curve_2D_PI_Na_Dva.GetHead();
                while (number != 0)
                {
                    number--;
                    Curr1 = Curr1->GetNext();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    if (GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2()) < promejutochnaya_tochka->second)
                    {
                        promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                        flag = 1;
                    }
                }
                else
                {
                    promejutochnaya_tochka_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                }
            }
            else
            {
                Curr1 = Curve_2D_PI_Na_Dva.GetEnd();
                while (pow(2.0, m * N) - number - 1 != 0)
                {
                    number++;
                    Curr1 = Curr1->GetPrevious();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    if (GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2()) < promejutochnaya_tochka->second)
                    {
                        promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                        flag = 1;
                    }
                }
                else
                {
                    promejutochnaya_tochka_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                }
            }
            number = (pred_i_sled_shag_Minus_PI_Na_Dva.second) / ((b - a) / pow(2.0, m * N));
            pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
            if (number <= pow(2.0, m * N) / 2.0)
            {
                Curr2 = Curve_2D_Minus_PI_Na_Dva.GetHead();
                while (number != 0)
                {
                    number--;
                    Curr2 = Curr2->GetNext();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    if (GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2()) < promejutochnaya_tochka->second)
                    {
                        promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                        flag = 2;
                    }
                }
                else
                {
                    promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                }
            }
            else
            {
                Curr2 = Curve_2D_Minus_PI_Na_Dva.GetEnd();
                while (pow(2.0, m * N) - number - 1 != 0)
                {
                    number++;
                    Curr2 = Curr2->GetPrevious();
                }
                if (schetchick % 10 == 0 && schetchick != 0)
                {
                    if (GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2()) < promejutochnaya_tochka->second)
                    {
                        promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                        flag = 2;
                    }
                }
                else
                {
                    promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                }
            }
            if (schetchick % 10 == 0 && schetchick != 0)
            {
                if (flag == 0)
                {
                    promejutochnaya_tochka->first = pred_i_sled_shag.second;
                    number = 0;
                    Curr1 = Curve_2D_PI_Na_Dva.GetHead();
                    while (Curr->Get_x1() != Curr1->Get_x1() && Curr->Get_x2() != Curr1->Get_x2())
                    {
                        number++;
                        Curr1 = Curr1->GetNext();
                    }
                    promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                    number = 0;
                    Curr2 = Curve_2D_Minus_PI_Na_Dva.GetHead();
                    while (Curr->Get_x1() != Curr2->Get_x1() && Curr->Get_x2() != Curr2->Get_x2())
                    {
                        number++;
                        Curr2 = Curr2->GetNext();
                    }
                    promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                }
                if (flag == 1)
                {
                    promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second;
                    number = 0;
                    Curr = Curve_2D.GetHead();
                    while (Curr1->Get_x1() != Curr->Get_x1() && Curr1->Get_x2() != Curr->Get_x2())
                    {
                        number++;
                        Curr = Curr->GetNext();
                    }
                    promejutochnaya_tochka->first = pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
                    number = 0;
                    Curr2 = Curve_2D_Minus_PI_Na_Dva.GetHead();
                    while (Curr1->Get_x1() != Curr2->Get_x1() && Curr1->Get_x2() != Curr2->Get_x2())
                    {
                        number++;
                        Curr2 = Curr2->GetNext();
                    }
                    promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                }
                else
                {
                    promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
                    number = 0;
                    Curr = Curve_2D.GetHead();
                    while (Curr2->Get_x1() != Curr->Get_x1() && Curr2->Get_x2() != Curr->Get_x2())
                    {
                        number++;
                        Curr = Curr->GetNext();
                    }
                    promejutochnaya_tochka->first = pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
                    number = 0;
                    Curr1 = Curve_2D_PI_Na_Dva.GetHead();
                    while (Curr2->Get_x1() != Curr1->Get_x1() && Curr2->Get_x2() != Curr1->Get_x2())
                    {
                        number++;
                        Curr1 = Curr1->GetNext();
                    }
                    promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                }
            }
            interval->Add(promejutochnaya_tochka, N, -1.0, 0, -1, a, b, epsilon_granichnoe);
            interval_PI_Na_Dva->Add(promejutochnaya_tochka_PI_Na_Dva, N, -1.0, 0, -1, a, b, epsilon_granichnoe);
            interval_Minus_PI_Na_Dva->Add(promejutochnaya_tochka_Minus_PI_Na_Dva, N, -1.0, 0, -1, a, b, epsilon_granichnoe);
            schetchick++;
            chislo_iteraciy++;
        }
        return min_xy;
    }
}
//
int main()
{
    //pair <double, double> Extr = Base_1_2_Mer_Algoritm(0.0, 1.0, 0.000005, 4.0, 2, 11, pow(10.0, -6.0));
    pair <double, double> Extr = Base_1_2_Mer_Algoritm(0.0, 10.0, 0.00005, 2.0, 1, -1, pow(10.0, -6.0));
    cout << "Xmin = " << Extr.first << " " << "Ymin = " << Extr.second << endl;
    cout << chislo_iteraciy << endl;
}


