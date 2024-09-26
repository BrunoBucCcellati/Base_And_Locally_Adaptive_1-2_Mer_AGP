#include<algorithm>
#include " À¿——€.h"
#include <Windows.h>
//
void Otrezok::ChangeCharacteristic(double _m, unsigned short _N)
{
    HINSTANCE load_function = LoadLibrary(L"‘”Õ ÷»».dll");
    typedef double (*characteristic) (double, double, double, double, double, unsigned short);
    characteristic Characteristic = (characteristic)GetProcAddress(load_function, "Characteristic");
    R = Characteristic(_m, start->first, end->first, start->second, end->second, _N);
    FreeLibrary(load_function);
}
double Otrezok::GetCharacteristic()
{
    return R;
}
void Otrezok::SetEnd(pair<double, double>* _end, unsigned short _N)
{
    end = _end;
    M = abs((end->second - start->second) / pow(abs(end->first - start->first), (1 / double(_N))));
}
pair<double, double>* Otrezok::GetEnd()
{
    return end;
}
pair<double, double>* Otrezok::GetStart()
{
    return start;
}
Otrezok* Otrezok::GetNext()
{
    return Next;
}
void Otrezok::SetNext(Otrezok* _Next)
{
    Next = _Next;
}
Otrezok* Otrezok::GetPrevious()
{
    return Previous;
}
void Otrezok::SetPrevious(Otrezok* _Previous)
{
    Previous = _Previous;
}
double Otrezok::GetM()
{
    return M;
}
void Otrezok::SetM(double _M)
{
    M = _M;
}
//
void Otrezki::Add(pair<double, double>* tmp, unsigned short _N, double _eta_0, unsigned short _global_local_iterations, unsigned short _K, double _a, double _b, double _epsilon_granichnoe)
{
    Otrezok* curr = Head;
    while (curr->GetEnd()->first < tmp->first)
    {
        curr = curr->GetNext();
    }
    Otrezok* curr1 = new Otrezok(tmp, curr->GetEnd(), _N);
    curr->SetEnd(tmp, _N);
    curr1->SetPrevious(curr);
    curr->GetNext()->SetPrevious(curr1);
    curr1->SetNext(curr->GetNext());
    curr->SetNext(curr1);
    //
    if (_K >= _global_local_iterations && _K % 2 == 0)
    {
        if (curr->GetM() >= curr1->GetM())
        {
            double eta_shtrih = max(curr->GetM(), max(curr->GetPrevious()->GetM(), curr->GetNext()->GetM()));
            double eta_2shtrih = Mmax * pow((curr->GetEnd()->first - curr->GetStart()->first), double(1 / _N)) / dmax;
            m = max((1.0 / r) * eta_shtrih + ((r - 1.0) / r) * eta_2shtrih, _eta_0) * r;
            x_Rmax.first = curr->GetStart()->first;
            x_Rmax.second = curr->GetEnd()->first;
            y_Rmax.first = curr->GetStart()->second;
            y_Rmax.second = curr->GetEnd()->second;
        }
        else
        {
            double eta_shtrih = max(curr1->GetM(), max(curr1->GetPrevious()->GetM(), curr1->GetNext()->GetM()));
            double eta_2shtrih = Mmax * pow((curr1->GetEnd()->first - curr1->GetStart()->first), double(1 / _N)) / dmax;
            m = max((1.0 / r) * eta_shtrih + ((r - 1.0) / r) * eta_2shtrih, _eta_0) * r;
            x_Rmax.first = curr1->GetStart()->first;
            x_Rmax.second = curr1->GetEnd()->first;
            y_Rmax.first = curr1->GetStart()->second;
            y_Rmax.second = curr1->GetEnd()->second;
        }
    }
    else
    {
        curr = Head;
        Otrezok* Otrezok_Rmax = curr;
        Mmax = curr->GetM();
        curr = curr->GetNext();
        while (curr != Head)
        {
            if (curr->GetM() > Mmax)
            {
                Mmax = curr->GetM();
                Otrezok_Rmax = curr;
            }
            curr = curr->GetNext();
        }
        Otrezok_Rmax->SetM(0.0);
        if (Mmax != 0)
        {
            m = r * Mmax;
        }
        else
        {
            m = 1;
        }
        //
        curr = Head;
        while (curr->GetNext() != Head)
        {
            curr = curr->GetNext();
            curr->ChangeCharacteristic(m, _N);
        }
        curr = curr->GetNext();
        curr->ChangeCharacteristic(m, _N);
        //
        double Rmax = -DBL_MAX;
        HINSTANCE load_function = LoadLibrary(L"‘”Õ ÷»».dll");
        typedef double (*shag) (double, double, double, double, double, unsigned short);
        shag Shag = (shag)GetProcAddress(load_function, "Shag");
        if (Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second, _N) <= _b - _epsilon_granichnoe && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second, _N) >= _a + _epsilon_granichnoe)
        {
            Rmax = curr->GetCharacteristic();
            x_Rmax.first = curr->GetStart()->first;
            x_Rmax.second = curr->GetEnd()->first;
            y_Rmax.first = curr->GetStart()->second;
            y_Rmax.second = curr->GetEnd()->second;
        }
        dmax = pow(curr->GetEnd()->first - curr->GetStart()->first, double(1 / _N));
        curr = curr->GetNext();
        while (curr != Head)
        {
            if (curr->GetCharacteristic() > Rmax && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second, _N) <= _b - _epsilon_granichnoe && Shag(m, curr->GetStart()->first, curr->GetEnd()->first, curr->GetStart()->second, curr->GetEnd()->second, _N) >= _a + _epsilon_granichnoe)
            {
                Rmax = curr->GetCharacteristic();
                x_Rmax.first = curr->GetStart()->first;
                x_Rmax.second = curr->GetEnd()->first;
                y_Rmax.first = curr->GetStart()->second;
                y_Rmax.second = curr->GetEnd()->second;
            }
            if (pow(curr->GetEnd()->first - curr->GetStart()->first, double(1 / _N)) > dmax)
            {
                dmax = pow(curr->GetEnd()->first - curr->GetStart()->first, double(1 / _N));
            }
            curr = curr->GetNext();
        }
        FreeLibrary(load_function);
    }
}
double Otrezki::Getm()
{
    return m;
}
pair <double, double> Otrezki::GetX_Rmax()
{
    return x_Rmax;
}
pair <double, double> Otrezki::GetY_Rmax()
{
    return y_Rmax;
}
//
double CenterSquare::Get_x1()
{
    return x1x2.first;
}
double CenterSquare::Get_x2()
{
    return x1x2.second;
}
List CenterSquare::GetType()
{
    return Type;
}
CenterSquare* CenterSquare::GetStart()
{
    return Start;
}
CenterSquare* CenterSquare::GetEnd()
{
    return End;
}
CenterSquare* CenterSquare::GetNext()
{
    return Next;
}
CenterSquare* CenterSquare::GetPrevious()
{
    return Previous;
}
void CenterSquare::SetNext(CenterSquare* _Next)
{
    Next = _Next;
}
void CenterSquare::SetPrevious(CenterSquare* _Previous)
{
    Previous = _Previous;
}
//
CenterSquare* Hilbert_Curve_2D::GetHead()
{
    return Head;
}
CenterSquare* Hilbert_Curve_2D::GetEnd()
{
    return End;
}