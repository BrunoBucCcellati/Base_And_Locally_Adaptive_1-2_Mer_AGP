using namespace std;
//
class Otrezok
{
protected:
	pair<double, double>* start;
	pair<double, double>* end;
	double M;
	double R;
	Otrezok* Next;
	Otrezok* Previous;
	double m;
public:
	Otrezok(pair<double, double>* _start, pair<double, double>* _end, unsigned short _N)
	{
		start = _start;
		end = _end;
		M = abs((end->second - start->second) / pow(abs(end->first - start->first), (1 / double(_N))));
		R = 0.0;
		Next = this;
		Previous = this;
	};
	void ChangeCharacteristic(double _m, unsigned short _N);
	double GetCharacteristic();
	void SetEnd(pair<double, double>* _end, unsigned short _N);
	pair<double, double>* GetEnd();
	pair<double, double>* GetStart();
	Otrezok* GetNext();
	void SetNext(Otrezok* _Next);
	Otrezok* GetPrevious();
	void SetPrevious(Otrezok* _Previous);
	double GetM();
	void SetM(double _M);
};
//
enum List { Top, Dawn, Left, Right };
class CenterSquare
{
protected:
	pair<double, double> x1x2;
	CenterSquare* DawnLeft;
	CenterSquare* TopLeft;
	CenterSquare* TopRight;
	CenterSquare* DawnRight;
	List Type;
	CenterSquare* Start;
	CenterSquare* End;
	CenterSquare* Next;
	CenterSquare* Previous;
public:
	CenterSquare(double _x1, double _x2, List _Type, CenterSquare* _Next, CenterSquare* _Previous, unsigned short i, double h1, double h2)
	{
		x1x2.first = _x1;
		x1x2.second = _x2;
		Type = _Type;
		if (_Next == nullptr)
		{
			Next = this;
		}
		else
		{
			Next = _Next;
		}
		if (_Previous == nullptr)
		{
			Previous = this;
		}
		else
		{
			Previous = _Previous;
		}
		if (i != 0)
		{
			i--;
			if (Type == Right)
			{
				TopLeft = new CenterSquare(_x1 - h1, _x2 + h2, Dawn, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopRight = new CenterSquare(_x1 + h1, _x2 + h2, Right, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnRight = new CenterSquare(_x1 + h1, _x2 - h2, Right, TopRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnLeft = new CenterSquare(_x1 - h1, _x2 - h2, Top, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopLeft->SetPrevious(TopRight);
				TopRight->SetPrevious(DawnRight);
				DawnRight->SetPrevious(DawnLeft);
				Start = DawnLeft;
				End = TopLeft;
			}
			if (Type == Top)
			{
				DawnRight = new CenterSquare(_x1 + h1, _x2 - h2, Left, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopRight = new CenterSquare(_x1 + h1, _x2 + h2, Top, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopLeft = new CenterSquare(_x1 - h1, _x2 + h2, Top, TopRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnLeft = new CenterSquare(_x1 - h1, _x2 - h2, Right, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnRight->SetPrevious(TopRight);
				TopRight->SetPrevious(TopLeft);
				TopLeft->SetPrevious(DawnLeft);
				Start = DawnLeft;
				End = DawnRight;
			}
			if (Type == Left)
			{
				DawnRight = new CenterSquare(_x1 + h1, _x2 - h2, Top, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnLeft = new CenterSquare(_x1 - h1, _x2 - h2, Left, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopLeft = new CenterSquare(_x1 - h1, _x2 + h2, Left, DawnLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopRight = new CenterSquare(_x1 + h1, _x2 + h2, Dawn, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnRight->SetPrevious(DawnLeft);
				DawnLeft->SetPrevious(TopLeft);
				TopLeft->SetPrevious(TopRight);
				Start = TopRight;
				End = DawnRight;
			}
			if (Type == Dawn)
			{
				TopLeft = new CenterSquare(_x1 - h1, _x2 + h2, Right, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnLeft = new CenterSquare(_x1 - h1, _x2 - h2, Dawn, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				DawnRight = new CenterSquare(_x1 + h1, _x2 - h2, Dawn, DawnLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopRight = new CenterSquare(_x1 + h1, _x2 + h2, Left, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
				TopLeft->SetPrevious(DawnLeft);
				DawnLeft->SetPrevious(DawnRight);
				DawnRight->SetPrevious(TopRight);
				Start = TopRight;
				End = TopLeft;
			}
		}
	}
	double Get_x1();
	double Get_x2();
	List GetType();
	CenterSquare* GetStart();
	CenterSquare* GetEnd();
	CenterSquare* GetNext();
	CenterSquare* GetPrevious();
	void SetNext(CenterSquare* _Next);
	void SetPrevious(CenterSquare* _Previous);
};
//
class Hilbert_Curve_2D
{
protected:
	CenterSquare* Head;
	CenterSquare* End;
public:
	Hilbert_Curve_2D(double a = 0.0, double b = 0.0, double c = 0.0, double d = 0.0, unsigned short m = 1, List Type = Top)
	{
		Head = new CenterSquare((b - a) / 2.0, (d - c) / 2.0, Type, nullptr, nullptr, 1, (b - a) / 4.0, (d - c) / 4.0);
		End = Head->GetEnd();
		Head = Head->GetStart();
		m--;
		double h1 = (b - a) / 4.0;
		double h2 = (d - c) / 4.0;
		while (m != 0)
		{
			m--;
			h1 = h1 / 2.0;
			h2 = h2 / 2.0;
			CenterSquare* Curr = Head;
			CenterSquare* Current = new CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
			Curr = Curr->GetNext();
			Head = Current->GetStart();
			while (Curr->GetNext() != Curr)
			{
				CenterSquare* Current1 = new CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
				Current1->GetStart()->SetPrevious(Current->GetEnd());
				Current->GetEnd()->SetNext(Current1->GetStart());
				Current = Current1;
				Curr = Curr->GetNext();
			}
			CenterSquare* Current1 = new CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
			Current1->GetStart()->SetPrevious(Current->GetEnd());
			Current->GetEnd()->SetNext(Current1->GetStart());
			Current = Current1;
			End = Current->GetEnd();
		}
	}
	CenterSquare* GetHead();
	CenterSquare* GetEnd();
};
//
class Otrezki
{
protected:
	Otrezok* Head;
	double Mmax;
	double m;
	double r;
	pair<double, double> x_Rmax;
	pair<double, double> y_Rmax;
	double dmax;
public:
	Otrezki(Otrezok* _Head, double _r, unsigned short _N)
	{
		Head = _Head;
		r = _r;
		if (Head->GetM() != 0)
		{
			Mmax = Head->GetM();
			m = r * Mmax;
		}
		else
		{
			Mmax = 0;
			m = 1;
		}
		Head->ChangeCharacteristic(m, _N);
		x_Rmax.first = Head->GetStart()->first;
		x_Rmax.second = Head->GetEnd()->first;
		y_Rmax.first = Head->GetStart()->second;
		y_Rmax.second = Head->GetEnd()->second;
	};
	void Add(pair<double, double>* _tmp, unsigned short _N, double _eta_0, unsigned short _global_local_iterations, unsigned short _K, double _a, double _b, double _epsilon_granichnoe);
	double Getm();
	pair <double, double> GetX_Rmax();
	pair <double, double> GetY_Rmax();
};