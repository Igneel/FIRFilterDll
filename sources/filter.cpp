#include "filter.h"


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
// ВНИМАНИЕ!!! Напрямую не вызывать!!! Пользоваться трамплином!!!---------------
double Filter (const double in[], double out[], int sizeIn, int length, double Fdisk, double Fpropysk,double Fzatyh)
{
int N = length; //Длина фильтра 20
long double Fd = Fdisk; //Частота дискретизации входных данных 2000
long double Fs = Fpropysk; //Частота конца полосы пропускания  10
long double Fx = Fzatyh; //Частота начала полосы затухания    20

long double *H= new long double [N] ; //Импульсная характеристика фильтра
long double *H_id= new long double [N]; //Идеальная импульсная характеристика
long double *W= new long double [N]; //Весовая функция

/*long double H [N] = {0}; //Импульсная характеристика фильтра
long double H_id [N] = {0}; //Идеальная импульсная характеристика
long double W [N] = {0}; //Весовая функция  */

//Расчет импульсной характеристики фильтра
long double Fc = (Fs + Fx) / (2 * Fd);

for (int i=0;i<N;i++)
{
if (i==0) H_id[i] = 2*M_PI*Fc;
else H_id[i] = sinl(2*M_PI*Fc*i )/(M_PI*i);
// весовая функция Блекмена
W [i] = 0.42 - 0.5 * cosl((2*M_PI*i) /( N-1)) + 0.08 * cosl((4*M_PI*i) /( N-1));
H [i] = H_id[i] * W[i];
}

//Нормировка импульсной характеристики
long double SUM=0;
for (int i=0; i<N; i++) SUM +=H[i];
for (int i=0; i<N; i++) H[i]/=SUM; //сумма коэффициентов равна 1

//----------------------------------------------------------------
// печать коэффициентов фильтра
//for(int i=0;i<N;i++)
//Form1->Memo2->Lines->Add(FloatToStr(H[i]));
//----------------------------------------------------------------

//Фильтрация входных данных
for (int i=0; i<sizeIn; i++)
{
out[i]=0.;
for (int j=0; j<(i>N-1?N-1:i); j++)// та самая формула фильтра
out[i]+= H[j]*in[i-j];
}
delete [] H;
delete [] H_id;
delete [] W;
return (N-1)/2.0;
}

// ВНИМАНИЕ!!! Напрямую не вызывать!!! Пользоваться трамплином!!!---------------
double Filter (const std::vector<long double> &in, std::vector<long double> & out, int length, double Fdisk, double Fpropysk,double Fzatyh)
{
int N = length; //Длина фильтра
long double Fd = Fdisk; //Частота дискретизации входных данных 2000
long double Fs = Fpropysk; //Частота конца полосы пропускания  10
long double Fx = Fzatyh; //Частота начала полосы затухания    20

std::vector<long double> H(N);  //Импульсная характеристика фильтра
std::vector<long double> H_id(N); //Идеальная импульсная характеристика
std::vector<long double> W(N);   //Весовая функция

//Расчет импульсной характеристики фильтра
long double Fc = (Fs + Fx) / (2.0 * Fd);

for (int i=0;i<N;i++)
{
if (i==0) H_id[i] = 2.0*M_PI*Fc;
else H_id[i] = sinl(2.0*M_PI*Fc*i )/(M_PI*i);
// весовая функция Блекмена
W [i] = 0.42 - 0.5 * cosl((2*M_PI*i) /( N-1)) + 0.08 * cosl((4*M_PI*i) /( N-1.0));
H [i] = H_id[i] * W[i];
}

//Нормировка импульсной характеристики
long double SUM=0;
for (int i=0; i<N; i++) SUM +=H[i];
for (int i=0; i<N; i++) H[i]/=SUM; //сумма коэффициентов равна 1

//Фильтрация входных данных
for (int i=0; i<in.size(); i++)
{
out[i]=0.0;
for (int j=0; j<(i>N-1?N-1:i); j++)// та самая формула фильтра
out[i]+= H[j]*in[i-j];
}

return (N-1)/2.0;
}


// функция трамплин для фильтра.  ----------------------------------------------
double TrForMassiveFilter(long double *inB,long double *inY,long double* outB,long double *outY,
int lengthMassive,int lengthFilter,double Fdisk, double Fpropysk,double Fzatyh)
{

// перед фильтрацией делаем массив нечетной функцией.


//int size=inS->YValues->Count;  // получаем размер
double *in=new double[lengthMassive];  // выделяем память
for(int i=0;i<lengthMassive;i++)       // копируем
{
    in[i]=inY[i];
}
double *out=new double[lengthMassive]; // выделяем память для выходных значений
double k=Filter(in,out,lengthMassive,lengthFilter,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
 k*=(inB[lengthMassive-1]-inB[0])/(double)lengthMassive;// вычисляем сдвиг фаз
 // если что тут максимум и минимум надо бы вычислять.

//----------------------------------------------
//---------добавление для фильтрации магнитного поля
// внимание - она пока работает не корректно ибо фильтр наполняется только за
// его длину.
/*double *in2=new double[size];  // выделяем память
for(int i=0;i<size;i++)       // копируем
in2[i]=inS->XValues->Value[i];
double *out2=new double[size]; // выделяем память для выходных значений
double k2=Filter(in2,out2,size,length,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
k2*=(inS->XValues->MaxValue-inS->XValues->MinValue)/(double)inS->XValues->Count;// вычисляем сдвиг фаз
*/

//----------------------------------------------
for(int i=0;i<lengthMassive;i++) // выводим
{
    outB[i]=inB[i]-k;
    outY[i]=out[i];
}

delete[] in;  // прибираемся
//delete[] in2;  // прибираемся
delete[] out;
//delete[] out2;
return k;
}
//-----------------------------------------------
// функция трамплин для фильтра.  ----------------------------------------------
long double TrForMassiveFilter(std::vector<long double> & inB,
std::vector<long double> & inY,std::vector<long double> & outB,
std::vector<long double> & outY,int lengthFilter,long double Fdisk,
long double Fpropysk,long double Fzatyh)
{
    int lengthMassive=inY.size();
    if(lengthMassive==0)
    {
    return 0;
    }
    std::vector<long double> in(inY);
    std::vector<long double> out(lengthMassive);

    double k=Filter(in,out,lengthFilter,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
    k*=(max_elem(inB)-min_elem(inB))/(double)lengthMassive;// вычисляем сдвиг фаз
    // разность максимума и минимума на длину массива

//----------------------------------------------
//---------добавление для фильтрации магнитного поля
// внимание - она пока работает не корректно ибо фильтр наполняется только за
// его длину.
/*double *in2=new double[size];  // выделяем память
for(int i=0;i<size;i++)       // копируем
in2[i]=inS->XValues->Value[i];
double *out2=new double[size]; // выделяем память для выходных значений
double k2=Filter(in2,out2,size,length,Fdisk,Fpropysk,Fzatyh); // вызываем фильтр
k2*=(inS->XValues->MaxValue-inS->XValues->MinValue)/(double)inS->XValues->Count;// вычисляем сдвиг фаз
*/

//----------------------------------------------
    outB.resize(lengthMassive);
    outY=out;
    for(int i=0;i<lengthMassive;i++) // выводим
    {
        outB[i]=inB[i]-k;
    }
    return k;
}


//----------------------------------------------
//----------------------------------------------

// Попытка фильтрации оптимальным методом.

long double D(long double deltaP, long double deltaS)
{
/*
Нужна для расчета длины фильтра.
*/
    long double a1=5.309E-3;
    long double a2=7.114E-2;
    long double a3=-4.761E-1;
    long double a4=-2.66E-3;
    long double a5=-5.941E-1;
    long double a6=-4.278E-1;

return log10(deltaS*(a1*log10(deltaP)*log10(deltaP)+a2*
    log10(deltaP)+a3)+a4*log10(deltaP)*log10(deltaP)+a5*log10(deltaP)+a6);
}

long double f(long double deltaP, long double deltaS)
{
/*
Нужна для расчета длины фильтра.
*/

return 11.01217+0.51244*(log10(deltaP)-log10(deltaS));

}

int calcutaleFilterLength(long double deltaP, long double deltaS, long double deltaF)
{
/*
Рассчитывает длину фильтра для метода оптимальных коэффициентов.
Согласно формуле из Айфичера Джервиса для фильтра нижних частот, страница 410.

deltaP - неравномерность в полосе пропускания.

deltaS - неравномерность в полосе подавления.

deltaF - ширина полосы пропускания, нормированная на частоту дискретизации.

N=D(deltaP,deltaS)/deltaF-f(deltaP,deltaS)*deltaF+1;

*/

return D(deltaP,deltaS)/deltaF-f(deltaP,deltaS)*deltaF+1;
}

long double fapprox(long double x)
{
return x*x;
}

double OptimFilter(std::vector<double> in, std::vector<double> out,
        long double deltaP, long double deltaS, long double FPropusk, long double Fpodavl, long double FDiskr)
{
long double FPropuskN=FPropusk/FDiskr;
long double FPodavlN=Fpodavl/FDiskr;
long double deltaF=FPodavlN-FPropuskN;

int filterLength=calcutaleFilterLength(deltaP,deltaS,deltaF);

double weigthP = 1;
double weigthS = deltaS/deltaP;

boost::math::tools::remez_minimax<double> *remez= new boost::math::tools::remez_minimax<double>
    ( fapprox, weigthP, weigthS, FPropuskN, FPodavlN);
    // функция для аппроксимации.
    // Order of the numerator polynomial.
    //Order of the denominator polynomial.
    // Крайние точки диапазона оптимизации. End points of the range to optimise over.
boost::numeric::ublas::vector<double> cheb_points;
cheb_points=remez->chebyshev_points();


//Фильтрация входных данных
for (int i=0; i<cheb_points.size(); i++)
{
out[i]=0.;
for (int j=0; j<(i>filterLength-1?filterLength-1:i); j++)// та самая формула фильтра
    out[i]+= cheb_points[j]*in[i-j];
}

return (cheb_points.size()-1) /2.0;
//filter();

}
