#include "TVector3.h"
#include <vector>
#include "TVector3.h"
#include "TVector2.h"
#include "TMath.h"
#include <cstring>
#include <fstream>
using namespace std;
using namespace ROOT;
double Pi=acos(-1.0);
double Angle(TVector2 A, TVector2 B); 

void Test()
{
TVector2 A(1.0,0.0);
TVector2 B(1.0,1.0);
//TVector2 B(0.0,1.0);

cout<<"Angle = "<<Angle(A,B)<<endl;
}

double Angle(TVector2 A,TVector2 B)
{
double Mag_A, Mag_B, Dot_Prod, angle;
Mag_A= sqrt( ( A.X()*A.X() ) + ( A.Y()*A.Y() ));
Mag_B= sqrt( ( B.X()*B.X() ) + ( B.Y()*B.Y() ));
Dot_Prod = ( A.X()*B.X() ) + ( A.Y()*B.Y() );

angle = Dot_Prod/(Mag_A*Mag_B); 
angle= acos(angle);
angle= angle*(180/Pi);

return angle;

}
