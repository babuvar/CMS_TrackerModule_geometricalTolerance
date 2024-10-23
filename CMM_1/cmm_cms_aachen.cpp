//******************************************
//*Original Author: Varghese Babu***********
//*Tata Institute of Fundamental Research***
//******************************************
#include "TVector3.h"
#include "TVector2.h"
#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include <cstring>
#include <fstream>
using namespace std;
using namespace ROOT;
double Pi=acos(-1.0);
double Angle(TVector2 A,TVector2 B); 

void cmm_cms_aachen()
{
TVector3 UnitXt, UnitYt, UnitZt, UnitXb, UnitYb, UnitZb, Temp, L1t, L2t, L3t, L4t, L1b, L2b, L3b, L4b;
TVector2  DisT;




//X-Y Projections
//Measurement 1 of Aachen
TVector2 L1t_xy(94.222,102.753);
TVector2 L2t_xy(94.222,0.0);
TVector2 L3t_xy(0.0,0.0);
TVector2 L4t_xy(0.0,102.753);

TVector2 L1b_xy(94.23827893,102.80387716);
TVector2 L2b_xy(94.23546288,0.0508772);
TVector2 L3b_xy(0.01346291,0.05345946);
TVector2 L4b_xy(0.01627897,102.80645942);



/*
//Measurement 2 of Aachen
TVector2 L1t_xy(94.222,102.753);
TVector2 L2t_xy(94.222,0.0);
TVector2 L3t_xy(0.0,0.0);
TVector2 L4t_xy(0.0,102.753);

TVector2 L1b_xy(94.2360086,102.80471385);
TVector2 L2b_xy(94.23259649,0.0517139);
TVector2 L3b_xy(0.01059654,0.05484273);
TVector2 L4b_xy(0.01400865,102.80784267);
*/




//Displacements of L-marks
TVector2 Dis1=L1t_xy-L1b_xy;
TVector2 Dis2=L2t_xy-L2b_xy;
TVector2 Dis3=L3t_xy-L3b_xy;
TVector2 Dis4=L4t_xy-L4b_xy;


cout<<"Displacements of bottom L-marks w.r.t. top Lmarks, i.e., L_top - L_bot"<<endl;
cout<<"\tDel-x(μm)\tDel-y(μm)"<<endl;
cout<<"L-1:\t"<<Dis1.X()*1000<<"\t"<<Dis1.Y()*1000<<endl;
cout<<"L-2:\t"<<Dis2.X()*1000<<"\t"<<Dis2.Y()*1000<<endl;
cout<<"L-3:\t"<<Dis3.X()*1000<<"\t"<<Dis3.Y()*1000<<endl;
cout<<"L-4:\t"<<Dis4.X()*1000<<"\t"<<Dis4.Y()*1000<<endl<<endl<<endl;

//Full sensor Displacement
DisT=(Dis1+Dis2+Dis3+Dis4)*(0.25);
cout<<"Displacement of center of gravity of bottom sensor w.r.t. top sensor"<<endl;
cout<<"\tDel-x(μm)\tDel-y(μm)"<<endl;
cout<<"C.G.:\t"<<DisT.X()*1000<<"\t"<<DisT.Y()*1000<<endl;
//cout<<"Angle between sensors = "<<Angle((L1t_xy-L4t_xy),(L2b_xy-L3b_xy))*(180/Pi)<<" degrees (Rotation about sensor normal)"<<endl;//Degrees
//cout<<"Angle2 between sensors = "<<Angle((L1t_xy-L2t_xy),(L2b_xy-L1b_xy))*(180/Pi)<<" degrees (Rotation about sensor normal)"<<endl;//Degrees
//cout<<"Angle3 between sensors = "<<Angle((L1t_xy-L3t_xy),(L2b_xy-L4b_xy))*(180/Pi)<<" degrees (Rotation about sensor normal)"<<endl;//Degrees



//Rotations
cout<<"Angle between sensors = "<<Angle((L1t_xy-L4t_xy),(L1b_xy-L4b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L4"<<endl;//Radians
cout<<"Angle2 between sensors = "<<Angle((L1t_xy-L2t_xy),(L1b_xy-L2b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L2"<<endl;//Radians
cout<<"Angle3 between sensors = "<<Angle((L1t_xy-L3t_xy),(L1b_xy-L3b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L3"<<endl;//Radians

cout<<"--------------------------------------------------------------------------"<<endl;



}



double Angle(TVector2 A,TVector2 B)
{
double Mag_A, Mag_B, Dot_Prod, angle;
Mag_A= sqrt( ( A.X()*A.X() ) + ( A.Y()*A.Y() ));
Mag_B= sqrt( ( B.X()*B.X() ) + ( B.Y()*B.Y() ));
Dot_Prod = ( A.X()*B.X() ) + ( A.Y()*B.Y() );

angle = Dot_Prod/(Mag_A*Mag_B);
angle= acos(angle);
//angle= angle*(180/Pi);

return angle;

}








