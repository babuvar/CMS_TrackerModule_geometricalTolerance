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

void cmm_cms_mod()
{
TVector3 UnitXt, UnitYt, UnitZt, UnitXb, UnitYb, UnitZb, Temp, L1t, L2t, L3t, L4t, L1b, L2b, L3b, L4b, Dis1, Dis2, Dis3, Dis4, DisT;
vector<TVector3> abc_top, abc_bot, rand_top, rand_bot;
Double_t x,y,z,one_third=(1.0/3.0); 
cout<<"one_third = "<<one_third<<endl;
//ABC Locations
//Top
TVector3 At(447.011,581.731,20.188);
TVector3 Bt(494.136,478.994,20.037);
TVector3 Ct(540.649,479.282,20.183);

//Bottom
TVector3 Ab(523.597,598.301,21.587);
TVector3 Bb(477.695,495.026,21.770);
TVector3 Cb(431.184,494.749,21.657);

//Making Local Axes (Top)
Temp=Ct-Bt;
UnitXt= Temp.Unit();
Temp=At-Bt;
Temp=Temp-((Temp.Dot(UnitXt))*UnitXt);
UnitYt= Temp.Unit();
UnitZt= UnitXt.Cross(UnitYt);
UnitZt= UnitZt.Unit();


//Making Local Axes (Bottom)
Temp=Cb-Bb;
UnitXb= Temp.Unit();
Temp=Ab-Bb;
Temp=Temp-((Temp.Dot(UnitXb))*UnitXb);
UnitYb= Temp.Unit();
UnitZb = UnitXb.Cross(UnitYb);
UnitZb= UnitZb.Unit();

/*
//Check 
UnitXt.Print();
UnitYt.Print();
UnitZt.Print();
UnitXb.Print();
UnitYb.Print();
UnitZb.Print(); 
*/



//------------------------------------------------------------------------
ifstream fin;
fin.open("abc_top.txt");
while(!fin.eof())
{
fin>>x>>y>>z; Temp.SetXYZ(x,y,z);
abc_top.push_back(Temp);
}
fin.close();
//for(int i=0;i<12;i++){cout<<"i = "<<i<<"\t";abc_top.at(i).Print();}

//------------------------------------------------------------------------
fin.open("abc_bot.txt");
while(!fin.eof())
{
fin>>x>>y>>z; Temp.SetXYZ(x,y,z);
abc_bot.push_back(Temp);
}
fin.close();
//for(int i=0;i<12;i++){cout<<"i = "<<i<<"\t";abc_bot.at(i).Print();}

//------------------------------------------------------------------------

int numrand_top=-1,numrand_bot=-1;

ifstream fin;
fin.open("rand_top.txt");
//fin.open("abc_top.txt");
while(!fin.eof())
{numrand_top++;
fin>>x>>y>>z; Temp.SetXYZ(x,y,z);
rand_top.push_back(Temp);
}
fin.close();
//for(int i=0;i<numrand_top;i++){cout<<"i = "<<i<<"\t";rand_top.at(i).Print();}

//------------------------------------------------------------------------
ifstream fin;
fin.open("rand_bot.txt");
//fin.open("abc_bot.txt");

while(!fin.eof())
{numrand_bot++;
fin>>x>>y>>z; Temp.SetXYZ(x,y,z);
rand_bot.push_back(Temp);
}
fin.close();
//for(int i=0;i<numrand_bot;i++){cout<<"i = "<<i<<"\t";rand_bot.at(i).Print();}


cout<<"numrand_top = "<<numrand_top<<endl;
cout<<"numrand_bot = "<<numrand_bot<<endl;

//------------------------------------------------------------------------


//Transforming Coordinates
for(int i=0;i<12;i++){
//top
abc_top.at(i)=abc_top.at(i)-Bt;
x=abc_top.at(i).Dot(UnitXt);
y=abc_top.at(i).Dot(UnitYt);
z=abc_top.at(i).Dot(UnitZt);
abc_top.at(i).SetXYZ(x,y,z);


//bottom
abc_bot.at(i)=abc_bot.at(i)-Bb;
x=abc_bot.at(i).Dot(UnitXb);
y=abc_bot.at(i).Dot(UnitYb);
z=abc_bot.at(i).Dot(UnitZb);
abc_bot.at(i).SetXYZ(x,y,z);
}

for(int i=0;i<numrand_top;i++){
rand_top.at(i)=rand_top.at(i)-Bt;
x=rand_top.at(i).Dot(UnitXt);
y=rand_top.at(i).Dot(UnitYt);
z=rand_top.at(i).Dot(UnitZt);
rand_top.at(i).SetXYZ(x,y,z);
}
for(int i=0;i<numrand_bot;i++){
rand_bot.at(i)=rand_bot.at(i)-Bb;
x=rand_bot.at(i).Dot(UnitXb);
y=rand_bot.at(i).Dot(UnitYb);
z=rand_bot.at(i).Dot(UnitZb);
rand_bot.at(i).SetXYZ(x,y,z);
}


//------------------------------------------------------------------------

//L-mark displacements
//top
L1t=(abc_top.at(0)+abc_top.at(1)+abc_top.at(2))*(one_third);
L2t=(abc_top.at(3)+abc_top.at(4)+abc_top.at(5))*(one_third);
L3t=(abc_top.at(6)+abc_top.at(7)+abc_top.at(8))*(one_third);
L4t=(abc_top.at(9)+abc_top.at(10)+abc_top.at(11))*(one_third);
//L1t.Print();
//L2t.Print();
//L3t.Print();
//L4t.Print();
cout<<"--------------------------------------------------------------------------"<<endl;

//bottom
L1b=(abc_bot.at(0)+abc_bot.at(1)+abc_bot.at(2))*(one_third);
L2b=(abc_bot.at(3)+abc_bot.at(4)+abc_bot.at(5))*(one_third);
L3b=(abc_bot.at(6)+abc_bot.at(7)+abc_bot.at(8))*(one_third);
L4b=(abc_bot.at(9)+abc_bot.at(10)+abc_bot.at(11))*(one_third);
//L1b.Print();
//L2b.Print();
//L3b.Print();
//L4b.Print();

//X-Y Projections
TVector2 L1t_xy(L1t.X(),L1t.Y());
TVector2 L2t_xy(L2t.X(),L2t.Y());
TVector2 L3t_xy(L3t.X(),L3t.Y());
TVector2 L4t_xy(L4t.X(),L4t.Y());

TVector2 L1b_xy(L1b.X(),L1b.Y());
TVector2 L2b_xy(L2b.X(),L2b.Y());
TVector2 L3b_xy(L3b.X(),L3b.Y());
TVector2 L4b_xy(L4b.X(),L4b.Y());





//Displacements of L-marks
Dis1=L1t-L2b;
Dis2=L2t-L1b;
Dis3=L3t-L4b;
Dis4=L4t-L3b;
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
//cout<<"Angle between sensors = "<<Angle((L1t_xy-L4t_xy),(L2b_xy-L3b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L4"<<endl;//Radians
//cout<<"Angle2 between sensors = "<<Angle((L1t_xy-L2t_xy),(L2b_xy-L1b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L2"<<endl;//Radians
//cout<<"Angle3 between sensors = "<<Angle((L1t_xy-L3t_xy),(L2b_xy-L4b_xy))*1000000<<" microRadians (Rotation about sensor normal) Line from L1 to L3"<<endl;//Radians
//cout<<"--------------------------------------------------------------------------"<<endl;

//Bottom sensor coordinates
TVector3 senV1,senV2,sen_unitX,sen_unit_Y,sen_unit_Z;
TVector2 Final_vect;

senV2=L2b-L3b; //Y-axis
senV1=L4b-L3b; //X-axis

sen_unit_Z=senV1.Cross(senV2);
sen_unit_Z=sen_unit_Z.Unit();
sen_unit_Y=sen_unit_Z.Cross(senV1);
sen_unit_Y=sen_unit_Y.Unit();
sen_unit_X=sen_unit_Y.Cross(sen_unit_Z);
sen_unit_X=sen_unit_X.Unit();

double x,y,z;
x=(L3t-L4t).Dot(sen_unit_X);
y=(L3t-L4t).Dot(sen_unit_Y)*(-1);



TVector2 Final_vect(x,y);
cout<<"Rotation is "<<Final_vect.Phi()<<endl;



//Plotting Range
Double_t tempmaxz_b=-1000000,tempminz_b=1000000;
for(Int_t  i=0;i<numrand_bot;i++){
if(rand_bot.at(i).Z()>=tempmaxz_b){tempmaxz_b=rand_bot.at(i).Z();}
if(rand_bot.at(i).Z()<=tempminz_b){tempminz_b=rand_bot.at(i).Z();}
}

Double_t tempmaxz_t=-1000000,tempminz_t=1000000;
for(Int_t  i=0;i<numrand_top;i++){
if(rand_top.at(i).Z()>=tempmaxz_t){tempmaxz_t=rand_top.at(i).Z();}
if(rand_top.at(i).Z()<=tempminz_t){tempminz_t=rand_top.at(i).Z();}
}


TCanvas * c1 =new TCanvas("c1","c1");
TCanvas * c2 =new TCanvas("c2","c2");
c1->cd();
//top
 TGraph2D *g_t = new TGraph2D();
g_t->SetTitle("Top Sensor");
for(i=0;i<numrand_top;i++)
{x=rand_top.at(i).X(); y= rand_top.at(i).Y(); z= rand_top.at(i).Z();
g_t->SetPoint(i,x,y,z);}
  g_t->Draw("TRI1");
  g_t->Draw("SAME P0");
  g_t->GetZaxis()->SetRangeUser(tempminz_t,tempmaxz_t);


c2->cd();
//bottom
 TGraph2D *g_b = new TGraph2D();
g_b->SetTitle("Bottom Sensor");
for(i=0;i<numrand_bot;i++)
{x=rand_bot.at(i).X(); y= rand_bot.at(i).Y(); z= rand_bot.at(i).Z();
g_b->SetPoint(i,x,y,z);}
  g_b->Draw("TRI1");
  g_b->Draw("SAME P0");
  g_b->GetZaxis()->SetRangeUser(tempminz_b,tempmaxz_b);






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








