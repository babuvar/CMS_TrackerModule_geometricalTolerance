#include "includefiles.h"

#define _CRT_SECURE_NO_WARNINGS

#pragma comment(lib,"libGpad.lib")
#pragma comment(lib,"libPhysics.lib")
#pragma comment(lib,"libTree.lib")
#pragma comment(lib,"libMathCore.lib")
#pragma comment(lib,"libGeom.lib")
#pragma comment(lib,"libGraf.lib")
#pragma comment(lib,"libRIO.lib")
#pragma comment(lib,"libCore.lib")
#pragma comment(lib,"libHist.lib")
#pragma comment(lib,"libMinuit")


#include <TTree.h>
#include <TRandom3.h>
#include <time.h>
#include <iostream>
#include <string>
#include <TNtuple.h>
#include <math.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iomanip>
#include <TAxis.h>
#include <thread>
#include <sstream>
#include <TLorentzVector.h>
#include <TMath.h>

const int NDSSD = 3;
const double pi = TMath::Pi();

double mesDSSDpoints_X[NDSSD][4];
double mesDSSDpoints_Y[NDSSD][4];
double mesDSSDpoints_Z[NDSSD][4];

const double desr = 68.66267181;
const double desrs1 = 65.86462567; // fw side
const double desrs2 = 69.51885291; // bw side
const double des_s0eta = 0.314330265;
const double des_s1eta = 0.448693582;
const double des_eta = 0.442478133;
const double slantang = 11.9*pi / 180;


double desX0[NDSSD];
double desY0[NDSSD];
double desZ0[NDSSD];

double X0, Y0, Z0;
int i_DSSD;

//#define twoDfit

/*
minuit->mnparm(0, "delx", vstart[0], step[0], 0, 0, ierflg[0]);
minuit->mnparm(1, "dely", vstart[1], step[1], 0, 0, ierflg[0]);
minuit->mnparm(2, "delz", vstart[2], step[2], 0, 0, ierflg[0]);
minuit->mnparm(3, "deleta", vstart[3], step[3], 0, 0, ierflg[0]);
minuit->mnparm(4, "delth", vstart[4], step[4], 0, 0, ierflg[0]);
minuit->mnparm(5, "ph", vstart[5], step[5], 0, 0, ierflg[0]);
*/

double getchisq(double delx,double dely, double delz
	, double deleta, double delth, double ph){ 

	double chisq = 0;

	double Desr[4];
	double th0[4];
	double th;
	if (i_DSSD != 4){
		Desr[0] = desr;
		Desr[1] = desr;
		Desr[2] = desr;
		Desr[3] = desr;
		th0[0] = des_eta;
		th0[1] = pi - des_eta;
		th0[2] = des_eta + pi;
		th0[3] = -des_eta;
		th = delth;
	}else{
		Desr[0] = desrs1;
		Desr[1] = desrs2;
		Desr[2] = desrs2;
		Desr[3] = desrs1;
		th0[0] = des_s0eta;
		th0[1] = pi - des_s1eta;
		th0[2] = des_s1eta + pi;
		th0[3] = -des_s0eta;
		th = delth;
	}
		
	
	TVector3 vert_ori(X0  , Y0  , Z0  );
	TVector3  del_ori(delx, dely, delz);
	for (int j = 0; j < 4; j++){
		TVector3 v(Desr[j]*cos(th0[j]+deleta), Desr[j]*sin(th0[j]+deleta), 0);
		v.RotateY(delth);
		v.RotateZ(ph);
		v += del_ori;
	
		TVector3 mes(mesDSSDpoints_X[i_DSSD][j], mesDSSDpoints_Y[i_DSSD][j]
			, mesDSSDpoints_Z[i_DSSD][j]);

		mes -= vert_ori;
		TVector3 dif = v - mes;
#ifdef twoDfit
		dif.SetZ(0);
#endif

		chisq += dif.Mag2()/(0.01*0.01); //
	}

	return chisq;
}


void fcn(int& nDim, double* gout, double& result, double par[], int flg){
	result =  getchisq(par[0],par[1],par[2],par[3],par[4], par[5]);
}


using namespace std;

int main(int argc, char ** argv){



	ifstream fin("mesDSSD.txt");
	
	for (int i = 0; i < NDSSD; i++){
		for (int j = 0; j < 4; j++){
			double x, y, z;
			fin >> x >> y >> z;
			//if (!x&&!y&&!z)break;
			mesDSSDpoints_X[i][j] = x;
			mesDSSDpoints_Y[i][j] = y;
			mesDSSDpoints_Z[i][j] = z - 106.2710872;
		}
	}
	fin.close();
	fin.open("desgrab.txt");
	for (int i = 0; i < NDSSD; i++){
		double x, y, z;
		fin >> x >> y >> z;
		desX0[i] = x;
		desY0[i] = y;
		desZ0[i] = z;
	}
	
	TMinuit * minuit  = new TMinuit(6);
	minuit  -> SetFCN(fcn);
	
	//minuit  -> SetPrintLevel(-1);
	
	int ierflg[2]={0,0};
	double arglist[10];

	arglist[0] = 1;
	minuit   -> mnexcm("SET ERR",arglist,1,ierflg[0]);
	
	static double vstart[6] = {0,0,0,0,0,0};
	static double step[6] = {0.001,0.001,0.001,0.001,0.001,0.001};
	static double min[6] = {0,0,0,-0.01,0 ,-pi/2};
	static double max[6] = {0,0,0, 0.01,pi,pi/2};
	minuit->mnparm(0,"delx",   vstart[0],step[0],min[0],max[0],ierflg[0]);
	minuit->mnparm(1, "dely",  vstart[1], step[1], min[1], max[1], ierflg[0]);
	minuit->mnparm(2, "delz",  vstart[2], step[2], min[2], max[2], ierflg[0]);
	minuit->mnparm(3, "deleta",vstart[3], step[3], min[3], max[3], ierflg[0]);
	minuit->mnparm(4, "delth", vstart[4], step[4], min[4], max[4], ierflg[0]);
	minuit->mnparm(5, "ph",    vstart[5], step[5], min[5], max[5], ierflg[0]);

#ifdef twoDfit
	minuit->FixParameter(2);
	minuit->FixParameter(4);
	minuit->FixParameter(5);
#endif

	arglist[0] = 5000;
	arglist[1] = 1.0;
	
	double fit_delx[NDSSD];
	double fit_dely[NDSSD];
	double fit_delz[NDSSD];
	double fit_delth[NDSSD];
	double fit_deleta[NDSSD];
	double fit_ph[NDSSD];
	
	double fit_delx_er[NDSSD];
	double fit_dely_er[NDSSD];
	double fit_delz_er[NDSSD];
	double fit_delth_er[NDSSD];
	double fit_deleta_er[NDSSD];
	double fit_ph_er[NDSSD];

	double amin[NDSSD];
	double edm[NDSSD];
	double errdef[NDSSD];
	int nvpar[NDSSD];
	int nparx[NDSSD];
	int icstat[NDSSD];

	for(int i=0; i<NDSSD;i++){
		
		X0 = desX0[i];
		Y0 = desY0[i];
		Z0 = desZ0[i];

		for (int j = 0; j < 6; j++){
			vstart[j] = 0;
		}

		i_DSSD = i;
		if (i == 4){
			vstart[4] = slantang;
			minuit->mnparm(4, "delth", vstart[4], step[4], 0, 0, ierflg[0]);
		}

		minuit -> mnexcm("MIGRAD",arglist,2,ierflg[0]);
		//minuit -> mnstat(amin[i],edm[i],errdef[i],nvpar[i],nparx[i],icstat[i]);
		//minuit -> mnprin(1,amin[i]);
		
		minuit -> GetParameter(0, fit_delx[i]  ,fit_delx_er[i]);
		minuit -> GetParameter(1, fit_dely[i]  ,fit_dely_er[i]);
		minuit -> GetParameter(2, fit_delz[i]  ,fit_delz_er[i]);
		minuit -> GetParameter(3, fit_deleta[i], fit_deleta_er[i]);
		minuit -> GetParameter(4, fit_delth[i] , fit_delth_er[i]);
		minuit -> GetParameter(5, fit_ph[i]    , fit_ph_er[i]);
		minuit->mnstat(amin[i], edm[i], errdef[i], nvpar[i], nparx[i], icstat[i]);
	}

	/*
	g_fitdssdnum = numberofdssd-1;
	sminuit -> mnexcm("MIGRAD",arglist,2,ierflg[1]);
	sminuit -> GetParameter(0,fitx[4],fitxerr[4]);
	sminuit -> GetParameter(1,fity[4],fityerr[4]);
	sminuit -> GetParameter(2,fittheta[4],fitthetaerr[4]);
	sminuit -> mnstat(amin[4],edm[4],errdef[4],nvpar[4],nparx[4],icstat[4]);
	//sminuit -> mnprin(1,amin[4]);
	*/
	
	
	
	string outlogname = "result.txt";
	ofstream fout(outlogname.c_str());


	cout << "************************************ fit result***************************************\n\n";
	for(int i=0;i<NDSSD;i++){

		double nx = sin(fit_delth[i])*cos(fit_ph[i]);
		double ny = sin(fit_delth[i])*sin(fit_ph[i]);
		double nz = cos(fit_delth[i]);



		double tilt = asin(ny / sqrt(nz*nz + ny*ny));
		double norm = asin(nx / sqrt(nz*nz + nx*nx));
		double dtilt = 180 * tilt / pi;
		double dnorm = 180 * norm / pi;

		TVector3 unit(cos(fit_deleta[i]), sin(fit_deleta[i]), 0);
		unit.RotateY(fit_delth[i]);
		unit.RotateZ(fit_ph[i]);
		
		TVector3 u = unit.Unit();
		double proj_rotZ = atan(u.Y() / u.X());
		double dproj_rotZ = 180*proj_rotZ/pi;

		cout << "DSSD-" << i+1 << "\n";
		cout << "reduced-Chi^2 = " << amin[i]/(8-1-3) << "\n";
		cout << "del_x = " << fit_delx[i] << " +/- " << fit_delx_er[i] << "\n";
		cout << "del_y = " << fit_dely[i] << " +/- " << fit_dely_er[i] << "\n";
		cout << "del_z = " << fit_delz[i] << " +/- " << fit_delz_er[i] << "\n";
		cout << "del_eta = " << fit_deleta[i] << " +/- " << fit_deleta_er[i] << "\n";
		if (i != 4){
			cout << "delth = " << fit_delth[i] << " +/- " << fit_delth_er[i] << "\n";
		}
		else{
			cout << "th = " << fit_delth[i] << " +/- " << fit_delth_er[i] << "\n";
		}
		cout << "ph = " << fit_ph[i] << " +/- " << fit_ph_er[i] << "\n";
		cout << "tilt   angle = " << tilt << " ( " << dtilt << " deg)\n";
		cout << "normal angle = " << norm << " ( " << dnorm << " deg)\n";
		cout << "xy-projected rotation = " << proj_rotZ << " ( " << dproj_rotZ << " deg)\n\n";
		

		fout << setprecision(15);
		fout << "DSSD-" << i + 1 << "\n\n";
		fout << "reduced-Chi^2 = " << amin[i] / (12 - 6) << "\n";
		fout << "del_x = " << fit_delx[i] << " +/- " << fit_delx_er[i] << "\n";
		fout << "del_y = " << fit_dely[i] << " +/- " << fit_dely_er[i] << "\n";
		fout << "del_z = " << fit_delz[i] << " +/- " << fit_delz_er[i] << "\n";
		fout << "del_eta = " << fit_deleta[i] << " +/- " << fit_deleta_er[i] << "\n";
		if (i != 4){
			fout << "delth = " << fit_delth[i] << " +/- " << fit_delth_er[i] << "\n";
		}
		else{
			fout << "th = " << fit_delth[i] << " +/- " << fit_delth_er[i] << "\n";
		}
		fout << "ph = " << fit_ph[i] << " +/- " << fit_ph_er[i] << "\n\n";
		fout << "tilt   angle = " << tilt << " ( " << dtilt << " deg)\n";
		fout << "normal angle = " << norm << " ( " << dnorm << " deg)\n\n";
		fout << "xy-projected rotation = " << proj_rotZ << " ( " << dproj_rotZ << " deg)\n\n";
	}

	fout.close();
	double r = 999;

}