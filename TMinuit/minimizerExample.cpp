// The function I want to minimize. Remember C++ rules,
// about function order (put it before minimizerExample)
double myFunction(double x, double y) {
double ans = ( (x-3.5)*(x-3.5) + (y+1)*(y+1) );
return ans;
}

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
result = myFunction(par[0], par[1]);
}

// Main function in minimizerExample.C
void minimizerExample() {
TFitter* minimizer = new TFitter(2);
// MAKE IT QUIET!!
{
double p1 = -1;
minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);
}
// Tell the minimizer about the function to be minimzed
minimizer->SetFCN(minuitFunction);
// Define the parameters
// arg1 – parameter number
// arg2 – parameter name
// arg3 – first guess at parameter value
// arg4 – estimated distance to minimum
// arg5, arg6 – ignore for now
minimizer->SetParameter(0,"X",1,0.5,0,0);
minimizer->SetParameter(1,"Y",1,0.5,0,0);

// Run the simplex minimizer to get close to the minimum
minimizer->ExecuteCommand("SIMPLEX",0,0);
// Run the migrad minimizer (an extended Powell's method) to improve the
// fit.
minimizer->ExecuteCommand("MIGRAD",0,0);
// Get the best fit values
double bestX = minimizer->GetParameter(0);
double bestY = minimizer->GetParameter(1);
// Get the function value at the best fit.
double minimum = myFunction(bestX, bestY);

cout<<"bestX = "<<bestX<<endl;
cout<<"bestY = "<<bestY<<endl;


// Get the best fit errors
double bestX_e = minimizer->GetParError(0);
double bestY_e = minimizer->GetParError(1);



cout<<"bestX_e = "<<bestX_e<<endl;
cout<<"bestY_e = "<<bestY_e<<endl;




/*
// Scan the X parameter to find it's uncertainty
double errX;
for (errX=0; errX<10; errX = errX + 0.001) {
minimizer->SetParameter(0,"X",errX,0,0,0);
minimizer->SetParameter(1,"Y",bestY,1,0,0);
minimizer->ExecuteCommand("MIGRAD",0,0);
double t = myFunction(minimizer->GetParameter(0),
minimizer->GetParameter(1));
if (t-minimum > 1.0) break;
}

cout<<"errX = "<<errX<<endl;
*/



}








